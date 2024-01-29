from pathlib import PurePath

rule fastq_dl:
    output:
        temp(expand('publicSamples/fastq/{accession}_R{N}.fastq.gz', N=(1,2), allow_missing=True)),
    params:
        _dir = lambda wildcards, output: PurePath(output[0]).parent
    threads: 1
    resources:
        mem_mb = 5000,
        proxy_load = 1
    envmodules:
        'eth_proxy'
    conda: 'fastq-dl'
    shell:
        '''
        fastq-dl --accession {wildcards.accession} --cpus {threads} --silent --outdir {params._dir} --prefix {wildcards.accession} --group-by-sample && [[ -s {output[0]} ]] && [[ -s {output[1]} ]]
        '''

rule fastq_dl_fixed:
    input:
        fastq = 'publicSamples/fastq/{accession}_R{N}.fastq.gz'
    output:
        temp('publicSamples/fastq/{accession}_R{N}.fixed.fastq.gz')
    threads: 4
    resources:
        mem_mb = 1500
    shell:
        '''
        pigz -dc -p {threads} {input.fastq} | awk '{{print (NR%4==1 && NF>1) ? "@"$2  : $0}}' | pigz -c -p {threads} > {output}
        '''

rule fastp_filter:
    input:
        lambda wildcards: expand(rules.fastq_dl_fixed.output,accession=SRA_samples[wildcards.sample],N=(1,2),allow_missing=True) if wildcards.sample in SRA_samples else expand(config['local_bams']+'{{sample}}_R{N}.fastq.gz',N=(1,2))
    output:
        fastq = temp(expand('publicSamples/fastq/{sample}.R{N}.fastq.gz',N=(1,2),allow_missing=True))
    params:
        min_quality = config.get('fastp',{}).get('min_quality',15),
        unqualified = config.get('fastp',{}).get('unqualified',40),
        min_length  = config.get('fastp',{}).get('min_length',15)
    threads: 4
    resources:
        mem_mb = 2500
    shell:
        '''
        fastp -q {params.min_quality} -u {params.unqualified} -g --length_required {params.min_length} --thread {threads} -i {input[0]} -o {output.fastq[0]} -I {input[1]} -O {output.fastq[1]} --json /dev/null --html /dev/null
        '''

## We don't index the reference, because it depends on read size which is variable. Do it on the fly
rule strobealign:
    input:
        fastq = expand(rules.fastp_filter.output,allow_missing=True),
        reference = lambda wildcards: config['references'][wildcards.reference]
    output:
        bam = temp(multiext('publicSamples/{sample}.{reference}.cram','','.crai')),
        dedup_stats = 'publicSamples/{sample}.{reference}.dedup.stats'
    params:
        rg = '"@RG\\tID:{sample}\\tCN:UNK\\tLB:{sample}\\tPL:illumina\\tSM:{sample}"',
    threads: 16
    resources:
        mem_mb = 2000,
        scratch = '50g',
        walltime = '4h'
    shell:
        '''
        strobealign {input.reference} {input.fastq} -t {threads} --rg-id {wildcards.sample} |\
        samtools collate -u -O -@ {threads} - |\
        samtools fixmate -m -u -@ {threads} - - |\
        samtools sort -T $TMPDIR -u -@ {threads} |\
        samtools markdup -T $TMPDIR -S -@ {threads} --write-index -f {output.dedup_stats} --reference {input.reference} --output-fmt-option version=3.1 - {output.bam[0]}
        '''

rule samtools_stats:
    input:
        rules.strobealign.output['bam']
    output:
        'coverage/{sample}.{reference}.cram.stats'
    threads: 2
    resources:
        mem_mb = 1500
    shell:
        '''
        samtools stats -d -f 3 -F 1284 --threads {threads} {input[0]} {{1..29}} > {output[0]} 
        '''

#hardcoded filenames in the awk loop
rule make_fancy_csv:
    input:
        csv = config['samples'],
        stats = expand(rules.samtools_stats.output,sample=samples,allow_missing=True),
        dedup = expand(rules.strobealign.output['dedup_stats'],sample=samples,allow_missing=True)
    output:
        'sample_information.{reference}.csv'
    localrule: True
    shell:
        '''
        echo "sample ID breed coverage read_length duplication_rate" > {output}
        while read -r a b c d
        do
          if [[ ${{a::1}} == "#" ]] ; then
            continue
          fi
          echo "$a $b $c $(awk '/total length/ {{ print $4/2489385779 }}' coverage/$a.{wildcards.reference}.cram.stats) $(awk '$1=="RL" {{L+=($2*$3);n+=$3}} END {{printf L/n}}' coverage/$a.{wildcards.reference}.cram.stats) $(awk '$1=="PAIRED:" {{P=$2}} {{if ($1=="DUPLICATE"&&$2=="PAIR:") {{ D=$3 }} }} END {{print D/P}}' publicSamples/$a.{wildcards.reference}.dedup.stats)"
        done < {input.csv} >> {output}
        '''
