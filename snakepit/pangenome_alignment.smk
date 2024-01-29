from pathlib import PurePath

wildcard_constraints:
    EXT = r'gaf.gz|gam|bam|gaf'

rule samtools_fastq:
    input:
        bam = rules.strobealign.output,
        reference = lambda wildcards: config['references'][wildcards.reference]
    output:
        '{pangenome}/extractedFastq_{reference}/{sample}.fastq.gz'
    params:
        region = lambda wildcards: config['graphs'][wildcards.pangenome]['region'][wildcards.reference]
    threads: 2
    resources:
        mem_mb = 2000
    shell:
        '''
        samtools view --reference {input.reference} -@ {threads} -u {input.bam[0]} {params.region} 197bp 6kb | samtools collate -O -@ {threads} -u -T $TMPDIR - | samtools fastq --threads {threads} -f 3 -F 1284 -N -o {output} -s /dev/null
        '''

rule vg_autoindex:
    input:
        gfa = expand(rules.odgi_extract.output,renamed='renamed',allow_missing=True),
        fasta = lambda wildcards: expand(rules.panSN_spec.output,chromosome=convert_region_to_chromosome(wildcards))
    output:
        multiext('{pangenome}/index','.giraffe.gbz','.min','.dist')
    params:
        prefix = lambda wildcards, output: PurePath(output[1]).with_suffix('')
    threads: 4
    resources: 
        mem_mb = 4000
    shell:
        '''
        vg autoindex -t {threads} --workflow giraffe -g {input.gfa} -r <(zcat {input.fasta[0]} | head -n 2)  -p {params.prefix}
        chmod 444 {output[2]}
        '''

rule vg_giraffe:
    input:
        indexed = rules.vg_autoindex.output,
        reads = rules.samtools_fastq.output
    output:
        '{pangenome}/{sample}.{reference}.giraffe.gaf'
    threads: 4
    resources: 
        mem_mb = 1500,
        walltime = '4h'
    shell:
        '''
        vg giraffe -t {threads} -Z {input.indexed[0]} -m {input.indexed[1]} -d {input.indexed[2]} --named-coordinates -o gaf --interleaved -f {input.reads} > {output}
        '''

rule gafpack_count:
    input:
        gfa = expand(rules.odgi_extract.output,renamed='renamed',allow_missing=True),
        gaf = rules.vg_giraffe.output[0]
    output:
        '{pangenome}/{sample}.{reference}.giraffe.count'
    resources:
        walltime = '30m'
    shell:
        '''
        gafpack -g {input.gfa} -a {input.gaf} -l -c | awk '!/#/ {{print "{wildcards.sample}",$1,$2}}' > {output}
        '''

rule make_node_coverage_csv:
    input:
        expand(rules.gafpack_count.output,sample=samples,allow_missing=True)
    output:
        '{pangenome}/node_coverage.{reference}.csv.gz'
    localrule: True
    shell:
        '''
        {{ echo "sample node coverage" ; cat {input} ; }} | pigz -p 2 -c > {output}
        '''
