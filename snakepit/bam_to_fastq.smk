def get_files():
    targets = []
    for samples in config['samples']:
        targets.extend(expand('fastq_F/{sample}.chr6.R{id}.fastq.gz',sample=config['samples'], id=["1","2"]))
    # print(targets)
    return targets

get_files()
rule all:
    input:
        get_files()

rule samtools_view:
    input:
        bam = lambda wildcards: config["samples"][wildcards.sample]
    output:
        out = 'bam/{sample}.chr6.bam'
    threads: 2
    params: 6
    resources:
        mem_mb = 2000
    shell:
        '''
        samtools view -b {input.bam} {params} > {output.out}
        '''

rule samtools_sort:
    input:
        bam = lambda wildcards: config["samples"][wildcards.sample],
        bam_chr = rules.samtools_view.output
    output:
        out = 'bam/{sample}.chr6.qsort.bam'
    threads: 2
    resources:
        mem_mb = 2000
    shell:
        '''
        samtools sort -m 30M -@ {threads} -T $TMPDIR -o {output.out} {input.bam_chr} --write-index
        '''

rule samtools_fastq:
    input:
        bam = rules.samtools_sort.output
    output:
        fq1 = 'fastq_F/{sample}.chr6.R1.fastq.gz',
        fq2 = 'fastq_F/{sample}.chr6.R2.fastq.gz'
    threads: 2
    resources:
        mem_mb = 2000
    shell:
        'samtools fastq --threads {threads} -1 {output.fq1} -2 {output.fq2} -F 0 {input.bam}'
