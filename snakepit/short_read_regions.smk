
rule all:
    input:
        expand('bam_slices/{region}.coverage.Q{Q}.tsv',region=config['regions'],Q=config.get('Q',10))

rule samtools_view:
    input:
        '{sample}.bam'
    output:
        'bam_slices/{sample}.{region}.bam'
    threads: 2
    resources:
        mem_mb = 3000,
        walltime = '30'
    params:
        region = '6:70063972-70155205' #lambda wildcards: config['regions'][wildcards.region]
    shell:
        '''
        samtools view --threads {threads} --write-index -o {output} {input} {params.region}
        '''

rule samtools_bedcov:
    input:
        bams = expand('{sample}.bam',sample=config['samples']),
        bed = 'bam_slices/{region}.bed'
    output:
        'bam_slices/{region}.coverage.Q{Q}.tsv'
    params:
        breeds = '\\n'.join(['ID\\tbreeds'] + [f'{K}\\t{V}' for K,V in config['samples'].items()])
    threads: 1
    resources:
        mem_mb = 25000
    shell:
        '''
        paste <(echo -e "{params.breeds}") <(samtools bedcov -Q {wildcards.Q} {input.bed} {input.bams} |\
        cut -f 4- |\
        awk '{{for(i=1;i<=NF;i++)a[i][NR]=$i}}END{{for(i in a)for(j in a[i])printf"%s"(j==NR?RS:FS),a[i][j]}}' |\
        tr ' ' '\\t') > {output}
        '''
        
