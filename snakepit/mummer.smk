def get_files():
    targets = []
    for samples in config['query']:
        targets.extend(expand('{sample}_vs_{reference}.mumplot',sample=samples,reference=config['references']))
    print(targets)
    return targets

get_files()
rule all:
    input:
        get_files()
        
rule run_nucmer:
    input:
        ref = lambda wildcards: config["references"][wildcards.ref],
        subject = lambda wildcards: config["query"][wildcards.query]
    output:
        '{query}_vs_{ref}.delta'
    params: 
        prefix = '{query}_vs_{ref}'
    threads: 2
    resources:
        mem_mb = 2000
    shell:
        '''
        nucmer -p={params.prefix} {input.ref} {input.subject}
        samtools faidx {input.ref}
        samtools faidx {input.subject}
        awk '{{print $1,$2,"+"}}' {input.ref}.fai > {input.ref}.tsv
        awk '{{print $1,$2,"+"}}' {input.subject}.fai > {input.subject}.tsv
        '''


rule run_mummerplot:
    input:
        deltafile = rules.run_nucmer.output,
        ref = lambda wildcards: config["references"][wildcards.ref],
        subject = lambda wildcards: config["query"][wildcards.query]
    output:
        out = '{query}_vs_{ref}.mumplot'
    params: 
        settings = '-l -f',
        prefix = '{query}_vs_{ref}'
    threads: 2
    resources:
        mem_mb = 2000
    shell:
        'mummerplot {params.settings} -R {input.ref}.tsv -Q {input.subject}.tsv {input.deltafile} -p={params.prefix}'