def get_files():
    targets = []
    for samples in config['samples']:
        targets.extend(expand('{sample}.{reference}.blastn.out',sample=samples,reference=config['references']))
    # print(targets)
    return targets

get_files()
rule all:
    input:
        get_files()
        
rule run_blast:
    input:
        ref = lambda wildcards: config["references"][wildcards.reference],
        subject = lambda wildcards: config["samples"][wildcards.sample]
    output:
        out = '{sample}.{reference}.blastn.out'
    params: 6
    threads: 2
    resources:
        mem_mb = 2000
    shell:
        'blastn -query {input.ref} -subject {input.subject} -outfmt {params} -out {output.out}'