def get_files():
    targets = []
    for read_type, samples in config['samples'].items():
        targets.extend(expand('{sample}.{reference}.{read_type}.minimap.bam',read_type=read_type,sample=samples,reference=config['references']))
    
    # print(targets)
    return targets

get_files()
rule all:
    input:
        get_files()

presets = {'ONT':'map-ont','HiFi':'map-hifi','CLR':'map-pb','SR_w':'sr','SR':'sr'}

rule map_minimap2:
    input:
        ref = lambda wildcards: config["references"][wildcards.reference],
        fastq = lambda wildcards: config["samples"][wildcards.read][wildcards.sample]
    output:
        multiext('{sample}.{reference}.{read}.minimap.bam','','.csi')
    threads: 12
    params: 
        preset = lambda wildcards: presets[wildcards.read]
    resources:
        mem_mb = 6000,
        walltime = "4:00",
        disk_scratch = 50
    shell:
        'minimap2 -ax {params.preset} -t {threads} {input.ref} {input.fastq} | samtools sort - -m 3000M -@ {threads} -T $TMPDIR -o {output[0]} --write-index'