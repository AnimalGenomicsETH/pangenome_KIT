rule all:
    input:
        expand('forcing_parameters/{samples}.{presets}.{references}.mm2.bam',**config)#sample=config['samples'],preset=config['presets'],reference=config['references'])

rule minimap2_align:
    input:
        reference = lambda wildcards: config['references'][wildcards.reference],
        sample = lambda wildcards: config['samples'][wildcards.sample]
    output:
        'forcing_parameters/{sample}.{preset}.{reference}.mm2.bam'
    params:
        preset = lambda wildcards: config['presets'][wildcards.preset]
    threads: 2
    resources:
        mem_mb = 5000,
        disk_scratch = 50
    shell:
        '''
        minimap2 -a {params.preset} -t {threads} {input.reference} {input.sample} | samtools sort - -m 3000M -@ 4 -T $TMPDIR --write-index -o {output}
        '''
