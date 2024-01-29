
rule DeepVariant:
    input:
        expand(rules.strobealign.output,sample=samples,allow_missing=True,reference='ARS')
    output:
        'PCA/{region}.Unrevised.vcf.gz'
    params:
        samples = '[' + ','.join(samples) + ']',
        config = 'config/DV.yaml' #or directly call variants and add here
    localrule: True
    shell:
        '''
        #other workflow available from here https://github.com/AnimalGenomicsETH/BSW_analysis/blob/main/snakepit/deepvariant.smk
        snakemake -s /snakepit/deepvariant.smk --configfile {params.config} \
                --config Run_name="PCA" bam_path="publicSamples/" bam_index=".csi" bam_name="{{sample}}.ARS.bam" model="WGS" samples={params.samples} \
                --profile "slurm/fullNT" --nolock {output}
        '''

rule plink_PCA:
    input:
        rules.DeepVariant.output
    output:
        multiext('PCA/{region}','.prune.in','.eigenval','.eigenvec')
    params:
        prefix = lambda wildcards, output: PurePath(output[1]).with_suffix('')
    threads: 2
    resources:
        mem_mb = 5000
    shell:
        '''
        plink2 --threads {threads} --vcf {input} --indep-pairwise 100kb 0.8 --maf 0.1 --out {params.prefix} --vcf-half-call m --snps-only --max-alleles 2
        plink2 --threads {threads} --vcf {input} --pca --out {params.prefix} --vcf-half-call m --extract {output[0]}
        '''
