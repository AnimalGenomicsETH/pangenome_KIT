def get_all_samples():
    return list(filter(None,[asm for sample,asms in config['assemblies'].items() for asm in asms if asm !='~']))

def panSN_naming():
    mapping = []
    for n, (sample_ID, sample_assemblies) in enumerate(config['assemblies'].items()):
        for haplotype, path in enumerate(sample_assemblies,len(sample_assemblies)!=1):
            if n == 0: #is reference
                mapping.append(f'"{sample_ID}#{haplotype}#  {path}"')
            elif path: #is generic/haplotype, #NOTE this is still unclear on vg path metadata model, as potentially needs a further #phaseblock after the chromosome, which we can't prefix here
                mapping.append(f'"{sample_ID}#{haplotype}#  {path}"')
    return mapping

localrules: pangenome_tree

rule panSN_spec:
    input:
        lambda wildcards: get_all_samples()
    output:
        multiext('graphs/{chromosome}.fa.gz','','.fai','.gzi')
    params:
        mapping = lambda wildcards: panSN_naming()
    threads: 4
    resources:
        mem_mb = 1500,
        walltime = '30m'
    shell:
        '''
        for i in {params.mapping}
        do
          set -- $i
          samtools faidx --length 0 $2 {wildcards.chromosome} | fastix -p $1 -
        done | bgzip -@ {threads} -c > {output[0]}
        samtools faidx {output[0]}
        '''

rule mash_triangle:
    input:
        rules.panSN_spec.output[0]
    output:
        'graphs/{chromosome}.mash'
    threads: 4
    resources:
        mem_mb = 1000,
        walltime = '1h'
    shell:
        '''
        mash triangle -s 10000 -k 25 -p {threads} -i {input} | awk 'NR>1' > {output}
        '''

rule mash_sketch:
    input:
        lambda wildcards: config['assemblies'][wildcards.sample][max(int(wildcards.haplotype)-1,0)]
    output:
        'graphs/{sample}.{haplotype}.msh'
    threads: 4
    resources:
        mem_mb = 1500
    shell:
        '''
        samtools faidx {input} {{1..29}} | mash sketch -s 10000 -k 25 -p {threads} -I {wildcards.sample}.{wildcards.haplotype} -o {output} -
        '''

def find_assembly_haplotypes():
    panSN = panSN_naming()
    samples, haplotypes = [], []
    for i in panSN:
        sample, haplotype, *_ = i.split('#')
        samples.append(sample[1:])
        haplotypes.append(int(haplotype))
    return {'sample':samples,'haplotype':haplotypes}

rule mash_big_triangle:
    input:
        expand(rules.mash_sketch.output,zip,**find_assembly_haplotypes())
    output:
        'mash.big_triangle'
    localrule: True
    shell:
        '''
        mash triangle -s 10000 -k 25 {input} | tail -n +2 > {output}
        '''

import numpy as np
def read_mash_triangle(mash_triangle,estimate_divegence=False):
    names, vals = [], []
    with open(mash_triangle,'r') as fin:
        for i,line in enumerate(fin):
            parts = line.rstrip().split()
            names.append(parts[0])
            vals.append(parts[1:]+[0])
    Q = np.asarray([np.pad(a, (0, len(vals) - len(a)), 'constant', constant_values=0) for a in vals],dtype=float)
    if estimate_divegence:
        return round((1-Q.max()*2.5)*100,1) #adjust max divergence by 2.5x factor
    return names, (Q+Q.T)

rule pangenome_tree:
    input:
        mash = rules.mash_triangle.output
    output:
        'graphs/{chromosome}.tree.svg'
    # localrule: True
    run:
        from scipy.spatial.distance import squareform
        from scipy.cluster import hierarchy
        import matplotlib.pyplot as plt
        names, dists = read_mash_triangle(input.mash[0])
        Z = hierarchy.linkage(squareform(dists),method='average',optimal_ordering=True)
        f, ax = plt.subplots()
        dn = hierarchy.dendrogram(Z,labels=names,orientation='left',get_leaves=True,ax=ax)
        f.savefig(output[0])

#Do we want to implement split_approx_mappings_in_chunks.py?
#implemented in nfcore-pangenome

rule pggb_construct:
    input:
        fasta = rules.panSN_spec.output,
        mash = rules.mash_triangle.output
    output:
        gfa = 'graphs/{chromosome}.pggb.gfa',
        og = 'graphs/{chromosome}.pggb.og'
    threads: 12
    resources:
        mem_mb = 3500,
        walltime = '24h',
        scratch = '50G'
    params:
        _dir = lambda wildcards, output: Path(output[0]).parent,
        divergence = lambda wildcards, input: read_mash_triangle(input.mash[0],True),
        n_haplotypes = lambda wildcards: len(get_all_samples()),
        min_match = 31, #adjusted from HPRC pipeline
        segment_length = 75000
    shell:
        '''
        pggb -i {input.fasta[0]} -t {threads} \
        -s {params.segment_length} -p {params.divergence} -n {params.n_haplotypes} \
        -k {params.min_match} \
        --skip-viz --temp-dir $TMPDIR \
        -o {params._dir}

        mv {params._dir}/{wildcards.chromosome}.*.smooth.final.gfa {output.gfa}
        mv {params._dir}/{wildcards.chromosome}.*.smooth.final.og {output.og}
        '''
