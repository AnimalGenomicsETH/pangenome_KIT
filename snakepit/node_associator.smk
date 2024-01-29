
rule all:
    input:
        'association.bed'
import glob
rule mgutils_merge:
    input:
        glob.glob('/cattle/*.bed')
        
    output:
        'association.vcf'
    shell:
        '''
        paste {input} | mgutils.js merge -s <(ls {input} | xargs -L1 -I{{}} basename {{}} .bed) - | mgutils.js merge2vcf | bcftools query -H -f '%INFO/VS %INFO/VE %INFO/AWALK [%GT0 ]\\n' > {output}
        '''


from collections import defaultdict
import regex
import numpy as np
import matplotlib.pyplot as plt

arrows = regex.compile('>|<')

def get_allele_nodes(start,end,nodes):
    alleles = defaultdict(list)
    for i,node in enumerate(nodes.split(',')):
        if node == '*':
            alleles[i] = [f'{start}-{end}']
        else:
            for n in arrows.split(node)[1:]:
                alleles[i].append(n)

    return dict(alleles)

def generate_matrix(fname):
    MATRIX = {}
    samples = None
    for line in open(fname):
        if not samples:
            samples = [i.split(']')[1].split(':')[0] for i in line.rstrip().split()[4:]]
            continue
        start, end, nodes, *GTs = line.rstrip().split()
        alleles = get_allele_nodes(start,end,nodes)
        for i, GT in enumerate(GTs):
            if GT == '.':
                continue
            GT = int(GT)
            for node in alleles[int(GT)]:
                if node not in MATRIX:
                    MATRIX[node] = [False]*len(GTs)

                MATRIX[node][i] = True
    IDs = np.array(list(MATRIX.keys()))
    vals = np.array(list(MATRIX.values()))
    return samples,IDs,vals

rule generate_assoication_bed:
    input:
        'association.vcf'
    output:
        'association.bed'
    run:
        samples,nodes,matrix = generate_matrix(input[0])
        #!! HARDCODED SAMPLE ORDER !!
        white_headed = np.array([False, True, False, False, False, False, False, False, False, False, False, False, False, False, False, True, False, True, False])
        #!! HARDCODED SAMPLE ORDER !!
        _sum = np.isclose(matrix,white_headed,equal_nan=True).sum(axis=1)
        cmap = plt.get_cmap('cividis',len(samples))
        with open(output[0],'w') as fout:
            for i,n in enumerate(nodes):
                if '-' in n:
                    for e in n.split('-')[1:]:
                        fout.write(f'{e[2:]}\t0\t1000\tX\t{_sum[i]}\t+\t0\t1000\t255,40,77\n')
                else:
                    col = ','.join([f'{255*i:.0f}' for i in cmap(_sum[i]/len(samples))[:3]])
                    fout.write(f'{n[1:]}\t0\t100000\tX\t{_sum[i]}\t+\t0\t0\t{col}\n')
