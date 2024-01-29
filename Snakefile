
## split out our samples into local and downloaded
local_samples, SRA_samples = [], {}
with open(config['samples'],'r') as fin:
    for line in fin:
        name, ID, *_ = line.split()
        if ID == 'local':
            local_samples.append(name)
        elif name[0] == '#':
            continue
        else:
            SRA_samples[name] = ID

samples = local_samples + list(SRA_samples.keys())

workflow._singularity_args = f'-B $TMPDIR'

include: 'snakepit/public_downloading.smk'
include: 'snakepit/pangenome_construction.smk'
include: 'snakepit/pangenome_analysis.smk'
include: 'snakepit/pangenome_alignment.smk'
include: 'snakepit/PCA.smk'

def make_targets():
    targets = []
    for graph, config_items in config['graphs'].items():
        for reference in config['references'].keys():
            targets.append(f'{graph}/node_coverage.{reference}.csv.gz')
            targets.append(f'sample_information.{reference}.csv')
            targets.append(f'PCA/control.eigenvec')

    for graph, config_items in config['gggenes'].items():
        targets.append(f'gggenes/{graph}.repeats.tsv')

    return targets

rule all:
    input:
        make_targets()
