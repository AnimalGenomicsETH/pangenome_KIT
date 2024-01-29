[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![snakemaker](https://github.com/AnimalGenomicsETH/KITKAT/actions/workflows/snakemaker.yaml/badge.svg?branch=main)](https://github.com/AnimalGenomicsETH/KITKAT/actions/workflows/snakemaker.yaml)

## Taurine pangenome uncovers a segmental duplication upstream of _KIT_ associated with breed-defining depigmentation in white-headed cattle


Our project is focused on associating the white-headed coat color phenotype observed in multiple cattle breeds including 
  - Hereford
  - Simmental

and identifying a large structural variant upstream _KIT_.

With additional public short read data, covering other breeds with both the white- or color-headed phenotype, mapped directly to the graph, we could further validate that the alleles of the structural variant segregate with head depigmentation.


### Main workflow steps

The main steps include:
- Pangenome graph construction per chromosome with [pggb](https://github.com/pangenome/pggb)
- Download short reads data from public databases ([fastq_dl](https://github.com/rpetit3/fastq-dl))
- Align the short reads to the ARS-UCD1.2 bovine reference genome ([strobealign](https://github.com/ksahlin/strobealign)) and extract the ones mapped in the region of interest ([samtools](https://github.com/samtools/samtools))
- Align the short reads on the graph with vg giraffe and calculate the coverage per node in the graph ([gafpack](https://github.com/ekg/gafpack))

The coverages can then be visualised in [Bandage](https://github.com/asl/BandageNG) to see sample coverage across the pangenome.

#### Rule graph


![workflow](https://github.com/AnimalGenomicsETH/KITKAT/assets/29678761/a695abff-f6a2-4e36-8a03-b37ef8531dca)


This can be generated with 

```
snakemake -s Snakefile --rulegraph --configfile .test/config/test.yaml | dot -Tsvg > workflow.svg
```

### Citation

