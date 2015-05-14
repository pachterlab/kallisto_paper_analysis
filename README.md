# kallisto paper analysis

This repo contains all the analysis to reproduce the results in the kallisto paper.

# Preliminaries

- The annotation used is gzipped in the annotation folder, unzip it before running the analysis
- Install kallisto
- Install [snakemake](https://bitbucket.org/johanneskoester/snakemake)


- Download the GEUVADIS sample http://www.ebi.ac.uk/ena/data/view/ERR188140
- Place it into simulations/NA12716_7/NA12716_7_1.fastq.gz, simulations/NA12716_7/NA12716_7_2.fastq.gz
- Make a symlink from those files to personalized_simulation/NA12716_7/ or copy them over again

# Directories

- `personalized_simulation` - contains the allele-specific expression
  simulations
- `seqc` - contails the analysis for the SEQC dataset used in the bootstrap
  analysis
- `simulations` - contains the simulations
