# kallisto paper analysis

This repo contains all the analysis to reproduce the results in the kallisto paper.

The latest version of this analysis is accessible via [GitHub](https://github.com/pachterlab/kallisto_paper_analysis).

# Preliminaries

- Install [snakemake](https://bitbucket.org/johanneskoester/snakemake)
- Download software listed in `config.py` and place it in the appropriate directory in `software`
- Download and install `R` along with dependencies listed below (R dependencies section)

# Running the scripts

Everything is run using `snakemake`.
Start with the top level script as there are some dependencies.
The data and dependencies should get pulled automatically (assuming software has been appropriately installed).
If you wish to run everything, you can simply run the `run_all.sh` script in the top level.

# R dependencies

The version of `R` we used to generate results is 3.2.3.

### from CRAN

Install using `install.packages()`

- `cowplot`
- `devtools`
- `dplyr`
- `data.table`
- `ggplot2`
- `jsonlite`
- `reshape2`
- `rjson`
- `scales`

### from Bioconductor

First, install Bioconductor:

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
```

Then, you should be able to install packages using the `biocLite()` function.

- `biomaRt`
- `GEOquery`
- `SRAdb`

### from GitHub

Install these packages using `devtools` (`install.packages('devtools')`)

- `sleuth` v0.27.3: `devtools::install_github('pachterlab/sleuth', ref = 'v0.27.3')`
- `mamabear` v0.1: `devtools::install_github('pimentel/mamabear', ref = 'v0.1')`
