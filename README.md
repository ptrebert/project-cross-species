# Project repository: epigenome-based prediction of gene expression across species


## Unique identifier

This repository: DOI tba

Publication: DOI tba

## Steps to replicate the results

0) It is strongly recommended to run the analysis on a compute cluster with a DRMAA-compatible job scheduler.
1) Create the software environment as specified in the Conda environment file (`/environment/conda_env_pplcs.yml`)
2) All pipelines of this project are implemented as [Ruffus pipelines](http://www.ruffus.org.uk/). It is recommended to install Ruffus to execute the pipelines. All (unformatted) command lines are also available in the `.ini` configuration files for each pipeline (see `/pipelines/main` and `/pipelines/prep`). Furthermore, you may want to use [PiedPiper](http://piedpiper.readthedocs.io) for convenience.
3) Download all raw data from DEEP, ENCODE and BLUEPRINT (full list contained in Additional file X of the publication). There are some helper scripts/notebooks available to bulk-download data, see `/scripts/enc_download` and `/notebooks/utils/load_sra_fastq.ipynb`; bulk downloads require some manual cleanup afterwards.
4) Build or download all necessary reference data (assembly sequence files, chromosome size files, gene models, chain files, gene ortholog annotation, gene age annotation, LOLA region databases, phyloP scores). For convenience, you can run (the relevant parts of) the pipeline in the [reference data repository](https://github.molgen.mpg.de/pebert/refdata). Note that this pipeline builds reference and annotation data for various projects, i.e. not all outputs of that pipeline are needed.
