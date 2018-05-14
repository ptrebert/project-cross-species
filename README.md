# Project repository: epigenome-based prediction of gene expression across species


## Unique identifier

This repository: [DOI: 10.17617/1.69](https://doi.org/10.17617/1.69)

Preprint: tba

Publication: tba

## Recreating figures

All figures in the publication can be recreated by running the relevant Jupyter notebook in `notebooks/plotting` with the following exceptions:

- figure 1: manually created
- figure S8A/B: downloaded from enrichR web service; manually added figure key

### Notebooks for figures

- figure 2, S2: `plot_transfer_epigenome_errbar.ipynb`
- figure 3: `plot_sigcorr_matrix.ipynb`
- figure S3A/B: `plot_region_signal.ipynb`
- figure 4A-D: `plot_lola_marked_bubble_inv.ipynb`
- figure 5, S4: `plot_exp_ortho_clust.ipynb`
- figure 6A/B, 7A/B: `plot_model_perf_comp_roc-joint.ipynb`
- figure 8A/B: `plot_model_perf_comp_auc-bars.ipynb`
- figure 9, S5, S6A/B: `plot_model_perfcomp_tss.ipynb`
- figure 10, S7: `plot_model_perfcomp_prob.ipynb`
- figure 11, S9: `plot_unaln_ortho_gene_bars.ipynb`
- figure 12, S10A/B: `plot_gene_ages_relbar.ipynb`

## Steps to replicate all results

0) It is strongly recommended to run the analysis on a compute cluster with a DRMAA-compatible job scheduler. The expected runtime for the entire project is 3-5 days (assuming that a handful of servers is available to run the pipeline in parallel).
1) Create the software environment as specified in the Conda environment file (`/environment/conda_env_pplcs.yml`)
2) All pipelines of this project are implemented as [Ruffus pipelines](http://www.ruffus.org.uk/). It is recommended to install Ruffus to execute the pipelines. All (unformatted) command lines are also available in the `.ini` configuration files for each pipeline (see `/pipelines/main` and `/pipelines/prep`). Furthermore, you may want to use [PiedPiper](http://piedpiper.readthedocs.io) for convenience.
3) Download all raw data from DEEP, ENCODE and BLUEPRINT (full list contained in Additional file X of the publication). There are some helper scripts/notebooks available to bulk-download data, see `/scripts/enc_download` and `/notebooks/utils/load_sra_fastq.ipynb`; bulk downloads require some manual cleanup afterwards.
4) Build or download all necessary reference data (assembly sequence files, chromosome size files, gene models, chain files, gene ortholog annotation, gene age annotation, LOLA region databases, phyloP scores). For convenience, you can run (the relevant parts of) the pipeline in the [reference data repository](https://github.molgen.mpg.de/pebert/refdata). Note that the reference data pipeline builds annotation data for various projects, i.e., not all outputs of that pipeline are needed for this project.
5) Run complete pipeline `/pipelines/prep/ppl_prep_epigenomes`
6) Run complete pipeline `/pipelines/prep/ppl_prep_transcriptomes`
7) Run pipeline `/pipelines/main/ppl_cs_main` up to the task `summ_perf_status`. This task aggregates the output of all model predictions into one single file that is the input for most of the notebooks generating the plots.
