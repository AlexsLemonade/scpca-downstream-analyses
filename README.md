# scPCA downstream analyses

This repository is intended to store our pipeline for scPCA downstream analyses.
This pipeline has been created in a modular format with the intent for users to begin processing at the module or script that is most appropriate for their data.
This pipeline also includes an optional genes of interest analysis, when a genes of interest list is provided.

Note that R 4.1 is required for running our pipeling, along with Bioconductor 3.14.

## How to run the analysis

### How to run the main scPCA downstream analysis pipeline

The core downstream scPCA analysis pipeline, which includes filtering, normalization, and dimension reduction, is implemented using a Snakemake workflow.
Therefore, you will first need to install Snakemake before running the pipeline.

#### 1) Install Snakemake
You can install Snakemake by following the [instructions provided in Snakemake's docs](https://snakemake.readthedocs.io/en/v7.3.8/getting_started/installation.html#installation-via-conda-mamba).

To reiterate the instructions there, you will want to enter the following commands to install Snakemake via conda and activate Snakemake:

```
conda install -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
```

#### 2) Restore `renv.lock` packages

Packages that are required for the pipeline are included in the `renv.lock` file as they are installed and used. 

When checking out a branch or implementing things in your local enviroment for the first time, `renv::restore()` can be used to sync local package installations with those required to run the pipeline. To keep this file up to date, `renv::snapshot()` should be run periodically during development, as to add any packages that are used in added/modified scripts and notebooks to the `renv.lock` file.

#### 3) Run Snakemake

Now the environment should be all set to implement the Snakemake workflow. Note that you will want to modify the `config.yaml` file to point to the files relevant to your project.

The default `config.yaml` variables that will need to be relevant to your project include:

- `data_dir`: path to the main data directory that holds your samples' files
- `results_dir`: path to a results directory to hold your project's output files
- `project_metadata`: path to your specific project metadata TSV file with the columns `sample_id`, `library_id`, and `filtering_method`, where `filtering_method` can be "manual" or "miQC".


You can tell the config file to point to your specific project variables by running Snakemake using the `snakemake --cores 2` command and modifying the relevant parameters using the `--config` flag as in the following example:

```
snakemake --cores 2 --config data_dir="path/to/main/data/directory" \
results_dir="path/to/relevant/results/directory" \
project_metadata="project-metadata/your-project-metadata.TSV"
```

### How to run the optional genes of interest analysis pipeline

There is an optional genes of interest analysis pipeline in the `optional-goi-analysis` subdirectory of this repository.

To run this analysis, you can run the following command from the main directory:

```
bash optional-goi-analysis/run-provided-goi-analysis.sh \
 --output_dir "path/to/output-results" \
 --sample_name "sample"  \
 --sample_matrix "path/to/sample/matrix" \
 --sample_metadata "path/to/sample/metadata" \
 --goi_list "path/to/goi-list"
```

Where `goi_list` is the path to the genes of interest TSV file relevant to your dataset.
