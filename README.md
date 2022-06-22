# ScPCA downstream analyses

This repository stores our pipeline for ScPCA downstream analyses.
This pipeline has been created in a modular format with the intent for users to begin processing at the module or script that is most appropriate for their data.
This pipeline also includes an optional genes of interest analysis, when a genes of interest list is provided.

Note that R 4.1.2 is required for running our pipeline, along with Bioconductor 3.14.

This analysis is meant to be performed on single cell gene expression data. 
More specifically, the expected input for each sample is a RDS file containing a `SingleCellExperiment` object that has been pre-processed using the [scpca-nf workflow](https://github.com/AlexsLemonade/scpca-nf) where the data is selectively aligned using [alevin-fry](https://alevin-fry.readthedocs.io/en/latest/). For more information on the this pre-processing, please see the [ScPCA Portal docs](https://scpca.readthedocs.io/en/latest/).

Package dependencies for the analysis workflows in this repository are managed using [`renv`](https://rstudio.github.io/renv/index.html). 

If you are making any changes to the workflow, `renv::snapshot()` should be run periodically during development, as to add any packages that are used in added/modified scripts and notebooks to the `renv.lock` file.

# Running the analysis workflows

## Running the core ScPCA downstream analysis pipeline

The core downstream scPCA analysis pipeline, which includes filtering, normalization, and dimension reduction, is implemented using a Snakemake workflow.
Therefore, you will first need to install Snakemake before running the pipeline.

### 1) Install Snakemake

You can install Snakemake by following the [instructions provided in Snakemake's docs](https://snakemake.readthedocs.io/en/v7.3.8/getting_started/installation.html#installation-via-conda-mamba).

As described in the Snakemake instructions, the recommended way to install snakemake is using the conda package manager. 
After installing [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html), you can follow the below steps to install snakemake:

```
conda install -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
```

### 2) Run Snakemake

Now the environment should be all set to implement the Snakemake workflow. 
Note that we have also provided a [configuration file](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html), `config.yaml` which sets the default values for variables needed to run the workflow.

The default `config.yaml` variables that are relevant to your project include the following:

- `results_dir`: path to a results directory to hold your project's output files
- `project_metadata`: path to your specific project metadata TSV file with the columns as follows:
    -  `sample_id`, unique ID for each piece of tissue or sample that cells were obtained from,  all libraries that were sampled from the same piece of tissue should have the same `sample_id`
    - `library_id`, unique ID used for each set of cells that has been prepped and sequenced separately
    - `filtering_method`, whose values should be one of "manual" or "miQC"
    - `filepath`, the relative path to the RDS file containing the pre-processed `SingleCellExperiment` object, each library ID should have a unique `filepath`


You can tell the config file to point to your specific project variables by running Snakemake using the `snakemake --cores 2` command and modifying the relevant parameters using the `--config` flag as in the following example:

```
snakemake --cores 2 --config data_dir="path/to/main/data/directory" \
results_dir="path/to/relevant/results/directory" \
project_metadata="project-metadata/your-project-metadata.TSV"
```

You can also use `snakemake --cores 1` to run the workflow as is, using the default values for the variables.

## Running the optional genes of interest analysis pipeline

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
