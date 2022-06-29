# ScPCA downstream analyses

This repository stores workflows used for performing downstream analyses on quantified single-cell and single-nuclei gene expression data available on the [Single-cell Pediatric Cancer Atlas portal](https://scpca.alexslemonade.org/).
More specifically, the repository currently contains a core workflow that performs initial pre-processing of gene expression data. 
Future development will include addition of optional workflows for extended data analyses to be applied to datasets after running the core workflow. 

The core workflow takes as input the gene expression data for each library being processed and performs the following steps: 

1. Filtering: Each library is filtered to remove any low quality cells.
Here filtering and removal of low quality cells can be performed using [`miQC::filterCells()`](https://rdrr.io/github/greenelab/miQC/man/filterCells.html) or through setting a series of manual thresholds (e.g. minimum number of UMI counts).
In addition to removing low quality cells, genes found in a low percentage of cells in a library are removed.
2. Normalization and dimensionality reduction: Cells are normalized using the [deconvolution method](https://doi.org/10.1186/s13059-016-0947-7) and reduced dimensions are calculated using both principal component analysis (PCA) and uniform manifold approximation and projection (UMAP). 
Normalized log counts and embeddings from PCA and UMAP are stored in the `SingleCellExperiment` object returned by the workflow. 
3. Clustering: Cells are assigned to cell clusters using graph-based clustering. 
The default clustering here is Louvain clustering with a nearest neighbors parameter of 10. 
Cluster assignments are stored in the `SingleCellExperiment` object returned by the workflow.

**Note** that R 4.1 is required for running our pipeline, along with Bioconductor 3.14.
Package dependencies for the analysis workflows in this repository are managed using [`renv`](https://rstudio.github.io/renv/index.html), and `renv` must be installed locally prior to running the workflow. 

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

- `results_dir`: relative path to a results directory to hold your project's output files
- `project_metadata`: relative path to your specific project metadata TSV file with the columns as follows:
    -  `sample_id`, unique ID for each piece of tissue or sample that cells were obtained from,  all libraries that were sampled from the same piece of tissue should have the same `sample_id`
    - `library_id`, unique ID used for each set of cells that has been prepped and sequenced separately
    - `filtering_method`, whose values should be one of "manual" or "miQC"
    - `filepath`, the relative path to the RDS file containing the pre-processed `SingleCellExperiment` object, each library ID should have a unique `filepath`
- `mito_file`: full path to a file containing a list of mitochondrial genes specific to the genome or transcriptome version used for alignment. 
By default, the workflow will use the mitochondrial gene list obtained from Ensembl version 104 which can be found in the `reference-files` directory. 

You can tell the config file to point to your specific project variables by running Snakemake using the `snakemake --cores 2` command and modifying the relevant parameters using the `--config` flag as in the following example:

```
snakemake --cores 2 --config data_dir="path/to/main/data/directory" \
results_dir="path/to/relevant/results/directory" \
project_metadata="project-metadata/your-project-metadata.TSV"
```

You can also use `snakemake --cores 1` to run the workflow as is, using the default values for the variables.

**Note:** To run the workflow while located outside of this directory, you will need to provide the full path to the Snakefile in this directory at the command line using `snakemake -s <full path to scpca-downstream-analyses/Snakefile>`.

## The optional genes of interest analysis pipeline (In development)

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
