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
Louvain clustering is performed using the [`bluster::NNGraphParam()`](https://rdrr.io/github/LTLA/bluster/man/NNGraphParam-class.html) function using the default nearest neighbors parameter of 10. 
The type of graph based clustering and number of nearest neighbors parameter can be altered if desired.
Cluster assignments are stored in the `SingleCellExperiment` object returned by the workflow.

**Note** that R 4.1 is required for running our pipeline, along with Bioconductor 3.14.
Package dependencies for the analysis workflows in this repository are managed using [`renv`](https://rstudio.github.io/renv/index.html), and `renv` must be installed locally prior to running the workflow. 

# The core downstream analyses workflow

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Input data format](#input-data-format)
- [How to install the core downstream analyses workflow](#how-to-install-the-core-downstream-analyses-workflow)
- [Metadata file format](#metadata-file-format)
- [Running the workflow](#running-the-workflow)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Input data format

The expected input for our core single-cell downstream analysis pipeline is a [`SingleCellExperiment` object](https://rdrr.io/bioc/SingleCellExperiment/man/SingleCellExperiment.html) that has been stored as a RDS file.
This`SingleCellExperiment` object should contain non-normalized gene expression data with barcodes as the column names and gene identifiers as the row names.
All barcodes included in the `SingleCellExperiment` object should correspond to droplets likely to contain cells and should not contain empty droplets (e.g. droplets with FDR < 0.01 calculated with [`DropletUtils::emptyDropsCellRanger()`](https://rdrr.io/github/MarioniLab/DropletUtils/man/emptyDropsCellRanger.html).
The full path to each individual RDS file should be defined in the project metadata described in the following "How to run the pipeline" section.

The pipeline in this repository is setup to process data available on the [Single-cell Pediatric Cancer Atlas portal](https://scpca.alexslemonade.org/) and output from the [scpca-nf workflow](https://github.com/AlexsLemonade/scpca-nf) where single-cell/single-nuclei gene expression data is mapped and quantified using [alevin-fry](https://alevin-fry.readthedocs.io/en/latest/).
For more information on the this pre-processing, please see the [ScPCA Portal docs](https://scpca.readthedocs.io/en/latest/).
Note however that the input for this pipeline is **not required** to be scpca-nf processed output.

## How to install the core downstream analyses workflow

### 1) Clone the repository

First you will want to clone the [`scpca-downstream-analyses` repository](https://github.com/AlexsLemonade/scpca-downstream-analyses) from GitHub.

You can do this by navigating to the `Code` button at the top of the repository page and copying the URL.
Next, open a local `Terminal` window and use `cd` to navigate to the desired local directory for storing the repository.
We recommend cloning this repository into a separate folder specifically for git repositories.

You can then implement the following command to clone the repository:

`git clone https://github.com/AlexsLemonade/scpca-downstream-analyses.git`

More instructions on cloning a GitHub repository can be found [here](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository).

Once the repository is successfully cloned, a folder named `scpca-downstream-analyses` containing a local copy of the contents of the repository will be created.

### 2) Install Snakemake

The core downstream single-cell analysis pipeline, which includes filtering, normalization, dimensionality reduction, and clustering is implemented using a Snakemake workflow.
Therefore, you will also need to install Snakemake before running the pipeline.

You can install Snakemake by following the [instructions provided in Snakemake's docs](https://snakemake.readthedocs.io/en/v7.3.8/getting_started/installation.html#installation-via-conda-mamba).

As described in the Snakemake instructions, the recommended way to install snakemake is using the conda package manager. 
After installing [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html), you can follow the below steps to install snakemake:

```
conda install -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
```

Note that `pandoc` must also be installed and in your path to successfully run the `Snakefile`.
Therefore, if you are using snakemake in a separate environment then pandoc must also be in that environment.
See [pandoc's installation instructions](https://pandoc.org/installing.html) for more information.

## Metadata file format

Now the environment should be all set to implement the Snakemake workflow.
Before running the workflow, you will need to create a project metadata file as a tab-separated value (TSV) file that contains the relevant data for your input files needed to run the workflow.
The file should contain the following columns: 

- `sample_id`, unique ID for each piece of tissue or sample that cells were obtained from,  all libraries that were sampled from the same piece of tissue should have the same `sample_id`
- `library_id`, unique ID used for each set of cells that has been prepped and sequenced separately
- `filtering_method`, whose values should be one of "manual" or "miQC"
- `filepath`, the full path to the RDS file containing the pre-processed `SingleCellExperiment` object, each library ID should have a unique `filepath`

## Running the workflow

We have provided a [configuration file](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html), `config.yaml` which sets the defaults for all parameters needed to run the workflow.

There are a set of parameters included in the `config.yaml` file that will need to be specified when running the workflow. 
These parameters are specific to the project or dataset being processed.
These include the following parameters:

- `results_dir`: relative path to the directory where output files from running the core workflow will be stored
- `project_metadata`: relative path to your specific project metadata TSV file
- `mito_file`: full path to a file containing a list of mitochondrial genes specific to the genome or transcriptome version used for alignment. 
By default, the workflow will use the mitochondrial gene list obtained from Ensembl version 104 which can be found in the `reference-files` directory.

You can tell the config file to point to your specific project variables by running Snakemake using the `snakemake --cores 2` command and modifying the relevant parameters using the `--config` flag as in the following example:

```
snakemake --cores 2 --config results_dir="path/to/relevant/results/directory" \
project_metadata="project-metadata/your-project-metadata.TSV" \
mito_file="reference-files/your-mito-file.txt"
```

You can also modify the relevant parameters by manually updating the `config.yaml` file using a text editor of your choice.
The required parameters mentioned above can be found under the [`Project-specific parameters` section](https://github.com/AlexsLemonade/scpca-downstream-analyses/blob/9e82725fe12bcfb6179158aa03e8674f59a9a259/config.yaml#L3) of the config file, while the remaining parameters that can be optionally modified are found under the [`Processing parameters` section](https://github.com/AlexsLemonade/scpca-downstream-analyses/blob/9e82725fe12bcfb6179158aa03e8674f59a9a259/config.yaml#L11).

**Note:** To run the workflow while located outside of this directory, you will need to provide the full path to the Snakefile in this directory at the command line using the `-s` flag as in the following example: 

```
snakemake --cores 2 \
 -s "path to snakemake file" \
 --config project_metadata="path to project metadata"
```

**Note:** For Data Lab staff members working on development, the default `config.yaml` file as well as the project metadata file have been set up to use the shared data present on the Rstudio server at `/shared/scpca/gawad_data/scpca_processed_output`.
The workflow can still be run from inside the directory that holds this repository without modifying any parameters, just by specifying the number of cores as in `snakemake --cores 2`.

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
