# ScPCA downstream analyses

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Core analysis overview](#core-analysis-overview)
- [Quick Start Guide](#quick-start-guide)
- [1. How to install the core downstream analyses workflow](#1-how-to-install-the-core-downstream-analyses-workflow)
  - [a) Clone the repository](#a-clone-the-repository)
  - [b) Install Snakemake](#b-install-snakemake)
    - [What is Snakemake?](#what-is-snakemake)
    - [How to install Snakemake](#how-to-install-snakemake)
  - [c) Additional dependencies](#c-additional-dependencies)
    - [Snakemake/conda installation](#snakemakeconda-installation)
- [2. Verify input data format](#2-verify-input-data-format)
- [3. Create metadata file](#3-create-metadata-file)
- [4. Configure config file](#4-configure-config-file)
- [5. Running the workflow](#5-running-the-workflow)
- [6. Expected output](#6-expected-output)
  - [What to expect in the output `SingleCellExperiment` object](#what-to-expect-in-the-output-singlecellexperiment-object)
- [Additional analysis modules](#additional-analysis-modules)
  - [Clustering analysis](#clustering-analysis)
  - [Genes of interest analysis](#genes-of-interest-analysis)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Core analysis overview

This repository stores workflows used for performing downstream analyses on quantified single-cell and single-nuclei gene expression data available on the [Single-cell Pediatric Cancer Atlas portal](https://scpca.alexslemonade.org/).
More specifically, the repository currently contains a core workflow that performs initial pre-processing of gene expression data.
Future development will include addition of optional workflows for extended data analyses to be applied to datasets after running the core workflow.

The core workflow takes as input the gene expression data for each library being processed and performs the following steps:

1. [Filtering](./additional-docs/processing-information.md#filtering-low-quality-cells): Each library is filtered to remove any low quality cells.
Here filtering and removal of low quality cells can be performed using [`miQC::filterCells()`](https://rdrr.io/github/greenelab/miQC/man/filterCells.html) or through setting a series of manual thresholds (e.g., minimum number of UMI counts).
In addition to removing low quality cells, genes found in a low percentage of cells in a library can optionally be removed.
2. [Normalization](./additional-docs/processing-information.md#normalization) and [dimensionality reduction](./processing-information.md#dimensionality-reduction): Cells are normalized using the [deconvolution method from Lun, Bach, and Marioni (2016)](https://doi.org/10.1186/s13059-016-0947-7) and reduced dimensions are calculated using both principal component analysis (PCA) and uniform manifold approximation and projection (UMAP).
Normalized log counts and embeddings from PCA and UMAP are stored in the `SingleCellExperiment` object returned by the workflow.
3. [Clustering](./additional-docs/processing-information.md#clustering): Cells are assigned to cell clusters using graph-based clustering.
By default, louvain clustering is performed using the [`bluster::NNGraphParam()`](https://rdrr.io/github/LTLA/bluster/man/NNGraphParam-class.html) function using the default nearest neighbors parameter of 10.
Alternatively, walktrap graph-based clustering can be specified, and the number of nearest neighbors parameter can be altered if desired.
Cluster assignments are stored in the `SingleCellExperiment` object returned by the workflow.

You can read more details about the individual steps of the workflow in the processing documentation linked below:

|[View Processing Information Documentation](./additional-docs/processing-information.md)|
|---|

## Quick Start Guide

To run the core analysis workflow you will want to implement the following steps in order:

1. Clone the repository and install Snakemake using the [instructions provided in the Snakemake docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba).
2. [Install the packages and dependencies](#c-additional-dependencies) that are required to run the workflow.
3. Ensure that the input single-cell gene expression data are stored as `SingleCellExperiment` objects in RDS files (see more on this in the ["Input data format" section](#2-input-data-format)).
The workflow can directly take as input the `filtered` RDS files downloaded from the [Single-cell Pediatric Cancer Atlas portal](https://scpca.alexslemonade.org/) or the output from the [scpca-nf workflow](https://github.com/AlexsLemonade/scpca-nf), a workflow that can be used to quantify your own single-cell/single-nuclei gene expression data. If working with data from the ScPCA portal, see more information on preparing that data to run the core workflow [here](./additional-docs/working-with-scpca-portal-data.md).
4. [Create a metadata tab-separated value (TSV) file](#3-metadata-file-format) that defines the sample id, library id, and filepath associated with the pre-processed `SingleCellExperiment` files to be used as input for the workflow.
5. [Configure the config file](#4-configure-config-file) to adjust the `results_dir` and `project_metadata` parameters to point to the full path to your desired results directory and project metadata file that you created in step 4.
6. Open terminal to run the workflow using the following snakemake command:

```
snakemake --cores 2 --use-conda
```

**Note** that R 4.1 is required for running our pipeline, along with Bioconductor 3.14.
Package dependencies for the analysis workflows in this repository are managed using [`renv`](https://rstudio.github.io/renv/index.html), and `renv` must be installed locally prior to running the workflow.
If you are using conda, dependencies can be installed as [part of the setup mentioned in step 2 above](#snakemakeconda-installation).

If you did not install dependencies with [conda via snakemake](#snakemakeconda-installation) in step 2, you will need to remove the `--use-conda` flag from the command above.
See the section on [running the workflow](#4-running-the-workflow) for more information.

**Output Files**
There are two expected output files thay will be associated with each provided `SingleCellExperiment` object and `library_id`:

- A processed `SingleCellExperiment` object containing normalized data and clustering results
- A summary HTML report detailing:
    - Filtering of low quality cells
    - Dimensionality reduction
    - Clustering that was performed within the workflow

See the [expected output section](#5-expected-output) for more information on these output files.
You can also download a ZIP file with an example of the output from running the core workflow, including the summary HTML report and processed `SingleCellExperiment` objects stored as RDS files [here](https://scpca-references.s3.amazonaws.com/example-data/scpca-downstream-analyses/core_example_results.zip).

## 1. How to install the core downstream analyses workflow

### a) Clone the repository

First you will want to clone the [`scpca-downstream-analyses` repository](https://github.com/AlexsLemonade/scpca-downstream-analyses) from GitHub.

We recommend cloning this repository into a separate folder specifically for git repositories.
Open a local `Terminal` window and use `cd` to navigate to the desired local directory for storing the repository, and run the following command to clone the repository:

`git clone https://github.com/AlexsLemonade/scpca-downstream-analyses.git`

More instructions on cloning a GitHub repository can be found [here](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository).

Once the repository is successfully cloned, a folder named `scpca-downstream-analyses` containing a local copy of the contents of the repository will be created.

### b) Install Snakemake

#### What is Snakemake?

The core downstream single-cell analysis pipeline, which includes filtering, normalization, dimensionality reduction, and clustering is implemented using the [Snakemake](https://snakemake.github.io/) workflow manager.
With Snakemake, users can run the core downstream analysis pipeline with a single call, rather than having to run each [core downstream analysis script](core-analysis/) on its own.
Therefore, you will also need to install Snakemake before running the pipeline.

#### How to install Snakemake

Note that the **minimum** version of Snakemake you will need to have installed is version **7.20.0**.

You can install Snakemake by following the [instructions provided in Snakemake's docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba).

Snakemake recommends installing it using the conda package manager.
Here are the instructions to [install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).
We recommend the Miniconda installation.

After installing conda, you can follow the steps below to set up the bioconda and conda-forge channels and install Snakemake in an isolated environment:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda install -n base mamba
conda activate base
mamba create -n snakemake snakemake
conda activate snakemake
```

### c) Additional dependencies

To run the Snakemake workflow, you will need to have R version 4.2 installed, as well as the `optparse` and `renv` packages and pandoc.
This can be done independently if desired, but we recommend using Snakemake's conda integration to set up the R environment and all dependencies that the workflow will use, as described below.


#### Snakemake/conda installation

Snakemake can handle the required dependencies by creating its own conda environments, which we have provided as an option.
To create the necessary environment, which includes an isolated version of R, pandoc, and all required dependencies, run the following command from the base of the repository:

```
bash setup_envs.sh
```

This script will use Snakemake to install all necessary components for the workflow in an isolated environment.
If you are on an Apple Silicon (M1/M2/Arm) Mac, this should properly handle setting up R to use an Intel-based build for compatibility with Bioconductor packages.

To use the environment you have just created, you will need to run Snakemake with the `--use-conda` flag each time.

If you would like to perform installation without the conda environments as described above, see the [independent installation instructions document](https://github.com/AlexsLemonade/scpca-downstream-analyses/blob/main/additional-docs/independent-installation-instructions.md).

## 2. Verify input data format

The expected input for our core single-cell downstream analysis pipeline is a [`SingleCellExperiment` object](https://rdrr.io/bioc/SingleCellExperiment/man/SingleCellExperiment.html) that has been stored as a RDS file.
This `SingleCellExperiment` object should contain non-normalized gene expression data with barcodes as the column names and gene identifiers as the row names.
All barcodes included in the `SingleCellExperiment` object should correspond to droplets likely to contain cells and should not contain empty droplets (e.g., droplets with FDR < 0.01 calculated with [`DropletUtils::emptyDropsCellRanger()`](https://rdrr.io/github/MarioniLab/DropletUtils/man/emptyDropsCellRanger.html)).
The full path to each individual RDS file should be defined in the project metadata described in the following "How to run the pipeline" section.

The pipeline in this repository is setup to process data available on the [Single-cell Pediatric Cancer Atlas portal](https://scpca.alexslemonade.org/) and output from the [scpca-nf workflow](https://github.com/AlexsLemonade/scpca-nf) where single-cell/single-nuclei gene expression data is mapped and quantified using [alevin-fry](https://alevin-fry.readthedocs.io/en/latest/).
For more information on the this pre-processing, please see the [ScPCA Portal docs](https://scpca.readthedocs.io/en/latest/).
Note however that the input for this pipeline is **not required** to be scpca-nf processed output.

If you are working with data downloaded from the ScPCA portal, see our guide on preparing that data to run the core workflow [here](./additional-docs/working-with-scpca-portal-data.md).

## 3. Create metadata file

Now the environment should be all set to implement the Snakemake workflow.
Before running the workflow, you will need to create a project metadata file as a tab-separated value (TSV) file that contains the relevant data for your input files needed to run the workflow.
The file should contain the following columns:

- `sample_id`, unique ID for each piece of tissue or sample that cells were obtained from,  all libraries that were sampled from the same piece of tissue should have the same `sample_id`.
- `library_id`, unique ID used for each set of cells that has been prepped and sequenced separately.
- `filepath`, the full path to the RDS file containing the pre-processed `SingleCellExperiment` object.
Each library ID should have a unique `filepath`.

|[View Example Metadata File](https://github.com/AlexsLemonade/scpca-downstream-analyses/blob/main/example-data/project-metadata/example-library-metadata.tsv)|
|---|

## 4. Configure config file

We have provided an [example snakemake configuration file](config/config.yaml), `config/config.yaml` which defines all parameters needed to run the workflow.
Learn more about snakemake configuration files [here](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).

You can modify the relevant parameters by manually updating the `config/config.yaml` file using a text editor of your choice.
There are a set of parameters included in the `config/config.yaml` file that will **need to be specified** when running the workflow with your own data.
These parameters are specific to the project or dataset being processed.
These project-specific parameters can be found under the [`Project-specific parameters` section](./config/config.yaml#L3) of the config file, while the remaining parameters that can be optionally modified are found under the [`Processing parameters` section](./config/config.yaml#L11).

To run the workflow on your data, modify the following parameters in the `config/config.yaml` file:

| Parameter        | Description |
|------------------|-------------|
| `results_dir` | full path to the directory where output files from running the core workflow will be stored |
| `project_metadata` | full path to your specific project metadata TSV file |
| `mito_file` | full path to a file containing a list of mitochondrial genes specific to the reference genome or transcriptome version used for alignment. By default, the workflow will use the mitochondrial gene list obtained from Ensembl version 104 of the Human transcriptome which can be found in the [`reference-files` directory](./reference-files). |

By default, these parameters point to the [example data](./example-data). 
The two example `_filtered.rds` files were both processed using the [`scpca-nf` workflow](https://github.com/AlexsLemonade/scpca-nf/blob/main/examples/README.md).
Therefore, if you would like to test this workflow using the example data, you can continue to the next step, running the workflow, without modifying the config file.

The config file also contains processing parameters like cutoffs for minimum genes detected, minimum unique molecular identifiers (UMI) per cell, etc. 
We have set default values for these parameters. 
Learn more about the [processing parameters](./additional-docs/additional-parameters.md#core-analysis-parameters) and how to modify them.

See the [processing information documentation](./additional-docs/processing-information.md) for more information on the individual workflow steps and how the parameters are used in each of the steps.

|[View Config File](https://github.com/AlexsLemonade/scpca-downstream-analyses/blob/main/config/config.yaml)|
|---|

## 5. Running the workflow

After you have successfully modified the required parameters in the config file, you can run the snakemake workflow with just the `--cores` and `--use-conda` flags as in the following example: 

```
snakemake --cores 2 --use-conda
```

It is mandatory to specify the number of CPU cores for snakemake to use by using the [`--cores` flag](https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html?highlight=cores#step-1-specifying-the-number-of-used-threads).
If `--cores` is given without a number, all available cores are used to run the workflow.

**Note:** If you did not install dependencies [with conda via snakemake](#snakemakeconda-installation), you will need to remove the `--use-conda` flag.

You can also modify the config file parameters at the command line, rather than manually as recommended in step 4.
See our [command line options](./additional-docs/command-line-options.md) documentation for more information.

**Note:** For Data Lab staff members working on development, new changes should be merged through a pull request to the `development` branch.
Changes will be pushed to the `main` branch once changes are ready for a new release (per the [release checklist document](.github/ISSUE_TEMPLATE/release-checklist.md)).

## 6. Expected output

For each `SingleCellExperiment` and associated `library_id` used as input, the workflow will return two files: a processed `SingleCellExperiment` object containing normalized data and clustering results, and a summary HTML report detailing the filtering of low quality cells, dimensionality reduction, and clustering that was performed within the workflow.
These files can be found in the `example_results` folder, as defined in the `config.yaml` file.
Within the `example_results` folder, output files for each library will be nested within folders labeled with the provided `sample_id`.
Each output filename will be prefixed with the associated `library_id`.
Below is an example of the nested file structure you can expect.

```
example_results
└── sample_id
	 ├── <library_id>_core_analysis_report.html
	 └── <library_id>_processed.rds
```

The `<library_id>_core_analysis_report.html` file is the [html file](https://bookdown.org/yihui/rmarkdown/html-document.html#html-document) that contains the summary report of the filtering, dimensionality reduction, and clustering results associated with the processed `SingleCellExperiment` object.

The `<library_id>_processed.rds` file is the [RDS file](https://rstudio-education.github.io/hopr/dataio.html#saving-r-files) that contains the final processed `SingleCellExperiment` object (which contains the filtered, normalized data and clustering results).

You can also download a ZIP file with an example of the output from running the core workflow, including the summary HTML report and processed `SingleCellExperiment` objects stored as RDS files [here](https://scpca-references.s3.amazonaws.com/example-data/scpca-downstream-analyses/core_example_results.zip).

### What to expect in the output `SingleCellExperiment` object

As a result of the normalization step of the workflow, a log-transformed normalized expression matrix can be accessed using [`logcounts(sce)`](https://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html#adding-more-assays).

In the [`colData`](https://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html#handling-metadata) of the output `SingleCellExperiment` object, you can find the following:

- Clustering results stored in a metadata column named using the associated clustering type and nearest neighbours values.
For example, if using the default values of Louvain clustering with a nearest neighbors parameter of 10, the column name would be `louvain_10` and can be accessed using `colData(sce)$louvain_10`.

In the [`metadata`](https://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html#other-metadata) of the output `SingleCellExperiment` object, which can be accessed using `metadata(sce)`,  you can find the following information:

| Metadata Key       | Description |
|----------------------------|-------------|
| `scpca_filter_method` | The type of filtering performed ([`miQC`](https://bioconductor.org/packages/release/bioc/html/miQC.html) or `manual`) on the expression data. |
| `prob_compromised_cutoff` | The maximum probability of cells being compromised, which is only present when the `filtering_method` is set to `miQC`. |
| `miQC_model` | The linear mixture model calculated by `miQC` and therefore is only present when `filtering_method` is set to `miQC`. |
| `mito_percent_cutoff` | Maximum percent mitochondrial reads per cell threshold, which is only present when `filtering_method` is set to `manual`. |
| `genes_filtered` | Indicates whether or not genes have been filtered. |
| `detected_gene_cutoff` | Minimum number of genes detected per cell, which is only present when `filtering_method` is set to `manual`. |
| `umi_count_cutoff` | Minimum unique molecular identifiers (UMI) per cell, which is only present when `filtering_method` is set to `manual`. |
| `num_filtered_cells_retained` | The number of cells retained after filtering using the specified filtering method, either `miQC` or `manual`. |
| `normalization` | The type of normalization performed (`log-normalization` or `deconvolution` when clustering of similar cells using [`scran::quickCluster()`](https://rdrr.io/bioc/scran/man/quickCluster.html) prior to normalization with `scater::logNormCounts()` is successful). |
| `highly_variable_genes` | The subset of the most variable genes, determined using [`scran::getTopHVGs()`](https://rdrr.io/bioc/scran/man/getTopHVGs.html). |

You can find more information on the above in the [processing information documentation](./additional-docs/processing-information.md).

## Additional analysis modules

### Clustering analysis

There is an optional clustering analysis workflow stored in the `optional-clustering-analysis` subdirectory of this repository.
This workflow can help users identify the optimal clustering method and parameters for each library in their dataset.
Libraries are unique, which means that the optimal clustering is likely to be library-dependent.
The clustering analysis workflow provided can be used to explore different methods of clustering and test a range of parameter values for the given clustering methods in order to identify the optimal clustering for each library.

For more on what's in the clustering analysis workflow and how to run the workflow, see the [`README.md`](optional-clustering-analysis/README.md#optional-clustering-analysis) file in the clustering analysis subdirectory.

### Genes of interest analysis

There is an optional genes of interest analysis pipeline in the `optional-goi-analysis` subdirectory of this repository.
This workflow can help users evaluate expression of a specific list of genes in their sample dataset and how the expression values compare to the remaining genes in the dataset.

For more on what's in the genes of interest analysis workflow and how to run the workflow, see the [`README.md`](optional-goi-analysis/README.md) file in the genes of interest analysis subdirectory.
