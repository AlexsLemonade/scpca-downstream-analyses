# Optional Integration Analysis

This directory includes a data integration workflow that can be used to combine data from multiple individual single-cell/single-nuclei libraries to obtain a single merged object. 
The merged objects contain the raw and normalized counts for all included libraries as well as batch corrected PCA and UMAP results.

**The data integration workflow cannot be implemented until after users have successfully run the main downstream analysis core workflow as described in this repository's main [README.md](../README.md) file or have downloaded data from the [ScPCA portal](https://scpca.alexslemonade.org/).**

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Analysis overview](#analysis-overview)
- [Expected input](#expected-input)
- [Create metadata file](#create-metadata-file)
- [Configure config file](#configure-config-file)
- [Running the workflow](#running-the-workflow)
- [Expected output](#expected-output)
  - [What to expect in the integrated `SingleCellExperiment` object](#what-to-expect-in-the-integrated-singlecellexperiment-object)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Analysis overview

The data integration workflow provided here can be used to create a single object containing data from multiple libraries.
First, the raw and normalized counts data from the original libraries are merged into a single object containing all cells and only shared genes across all libraries.
Then, batch correction is performed using either [`fastMNN`](https://rdrr.io/github/LTLA/batchelor/man/fastMNN.html) and/or [`harmony`](https://portals.broadinstitute.org/harmony/articles/quickstart.html).
Batch-corrected PCA and UMAP results are added to the merged object.
The merged objects are useful for exploring and comparing data in a group of single-cell/single-nuclei libraries together.

There are three main steps of this data integration workflow:

1. **Merge**: The counts data from each of the libraries are first merged into one `SingleCellExperiment` object.
2. **Integrate**: The merged data is integrated using [`fastMNN`](https://rdrr.io/github/LTLA/batchelor/man/fastMNN.html) and/or [`harmony`](https://portals.broadinstitute.org/harmony/articles/quickstart.html) to obtain batch corrected reduced dimensionality embeddings.
3. **Generate integration summary report**: Once the data from the multiple libraries are integrated and stored in the `SingleCellExperiment` object, the results from batch correction are evaluated to create a summary report.
This report includes UMAPs showing data from all libraries before and after batch correction and summaries of a set of integration metrics.
The plots are displayed in a html report for ease of reference.

The metrics we use to evaluate integration, as mentioned in step 3 above, include:

- Batch average silhouette width (ASW), which assesses consistency within clusters and measures how close a given data point adheres to cells from the same cluster vs. other clusters.
- Within-batch ARI, which measures how well a given integration method preserves biological heterogeneity by comparing pre- and post-integration clustering assignments for each library.

**Note** that the same [software requirements for the core workflow](../README.md#3-additional-dependencies) are also required for this clustering workflow.
R 4.2 is required for running our pipeline, along with Bioconductor 3.15.
If you are using conda, dependencies can be installed as [part of the initial setup](../README.md#snakemakeconda-installation).
Package dependencies for the analysis workflows in this repository are managed using [`renv`](https://rstudio.github.io/renv/index.html), which can be installed independently if desired, but we recommend using Snakemake's conda integration to set up the R environment and all dependencies that the workflow will use.

## Expected input

To run this workflow, you will need to provide:

1. Two or more RDS files containing normalized [`SingleCellExperiment` objects](../README.md#expected-output).
These can be objects output from the core dowstream analyses workflow or the processed `SingleCellExperiment` objects downloaded from the ScPCA portal (found in the `_processed.rds` file).
Each `SingleCellExperiment` object must contain a log-normalized counts matrix in an assay named `logcounts`.
2. A project metadata tab-separated value (TSV) file containing relevant information about your data necessary for processing, including `sample_id`, `library_id`, `processed_sce_filepath`, and `integration_group` (see more on this in the ["Create metadata file" section](#create-metadata-file) and an example of this metadata file [here](../project-metadata/example-integration-library-metadata.tsv)).

If working with data from the ScPCA portal, see our guide on preparing that data to run the data integration workflow [here](../additional-docs/working-with-scpca-portal-data.md).

## Create metadata file

Before running the workflow, you must create a project metadata file as a tab-separated value (TSV) file.
Each row in the metadata should correspond to a single-cell/single-nuclei library to include in integration.
The file should contain the following columns:

| column name | description |
| ----------- | ----------- |
| `sample_id` | Unique ID for each piece of tissue or sample that cells were obtained from,  all libraries that were sampled from the same piece of tissue should have the same `sample_id`. |
| `library_id` | Unique ID used for each set of cells that has been prepped and sequenced separately. |
| `processed_sce_filepath` | The full path to the RDS file containing the processed `SingleCellExperiment` object, each library ID should have a unique `processed_sce_filepath`; note that this RDS file is the output RDS file from the [core analysis workflow](../README.md#6-expected-output) OR the processed file from the [Single-cell Pediatric Cancer Atlas portal](https://scpca.alexslemonade.org/). The filename should contain the following format, `<library_id>_processed.rds`. |
| `integration_group` | Name used to specify which libraries should be merged. All libraries with the same `integration_group` will be merged and integrated into a single `SingleCellExperiment` object. |

|[View Example Metadata File](https://github.com/AlexsLemonade/scpca-downstream-analyses/blob/main/example-data/project-metadata/example-integration-library-metadata.tsv)|
|---|


## Configure config file

As in the main core workflow, we have provided an [example snakemake configuration file, `config/integration_config.yaml`](../config/integration_config.yaml), which defines all parameters needed to run the workflow.
Learn more about snakemake configuration files [here](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).

The config file contains two sets of parameters:

- **[Project-specific Parameters](../config/integration_config.yaml#L3)**: This set of parameters are for specifying dataset or project related details. 
These parameters are **required** to run the workflow on your data.
- **[Processing Parameters](../config/integration_config.yaml#L7)**: This set of parameters specify configurations for the integration method(s) to be performed.
You can change them to explore your data but it is optional.
You can modify the relevant parameters by manually updating the `config/integration_config.yaml` file using a text editor of your choice.

To run the workflow on your data, modify the following parameters in the `config/integration_config.yaml` file:

| Parameter        | Description |
|------------------|-------------|
| `results_dir` | full path to the directory where output files will be stored |
| `integration_project_metadata` | full path to your integration-specific project metadata TSV file |

|[View Integration Config File](../config/integration_config.yaml)|
|---|

The [`config/integration_config.yaml`](../config/integration_config.yaml) file also contains additional processing parameters like the integration method(s) that should be used and the number of mulit-processing threads to use when merging the data.
We have set default values for these parameters. 
Learn more about the [processing parameters](../additional-docs/additional-parameters.md#integration-analysis-parameters) and how to modify them.

## Running the workflow

The execution file with the data integration Snakemake workflow is named `integration.snakefile` and can be found in the root directory. 
To tell snakemake to run the specific clustering workflow be sure to use the `--snakefile` or `-s` option followed by the name of the snakefile, `integration.snakefile`.

After you have successfully modified the required project-specific parameters in the `integration_config.yaml` file and navigated to within the root directory of the `scpca-downstream-analyses` repository, you can run the clustering Snakemake workflow with just the `--cores` flag as in the following example:  

```
snakemake --snakefile integration.snakefile --cores 2
```

It is mandatory to specify the number of CPU cores for snakemake to use by using the [`--cores` flag](https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html?highlight=cores#step-1-specifying-the-number-of-used-threads).
If `--cores` is given without a number, all available cores are used to run the workflow.

You can also modify the config file parameters at the command line, rather than manually as recommended in the configure config file section above.
See our [command line options](../additional-docs/command-line-options.md) documentation for more information.

## Expected output

For each provided `integration_group`, the workflow will return two files saved in the `results_dir` specified in the config file:

1. The `<integration_group>_integrated_sce.rds` file containing the `SingleCellExperiment` object containing batch-corrected embeddings for each `integration_method` specified in the project metadata.
2. The `<integration_group>_integration_report.html` file is a summary html report containing plots that summarize the data integration results.

Below is an example of the nested file structure you can expect.

```
example_results
├── <integration_group>_integrated_sce.rds
└── <integration_group>_integration_report.html
```

You can also download a ZIP file with an example of the output from running the data integration workflow, including the summary HTML report and the integrated `SingleCellExperiment` object stored as an RDS file, [here](https://scpca-references.s3.amazonaws.com/example-data/scpca-downstream-analyses/integration_example_results.zip).

### What to expect in the integrated `SingleCellExperiment` object

In the [`reducedDim`](https://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html#dimensionality-reduction-results) of the integrated `SingleCellExperiment` object, you can find the following:

Batch-corrected embeddings stored in the reduced dimensions and named using the associated type of integration performed.
For example, if `fastMNN` is the type of integration performed, the PCA integration results are stored in `fastMNN_PCA` and can be accessed using `reducedDim(sce, "fastMNN_PCA")`.
The UMAP results can similary be accessed using `reducedDim(sce, "fastMNN_UMAP")`.

You can find more information on the above in the [processing information documentation](../additional-docs/processing-information.md).

