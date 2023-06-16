# Optional Integration Analysis

This directory includes a data integration analysis workflow that can help users integrate data from multiple libraries. 

**The data integration analysis workflow cannot be implemented until after users have successfully run the main downstream analysis core workflow as described in this repository's main [README.md](../README.md) file or have downloaded data from the [ScPCA portal](https://scpca.alexslemonade.org/).**

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Analysis overview](#analysis-overview)
- [Create metadata file](#create-metadata-file)
- [Configure config file](#configure-config-file)
- [Running the workflow](#running-the-workflow)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Analysis overview

The data integration analysis workflow provided here can be used to explore multiple libraries simultaneously by integrating the data associated with specified libraries.

There are three main steps of this data integration analysis workflow:

1. **Merge**: The data from each of the libraries are first merged into one `SingleCellExperiment` object.
2. **Integrate**: The merged data is integrated using [`fastMNN`](https://rdrr.io/github/LTLA/batchelor/man/fastMNN.html) and/or [`harmony`](https://portals.broadinstitute.org/harmony/articles/quickstart.html) as specified by the user.
3. **Plot**: Once the data from the multiple libraries are integrated and stored in the `SingleCellExperiment` object, the results from each of the integration methods used are displayed in a UMAP plot.
Additionally, integration performance metrics are shown for each of the integration methods.
These metrics are displayed in plots and include iLISI, batch ARI, within-batch ARI, and batch average silhouette width (ASW).
The plots are displayed in a html report for ease of reference.

**Note** that the same [software requirements for the core workflow](../README.md#3-additional-dependencies) are also required for this clustering workflow.
R 4.2 is required for running our pipeline, along with Bioconductor 3.15.
Package dependencies for the analysis workflows in this repository are managed using [`renv`](https://rstudio.github.io/renv/index.html), which must be installed locally prior to running the workflow.
If you are using conda, dependencies can be installed as [part of the initial setup](../README.md#snakemakeconda-installation).

## Create metadata file

Before running the workflow, you will need to create a project metadata file as a tab-separated value (TSV) file that contains the relevant data for your input files needed to run the data integration workflow.
The file should contain the following columns:

- `sample_id`, unique ID for each piece of tissue or sample that cells were obtained from,  all libraries that were sampled from the same piece of tissue should have the same `sample_id`.
- `library_id`, unique ID used for each set of cells that has been prepped and sequenced separately.
- `processed_sce_filepath`, the full path to the RDS file containing the processed `SingleCellExperiment` object.
Each library ID should have a unique `processed_sce_filepath`.
- `integration_group`, the variable specifying which libraries should be grouped and integrated.

|[View Example Metadata File](https://github.com/AlexsLemonade/scpca-downstream-analyses/blob/main/example-data/project-metadata/example-integration-library-metadata.tsv)|
|---|


## Configure config file

As in the main core workflow, we have provided an [example snakemake configuration file](../config/config.yaml), `config/config.yaml`, which defines all parameters needed to run the workflow.
Learn more about snakemake configuration files [here](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).

The config file contains two sets of parameters:

- **[Project-specific Parameters](../config/config.yaml#L3)**: This set of parameters are for specifying dataset or project related details. 
These parameters are **required** to run the workflow on your data.
- **[Processing Parameters](../config/integration_config.yaml)**: This set of parameters specify configurations for the integration method(s) to be performed.
You can change them to explore your data but it is optional.
You can modify the relevant parameters by manually updating the `config/integration_config.yaml` file using a text editor of your choice.

To run the workflow on your data, modify the following parameters in the `config/config.yaml` and `config/integration_config.yaml` files:

| Parameter        | Description |
|------------------|-------------|
| `results_dir` | full path to the directory where output files will be stored |
| `integration_project_metadata` | full path to your integration-specific project metadata TSV file |

|[View Config File](../config/config.yaml)|
|---|

|[View Integration Config File](../config/integration_config.yaml)|
|---|

The [`config/integration_config.yaml`](../config/integration_config.yaml) file also contains additional processing parameters like the integration method(s) that should be used and the number of mulit-processing threads to use when merging the data.
We have set default values for these parameters. 
Learn more about the [processing parameters](../additional-docs/processing-parameters.md#integration-analysis-parameters) and how to modify them.

## Running the workflow

The execution file with the data integration Snakemake workflow is named `integration.snakefile` and can be found in the root directory. To tell snakemake to run the specific clustering workflow be sure to use the `--snakefile` or `-s` option followed by the name of the snakefile, `integration.snakefile`.

After you have successfully modified the required project-specific parameters in the config file and navigated to within the root directory of the `scpca-downstream-analyses` repository, you can run the clustering Snakemake workflow with just the `--cores` and `--use-conda` flags as in the following example: 

```
snakemake --snakefile integration.snakefile --cores 2 --use-conda
```

It is mandatory to specify the number of CPU cores for snakemake to use by using the [`--cores` flag](https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html?highlight=cores#step-1-specifying-the-number-of-used-threads).
If `--cores` is given without a number, all available cores are used to run the workflow.

**Note:** If you did not install dependencies [with conda via snakemake](#snakemakeconda-installation), you will need to remove the `--use-conda` flag.

You can also modify the config file parameters at the command line, rather than manually as recommended in the configure config file section above.
See our [command line options](../additional-docs/command-line-options.md) documentation for more information.
