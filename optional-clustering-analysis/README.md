# Optional Clustering Analysis

This directory includes a clustering analysis workflow that can help users identify the optimal clustering method and parameters for each library in their dataset. 
Libraries are unique, which means that the optimal clustering is likely to be library-dependent.
The clustering analysis workflow provided here can be used to explore different methods of clustering and test a range of parameter values for a given clustering method in order to identify the optimal clustering for each library.

**The clustering analysis workflow cannot be implemented until after users have successfully run the main downstream analysis core workflow as described in this repository's main [README.md](../README.md) file.**

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Analysis overview](#analysis-overview)
- [Expected input](#expected-input)
- [Expected output](#expected-output)
- [Running the workflow](#running-the-workflow)
  - [Project-specific parameters](#project-specific-parameters)
  - [Processing parameters](#processing-parameters)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Analysis overview

There are two main steps of this clustering analysis workflow:

1. **Cluster Calculations**: clustering results are calculated for the provided type of graph-based clustering, and range of nearest neighbor values.
Cluster validity and stability results are also calculated.
2. **Cluster Plots**: once clustering results are calculated and stored in the `SingleCellExperiment` object, the results from each of the clustering methods tested are displayed in a UMAP plot.
Additionally, metrics associated with each of the clustering results such as silhouette width, cluster purity, and cluster stability are calculated and plotted.
The plots are displayed in a html report for ease of reference.

**Note** that the same [software requirements for the core workflow](https://github.com/AlexsLemonade/scpca-downstream-analyses/tree/development#3-additional-dependencies) are also required for this clustering workflow.
R 4.2 is required for running our pipeline, along with Bioconductor 3.15.
Package dependencies for the analysis workflows in this repository are managed using [`renv`](https://rstudio.github.io/renv/index.html), which must be installed locally prior to running the workflow.
If you are using conda, dependencies can be installed as [part of the initial setup](../README.md#snakemakeconda-installation).

## Expected input

To run this workflow, you will need to provide:

1. The RDS file containing the [output `SingleCellExperiment` object](../README.md#expected-output) from the core dowstream analyses workflow.
2. The [`library_id`](../README.md#metadata-file-format) associated with the `SingleCellExperiment` object.
2. The desired type of graph-based clustering (can be "louvain" or "walktrap"), along with the minimum, maximum, and incremental values that will be used to define the range of nearest neighbor values to be tested.
For example, if you would like a range of `5:25` to be tested in increments of 5 (as in, `5, 10, 15, 20, 25`), you will provide `nearest_neighbors_min` = 5, `nearest_neighbors_max` = 25, and `nearest_neighbors_increment` = 5.

## Expected output

For each provided `SingleCellExperiment` RDS file and associated `library_id`, the workflow will return five files in the same directory where the inputted RDS file is stored:

1. The `SingleCellExperiment` RDS file containing the added clustering results that were calculated in the first step of the worflow.
2. The `_optional_clustering_report.html` file, which is the summary html report with plots containing the clustering results.
3. A `_clustering_all_validity_stats.tsv` file with all of the calculated cluster validity statistics.
4. A `_clustering_summary_validity_stats.tsv` file with the cluster validity stats averaged across each cluster assignment.
5. A `_clustering_summary_stability_stats.tsv` file with the average adjusted Rand Index (ARI) values across each cluster assignment.

The TSV files will be saved in a subdirectory called `clustering_stats/`.

## Running the workflow

As in the main core workflow, we have provided a [configuration file](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html), `config.yaml` which sets the defaults for all parameters needed to run the clustering workflow.

### Project-specific parameters

There are a set of parameters included in the `config.yaml` file that will need to be specified when running the workflow.
These parameters are specific to the project or dataset being processed.
These include the following parameters:

| Parameter        | Description |
|------------------|-------------|
| `results_dir` | relative path to the directory where output files will be stored (use the same `results_dir` used in the prerequisite core workflow) |
| `project_metadata` | relative path to your specific project metadata TSV file (use the same `project_metadata` used in the prerequisite core workflow) |

The above parameters can be modified at the command line by using the [`--config` flag](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).
You must also specify the number of CPU cores for snakemake to use by using the [`--cores` flag](https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html?highlight=cores#step-1-specifying-the-number-of-used-threads).
If `--cores` is given without a number, all available cores are used to run the workflow.
If you installed dependencies for the workflow [with conda via snakemake](../README.md#snakemakeconda-installation), you will need to provide the `--use_conda` flag as well.

The execution file with the clustering Snakemake workflow is named `cluster.snakefile` and can be found in the root directory.
To run tell snakemake to run the specific clustering workflow be sure to use the `--snakefile` or `-s` option followed by the name of the snakefile, `cluster.snakefile`.
The below code is an example of running the clustering workflow using the project-specific parameters.

```
snakemake --snakefile cluster.snakefile \ 
  --cores 2 \
  --use-conda \
  --config results_dir="relative path to relevant results directory" \
  project_metadata="relative path to your-project-metadata.TSV"
```

**Note:**  If you did not install dependencies [with conda via snakemake](../README.md#snakemakeconda-installation), you will need to remove the `--use-conda` flag.

### Processing parameters

The parameters found under the `Processing parameters` section of the config file can be optionally modified.
Those that are relevant to the clustering workflow are as follows:

| Parameter        | Description | Default value |
|------------------|-------------|---------------|
| `nearest_neighbors_min` | the minimum value to use for a range of neareast neighbors values for exploration | 5 |
| `nearest_neighbors_max` | the maximum value to use for a range of neareast neighbors values for exploration | 25 |
| `nearest_neighbors_increment` | the increment to use when implementing the range number of nearest neighbors for cluster stats (e.g. a value of 5 with min of 5 and max of 25 will test the nearest neighbors values of 5, 10, 15, 20, and 25) | 5 |
| `overwrite_results` | a binary value indicating whether or not to overwrite any existing clustering results | `TRUE` |

