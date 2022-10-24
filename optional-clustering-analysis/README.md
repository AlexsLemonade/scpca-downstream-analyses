# Optional Clustering Analysis

This directory includes a clustering analysis workflow that can help users identify the optimal clustering method and parameters for each library in their dataset. 
Libraries are unique, which means that the optimal clustering is likely to be library-dependent.
The clustering analysis workflow provided here can be used to explore different methods of clustering and test a range of parameter values for given graph-based clustering methods in order to identify the optimal clustering for each library.

We evaluate the clustering results across each of the parameters tested by looking at the following measures:

1. **Cluster purity** - The cluster purity is calculated for each cell as the proportion of neighboring cells that are assigned to the same cluster, with a range of 0 to 1. Well-separated clusters should show little overlap between member and neighboring cells, and therefore high purity values for all member cells.
2. **Silhouette width** - Silhouette width is calculated for each cell across clusters and measures how well-separated each of the clusters are, with a range of -1 to 1. For each cell, the average distance to all cells in the same cluster and the average distance to all cells in another cluster, are calculated. The silhouette width for each cell is defined as the difference between these two values divided by their maximum distance. Therefore, cells will ideally have large positive silhouette widths, as this would mean that the cells of one cluster are well-separated from other clusters.
3. **Cluster stability** - The stability of the clustering results associated with each of the clustering parameters tested is also reported. Here, cells within each dataset are sampled using bootstrapping and the sampled cells are re-clustered. The new clustering assignments are compared to the original assignment by obtaining an adjusted rand index, and this process is repeated 20 times. Clustering results with high stability would reveal that clustering after each bootstrap replicate is consistent with the original clustering results. We use summary adjusted Rand Index (ARI) values in the plots shown in the clustering report output by the workflow to represent the average calculated cluster stability values. The range here is 0 to 1, where stable clusters have values closer to 1.

For a more in depth discussion on these metrics and how they can be used to identify the optimal clustering results, see the [advanced clustering chapter in Orchestrating Single Cell Analysis](http://bioconductor.org/books/3.15/OSCA.advanced/clustering-redux.html#motivation).

**The clustering analysis workflow cannot be implemented until after users have successfully run the main downstream analysis core workflow as described in this repository's main [README.md](../README.md) file.**

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Analysis overview](#analysis-overview)
- [Expected input](#expected-input)
- [Expected output](#expected-output)
- [Running the workflow](#running-the-workflow)
  - [Project-specific parameters](#project-specific-parameters)
  - [Clustering parameters](#clustering-parameters)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Analysis overview

There are two main steps of this clustering analysis workflow:

1. **Cluster Calculations**: Clustering results are calculated for the provided type(s) of graph-based clustering, and range of nearest neighbor values.
Cluster validity and stability results are also calculated.
2. **Cluster Plots**: Once clustering results are calculated and stored in the `SingleCellExperiment` object, the results from each of the clustering methods tested are displayed in a UMAP plot.
Additionally, metrics associated with each of the clustering results such as silhouette width, cluster purity, and cluster stability (as described above) are calculated and plotted.
The plots are displayed in a html report for ease of reference.

**Note** that the same [software requirements for the core workflow](../README.md#3-additional-dependencies) are also required for this clustering workflow.
R 4.2 is required for running our pipeline, along with Bioconductor 3.15.
Package dependencies for the analysis workflows in this repository are managed using [`renv`](https://rstudio.github.io/renv/index.html), which must be installed locally prior to running the workflow.
If you are using conda, dependencies can be installed as [part of the initial setup](../README.md#snakemakeconda-installation).

## Expected input

To run this workflow, you will need to provide:

1. The RDS file containing the [output `SingleCellExperiment` object](../README.md#expected-output) from the core dowstream analyses workflow.
2. A project metadata tab-separated value (TSV) file containing relevant information about your data necessary for processing, the same metadata file used in the core workflow should be used here (see more on this in the ["Metadata file format" section](../README.md#metadata-file-format) and an example of this metadata file [here](../project-metadata/example-library-metadata.tsv)).
3. The desired type of graph-based clustering (can be "louvain" and/or "walktrap"), along with the minimum, maximum, and incremental values that will be used to define the range of nearest neighbor values to be tested.
For example, if you would like a range of `5:25` to be tested in increments of 5 (as in, `5, 10, 15, 20, 25`), you will provide `nearest_neighbors_min` = 5, `nearest_neighbors_max` = 25, and `nearest_neighbors_increment` = 5.

## Expected output

For each provided `SingleCellExperiment` RDS file and associated `library_id`, the workflow will return five files in the same directory where the inputted RDS file is stored:

1. The `SingleCellExperiment` RDS file containing the added clustering results that were calculated in the first step of the worflow.
2. The `_clustering_report.html` file, which is the summary html report with plots containing the clustering results.
3. A `_clustering_all_validity_stats.tsv` file with all of the calculated cluster validity statistics.
4. A `_clustering_summary_validity_stats.tsv` file with the cluster validity stats averaged across each cluster assignment.
5. A `_clustering_summary_stability_stats.tsv` file with the average adjusted Rand Index (ARI) values across each cluster assignment.

The TSV files will be saved in a subdirectory called `{library_id}_{filtering_method}_clustering_stats/`, where both the `{library_id}` and `{filtering_method}` are obtained from the input metadata file.

## Running the workflow

As in the main core workflow, we have provided a [configuration file](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html), `config/config.yaml` which sets the defaults for project-specific parameters needed to run the clustering workflow.

### Project-specific parameters

There are a set of parameters included in the `config/config.yaml` file that will need to be specified when running the workflow.
These parameters are specific to the project or dataset being processed.
These include the following parameters:

| Parameter        | Description |
|------------------|-------------|
| `results_dir` | relative path to the directory where output files will be stored (use the same `results_dir` used in the prerequisite core workflow) |
| `project_metadata` | relative path to your specific project metadata TSV file (use the same `project_metadata` used in the prerequisite core workflow) |

|[View Config File](../config/config.yaml)|
|---|

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

### Clustering parameters

We have provided an additional [configuration file](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html), `config/cluster_config.yaml` which sets the defaults for all parameters needed to run the clustering workflow.
It is **not required** to alter these parameters to run the workflow, but if you would like to change the type of clustering or range of nearest neighbor parameters, you can do so by changing these parameters. 

The parameters found in the `config/cluster_config.yaml` file can be optionally modified and are as follows:

| Parameter        | Description | Default value |
|------------------|-------------|---------------|
| `optional_cluster_types` | a comma-separated list of graph-based clustering type(s), options include "louvain" and/or "walktrap" (see more on these graph-based clustering methods in this [Community Detection Algorithms article](https://towardsdatascience.com/community-detection-algorithms-9bd8951e7dae)) | "louvain,walktrap" |
| `nearest_neighbors_min` | the minimum value to use for a range of neareast neighbors values for exploration | 5 |
| `nearest_neighbors_max` | the maximum value to use for a range of neareast neighbors values for exploration | 25 |
| `nearest_neighbors_increment` | the increment to use when implementing the range number of nearest neighbors for cluster stats (e.g. a value of 5 with min of 5 and max of 25 will test the nearest neighbors values of 5, 10, 15, 20, and 25) | 5 |
| `overwrite_results` | a binary value indicating whether or not to overwrite any existing clustering results | `TRUE` |

|[View Clustering Config File](../config/cluster_config.yaml)|
|---|