# Optional Clustering Analysis

This directory includes a clustering analysis workflow that can be implemented after users have successfully ran the main downstream analysis core workflow as described in this repository's main [README.md](../README.md) file.

There are two main steps of this clustering analysis workflow:

1. **Cluster Calculations**: clustering results are calculated for the provided type of graph-based clustering, and range of nearest neighbor values.
Cluster validity and stability results are also calculated.
2. **Cluster Plots**: once clustering results are calculated and stored in the `SingleCellExperiment` object, each of the results are associated with UMAP plots, as well as cluster validity and stability plots to represent the reliability of each set of clustering results.
The plots are displayed in a html report for ease of reference.

**Note** that the same software requirements for the core workflow are also required for this clustering workflow.
R 4.2 is required for running our pipeline, along with Bioconductor 3.15.
Package dependencies for the analysis workflows in this repository are managed using [`renv`](https://rstudio.github.io/renv/index.html), and `renv` must be installed locally prior to running the workflow.
If you are using conda, dependencies can be installed as [part of the initial setup](../README.md#snakemakeconda-installation).

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Expected input](#expected-input)
- [Expected output](#expected-output)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Expected input

To run this workflow, you will need to provide:

1. The [output `SingleCellExperiment` object](../README.md#expected-output) from the core dowstream analyses workflow.
2. The desired type of graph-based clustering (can be "louvain" or "walktrap"), along with the minimum, maximum, and incremental values that will be used to define the range of nearest neighbor values to be tested.
For example, if you would like a range of 5:25 to be tested in increments of 5 (5, 10, 15, 20, 25), you will provide `nearest_neighbors_min` = 5, `nearest_neighbors_max` = 25, and `nearest_neighbors_increment` = 5.

## Expected output

For each `SingleCellExperiment` RDS file and associated `library_id` used as input, the workflow will return five files:

1. The `SingleCellExperiment` RDS file containing the added clustering results that were calculated in the first step of the worflow.
2. The summary html report with the plots representing the clustering results.
3. A `_clustering_all_validity_stats.tsv` file with all of the calculated cluster validity statistics.
4. A `_clustering_summary_validity_stats.tsv` file with the cluster validity stats averaged across each cluster assignment.
5. A `_clustering_summary_stability_stats.tsv` file with the average adjusted Rand Index (ARI) values across each cluster assignment.

The latter three files described above will be found in a `clustering_stats` subdirectory stored in the same results directory where the `SingleCellExperiment` RDS file is stored.
 


