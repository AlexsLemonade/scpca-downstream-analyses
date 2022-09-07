# Optional Clustering Analysis

This directory includes a clustering analysis workflow that can be implemented after users have successfully ran the main downstream analysis core workflow as described in this repository's main [README.md](../README.md) file.

The clustering analysis workflow can be used to explore different types of clustering and test a range of parameter values for a given type of clustering.
Libraries are unique, which means that the optimal clustering is going to be library-dependent.
Exploring clustering types and parameters can be helpful in choosing the optimal clustering for each library.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Analysis overview](#analysis-overview)
- [Expected input](#expected-input)
- [Expected output](#expected-output)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Analysis overview

There are two main steps of this clustering analysis workflow:

1. **Cluster Calculations**: clustering results are calculated for the provided type of graph-based clustering, and range of nearest neighbor values.
Cluster validity and stability results are also calculated.
2. **Cluster Plots**: once clustering results are calculated and stored in the `SingleCellExperiment` object, each of the results are associated with UMAP plots, as well as cluster validity and stability plots to represent the reliability of each set of clustering results.
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
 


