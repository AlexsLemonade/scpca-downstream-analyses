# Core analysis workflow decisions

This documentation contains information on the decisions made throughout the core analysis workflow in this repository.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Filtering](#filtering)
- [Normalization](#normalization)
- [Dimensionality reduction](#dimensionality-reduction)
- [Clustering](#clustering)


## Filtering

The options for filtering methods that can be implemented in the core analysis workflow include:

1. miQC - this is the default filtering method, where [miQC](https://bioconductor.org/packages/release/bioc/html/miQC.html) is implemented using the user provided threshold for `prob_compromised_cutoff`.
This threshold sets the maximum probability of a cell being compromised as calculated by miQC.
Our default `prob_compromised_cutoff` value is **0.75**, which is also the default used by miQC.
2. manual - if chosen, the input data are filtered using the user provided thresholds for the `mito_percent_cutoff` (maximum percent mitochondrial reads per cell), `detected_gene_cutoff` (minimum number of genes detected per cell), and `umi_count_cutoff` (minimum unique molecular identifiers (UMI) per cell).
We chose default values for these thresholds based on what looked reasonable for the datasets we used for development, and what was most common in literature.

See the main [`README.md`](./README.md) file for more information on the default thresholds set throughout the workflow.

## Normalization

To normalize and log transform the filtered data, we implement the [`scater::logNormCounts()`](https://rdrr.io/github/LTLA/scuttle/man/logNormCounts.html) function.
Before doing so, we use [`scran::quickCluster()`](https://rdrr.io/bioc/scran/man/quickCluster.html) to cluster similar cells where possible, making a note in the metadata of the sce oject recording whether or not similar cells were successfully clustered.

## Dimensionality reduction

We calculate and add PCA and UMAP results to the normalized `SingleCellExperiment` object using the scater package functions [`runPCA()`](https://rdrr.io/bioc/scater/man/runPCA.html) and [`runUMAP()`](https://rdrr.io/bioc/scater/man/runUMAP.html), respectively.
These functions take in a subset of the most variable genes to specify the subset of features to use for the dimensionality reduction.
The default behavior at this step in the workflow is to subset the top **2000** high variance genes to use for calculating the dimensionality results.

## Clustering

We made the decision to implement [graph-based clustering](http://bioconductor.org/books/3.15/OSCA.basic/clustering.html#clustering-graph) options in our workflow. The clustering methods that can be implemented in the workflow include:

1. `louvain` - this is the default clustering method for our workflow (which we chose due to the fact that it is also the default clustering method implemented by [Seurat](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)).
In this case the weighting type will be `jaccard`, which means edges will be weighted according to the Jaccard index of the two sets of neighbors.
2. `walktrap` - this clustering method is based on random walks and uses the `rank` weighting type, which means edges will be weighted edge according to the highest average rank of the shared neighbors.
