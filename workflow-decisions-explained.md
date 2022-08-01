# Processing steps within the core analysis workflow

This documentation contains detailed information on the processing that occurs during each step of the core analysis workflow in this repository.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Filtering low quality cells](#filtering-low-quality-cells)
- [Normalization](#normalization)
- [Dimensionality reduction](#dimensionality-reduction)
- [Clustering](#clustering)


## Filtering low quality cells

The first step in the pipeline includes filtering and removal of low quality cells from each input library. 
There are two different methods that can be implemented to perform filtering of low quality cells in the core analysis workflow, which include:

1. miQC - this is the recommended filtering method, which implements [miQC](https://bioconductor.org/packages/release/bioc/html/miQC.html), a package that jointly models proportion of reads belonging to mitochondrial genes and number of unique genes detected to predict low quality cells.
If using the `miQC` option for filtering, cells above the maximum probability of a cell being compromised as calculated by miQC are removed. 
This value can be defined by altering the default `prob_compromised_cutoff` value of **0.75**, which is also the default used by miQC.
Note that if miQC fails, the manual filtering described below will be implemented as a default backup method.
2. manual - if chosen, low quality cells are removed that fall below or above the user provided thresholds for the `mito_percent_cutoff` (maximum percent mitochondrial reads per cell, default: 20), `detected_gene_cutoff` (minimum number of genes detected per cell, default: 500), and `umi_count_cutoff` (minimum unique molecular identifiers (UMI) per cell, default: 500).

See the main [`README.md`](https://github.com/AlexsLemonade/scpca-downstream-analyses#scpca-downstream-analyses) file for more information on the default thresholds set throughout the workflow.
The filtering method can be specified for each sample in the [metadata file](https://github.com/AlexsLemonade/scpca-downstream-analyses#metadata-file-format).
You can also find more on modifying additional filtering parameters under [processing parameters](https://github.com/AlexsLemonade/scpca-downstream-analyses#processing-parameters).

## Normalization

To normalize and log transform the filtered data, we implement the normalization by deconvolution method [(Lun *et al.* (2016))](https://doi.org/10.1186/s13059-016-0947-7), implemented in the `scran`/`scater` packages. 
We first use [`scran::quickCluster()`](https://rdrr.io/bioc/scran/man/quickCluster.html) to cluster similar cells where possible, making a note in the metadata of the sce oject recording whether or not similar cells were successfully clustered, followed by using the [`scater::logNormCounts()`](https://rdrr.io/github/LTLA/scuttle/man/logNormCounts.html) function.

## Dimensionality reduction

We calculate and add PCA and UMAP results to the normalized `SingleCellExperiment` object using the scater package functions [`runPCA()`](https://rdrr.io/bioc/scater/man/runPCA.html) and [`runUMAP()`](https://rdrr.io/bioc/scater/man/runUMAP.html), respectively.
These functions take in a subset of the most variable genes to specify the subset of features to use for the dimensionality reduction.
The default behavior at this step in the workflow is to subset the top **2000** high variance genes to use for calculating the dimensionality results.
See more on altering this top `n` value in the [processing parameters](https://github.com/AlexsLemonade/scpca-downstream-analyses#dimensionality-reduction-and-clustering-parameters) section of the main README file.

## Clustering


Here clustering is performed using [graph-based clustering](http://bioconductor.org/books/3.15/OSCA.basic/clustering.html#clustering-graph) options in our workflow. 
The clustering methods that can be implemented in the workflow include:

1. `louvain` - this is the default clustering method for the workflow for community detection after creation of the nearest neighbors graph.
In this case the weighting type will be `jaccard`, which means edges will be weighted according to the Jaccard index of the two sets of neighbors.
2. `walktrap` - this clustering method is based on random walks and uses the `rank` weighting type, which means edges will be weighted edge according to the highest average rank of the shared neighbors.

The default nearest neighbors value is set to 10. See more on modifying this value and the clustering method in the  [processing parameters](https://github.com/AlexsLemonade/scpca-downstream-analyses#dimensionality-reduction-and-clustering-parameters) section of the main README.
