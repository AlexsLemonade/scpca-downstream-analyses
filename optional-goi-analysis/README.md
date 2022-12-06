# Optional Genes of Interest (GOI) Analysis

This directory includes a genes of interest analysis workflow to help users evaluate expression of a specific list of genes in a sample dataset.


**The genes of interest (GOI) analysis workflow cannot be implemented until after users have successfully run the main downstream analysis core workflow as described in this repository's main [README.md](../README.md) file.**

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Analysis overview](#analysis-overview)
- [Expected input](#expected-input)
- [Expected output](#expected-output)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Analysis overview

There are three main steps of this genes of interest analysis workflow:

1. **Mapping gene identifiers (Optional)**: If input gene identifiers do not match the gene identifiers used in the `SingleCellExperiment` object used as input (e.g., input genes are gene symbols while the `SingleCellExperiment` object contains Ensembl identifiers), then input gene identifiers must first be mapped to a specified type.
2. **Heirarchical clustering of the dataset with genes of interest**: Hierarchical clustering is performed on the normalized data associated with the provided genes of interest.
These clustering results are used to generate a heatmap that displays genes on the x-axis, cells on the y-axis, and uses color to indicate gene expression scaled using a z-score.
This heatmap can be used to identify if any clear structure or groupings of cells with similar gene expression patterns are present.
3. **Visualization of genes of interest expression** (Sina, PCA, and UMAP plots):
The normalized and transformed expression for each of the provided genes of interest is compared to the mean expression of all genes present in the dataset and shown using a sina plot.
UMAP and PCA plots are also provided, where each dot represents a cell and the color indicates the individual gene of interest’s expression.


**Note** that the same [software requirements for the core workflow](../README.md#3-additional-dependencies) are also required for this clustering workflow.
R 4.2 is required for running our pipeline, along with Bioconductor 3.15.
Package dependencies for the analysis workflows in this repository are managed using [`renv`](https://rstudio.github.io/renv/index.html), which must be installed locally prior to running the workflow.
If you are using conda, dependencies can be installed as [part of the initial setup](../README.md#snakemakeconda-installation).

## Expected input

To run this workflow, you will need to provide:

1. The RDS file containing the normalized [output `SingleCellExperiment` object](../README.md#expected-output) from the core dowstream analyses workflow.
This `SingleCellExperiment` object must contain a log-normalized counts matrix in an assay named `logcounts`.
2. A project metadata tab-separated value (TSV) file containing relevant information about your data necessary for processing, the same metadata file used in the core workflow should be used here (see more on this in the ["Metadata file format" section](../README.md#metadata-file-format) and an example of this metadata file [here](../project-metadata/example-library-metadata.tsv)).
3. The genes of interest list, stored as a tab-separated value (TSV) file.
This file should contain at least one column named `gene_id` with the relevant gene identifiers, and can optionally contain an additional column named `gene_set` that denotes the gene set that each gene identifier belongs to (see an example of this genes of interest file [here](../example-data/goi-lists/sample01_goi_list.tsv)).

## Expected output

For each provided `SingleCellExperiment` RDS file and associated `library_id`, the workflow will return the following files in the specified output directory:

```
goi_stats
├── library01_mapped_genes.tsv
├── library01_normalized_zscores.mtx
├── library01_heatmap_annotation.rds
└── library01_goi-report.html
```

1. A `_mapped_genes.tsv` file with mapped genes of interest if the `--perform_mapping` is `TRUE` upon running the workflow. Otherwise this file will store just the provided genes of interest, and a new column named `sce_rownames_identifier` with the genes of interest identifiers copied over to the column.
2. A `_normalized_zscores.mtx` matrix file with the z-scored matrix calculated using the normalized data specific to the provided genes of interest.
3. A `_heatmap_annotation.rds` file with the annotations to be used when plotting the heatmap for the html report.
4. The `_goi_report.html` file, which is the summary html report with plots containing the GOI results.

