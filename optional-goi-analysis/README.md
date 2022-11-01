# Optional Genes of Interest (GOI) Analysis

This directory includes a genes of interest analysis workflow that can help users evaluate values specific to a list of genes against the rest of a sample dataset.


**The genes of interest (GOI) analysis workflow cannot be implemented until after users have successfully run the main downstream analysis core workflow as described in this repository's main [README.md](../README.md) file.**

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Analysis overview](#analysis-overview)
- [Expected input](#expected-input)
- [Expected output](#expected-output)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Analysis overview

There are two main steps of this genes of interest analysis workflow:

1. **Genes of Interest Calculations**: Provided gene identifiers are first optionally mapped to a specified keytype.
Then the normalized data are used to prepare heatmap matrix and annotation objects using only the provided genes of interest.
2. **Genes of Interest Plots**: Once GOI results are calculated and stored in separate RDS files, the heatmap matrix and annotation objects are used to plot a heatmap that displays genes on the x-axis, cells on the y-axis, and uses color to indicate gene expression scaled using a z-score.
This heatmap can used to identify if any clear structure or groupings of cells with similar gene expression patterns are present.
Additionally, the normalized and transformed expression for each of the provided genes of interest in comparison to the mean expression of all genes present in the dataset are plotted using `geom_sina()`.
Following are UMAP and PCA plots where each dot represents a cell and the color indicates the individual gene of interest’s expression.


**Note** that the same [software requirements for the core workflow](../README.md#3-additional-dependencies) are also required for this clustering workflow.
R 4.2 is required for running our pipeline, along with Bioconductor 3.15.
Package dependencies for the analysis workflows in this repository are managed using [`renv`](https://rstudio.github.io/renv/index.html), which must be installed locally prior to running the workflow.
If you are using conda, dependencies can be installed as [part of the initial setup](../README.md#snakemakeconda-installation).

## Expected input

To run this workflow, you will need to provide:

1. The RDS file containing the normalized [output `SingleCellExperiment` object](../README.md#expected-output) from the core dowstream analyses workflow.
2. A project metadata tab-separated value (TSV) file containing relevant information about your data necessary for processing, the same metadata file used in the core workflow should be used here (see more on this in the ["Metadata file format" section](../README.md#metadata-file-format) and an example of this metadata file [here](../project-metadata/example-library-metadata.tsv)).
3. The genes of interest list, stored as a tab-separated value (TSV) file.
This file should contain at least one column named `gene_id` with the relevant gene identifiers, and can optionally contain an additional column named `gene_set` that denotes the gene set that each gene identifier belongs to (see an example of this genes of interest file [here](../example-data/goi-lists/sample01_goi_list.tsv)).

## Expected output

For each provided `SingleCellExperiment` RDS file and associated `library_id`, the workflow will return five files in the same directory where the inputted RDS file is stored:

1. **Optionally**, a `_mapped_genes.tsv` file with mapped genes of interest if the `--perform_mapping` is `TRUE` upon running the workflow.
2. A `_normalized_sce_zscores.rds` file with the z-scored matrix calculated using the normalized data specific to the provided genes of interest.
3. A `_heatmap_annotation.rds` file with the annotations to be used when plotting the heatmap for the html report.
4. The `_goi_report.html` file, which is the summary html report with plots containing the GOI results.

