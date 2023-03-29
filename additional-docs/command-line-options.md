# Command Line Options for Running the Workflow(s)

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Running the workflow(s) at the command line](#running-the-workflows-at-the-command-line)
- [List of core workflow parameters that can be modified](#list-of-core-workflow-parameters-that-can-be-modified)
- [List of additional analysis module parameters that can be modified](#list-of-additional-analysis-module-parameters-that-can-be-modified)
  - [Clustering analysis parameters](#clustering-analysis-parameters)
  - [Genes of interest analysis parameters](#genes-of-interest-analysis-parameters)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Running the workflow(s) at the command line

In [step 4](../README.md#4-configure-config-file) of the main `README.md` file of this repository, we note that there are parameters that need to be modified in the provided config file to successfully run the workflow.
There we recommend modifying the `config/config.yaml` file manually.

All of the parameters in the config file can also be modified at the command line by using the [`--config` flag](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).
The below code is an example of running the Snakemake workflow using the project-specific parameters.

```
snakemake --cores 2 \
  --use-conda \
  --config results_dir="<FULL PATH TO RESULTS DIRECTORY>" \
  project_metadata="<FULL PATH TO YOUR PROJECT METADATA TSV>"
```

You will want to replace the paths for both `results_dir` and `project_metadata` to successfully run the workflow in the example above.
Where `results_dir` is the full path to the results directory where all results from running the workflow will be stored and `project_metadata` is the full path to the TSV file containing the relevant information about your input files.

## List of core workflow parameters that can be modified

Below are a list of the core workflow parameters that can be modified at the command line using the `--config` flag:

| Parameter        | Description | Default value |
|------------------|-------------| --------------|
| `results_dir` | relative path to the directory where output files from running the core workflow will be stored | - |
| `project_metadata` | relative path to your specific project metadata TSV file | - |
| `mito_file` | full path to a file containing a list of mitochondrial genes specific to the reference genome or transcriptome version used for alignment. By default, the workflow will use the mitochondrial gene list obtained from Ensembl version 104 of the Human transcriptome which can be found in the [`reference-files` directory](./reference-files). | - |
| `seed` | an integer to be used to set a seed for reproducibility when running the workflow | 2021 |
| `filtering_method` | `filtering_method`, the specified filtering method which can be one of "miQC" or "manual". For more information on choosing a filtering method, see [Filtering low quality cells](./additional-docs/processing-information.md#filtering-low-quality-cells) in the [processing information documentation](./additional-docs/processing-information.md) | "miQC" |
| `prob_compromised_cutoff` | the maximum probability of a cell being compromised as calculated by [miQC](https://bioconductor.org/packages/release/bioc/html/miQC.html), which is required when the `filtering_method` is set to `miQC` in the project metadata | 0.75 |
| `filter_genes` | a binary value indicating whether or not to perform gene fitering | `FALSE` |
| `gene_detected_row_cutoff` | the percent of cells a gene must be detected in; genes detected are only filtered if `filter_genes` is set to `TRUE` | 5 |
| `gene_means_cutoff` | mean gene expression minimum threshold; mean gene expression is only filtered if `filter_genes` is set to `TRUE` | 0.1 |
| `mito_percent_cutoff` | maximum percent mitochondrial reads per cell threshold, which is only required when `filtering_method` is set to `manual` | 20 |
| `min_gene_cutoff` | minimum number of genes detected per cell | 200 |
| `umi_count_cutoff` | minimum unique molecular identifiers (UMI) per cell, which is only required when `filtering_method` is set to `manual` | 500 |


## List of additional analysis module parameters that can be modified

### Clustering analysis parameters

For the clustering analysis module, the parameters found in the `config/cluster_config.yaml` file can be optionally modified at the command line using the `--config` flag and are as follows:

| Parameter        | Description | Default value |
|------------------|-------------|---------------|
| `input_data_dir` | relative path to the directory where the input data files can be found (default will be the `results_dir` used in the core workflow) | `"example_results"` |
| `optional_cluster_types` | a comma-separated list of graph-based clustering type(s), options include "louvain" and/or "walktrap" (see more on these graph-based clustering methods in this [Community Detection Algorithms article](https://towardsdatascience.com/community-detection-algorithms-9bd8951e7dae)) | "louvain,walktrap" |
| `nearest_neighbors_min` | the minimum value to use for a range of neareast neighbors values for exploration | 5 |
| `nearest_neighbors_max` | the maximum value to use for a range of neareast neighbors values for exploration | 25 |
| `nearest_neighbors_increment` | the increment to use when implementing the range number of nearest neighbors for cluster stats (e.g. a value of 5 with min of 5 and max of 25 will test the nearest neighbors values of 5, 10, 15, 20, and 25) | 5 |
| `overwrite_results` | a binary value indicating whether or not to overwrite any existing clustering results | `TRUE` |

### Genes of interest analysis parameters

For the genes of interest analysis module, the parameters found in the `config/goi_config.yaml` file can be optionally modified at the command line using the `--config` flag and are as follows:

| Parameter        | Description | Default value |
|------------------|-------------| --------------|
| `input_data_dir` | relative path to the directory where the input data files can be found (default will be the `results_dir` used in the core workflow) | `"example_results"` |
| `goi_list` | the file path to a tsv file containing the list of genes that are of interest | `"example-data/goi-lists/example_goi_list.tsv"` |
| `provided_identifier` | the type of gene identifiers used to populate the genes of interest list; example values that can implemented here include `"ENSEMBL"`, `"ENTREZID"`, `"SYMBOL"`; see more keytypes [here](https://jorainer.github.io/ensembldb/reference/EnsDb-AnnotationDbi.html) | `"SYMBOL"` |
| `overwrite` | a binary value indicating whether or not to overwrite existing output files | `TRUE` |
| `perform_mapping` | a binary value indicating whether or not to perform gene identifier mapping | `TRUE` |
| `organism` | the organism associated with the provided genes and data | `"Homo sapiens"` |
| `sce_rownames_identifier` | the type of gene identifiers found in the rownames of the `SingleCellExperiment` object; note that this parameter should not be changed unless the output of the core analysis workflow has been altered in some way prior to running this module | `"ENSEMBL"` |
| `multi_mappings` | how to handle multiple gene identifier mappings when `perform_mapping` is `TRUE` | `"list"` |