# Additional Workflow Parameters

Below is a list of all additional parameters that are needed to run the analysis workflows in this repository.
These parameters are all included in the config files and can optionally be altered when running the workflows, but are **not required to be changed**.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Core analysis parameters](#core-analysis-parameters)
  - [Filtering parameters](#filtering-parameters)
  - [Dimensionality reduction and clustering parameters](#dimensionality-reduction-and-clustering-parameters)
- [Clustering analysis parameters](#clustering-analysis-parameters)
- [Genes of interest analysis parameters](#genes-of-interest-analysis-parameters)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Core analysis parameters


The [configuration file](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html), `config/config.yaml` sets the defaults for all parameters needed to run the core analysis workflow.
It is **not required** to alter these parameters to run the workflow, but if you would like to change the filtering method or the minimum cutoff for the number of genes detected per cell, you can do so by changing these parameters via a text editor of your choice or at the command line per our documentation [here](./command-line-options.md).

### Filtering parameters

There are two types of filtering methods that can be specified in the project metadata file, [`miQC`](https://bioconductor.org/packages/release/bioc/html/miQC.html) or `manual` filtering.
For more information on choosing a filtering method, see [Filtering low quality cells](./processing-information.md#filtering-low-quality-cells) in the [processing information documentation](./processing-information.md).
Below are the parameters required to run either of the filtering methods.

| Parameter        | Description | Default value |
|------------------|-------------|---------------|
| `seed` | an integer to be used to set a seed for reproducibility when running the workflow | 2021 |
| `filtering_method` | `filtering_method`, the specified filtering method which can be one of "miQC" or "manual". For more information on choosing a filtering method, see [Filtering low quality cells](./processing-information.md#filtering-low-quality-cells) in the [processing information documentation](./processing-information.md) | "miQC" |
| `prob_compromised_cutoff` | the maximum probability of a cell being compromised as calculated by [miQC](https://bioconductor.org/packages/release/bioc/html/miQC.html), which is required when the `filtering_method` is set to `miQC` in the project metadata | 0.75 |
| `filter_genes` | a binary value indicating whether or not to perform gene fitering | `FALSE` |
| `gene_detected_row_cutoff` | the percent of cells a gene must be detected in; genes detected are only filtered if `filter_genes` is set to `TRUE` | 5 |
| `gene_means_cutoff` | mean gene expression minimum threshold; mean gene expression is only filtered if `filter_genes` is set to `TRUE` | 0.1 |
| `mito_percent_cutoff` | maximum percent mitochondrial reads per cell threshold, which is only required when `filtering_method` is set to `manual` | 20 |
| `min_gene_cutoff` | minimum number of genes detected per cell | 200 |
| `umi_count_cutoff` | minimum unique molecular identifiers (UMI) per cell, which is only required when `filtering_method` is set to `manual` | 500 |

### Dimensionality reduction and clustering parameters

In the core workflow, PCA and UMAP results are calculated and stored, and the PCA coordinates are used for graph-based clustering.
For more details on how the workflow performs [dimensionality reduction](./processing-information.md#dimensionality-reduction) and [clustering](./processing-information.md#clustering) see the documentation on [workflow processing information](./processing-information.md).
Below are the parameters required to run the dimensionality reduction and clustering steps of the workflow.

| Parameter        | Description | Default value |
|------------------|-------------|---------------|
| `n_genes_pca` | the `n` number of highly variable genes to select as input for PCA| 2000 |
| `cluster_type` | the type of clustering to be performed, values can be "louvain" or "walktrap" (see more on these graph-based clustering methods in this [Community Detection Algorithms article](https://towardsdatascience.com/community-detection-algorithms-9bd8951e7dae)) | "louvain" |
| `nearest_neighbors` | the `n` number of nearest neighbors when performing the chosen graph-based clustering method | 10 |

|[View Config File](../config/config.yaml)|
|---|

## Clustering analysis parameters

The [configuration file](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html), `config/cluster_config.yaml` sets the defaults for all parameters needed to run the clustering workflow.
It is **not required** to alter these parameters to run the workflow, but if you would like to change the type of clustering or range of nearest neighbor parameters, you can do so by changing these parameters via a text editor of your choice or at the command line per our documentation [here](./command-line-options.md). 

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

## Genes of interest analysis parameters

The [configuration file](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html), `config/goi_config.yaml`, sets the defaults for additional parameters required for mapping the gene identifiers in the genes of interest workflow.
These parameters will only be used if `perform_mapping` is set to `TRUE`.
It is **not required** to alter these parameters to run the workflow, but if you would like to modify the gene identifier mapping, you can do so by changing these parameters. 

The following gene mapping parameters found in the `config/goi_config.yaml` file can be optionally modified:

| Parameter        | Description | Default value |
|------------------|-------------|---------------|
| `perform_mapping` | a binary value indicating whether or not to perform gene identifier mapping | `TRUE` |
| `organism` | the organism associated with the provided genes and data | `"Homo sapiens"` |
| `sce_rownames_identifier` | the type of gene identifiers found in the rownames of the `SingleCellExperiment` object; note that this parameter should not be changed unless the output of the core analysis workflow has been altered in some way prior to running this module | `"ENSEMBL"` |
| `multi_mappings` | how to handle multiple gene identifier mappings when `perform_mapping` is `TRUE` | `"list"` |


|[View Genes of Interest Config File](../config/goi_config.yaml)|
|---|

## Integration analysis parameters

The [configuration file](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html), `config/integration_config.yaml` sets the defaults for all parameters needed to run the data integration workflow.
It is **not required** to alter these parameters to run the workflow, but if you would like to change the integration method(s) or the number of multi-processing threads to use, you can do so by changing these parameters via a text editor of your choice or at the command line per our documentation [here](./command-line-options.md). 

The parameters found in the `config/integration_config.yaml` file can be optionally modified and are as follows:

| Parameter        | Description | Default value |
|------------------|-------------|---------------|
| `threads` | the number of multiprocessing threads to use | 1 |
| `integration_method` | the method(s) to be used for integration | `"fastMNN,harmony"` |
| `batch_column` | the name of the `SingleCellExperiment` column that indicates the grouping of interest; this will generally be batches or cell types | `"library_id"` |
| `cell_id_column` | the name of the `SingleCellExperiment` column that contains the cell barcodes | `"cell_id"` |

|[View Integration Config File](../config/integration_config.yaml)|
|---|
