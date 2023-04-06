# Optional Genes of Interest (GOI) Analysis

This directory includes a genes of interest analysis workflow to help users evaluate expression of a specific list of genes in a sample dataset.


**The genes of interest (GOI) analysis workflow cannot be implemented until after users have successfully run the main downstream analysis core workflow as described in this repository's main [README.md](../README.md) file or have downloaded data from the [ScPCA portal](https://scpca.alexslemonade.org/).**

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Analysis overview](#analysis-overview)
- [Expected input](#expected-input)
- [Configure config file](#configure-config-file)
- [Gene mapping](#gene-mapping)
- [Running the workflow](#running-the-workflow)
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

1. The RDS file containing the normalized [output `SingleCellExperiment` object](../README.md#expected-output) from the core dowstream analyses workflow or the processed `SingleCellExperiment` object downloaded from the ScPCA portal (found in the `_processed.rds` file).
This `SingleCellExperiment` object must contain a log-normalized counts matrix in an assay named `logcounts`.
2. A project metadata tab-separated value (TSV) file containing relevant information about your data necessary for processing, the same metadata file used in the core workflow can be used here -- although only the `sample_id` and `library_id` columns are required for this workflow (see more on this in the ["Metadata file format" section](../README.md#metadata-file-format) and an example of this metadata file [here](../project-metadata/example-library-metadata.tsv)).
3. The genes of interest list, stored as a tab-separated value (TSV) file.
This file should contain at least one column named `gene_id` with the relevant gene identifiers, and can optionally contain an additional column named `gene_set` that denotes the gene set that each gene identifier belongs to (see an example of this genes of interest file [here](../example-data/goi-lists/sample01_goi_list.tsv)).

If working with data from the ScPCA portal, see our guide on preparing that data to run the genes of interest workflow [here](./additional-docs/working-with-scpca-portal-data.md).

## Configure config file

As in the main core workflow, we have provided an [example snakemake configuration file](../config/config.yaml), `config/config.yaml`, which defines all parameters needed to run the workflow.
Learn more about snakemake configuration files [here](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).

The config file contains two sets of parameters:

- **[Project-specific Parameters](../config/config.yaml#L3)**: This set of parameters are for specifying dataset or project related details. 
These parameters are **required** to run the workflow on your data.
- **[Processing Parameters](../config/goi_config.yaml)**: This set of parameters specify configurations for how to handle multiple gene identifier mappings, for example.
You can change them to explore your data but it is optional.
You can modify the relevant parameters by manually updating the `config/goi_config.yaml` file using a text editor of your choice.

To run the workflow on your data, modify the following parameters in the `config/config.yaml` and `config/goi_config.yaml` files:

| Parameter        | Description | Default value |
|------------------|-------------| --------------|
| `input_data_dir` | full path to the directory where the input data files can be found (default will be the `results_dir` used in the core workflow) | `"example_results"` |
| `results_dir` | full path to the directory where output files will be stored | `"example-results"` |
| `project_metadata` | full path to your specific project metadata TSV file (use the same `project_metadata` used in the prerequisite core workflow) | `"example-data/project-metadata/example-library-metadata.tsv"` |
| `goi_list` | full path to a tsv file containing the list of genes that are of interest | `"example-data/goi-lists/example_goi_list.tsv"` |
| `provided_identifier` | the type of gene identifiers used to populate the genes of interest list; example values that can implemented here include `"ENSEMBL"`, `"ENTREZID"`, `"SYMBOL"`; see more keytypes [here](https://jorainer.github.io/ensembldb/reference/EnsDb-AnnotationDbi.html) | `"SYMBOL"` |
| `overwrite` | a binary value indicating whether or not to overwrite existing output files | `TRUE` |

## Gene mapping

The first step in the workflow is to ensure that the input gene identifiers provided in the genes of interest list match the gene identifiers present in the `SingleCellExperiment` object. 
The `SingleCellExperiment` objects that are returned by the core workflow will contain gene names as Ensembl ids.
If the provided genes of interest list contains another identifier type (e.g., gene symbol, Entrez id) that does not match the gene names present in the `SingleCellExperiment` objects, then mapping must be performed.
The `perform_mapping` flag is used to indicate whether or not gene mapping will be performed and by default is set to `TRUE`.

To run the workflow without the gene mapping step, you will need to modify the `perform_mapping` parameter to be `FALSE` in the `config/goi_config.yaml` file.

The `config/goi_config.yaml` file also contains additional processing parameters like how to handle multiple gene identifier mappings. 
We have set default values for these parameters. 
Learn more about the [gene mapping parameters](../additional-docs/additional-parameters.md#genes-of-interest-analysis-parameters) and how to modify them.

## Running the workflow

The execution file with the genes of interest Snakemake workflow is named `goi.snakefile` and can be found in the root directory.
To tell snakemake to run the specific genes of interest workflow be sure to use the `--snakefile` or `-s` option followed by the name of the snakefile, `goi.snakefile`.

After you have successfully modified the required project-specific parameters in the config file and navigated to within the root directory of the `scpca-downstream-analyses` repository, you can run the clustering Snakemake workflow with just the `--cores` and `--use-conda` flags as in the following example: 

```
snakemake --snakefile goi.snakefile --cores 2 --use-conda
```

It is mandatory to specify the number of CPU cores for snakemake to use by using the [`--cores` flag](https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html?highlight=cores#step-1-specifying-the-number-of-used-threads).
If `--cores` is given without a number, all available cores are used to run the workflow.

**Note:** If you did not install dependencies [with conda via snakemake](#snakemakeconda-installation), you will need to remove the `--use-conda` flag.

You can also modify the config file parameters at the command line, rather than manually as recommended in the configure config file section above.
See our [command line options](../additional-docs/command-line-options.md) documentation for more information.

## Expected output

For each provided `SingleCellExperiment` RDS file and associated `library_id`, the workflow will return the following files in the specified output directory:

```
{library_id}_goi_report.html
{library_id}_goi_stats
├── {library_id}_mapped_genes.tsv
├── {library_id}_normalized_zscores.mtx
└── {library_id}_heatmap_annotation.rds
```

1. A `_mapped_genes.tsv` file with mapped genes of interest if the `--perform_mapping` is `TRUE` upon running the workflow. Otherwise this file will store just the provided genes of interest, and a new column named `sce_rownames_identifier` with the genes of interest identifiers copied over to the column.
2. A `_normalized_zscores.mtx` matrix file with the z-scored matrix calculated using the normalized data specific to the provided genes of interest.
3. A `_heatmap_annotation.rds` file with the annotations to be used when plotting the heatmap for the html report.
4. The `_goi_report.html` file, which is the summary html report with plots containing the GOI results.

You can download a ZIP file with an example of the output from running the genes of interest workflow, including the summary HTML report, the mapped genes TSV file, a `mtx` file containing z-scores, and the heatmap annotation RDS file [here](https://scpca-references.s3.amazonaws.com/example-data/scpca-downstream-analyses/goi_example_results.zip).