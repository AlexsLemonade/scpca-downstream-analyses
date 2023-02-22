# Optional Genes of Interest (GOI) Analysis

This directory includes a genes of interest analysis workflow to help users evaluate expression of a specific list of genes in a sample dataset.


**The genes of interest (GOI) analysis workflow cannot be implemented until after users have successfully run the main downstream analysis core workflow as described in this repository's main [README.md](../README.md) file.**

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Analysis overview](#analysis-overview)
- [Running the workflow](#running-the-workflow)
- [Expected input](#expected-input)
- [Parameters and config file](#parameters-and-config-file)
  - [Project-specific parameters](#project-specific-parameters)
- [Gene mapping](#gene-mapping)
  - [Gene mapping parameters](#gene-mapping-parameters)
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

## Running the workflow

The execution file with the genes of interest Snakemake workflow is named `goi.snakefile` and can be found in the root directory.
To tell snakemake to run the specific genes of interest workflow be sure to use the `--snakefile` or `-s` option followed by the name of the snakefile, `goi.snakefile`.
After navigating to within the root directory of the `scpca-downstream-analyses` repository, the below example command can be used to run the genes of interest workflow:

```
snakemake --snakefile goi.snakefile \ 
  --cores 2 \
  --use-conda \
  --config results_dir="<RELATIVE PATH TO RESULTS DIRECTORY>" \
  project_metadata="<RELATIVE PATH TO YOUR PROJECT METADATA TSV>"
```

**You will want to replace the paths for both `results_dir` and `project_metadata` to successfully run the workflow.** 
Where `results_dir` is the relative path to the directory where all results from running the workflow will be stored, and `project_metadata` is the relative path to the TSV file containing the relevant information about your input files.
See more information on project metadata in the [expected input section](#expected-input) below.

**Note:**  If you did not install dependencies [with conda via snakemake](../README.md#snakemakeconda-installation), you will need to remove the `--use-conda` flag.


## Expected input

To run this workflow, you will need to provide:

1. The RDS file containing the normalized [output `SingleCellExperiment` object](../README.md#expected-output) from the core dowstream analyses workflow.
This `SingleCellExperiment` object must contain a log-normalized counts matrix in an assay named `logcounts`.
2. A project metadata tab-separated value (TSV) file containing relevant information about your data necessary for processing, the same metadata file used in the core workflow should be used here (see more on this in the ["Metadata file format" section](../README.md#metadata-file-format) and an example of this metadata file [here](../project-metadata/example-library-metadata.tsv)).
3. The genes of interest list, stored as a tab-separated value (TSV) file.
This file should contain at least one column named `gene_id` with the relevant gene identifiers, and can optionally contain an additional column named `gene_set` that denotes the gene set that each gene identifier belongs to (see an example of this genes of interest file [here](../example-data/goi-lists/sample01_goi_list.tsv)).

## Parameters and config file

As in the main core workflow, we have provided a [configuration file](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html), `config/config.yaml` which sets the defaults for project-specific parameters needed to run the genes of interest workflow.

### Project-specific parameters

There are a set of parameters included in the provided configuration files (`config/config.yaml`, `config/goi_config.yaml`) that will always need to be specified when running the genes of interest workflow.
These parameters are specific to the project or dataset being processed.
These include the following parameters:

| Parameter        | Description | Default value |
|------------------|-------------| --------------|
| `results_dir` | relative path to the directory where output files will be stored (use the same `results_dir` used in the prerequisite core workflow) | `"example-results"` |
| `project_metadata` | relative path to your specific project metadata TSV file (use the same `project_metadata` used in the prerequisite core workflow) | `"example-data/project-metadata/example-library-metadata.tsv"` |
| `goi_list` | the file path to a tsv file containing the list of genes that are of interest | `"example-data/goi-lists/example_goi_list.tsv"` |
| `provided_identifier` | the type of gene identifiers used to populate the genes of interest list; example values that can implemented here include `"ENSEMBL"`, `"ENTREZID"`, `"SYMBOL"`; see more keytypes [here](https://jorainer.github.io/ensembldb/reference/EnsDb-AnnotationDbi.html) | `"SYMBOL"` |
| `overwrite` | a binary value indicating whether or not to overwrite existing output files | `TRUE` |

The above parameters can be modified at the command line by using the [`--config` flag](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).
You must also specify the number of CPU cores for snakemake to use by using the [`--cores` flag](https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html?highlight=cores#step-1-specifying-the-number-of-used-threads).
If `--cores` is given without a number, all available cores are used to run the workflow.
If you installed dependencies for the workflow [with conda via snakemake](../README.md#snakemakeconda-installation), you will need to provide the `--use_conda` flag as well.

The below code is an example of running the genes of interest workflow using the required parameters.

```
snakemake --snakefile goi.snakefile \ 
  --cores 2 \
  --use-conda \
  --config results_dir="<RELATIVE PATH TO RESULTS DIRECTORY>" \
  project_metadata="<RELATIVE PATH TO YOUR PROJECT METADATA TSV>" \
  goi_list="<RELATIVE PATH TO YOUR GENES OF INTEREST LIST>" \
  provided_identifier="<TYPE OF GENE IDENTIFIER USED IN GOI LIST>"
```

## Gene mapping

The first step in the workflow is to ensure that the input gene identifiers provided in the genes of interest list match the gene identifiers present in the `SingleCellExperiment` object. 
The `SingleCellExperiment` objects that are returned by the core workflow will contain gene names as Ensembl ids.
If the provided genes of interest list contains another identifier type (e.g., gene symbol, Entrez id) that does not match the gene names present in the `SingleCellExperiment` objects, then mapping must be performed.
The `perform_mapping` flag is used to indicate whether or not gene mapping will be performed and by default is set to `TRUE`.

To run the workflow with the gene mapping step, you will only need to provide the required parameters mentioned in the [project-specific parameters section](#project-specific-parameters).

```
snakemake --snakefile goi.snakefile \ 
  --cores 2 \
  --use-conda \
  --config results_dir="<RELATIVE PATH TO RESULTS DIRECTORY>" \
  project_metadata="<RELATIVE PATH TO YOUR PROJECT METADATA TSV>" \
  goi_list="<RELATIVE PATH TO YOUR GENES OF INTEREST LIST>" \
  provided_identifier="<TYPE OF GENE IDENTIFIER USED IN GOI LIST>"
```

To run the workflow without the gene mapping step, you can run the below command:

```
snakemake --snakefile goi.snakefile \ 
  --cores 2 \
  --use-conda \
  --config results_dir="<RELATIVE PATH TO RESULTS DIRECTORY>" \
  project_metadata="<RELATIVE PATH TO YOUR PROJECT METADATA TSV>" \
  goi_list="<RELATIVE PATH TO YOUR GENES OF INTEREST LIST>" \
  provided_identifier="<TYPE OF GENE IDENTIFIER USED IN GOI LIST>" \
  perform_mapping=FALSE
```

### Gene mapping parameters

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

You can download a ZIP file with an example of the output from running the genes of interest workflow, including the summary HTML report, the mapped genes TSV file, the z-scored matrix file, and the heatmap annotation RDS file [here](https://scpca-references.s3.amazonaws.com/example-data/scpca-downstream-analyses/goi_example_results.zip).