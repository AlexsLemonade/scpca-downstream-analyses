# Working with ScPCA portal data

Downloads from the ScPCA Portal include gene expression data, a QC report, and associated metadata for each processed sample.
These files are delivered as a zip file.
To find more information, visit the [ScPCA portal](https://scpca.alexslemonade.org/).

## Preparing the downloaded data for core workflow

You can run the core downstream analysis workflow using the filtered RDS file and the metadata file found within the zip file downloaded from the ScPCA portal.
Before running the workflow, you will need to create a project metadata file as a tab-separated value (TSV) file that contains the relevant data for your input files needed to run the workflow.
The metadata file should contain the following columns:

- `sample_id`, the unique ID for each piece of tissue or sample that cells were obtained from, which can be found in the `scpca_sample_id` column of the downloaded `single_cell_metadata.tsv` file.
- `library_id`, the unique ID used for each set of cells that has been prepped and sequenced separately, which can be found in the `scpca_library_id` columns of the downloaded `single_cell_metadata.tsv` file.
- `filepath`, the full path to the `_filtered.rds` file containing the filtered `SingleCellExperiment` object.
Each library ID should have a unique `filepath`.

An example of the project metadata file would be as follows:

| sample_id | library_id | filepath |
| --------- | ---------- | --------- |
| SCPCS000122 | SCPCL000141 | example-data/SCPCS000122/SCPCL000141_filtered.rds |
| SCPCS000123	 | SCPCL000142 | example-data/SCPCS000123/SCPCL000142_filtered.rds |

### Running the workflow

The below code is an example of running the Snakemake workflow using the project-specific parameters.

```
snakemake --cores 2 \
  --use-conda \
  --config results_dir="<FULL PATH TO RESULTS DIRECTORY>" \
  project_metadata="<FULL PATH TO YOUR PROJECT METADATA TSV>" \
  mito_file="<FULL PATH TO MITOCHONDRIAL GENES TXT FILE>"
```

## Preparing the downloaded data to run additional modules

You can run the downstream analysis clustering and genes of interest workflows using the processed RDS file and the metadata file found within the zip file downloaded from the ScPCA portal.
Before running any of these workflows, you will need to create a project metadata file as a tab-separated value (TSV) file that contains the relevant data for your input files needed to run the workflow.
The metadata file should contain the following columns:

- `sample_id`, the unique ID for each piece of tissue or sample that cells were obtained from, which can be found in the `scpca_sample_id` column of the downloaded `single_cell_metadata.tsv` file.
- `library_id`, the unique ID used for each set of cells that has been prepped and sequenced separately, which can be found in the `scpca_library_id` columns of the downloaded `single_cell_metadata.tsv` file.
- `filepath`, the full path to the `_processed.rds` file containing the processed `SingleCellExperiment` object.
Each library ID should have a unique `filepath`.

An example of the project metadata file would be as follows:

| sample_id | library_id | filepath |
| --------- | ---------- | --------- |
| SCPCS000122 | SCPCL000141 | example-data/SCPCS000122/SCPCL000141_processed.rds |
| SCPCS000123	 | SCPCL000142 | example-data/SCPCS000123/SCPCL000142_processed.rds |

### Running the workflow

The below code is an example of running the Snakemake clustering workflow using the project-specific parameters.

```
snakemake --snakefile cluster.snakefile \ 
  --cores 2 \
  --use-conda \
  --config input_data_dir="<FULL PATH TO INPUT DATA DIRECTORY>" \
  results_dir="<FULL PATH TO RESULTS DIRECTORY>" \
  project_metadata="<FULL PATH TO YOUR PROJECT METADATA TSV>"
```

The below code is an example of running the Snakemake genes of interest workflow using the project-specific parameters.

```
snakemake --snakefile goi.snakefile \ 
  --cores 2 \
  --use-conda \
  --config input_data_dir="<FULL PATH TO INPUT DATA DIRECTORY>" \
  results_dir="<FULL PATH TO RESULTS DIRECTORY>" \
  project_metadata="<FULL PATH TO YOUR PROJECT METADATA TSV>"
```
