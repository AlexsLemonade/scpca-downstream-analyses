# Working with ScPCA portal data

Downloads from the ScPCA Portal include gene expression data, a QC report, and associated metadata for each processed sample.
Each sample's zip file should include a `SingleCellExperiment` object stored as `_filtered.rds` that can be used as input to the core analysis workflow, or to skip directly to the additional analysis modules, use the processed `SingleCellExperiment` object stored as `_processed.rds`.
To find more information on data available in the portal, visit the [ScPCA portal](https://scpca.alexslemonade.org/).

## Preparing the downloaded data for core workflow

The core downstream analysis workflow can take as input the `SingleCellExperiment` objects that are available in the filtered RDS files on the portal, files ending in `_filtered.rds`.
You will also need information obtained in the metadata file found within the zip file downloaded from the ScPCA portal, `single_cell_metadata.tsv`.
Before running the workflow, you will need to create a project metadata file as a tab-separated value (TSV) file that contains the relevant data for your input files needed to run the workflow.
The metadata file should contain the following columns:

- `sample_id`, the unique ID for each piece of tissue or sample that cells were obtained from, which can be found in the `scpca_sample_id` column of the downloaded `single_cell_metadata.tsv` file.
- `library_id`, the unique ID used for each set of cells that has been prepped and sequenced separately, which can be found in the `scpca_library_id` columns of the downloaded `single_cell_metadata.tsv` file.
- `filepath`, the full path to the `_filtered.rds` file containing the filtered `SingleCellExperiment` object downloaded from the ScPCA portal.
Each library ID should have a unique `filepath` and will probably look something like this for each file: 
`<full path to directory storing data downloaded from ScPCA/scpca_sample_id/scpca_library_id_filtered.rds`.

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
  project_metadata="<FULL PATH TO YOUR PROJECT METADATA TSV>"
```

## Preparing the downloaded data to run additional modules

You can run any of the additional modules available (e.g., clustering or the genes of interest workflows) using the processed RDS file (ending in `_processed.rds`) and the metadata file (`single_cell_metadata.tsv`) found within the zip file downloaded from the ScPCA portal.
Before running any of these workflows, you will need to create a project metadata file as a tab-separated value (TSV) file that contains the relevant data for your input files needed to run the workflow.
The metadata file should contain the following columns:

- `sample_id`, the unique ID for each piece of tissue or sample that cells were obtained from, which can be found in the `scpca_sample_id` column of the downloaded `single_cell_metadata.tsv` file.
- `library_id`, the unique ID used for each set of cells that has been prepped and sequenced separately, which can be found in the `scpca_library_id` columns of the downloaded `single_cell_metadata.tsv` file.

An example of the project metadata file would be as follows:

| sample_id | library_id |
| --------- | ---------- |
| SCPCS000122 | SCPCL000141 |
| SCPCS000123	 | SCPCL000142 |

**Note:** Any of the additional modules will also need the parameter `input_data_dir` to be specified.
The required `input_data_dir` should contain the downloaded input files as follows:

```
input_data_dir
   ├── SCPCS000122
   │   ├── SCPCL000141_processed.rds
   │   └── single_cell_metadata.tsv
   ├── SCPCS000123
   │   ├── SCPCL000142_processed.rds
   │   └── single_cell_metadata.tsv
```

### Running the workflow

The below code is an example of running the additional workflow(s) using the project-specific parameters, where `<module.snakefile>` would be replaced with the name of the snakefile associated with the desired module.

```
snakemake --snakefile <module.snakefile> \ 
  --cores 2 \
  --use-conda \
  --config input_data_dir="<FULL PATH TO INPUT DATA DIRECTORY>" \
  results_dir="<FULL PATH TO RESULTS DIRECTORY>" \
  project_metadata="<FULL PATH TO YOUR PROJECT METADATA TSV>"
```

To run the clustering workflow, see ["running the workflow"](../optional-clustering-analysis/README.md#running-the-workflow) in the clustering module `README.md`.

To run the genes of interest workflow, see ["running the workflow"](../optional-goi-analysis/README.md#running-the-workflow) in the genes of interest module `README.md`.