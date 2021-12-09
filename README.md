# scPCA downstream analyses

This repository is intended to store our pipeline for scPCA downstream analyses.
This pipeline has been created in a modular format with the intent for users to begin processing at the module or script that is most appropriate for their data.
This pipeline also includes an optional genes of interest analysis, when a genes of interest list is provided.

Note that R 4.1 is required for running our pipeling, along with Bioconductor 3.14.

## How to run the analysis

### How to run the main scPCA downstream analysis pipeline

The core downstream scPCA analysis pipeline, which includes filtering, normalization, and dimension reduction, can be implemented by running the following command from the main directory:

`bash run-scpca-downstream-analysis.sh`

Note that you will want to provide the metadata TSV file related to your project (with `sample_id`, `library_id`, and `filtering_method`) to the `--sample_metadata` command line flag.

### How to run the optional genes of interest analysis pipeline

There is an optional genes of interest analysis pipeline in the `optional-goi-analysis` subdirectory of this repository.

To run this analysis, you can run the following command from the main directory:

```
bash optional-goi-analysis/run-provided-goi-analysis.sh \
 --output_dir "path/to/output-results" \
 --sample_name "sample"  \
 --sample_matrix "path/to/sample/matrix" \
 --sample_metadata "path/to/sample/metadata" \
 --goi_list "path/to/goi-list
```

Where `goi_list` is the path to the genes of interest TSV file relevant to your dataset.
