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

## Running the core workflow at the command line

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

You can find a list of parameters that can be modified at the command line for the core workflow under the [`project-specific parameters`](../README.md#project-specific-parameters) and [`processing parameters`](../README.md#processing-parameters) sections of the `README.md`.


## Running the additional modules at the command line

All of the parameters for each of the additional analysis modules can also be modified at the command line by using the [`--config` flag](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).
The below code is an example of running the additional workflow(s) using the project-specific parameters, where `<module.snakefile>` would be replaced with the name of the snakefile associated with the desired module.

```
snakemake --snakefile <module.snakefile> \ 
  --cores 2 \
  --use-conda \
  --config input_data_dir="<FULL PATH TO INPUT DATA DIRECTORY>" \
  results_dir="<FULL PATH TO RESULTS DIRECTORY>" \
  project_metadata="<FULL PATH TO YOUR PROJECT METADATA TSV>"
```

See the [clustering module `parameters and config file` documentation](../optional-clustering-analysis/README.md#parameters-and-config-file) and the [genes of interest module `parameters and config file` documentation](../optional-goi-analysis/README.md#parameters-and-config-file) docs for a list of parameters that can be modified at the command line for each of the respective modules.

