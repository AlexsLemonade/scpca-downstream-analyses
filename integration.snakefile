import pandas as pd

configfile: "config/config.yaml"
configfile: "config/integration_config.yaml"

# getting the samples information
if os.path.exists(config['integration_project_metadata']):
  samples_information = pd.read_csv(config['integration_project_metadata'], sep='\t', index_col=False)

  # get a list of the file paths and integration groups
  GROUP = list(samples_information['integration_group'])
else:
  # If the metadata file is missing, warn and fill with empty lists
  print(f"Warning: Project metadata file '{config['integration_project_metadata']}' is missing.")
  samples_information = None
  GROUP = list()

rule target:
    input:
        expand(os.path.join(config["results_dir"], "{group}_merged_sce.rds"),
               zip,
               group = GROUP)

rule merge_sces:
    input:
        config["integration_project_metadata"]
    output:
        os.path.join(config["results_dir"], "{group}_merged_sce.rds")
    log: os.path.join("logs", config["results_dir"], "{group}_integration.log")
    conda: "envs/scpca-renv.yaml"
    shell:
        " Rscript --vanilla 'optional-integration-analysis/merge-sce.R'"
        "  --input_metadata_tsv {input}"
        "  --integration_group {wildcards.group}"
        "  --output_sce_file {output}"
        "  --n_hvg {config[n_hvg]}"
        "  --threads {config[threads]}"
        "  --project_root $PWD"
        " &> {log}"
