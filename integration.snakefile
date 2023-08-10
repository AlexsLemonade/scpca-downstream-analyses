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
        expand(os.path.join(config["results_dir"], "{group}_integrated_sce.rds"),
               zip,
               group = GROUP),
        expand(os.path.join(config["results_dir"], "{group}_integration_report.html"),
               zip,
               group = GROUP)

rule merge_sces:
    input:
        config["integration_project_metadata"]
    output:
        temp(os.path.join(config["results_dir"], "{group}_merged_sce.rds"))
    log: os.path.join("logs", config["results_dir"], "{group}_merge_sce.log")
    conda: "envs/scpca-renv.yaml"
    shell:
        " Rscript --vanilla 'optional-integration-analysis/merge-sce.R'"
        "  --input_metadata_tsv {input}"
        "  --integration_group {wildcards.group}"
        "  --output_sce_file {output}"
        "  --n_hvg {config[n_genes_pca]}"
        "  --threads {config[threads]}"
        "  --project_root $PWD"
        " &> {log}"
        
rule perform_integration:
    input:
        "{basedir}/{group}_merged_sce.rds"
    output:
        "{basedir}/{group}_integrated_sce.rds"
    log: "logs/{basedir}/{group}_perform_integration.log"
    conda: "envs/scpca-renv.yaml"
    shell:
        " Rscript --vanilla 'optional-integration-analysis/perform-integration.R'"
        "  --merged_sce_file {input}"
        "  --integration_method {config[integration_method]}"
        "  --fastmnn_auto_merge"
        "  --output_sce_file {output}"
        "  --project_root $PWD"
        " &> {log}"
        
rule generate_integration_report:
    input:
        integrated_sce = "{basedir}/{group}_integrated_sce.rds"
    output:
        "{basedir}/{group}_integration_report.html"
    log: "logs/{basedir}/{group}/integration_report.log"
    conda: "envs/scpca-renv.yaml"
    shell:
        """
        Rscript --vanilla -e "
          rmarkdown::render('optional-integration-analysis/integration-report-template.Rmd',
                            clean = TRUE,
                            output_file = '{output}',
                            output_dir = dirname('{output}'),
                            params = list(integration_group = '{wildcards.group}',
                                          integrated_sce = '{input.integrated_sce}',
                                          integration_method = '{config[integration_method]}',
                                          batch_column = '{config[batch_column]}',
                                          cell_id_column = '{config[cell_id_column]}'),
                           envir = new.env())
        " &> {log}
        """
