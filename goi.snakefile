import pandas as pd

configfile: "config/config.yaml"
configfile: "config/goi_config.yaml"

# getting the samples information
if os.path.exists(config['project_metadata']):
  samples_information = pd.read_csv(config['project_metadata'], sep='\t', index_col=False)

  # get a list of the sample and library ids
  SAMPLES = list(samples_information['sample_id'])
  LIBRARY_ID = list(samples_information['library_id'])
else:
  # If the metadata file is missing, warn and fill with empty lists
  print(f"Warning: Project metadata file '{config['project_metadata']}' is missing.")
  samples_information = None
  SAMPLES = list()
  LIBRARY_ID = list()

rule target:
    input:
        expand(os.path.join(config["results_dir"], "{sample}/{library}_goi_stats"),
               zip,
               sample = SAMPLES,
               library = LIBRARY_ID),
        expand(os.path.join(config["results_dir"], "{sample}/{library}_goi_report.html"),
               zip,
               sample = SAMPLES,
               library = LIBRARY_ID)

rule calculate_goi:
    input:
        os.path.join(config["input_data_dir"], "{sample_id}/{library_id}_processed.rds")
    output:
        output_dir = directory(os.path.join(config["results_dir"], "{sample_id}/{library_id}_goi_stats"))
    log: os.path.join("logs", config["results_dir"], "{sample_id}/{library_id}/calculate_goi.log")
    conda: "envs/scpca-renv.yaml"
    shell:
        " Rscript --vanilla 'optional-goi-analysis/goi-calculations.R'"
        "  --sce {input}"
        "  --library_id {wildcards.library_id}"
        "  --input_goi_list {config[goi_list]}"
        "  --organism '{config[organism]}'"
        "  --output_directory {output.output_dir}"
        "  --seed {config[seed]}"
        "  --provided_identifier {config[provided_identifier]}"
        "  --sce_rownames_identifier {config[sce_rownames_identifier]}"
        "  --multi_mappings {config[multi_mappings]}"
        "  --perform_mapping {config[perform_mapping]}"
        "  --project_root $PWD"
        " &> {log}"

rule generate_goi_report:
    input:
        processed_sce = os.path.join(config["input_data_dir"], "{sample_id}/{library_id}_processed.rds"),
        goi_dir = os.path.join(config["results_dir"], "{sample_id}/{library_id}_goi_stats")
    output:
        os.path.join(config["results_dir"], "{sample_id}/{library_id}_goi_report.html")
    log: os.path.join("logs", config["results_dir"], "{sample_id}/{library_id}/goi_report.log")
    conda: "envs/scpca-renv.yaml"
    shell:
        """
        Rscript --vanilla -e "
          source(file.path('$PWD', 'utils', 'setup-functions.R'))
          setup_renv(project_filepath = '$PWD')
          rmarkdown::render('optional-goi-analysis/goi-report-template.Rmd',
                            clean = TRUE,
                            output_file = '{output}',
                            output_dir = dirname('{output}'),
                            params = list(library = '{wildcards.library_id}',
                                          normalized_sce = '{input.processed_sce}',
                                          goi_input_directory = '{input.goi_dir}',
                                          project_root = '$PWD'),
                            envir = new.env())
        " &> {log}
        """

