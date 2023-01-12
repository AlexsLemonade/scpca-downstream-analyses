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
               library = LIBRARY_ID)

rule calculate_goi:
    input:
        "{basedir}/{library_id}_processed_sce.rds"
    output:
        output_dir = directory("{basedir}/{library_id}_goi_stats")
    log: "logs/{basedir}/{library_id}/calculate_goi.log"
    conda: "envs/scpca-renv.yaml"
    shell:
        " Rscript 'optional-goi-analysis/goi-calculations.R'"
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
