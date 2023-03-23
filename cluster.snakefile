import pandas as pd

configfile: "config/config.yaml"
configfile: "config/cluster_config.yaml"

# getting the samples information
if os.path.exists(config['project_metadata']):
  samples_information = pd.read_csv(config['project_metadata'], sep='\t', index_col=False)

  # get a list of the sample and library ids
  SAMPLES = list(samples_information['sample_id'])
  LIBRARY_ID = list(samples_information['library_id'])
else:
  # If the metadata file is missing, warn and fill with empty lists
  print(f"Warning: Project metadata file '{config['clustering_project_metadata']}' is missing.")
  samples_information = None
  SAMPLES = list()
  LIBRARY_ID = list()

rule target:
    input:
        expand(os.path.join(config["results_dir"], "{sample}/{library}_clustered_sce.rds"),
               zip,
               sample = SAMPLES,
               library = LIBRARY_ID),
        expand(os.path.join(config["results_dir"], "{sample}/{library}_clustering_stats"),
               zip,
               sample = SAMPLES,
               library = LIBRARY_ID),
        expand(os.path.join(config["results_dir"], "{sample}/{library}_clustering_report.html"),
               zip,
               sample = SAMPLES,
               library = LIBRARY_ID)

rule calculate_clustering:
    input:
        os.path.join(config["input_data_dir"], "{sample_id}/{library_id}_processed.rds")
    output:
        sce = os.path.join(config["results_dir"], "{sample_id}/{library_id}_clustered_sce.rds"),
        stats_dir = directory(os.path.join(config["results_dir"], "{sample_id}/{library_id}_clustering_stats"))
    log: os.path.join("logs", config["results_dir"], "{sample_id}/{library_id}/calculate_clustering.log")
    conda: "envs/scpca-renv.yaml"
    shell:
        " Rscript --vanilla 'optional-clustering-analysis/clustering-calculations.R'"
        "  --sce {input}"
        "  --library_id {wildcards.library_id}"
        "  --cluster_types {config[optional_cluster_types]}"
        "  --output_directory {output.stats_dir}"
        "  --output_sce {output.sce}"
        "  --seed {config[seed]}"
        "  --nearest_neighbors_min {config[nearest_neighbors_min]}"
        "  --nearest_neighbors_max {config[nearest_neighbors_max]}"
        "  --nearest_neighbors_increment {config[nearest_neighbors_increment]}"
        "  --project_root $PWD"
        "  --overwrite {config[overwrite_results]}"
        " &> {log}"

rule generate_cluster_report:
    input:
        processed_sce = "{basedir}/{library_id}_clustered_sce.rds",
        stats_dir = "{basedir}/{library_id}_clustering_stats"
    output:
        "{basedir}/{library_id}_clustering_report.html"
    log: "{basedir}/{library_id}/cluster_report.log"
    conda: "envs/scpca-renv.yaml"
    shell:
        """
        Rscript --vanilla -e "
          source(file.path('$PWD', 'utils', 'setup-functions.R'))
          setup_renv(project_filepath = '$PWD')
          rmarkdown::render('optional-clustering-analysis/clustering-report-template.Rmd',
                            clean = TRUE,
                            output_file = '{output}',
                            output_dir = dirname('{output}'),
                            params = list(library = '{wildcards.library_id}',
                                          processed_sce = '{input.processed_sce}',
                                          stats_dir = '{input.stats_dir}',
                                          cluster_type = '{config[optional_cluster_types]}',
                                          nearest_neighbors_min = {config[nearest_neighbors_min]},
                                          nearest_neighbors_max = {config[nearest_neighbors_max]},
                                          nearest_neighbors_increment = {config[nearest_neighbors_increment]}),
                           envir = new.env())
        " &> {log}
        """
