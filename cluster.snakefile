import pandas as pd

configfile: "config.yaml"

# getting the samples information
if os.path.exists(config['project_metadata']):
  samples_information = pd.read_csv(config['project_metadata'], sep='\t', index_col=False)

  # get a list of the sample and library ids
  SAMPLES = list(samples_information['sample_id'])
  LIBRARY_ID = list(samples_information['library_id'])
  FILTERING_METHOD = list(samples_information['filtering_method'])
else:
  # If the metadata file is missing, warn and fill with empty lists
  print(f"Warning: Project metadata file '{config['clustering_project_metadata']}' is missing.")
  samples_information = None
  SAMPLES = list()
  LIBRARY_ID = list()
  FILTERING_METHOD = list()

rule target:
    input:
        expand(os.path.join(config["results_dir"], "{sample}/{library}_{filter_method}_clustered_sce.rds"),
               zip,
               sample = SAMPLES,
               library = LIBRARY_ID,
               filter_method = FILTERING_METHOD),
        expand(os.path.join(config["results_dir"], "{sample}/{library}_{filter_method}_clustering_stats"),
               zip,
               sample = SAMPLES,
               library = LIBRARY_ID,
               filter_method = FILTERING_METHOD),
        expand(os.path.join(config["results_dir"], "{sample}/{library}_{filter_method}_clustering_report.html"),
               zip,
               sample = SAMPLES,
               library = LIBRARY_ID,
               filter_method = FILTERING_METHOD)       

rule calculate_clustering:
    input:
        "{basedir}/{sample_id}/{library_id}_{filter_method}_processed_sce.rds"
    output:
        sce = "{basedir}/{sample_id}/{library_id}_{filter_method}_clustered_sce.rds",
        stats_dir = directory("{basedir}/{sample_id}/{library_id}_{filter_method}_clustering_stats")
    conda: "envs/scpca-renv.yaml"
    shell:
        "R_PROFILE_USER='{workflow.basedir}/.Rprofile'"
        " Rscript '{workflow.basedir}/optional-clustering-analysis/clustering-calculations.R'"
        "  --sce {input}"
        "  --library_id {wildcards.library_id}"
        "  --cluster_types {config[cluster_types]}"
        "  --output_directory {output.stats_dir}"
        "  --output_sce {output.sce}"
        "  --seed {config[seed]}"
        "  --nearest_neighbors_min {config[nearest_neighbors_min]}"
        "  --nearest_neighbors_max {config[nearest_neighbors_max]}"
        "  --nearest_neighbors_increment {config[nearest_neighbors_increment]}"
        "  --project_root {workflow.basedir}"
        "  --overwrite {config[overwrite_results]}"
        
rule generate_cluster_report:
    input:
        processed_sce = "{basedir}/{sample_id}/{library_id}_{filter_method}_clustered_sce.rds",
        stats_dir = "{basedir}/{sample_id}/{library_id}_{filter_method}_clustering_stats"
    output:
        "{basedir}/{sample_id}/{library_id}_{filter_method}_clustering_report.html"
    conda: "envs/scpca-renv.yaml"
    shell:
        """
        R_PROFILE_USER='{workflow.basedir}/.Rprofile' \
        Rscript -e \
        "rmarkdown::render('{workflow.basedir}/optional-clustering-analysis/clustering-report-template.Rmd', \
                           clean = TRUE, \
                           output_file = '{output}', \
                           output_dir = dirname('{output}'), \
                           params = list(library = '{wildcards.library_id}', \
                                         processed_sce = '{input.processed_sce}', \
                                         stats_dir = '{input.stats_dir}', \
                                         cluster_type = '{config[cluster_type]}', \
                                         nearest_neighbors_min = {config[nearest_neighbors_min]}, \
                                         nearest_neighbors_max = {config[nearest_neighbors_max]}, \
                                         nearest_neighbors_increment = {config[nearest_neighbors_increment]}), \
                           envir = new.env())"
        """
