import pandas as pd

configfile: "config.yaml"

# getting the samples information
if os.path.exists(config['clustering_project_metadata']):
  samples_information = pd.read_csv(config['clustering_project_metadata'], sep='\t', index_col=False)

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
        expand(os.path.join(config["results_dir"], "{sample}/{library}_processed_sce_clustering.rds"),
               zip,
               sample = SAMPLES,
               library = LIBRARY_ID),
        expand(os.path.join(config["results_dir"], "{sample}/clustering_stats/{library}_clustering_all_validity_stats.tsv"),
               zip,
               sample = SAMPLES,
               library = LIBRARY_ID),
        expand(os.path.join(config["results_dir"], "{sample}/clustering_stats/{library}_clustering_summary_validity_stats.tsv"),
               zip,
               sample = SAMPLES,
               library = LIBRARY_ID),
        expand(os.path.join(config["results_dir"], "{sample}/clustering_stats/{library}_clustering_summary_stability_stats.tsv"),
               zip,
               sample = SAMPLES,
               library = LIBRARY_ID),
        expand(os.path.join(config["results_dir"], "{sample}/{library}_optional_clustering_report.html"),
               zip,
               sample = SAMPLES,
               library = LIBRARY_ID)       


def get_input_rds_files(wildcards):
    lib_info = samples_information.set_index('library_id')
    return lib_info.loc[wildcards.library_id]['filepath']


rule calculate_clustering:
    input:
        get_input_rds_files
    output:
        os.path.join(config['results_dir'], "{sample_id}/{library_id}_processed_sce_clustering.rds"),
        os.path.join(config['results_dir'], "{sample_id}/clustering_stats/{library_id}_clustering_all_validity_stats.tsv"),
        os.path.join(config['results_dir'], "{sample_id}/clustering_stats/{library_id}_clustering_summary_validity_stats.tsv"),
        os.path.join(config['results_dir'], "{sample_id}/clustering_stats/{library_id}_clustering_summary_stability_stats.tsv")
    conda: "envs/scpca-renv.yaml"
    shell:
        "R_PROFILE_USER='{workflow.basedir}/.Rprofile'"
        " Rscript '{workflow.basedir}/optional-clustering-analysis/clustering-calculations.R'"
        "  --sce {input}"
        "  --library_id {wildcards.library_id}"
        "  --cluster_type {config[cluster_type]}"
        "  --output_directory {config[results_dir]}/{wildcards.sample_id}"
        "  --seed {config[seed]}"
        "  --nearest_neighbors_min {config[nearest_neighbors_min]}"
        "  --nearest_neighbors_max {config[nearest_neighbors_max]}"
        "  --nearest_neighbors_increment {config[nearest_neighbors_increment]}"
        "  --project_root {workflow.basedir}"
        "  --overwrite {config[overwrite_results]}"
        
rule generate_cluster_report:
    input:
        processed_sce = "{base_dir}/{sample_id}/{library_id}_processed_sce_clustering.rds"
    output:
        "{base_dir}/{sample_id}/{library_id}_optional_clustering_report.html"
    conda: "envs/scpca-renv.yaml"
    shell:
        """
        R_PROFILE_USER='.Rprofile' \
        Rscript -e \
        "rmarkdown::render('{workflow.basedir}/optional-clustering-analysis/clustering-report-template.Rmd', \
                           clean = TRUE, \
                           output_file = '{output}', \
                           output_dir = dirname('{output}'), \
                           params = list(library = '{wildcards.library_id}', \
                                         processed_sce = '{input.processed_sce}', \
                                         stats_dir = '{wildcards.base_dir}/{wildcards.sample_id}/clustering_stats', \
                                         cluster_type = '{config[cluster_type]}', \
                                         nearest_neighbors_min = {config[nearest_neighbors_min]}, \
                                         nearest_neighbors_max = {config[nearest_neighbors_max]}, \
                                         nearest_neighbors_increment = {config[nearest_neighbors_increment]}), \
                           envir = new.env())"
        """
