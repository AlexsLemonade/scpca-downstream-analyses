import pandas as pd

configfile: "config/config.yaml"

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
        expand(os.path.join(config["results_dir"], "{sample}/{library}_processed.rds"),
               zip,
               sample = SAMPLES,
               library = LIBRARY_ID),
        expand(os.path.join(config["results_dir"], "{sample}/{library}_core_analysis_report.html"),
               zip,
               sample = SAMPLES,
               library = LIBRARY_ID)


# Rule used for building conda & renv environment
rule build_renv:
    input: "renv.lock"
    output: "renv/.snakemake_timestamp"
    log: "logs/build_renv.log"
    conda: "envs/scpca-renv.yaml"
    shell:
      """
      Rscript -e "renv::restore(lockfile = '{input}')" &> {log}
      date -u -Iseconds  > {output}
      """


def get_input_rds_files(wildcards):
    lib_info = samples_information.set_index('library_id')
    return lib_info.loc[wildcards.library_id]['filepath']

rule filter_data:
    input:
        get_input_rds_files
    output:
        temp("{basedir}/{sample_id}/{library_id}_filtered.rds")
    log: "logs/{basedir}/{sample_id}/{library_id}/filter_data.log"
    conda: "envs/scpca-renv.yaml"
    shell:
        " Rscript 'core-analysis/01-filter-sce.R'"
        "  --sample_sce_filepath {input}"
        "  --sample_id {wildcards.sample_id}"
        "  --library_id {wildcards.library_id}"
        "  --mito_file {config[mito_file]}"
        "  --output_filepath {output}"
        "  --seed {config[seed]}"
        "  --gene_detected_row_cutoff {config[gene_detected_row_cutoff]}"
        "  --gene_means_cutoff {config[gene_means_cutoff]}"
        "  --mito_percent_cutoff {config[mito_percent_cutoff]}"
        "  --min_gene_cutoff {config[min_gene_cutoff]}"
        "  --umi_count_cutoff {config[umi_count_cutoff]}"
        "  --prob_compromised_cutoff {config[prob_compromised_cutoff]}"
        "  --filtering_method {config[filtering_method]}"
        "  --project_root $PWD"
        "  &> {log}"

rule normalize_data:
    input:
        "{basename}_filtered.rds"
    output:
        temp("{basename}_normalized.rds")
    log: "logs/{basename}/normalize_data.log"
    conda: "envs/scpca-renv.yaml"
    shell:
        " Rscript '{workflow.basedir}/core-analysis/02-normalize-sce.R'"
        "  --sce {input}"
        "  --seed {config[seed]}"
        "  --output_filepath {output}"
        "  --project_root $PWD"
        "  &> {log}"

rule dimensionality_reduction:
    input:
        "{basename}_normalized.rds"
    output:
        temp("{basename}_dimreduced.rds")
    log: "logs/{basename}/dimensionality_reduction.log"
    conda: "envs/scpca-renv.yaml"
    shell:
        " Rscript 'core-analysis/03-dimension-reduction.R'"
        "  --sce {input}"
        "  --seed {config[seed]}"
        "  --top_n {config[n_genes_pca]}"
        "  --output_filepath {output}"
        "  --overwrite"
        "  --project_root $PWD"
        "  &> {log}"

rule clustering:
    input:
        "{basename}_dimreduced.rds"
    output:
        "{basename}_processed.rds"
    log: "logs/{basename}/clustering.log"
    conda: "envs/scpca-renv.yaml"
    shell:
        " Rscript 'core-analysis/04-clustering.R'"
        "  --sce {input}"
        "  --seed {config[seed]}"
        "  --cluster_type {config[core_cluster_type]}"
        "  --nearest_neighbors {config[nearest_neighbors]}"
        "  --output_filepath {output}"
        "  --project_root $PWD"
        "  &> {log}"

rule generate_report:
    input:
        pre_processed_sce = get_input_rds_files,
        processed_sce =  "{basedir}/{sample_id}/{library_id}_processed.rds"
    output:
        "{basedir}/{sample_id}/{library_id}_core_analysis_report.html"
    log: "logs/{basedir}/{sample_id}/{library_id}/generate_report.log"
    conda: "envs/scpca-renv.yaml"
    shell:
        """
        Rscript -e \
        "rmarkdown::render('core-analysis/core-analysis-report-template.Rmd', \
                           clean = TRUE, \
                           output_dir = '{wildcards.basedir}/{wildcards.sample_id}', \
                           output_file = '{output}', \
                           params = list(library = '{wildcards.library_id}', \
                                         pre_processed_sce = '{input.pre_processed_sce}', \
                                         processed_sce = '{input.processed_sce}', \
                                         cluster_type = '{config[core_cluster_type]}', \
                                         nearest_neighbors = {config[nearest_neighbors]}, \
                                         mito_file = '{config[mito_file]}', \
                                         project_root = '$PWD'), \
                           envir = new.env())" \
        &> {log}
        """
