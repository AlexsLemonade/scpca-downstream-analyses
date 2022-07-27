import pandas as pd

configfile: "config.yaml"

# getting the samples information
samples_information = pd.read_csv(config["project_metadata"], sep='\t', index_col=False)

# get a list of the sample and library ids
SAMPLES = list(samples_information['sample_id'])
LIBRARY_ID = list(samples_information['library_id'])
FILTERING_METHOD = list(samples_information['filtering_method'])

rule target:
    input:
        expand(os.path.join(config["results_dir"], "{sample}/{library}_{filtering_method}_processed_sce.rds"),
               zip,
               sample = SAMPLES,
               library = LIBRARY_ID,
               filtering_method = FILTERING_METHOD),
        expand(os.path.join(config["results_dir"], "{sample}/{library}_{filtering_method}_core_analysis_report.html"),
               zip,
               sample = SAMPLES,
               library = LIBRARY_ID,
               filtering_method = FILTERING_METHOD)


def get_input_rds_files(wildcards):
    lib_info = samples_information.set_index('library_id')
    return lib_info.loc[wildcards.library_id]['filepath']

# Dummy rule used for building conda & renv environment
rule build_renv:
    input: "renv.lock"
    output: "renv/.snakemake_timestamp"
    conda: "envs/scpca-renv.yaml"
    shell: "date -u -Iseconds  > {output}"


rule filter_data:
    input:
        get_input_rds_files
    output:
        downstream_filtered_rds = temp(os.path.join(config["results_dir"], "{sample_id}/{library_id}_{filtering_method}_downstream_processed_sce.rds"))
    conda: "envs/scpca-renv.yaml"
    shell:
        "Rscript --vanilla {workflow.basedir}/01-filter-sce.R"
        "  --sample_sce_filepath {input}"
        "  --sample_id {wildcards.sample_id}"
        "  --library_id {wildcards.library_id}"
        "  --mito_file {config[mito_file]}"
        "  --output_filepath {output.downstream_filtered_rds}"
        "  --seed {config[seed]}"
        "  --gene_detected_row_cutoff {config[gene_detected_row_cutoff]}"
        "  --gene_means_cutoff {config[gene_means_cutoff]}"
        "  --mito_percent_cutoff {config[mito_percent_cutoff]}"
        "  --detected_gene_cutoff {config[detected_gene_cutoff]}"
        "  --umi_count_cutoff {config[umi_count_cutoff]}"
        "  --prob_compromised_cutoff {config[prob_compromised_cutoff]}"
        "  --filtering_method {wildcards.filtering_method}"
        "  --project_root {workflow.basedir}"

rule normalize_data:
    input:
        "{basename}_downstream_processed_sce.rds"
    output:
        temp("{basename}_downstream_processed_normalized_sce.rds")
    conda: "envs/scpca-renv.yaml"
    shell:
        "Rscript --vanilla {workflow.basedir}/02-normalize-sce.R"
        "  --sce {input}"
        "  --seed {config[seed]}"
        "  --output_filepath {output}"
        "  --project_root {workflow.basedir}"

rule dimensionality_reduction:
    input:
        "{basename}_downstream_processed_normalized_sce.rds"
    output:
        temp("{basename}_downstream_processed_normalized_reduced_sce.rds")
    conda: "envs/scpca-renv.yaml"
    shell:
        "Rscript --vanilla {workflow.basedir}/03-dimension-reduction.R"
        "  --sce {input}"
        "  --seed {config[seed]}"
        "  --top_n {config[n_genes_pca]}"
        "  --output_filepath {output}"
        "  --overwrite"
        "  --project_root {workflow.basedir}"

rule clustering:
    input:
        "{basename}_downstream_processed_normalized_reduced_sce.rds"
    output:
        "{basename}_processed_sce.rds"
    conda: "envs/scpca-renv.yaml"
    shell:
        "Rscript --vanilla {workflow.basedir}/04-clustering.R"
        "  --sce {input}"
        "  --seed {config[seed]}"
        "  --cluster_type {config[cluster_type]}"
        "  --nearest_neighbors {config[nearest_neighbors]}"
        "  --output_filepath {output}"
        "  --project_root {workflow.basedir}"

rule generate_report:
    input:
        pre_processed_sce = get_input_rds_files,
        processed_sce =  "{basedir}/{library_id}_{filtering_method}_processed_sce.rds"
    output:
        "{basedir}/{library_id}_{filtering_method}_core_analysis_report.html"
    conda: "envs/scpca-renv.yaml"
    shell:
        """
        Rscript -e \
        "rmarkdown::render('{workflow.basedir}/core-analysis-report-template.Rmd', \
                           clean = TRUE, \
                           output_file = '{output}', \
                           params = list(library = '{wildcards.library_id}', \
                                         pre_processed_sce = '{input.pre_processed_sce}', \
                                         processed_sce = '{input.processed_sce}', \
                                         cluster_type = '{config[cluster_type]}', \
                                         nearest_neighbors = {config[nearest_neighbors]}, \
                                         mito_file = '{config[mito_file]}', \
                                         project_root = '{workflow.basedir}'), \
                           envir = new.env())"
        """
