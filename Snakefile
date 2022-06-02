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
        expand(os.path.join(config["results_dir"], "{sample}/plots/{library}_{filtering_method}_cell_filtering.png"), 
               zip, 
               sample = SAMPLES, 
               library = LIBRARY_ID,
               filtering_method = FILTERING_METHOD)

def get_input_rds_files(wildcards):
    lib_info = samples_information.set_index('library_id')
    return lib_info.loc[wildcards.library_id]['filepath']
    
rule filter_data:
    input:
        get_input_rds_files
    output:
        downstream_filtered_rds = temp(os.path.join(config["results_dir"], "{sample_id}/{library_id}_{filtering_method}_downstream_processed_sce.rds")),
        plot = os.path.join(config["results_dir"], "{sample_id}/plots/{library_id}_{filtering_method}_cell_filtering.png")
    shell:
        "Rscript --vanilla 01-filter-sce.R"
        "  --sample_sce_filepath {input}"
        "  --sample_id {wildcards.sample_id}"
        "  --library_id {wildcards.library_id}"
        "  --mito_file {config[mito_file]}"
        "  --output_plots_directory $(dirname {output.plot})"
        "  --output_filepath {output.downstream_filtered_rds}"
        "  --seed 2021"
        "  --gene_detected_row_cutoff 5"
        "  --gene_means_cutoff 0.1"
        "  --mito_percent_cutoff 20"
        "  --detected_gene_cutoff 500"
        "  --umi_count_cutoff 500"
        "  --filtering_method {wildcards.filtering_method}"

rule normalize_data:
    input:
        "{basename}_downstream_processed_sce.rds"
    output:
        temp("{basename}_downstream_processed_normalized_sce.rds")
    shell:
        "Rscript --vanilla 02-normalize-sce.R"
        "  --sce {input}"
        "  --seed 2021"
        "  --output_filepath {output}"
        
rule dimensionality_reduction:
    input:
        "{basename}_downstream_processed_normalized_sce.rds"
    output:
        temp("{basename}_downstream_processed_normalized_reduced_sce.rds")
    shell:
        "Rscript --vanilla 03-dimension-reduction.R"
        "  --sce {input}"
        "  --seed 2021"
        "  --top_n {config[n_genes_pca]}"
        "  --output_filepath {output}"
        "  --overwrite"

rule clustering:
    input:
        "{basename}_downstream_processed_normalized_reduced_sce.rds"
    output:
        "{basename}_processed_sce.rds"
    params:
        cluster_type=config["cluster_type"],
        nearest_neighbors=config["nearest_neighbors"]
    shell:
        "Rscript --vanilla 04-clustering.R"
        "  --sce {input}"
        "  --seed 2021"
        "  --cluster_type {params.cluster_type}"
        "  --nearest_neighbors {params.nearest_neighbors}"
        "  --output_filepath {output}"
