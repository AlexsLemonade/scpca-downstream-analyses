import pandas as pd

configfile: os.path.join(workflow.basedir, "config.yaml")

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
               filtering_method = FILTERING_METHOD)

def get_input_rds_files(wildcards):
    lib_info = samples_information.set_index('library_id')
    return lib_info.loc[wildcards.library_id]['filepath']
    
rule filter_data:
    input:
        get_input_rds_files
    output:
        downstream_filtered_rds = temp(os.path.join(config["results_dir"], "{sample_id}/{library_id}_{filtering_method}_downstream_processed_sce.rds"))
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

rule normalize_data:
    input:
        "{basename}_downstream_processed_sce.rds"
    output:
        temp("{basename}_downstream_processed_normalized_sce.rds")
    shell:
        "Rscript --vanilla {workflow.basedir}/02-normalize-sce.R"
        "  --sce {input}"
        "  --seed {config[seed]}"
        "  --output_filepath {output}"
        
rule dimensionality_reduction:
    input:
        "{basename}_downstream_processed_normalized_sce.rds"
    output:
        temp("{basename}_downstream_processed_normalized_reduced_sce.rds")
    shell:
        "Rscript --vanilla {workflow.basedir}/03-dimension-reduction.R"
        "  --sce {input}"
        "  --seed {config[seed]}"
        "  --top_n {config[n_genes_pca]}"
        "  --output_filepath {output}"
        "  --overwrite"

rule clustering:
    input:
        "{basename}_downstream_processed_normalized_reduced_sce.rds"
    output:
        "{basename}_processed_sce.rds"
    shell:
        "Rscript --vanilla {workflow.basedir}/04-clustering.R"
        "  --sce {input}"
        "  --seed {config[seed]}"
        "  --cluster_type {config[cluster_type]}"
        "  --nearest_neighbors {config[nearest_neighbors]}"
        "  --output_filepath {output}"
