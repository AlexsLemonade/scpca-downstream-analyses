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
        expand(os.path.join(config["data_dir"], "results/{sample}/{library}_{filtering_method}_sce.rds"), 
               zip, 
               sample = SAMPLES, 
               library = LIBRARY_ID, 
               filtering_method = FILTERING_METHOD),
        expand(os.path.join(config["data_dir"], "results/{sample}/{library}_{filtering_method}_normalized_sce.rds"), 
               zip, 
               sample = SAMPLES, 
               library = LIBRARY_ID, 
               filtering_method = FILTERING_METHOD),
        expand(os.path.join(config["data_dir"], "results/{sample}/plots/{library}_{filtering_method}_cell_filtering.png"), 
               zip, 
               sample = SAMPLES, 
               library = LIBRARY_ID, 
               filtering_method = FILTERING_METHOD),
               

rule filter_data:
    input:
        "{base_dir}/{sample_id}/{library_id}_filtered.rds"
    output:
        filtered_rds = "{base_dir}/results/{sample_id}/{library_id}_{filtering_method}_sce.rds",
        plot = "{base_dir}/results/{sample_id}/plots/{library_id}_{filtering_method}_cell_filtering.png"
    shell:
        "Rscript --vanilla 01-filter-sce.R"
        "  --sample_sce_filepath {input}"
        "  --sample_name {wildcards.library_id}"
        "  --mito_file {config[mito_file]}"
        "  --output_plots_directory $(dirname {output.plot})"
        "  --output_filepath {output.filtered_rds}"
        "  --seed 2021"
        "  --gene_detected_row_cutoff 5"
        "  --gene_means_cutoff 0.1"
        "  --mito_percent_cutoff 20"
        "  --detected_gene_cutoff 500"
        "  --umi_count_cutoff 500"
        "  --filtering_method {wildcards.filtering_method}"

rule normalize_data:
    input:
        "{base_dir}/results/{sample_id}/{library_id}_{filtering_method}_sce.rds"
    output:
        normalized_rds = "{base_dir}/results/{sample_id}/{library_id}_{filtering_method}_normalized_sce.rds"
    shell:
        "Rscript --vanilla 02-normalize-sce.R"
        "  --sce {input}"
        "  --seed 2021"
        "  --output_filepath {output.normalized_rds}"
