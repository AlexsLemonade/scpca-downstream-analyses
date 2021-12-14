import pandas as pd

# getting the samples information
samples_information = pd.read_csv("project-metadata/gawad-library-metadata.tsv", sep='\t', index_col=False)

# get a list of the sample and library ids
SAMPLES = list(samples_information['sample_id'])
LIBRARY_ID = list(samples_information['library_id'])
FILTERING_METHOD = list(samples_information['filtering_method'])
          
rule target:
    input:
        expand("data/Gawad_processed_data/results/{ids[0]}/{ids[1]}_filtered_{ids[2]}_sce.rds", ids = zip(SAMPLES, LIBRARY_ID, FILTERING_METHOD)),
        expand("data/Gawad_processed_data/results/{ids[0]}/plots/{ids[1]}_{ids[2]}_cell_filtering.png", ids = zip(SAMPLES, LIBRARY_ID, FILTERING_METHOD))

rule filter_data:
    input:
        "{base_dir}/{sample_id}/{library_id}_filtered.rds"
    output:
        filtered_rds = "{base_dir}/results/{sample_id}/{library_id}_filtered_{filtering_method}_sce.rds",
        plot = "{base_dir}/results/{sample_id}/plots/{library_id}_{filtering_method}_cell_filtering.png"
    shell:
        "Rscript --vanilla 01-filter-sce.R"
        "  --sample_sce_filepath {input}"
        "  --sample_name {wildcards.library_id}"
        "  --mito_file data/Homo_sapiens.GRCh38.103.mitogenes.txt"
        "  --output_plots_directory $(dirname {output.plot})"
        "  --output_filepath {output.filtered_rds}"
        "  --seed 2021"
        "  --gene_detected_row_cutoff 5"
        "  --gene_means_cutoff 0.1"
        "  --mito_percent_cutoff 20"
        "  --detected_gene_cutoff 500"
        "  --umi_count_cutoff 500"
        "  --filtering_method {wildcards.filtering_method}"
        "  --project_directory ."
