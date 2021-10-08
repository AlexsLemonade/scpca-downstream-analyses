#!/bin/bash
set -euo pipefail
 
###############################################################################
# Run the marker genes analysis scripts including the optional marker genes
# mapping script in `utils`, and the prepare SingleCellExperiment object,
# filter, and normalize scripts (which are also optional based on the format of
# the data).

###############################################################################

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

RUN_MAPPING=${RUN_MAPPING:-FALSE}

main_data_dir="data/anderson-single-cell"
sample_name="GSM4186961"
sample_matrix="data/anderson-single-cell/GSE140819/GSM4186961_HTAPP-312-SMP-901_fresh-T1_channel1_raw_gene_bc_matrices_h5.h5"
sample_metadata="data/anderson-single-cell/GSE140819/GSM4186961_metadata_HTAPP-312-SMP-901_fresh-T1_channel1.csv.gz"
marker_genes="data/anderson-single-cell/marker-genes/nb_marker_genes.tsv"
mito_file="data/Homo_sapiens.GRCh38.103.mitogenes.txt"
input_identifiers="SYMBOL"
output_identifiers="ENSEMBL"
identifier_column_name="gene_symbol"
organism="Homo sapiens"
SEED=2021

# Run the gene identifier mapping script if `RUN_MAPPING = TRUE`
if [[ RUN_MAPPING == "TRUE" ]]
then
  Rscript --vanilla utils/map-marker-genes.R \
    --input_marker_gene_list ${marker_genes} \
    --input_identifiers ${input_identifiers} \
    --output_identifiers ${output_identifiers} \
    --identifier_column_name ${identifier_column_name} \
    --organism ${organism} \
    --multi_mappings "list" \
    --output_file "${main_data_dir}/marker-genes/mapped_marker_genes.tsv"
else
  echo "Skipping gene identifier mapping."
fi

# Run the 00 prepare SingleCellExperiment object script
Rscript --vanilla 00-prepare-sce.R \
 --sample_matrix_filepath ${sample_matrix} \
 --sample_metadata_filepath ${sample_metadata} \
 --input_file_type "h5" \
 --output_dir "${main_data_dir}/pre-filtered-sce" \
 --output_filename "${sample_name}_pre-filtered_sce.rds"
 
# Run the filtering script on pre-filtered SingleCellExperiment object (the
# implementation below incorporates `miQC` filtering)
Rscript --vanilla 01-filter-sce.R \
  --sample_sce_filepath "${main_data_dir}/pre-filtered-sce/${sample_name}_pre-filtered_sce.rds" \
  --sample_name ${sample_name} \
  --mito_file ${mito_file} \
  --output_data_directory "${main_data_dir}/filtered" \
  --output_plots_directory "plots" \
  --seed ${SEED} \
  --gene_detected_row_cutoff 5 \
  --gene_means_cutoff 0.1 \
  --filtering_method "miQC"

# Run the normalization script on filtered SingleCellExperiment object
Rscript --vanilla 02-normalize-sce.R \
  --sce "${main_data_dir}/filtered/filtered_${sample_name}_miQC_sce.rds" \
  --output_filepath "${main_data_dir}/normalized/normalized_${sample_name}_sce.rds" \
  --seed ${SEED}
