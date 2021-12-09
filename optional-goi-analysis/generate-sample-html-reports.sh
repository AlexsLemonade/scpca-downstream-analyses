#!/bin/bash
set -euo pipefail
 
###############################################################################
# Run the scPCA downstream analysis scripts including the optional genes of
# interest mapping script in `utils`, and the prepare SingleCellExperiment object,
# filter, and normalize scripts (which are also optional based on the format of
# the data).

# Usage
# Run the scPCA downstream analysis pipeline from start to finish for all samples in
# the GSE140819 public NB dataset
# bash generate-sample-html-reports.sh

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

sample_data_dir="../data/anderson-single-cell/GSE140819"
sample_names=(GSM4186961 GSM4186962 GSM4186963 GSM4186964 GSM4186965 GSM4186966 GSM4186967 GSM4186968 GSM4186969 GSM4186970)

declare -A sample_matrix=(
  [GSM4186961]="GSM4186961_HTAPP-312-SMP-901_fresh-T1_channel1_raw_gene_bc_matrices_h5.h5"
  [GSM4186962]="GSM4186962_HTAPP-312-SMP-902_fresh-C4-T2_channel1_raw_gene_bc_matrices_h5.h5"
  [GSM4186963]="GSM4186963_HTAPP-656-SMP-3481_fresh-T1_channel1_raw_gene_bc_matrices_h5.h5"
  [GSM4186964]="GSM4186964_HTAPP-StJude-SMP-PDX1_cell_channel1_raw_gene_bc_matrices_h5.h5"
  [GSM4186965]="GSM4186965_HTAPP-244-SMP-451_CST_channel1_raw_gene_bc_matrices_h5.h5"
  [GSM4186966]="GSM4186966_HTAPP-244-SMP-451_EZ_channel1_raw_gene_bc_matrices_h5.h5"
  [GSM4186967]="GSM4186967_HTAPP-244-SMP-451_NST_channel1_raw_gene_bc_matrices_h5.h5"
  [GSM4186968]="GSM4186968_HTAPP-244-SMP-451_TST_channel1_raw_gene_bc_matrices_h5.h5"
  [GSM4186969]="GSM4186969_HTAPP-656-SMP-3481_TST_channel1_raw_gene_bc_matrices_h5.h5"
  [GSM4186970]="GSM4186970_HTAPP-StJude-SMP-PDX1_nuclei_channel1_raw_gene_bc_matrices_h5.h5"
)

declare -A sample_metadata=(
  [GSM4186961]="GSM4186961_metadata_HTAPP-312-SMP-901_fresh-T1_channel1.csv.gz"
  [GSM4186962]="GSM4186962_metadata_HTAPP-312-SMP-902_fresh-C4-T2_channel1.csv.gz"
  [GSM4186963]="GSM4186963_metadata_HTAPP-656-SMP-3481_fresh-T1_channel1.csv.gz"
  [GSM4186964]="GSM4186964_metadata_HTAPP-StJude-SMP-PDX1_cell_channel1.csv.gz"
  [GSM4186965]="GSM4186965_metadata_HTAPP-244-SMP-451_CST_channel1.csv.gz"
  [GSM4186966]="GSM4186966_metadata_HTAPP-244-SMP-451_EZ_channel1.csv.gz"
  [GSM4186967]="GSM4186967_metadata_HTAPP-244-SMP-451_NST_channel1.csv.gz"
  [GSM4186968]="GSM4186968_metadata_HTAPP-244-SMP-451_TST_channel1.csv.gz"
  [GSM4186969]="GSM4186969_metadata_HTAPP-656-SMP-3481_TST_channel1.csv.gz"
  [GSM4186970]="GSM4186970_metadata_HTAPP-StJude-SMP-PDX1_nuclei_channel1.csv.gz"
)

for sample_name in "${sample_names[@]}"; do
  
  if [[ -e "${sample_data_dir}"/"${sample_matrix[$sample_name]}" ]]; then
  
    bash run-provided-goi-analysis.sh --output_dir "../data/anderson-single-cell/results" \
    --sample_name "${sample_name}" \
    --sample_matrix "${sample_data_dir}"/"${sample_matrix[$sample_name]}" \
    --sample_metadata "${sample_data_dir}"/"${sample_metadata[$sample_name]}" \
    --goi_list "../data/anderson-single-cell/goi-lists/nb_goi_list.tsv" \
    --input_file_type "h5" \
    --gene_set_column_name "gene_set"
  else
    echo "${sample_data_dir}"/"${sample_matrix[$sample_name]}" "does not exist."
  fi
    
done
