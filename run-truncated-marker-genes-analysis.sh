#!/bin/bash
set -euo pipefail
 
###############################################################################
# Run the `01-filter-sce`, `02-normalize-sce`, and `03-dimension-reduction` 
# marker genes analysis scripts on the already pre-processed SingleCellExperiment
# RDS files associated with the Gawad project analysis.

# Usage
# Run truncated marker gene analysis for RDS files using miQC filtering
# run-truncated-marker-genes-analysis.sh --output_dir "data/results" \
# --sample_name "GSM4186961" \
# --sample_matrix "path/to/sample-file-matrix.h5" \
# --sample_metadata "path/to/sample-file-metadata.csv" \
# --marker_genes "data/marker-genes/nb_marker_genes.tsv" \
# --input_file_type "h5" \
# --gene_set_column_name "gene_set" \
# --plotting_identifier_type "symbol"

###############################################################################

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

mito_file=${mito_file:-data/Homo_sapiens.GRCh38.103.mitogenes.txt}
filtering_method=${filtering_method:-miQC}
SEED=${SEED:-2021}
TOP_N=${TOP_N:-2000}
output_dir="data/Gawad_processed_data/results"
sample_name="SCPCS000216"
library_name="SCPCL000290"

# grab variables from command line
while [ $# -gt 0 ]; do
    if [[ $1 == *'--'* ]]; then
        v="${1/--/}"
        declare $v="$2"
    fi
    shift
done

# Run the filtering script on pre-filtered SingleCellExperiment object (the
# implementation below incorporates `miQC` filtering)
Rscript --vanilla 01-filter-sce.R \
  --sample_sce_filepath "data/Gawad_processed_data/${sample_name}/${library_name}_filtered.rds" \
  --sample_name ${library_name} \
  --mito_file ${mito_file} \
  --output_plots_directory "${output_dir}/${sample_name}/plots" \
  --output_filepath "${output_dir}/${sample_name}/${library_name}_filtered_${filtering_method}_sce.rds" \
  --seed ${SEED} \
  --gene_detected_row_cutoff 5 \
  --gene_means_cutoff 0.1 \
  --prob_compromised_cutoff 0.75 \
  --filtering_method ${filtering_method}
  
# Run the normalization script on filtered SingleCellExperiment object
Rscript --vanilla 02-normalize-sce.R \
  --sce "${output_dir}/${sample_name}/${library_name}_filtered_${filtering_method}_sce.rds" \
  --output_filepath "${output_dir}/${sample_name}/${library_name}_normalized_sce.rds" \
  --seed ${SEED}
  
# Run the dimension reduction script on the normalized SingleCellExperiment
# object
Rscript --vanilla 03-dimension-reduction.R \
  --sce "${output_dir}/${sample_name}/${library_name}_normalized_sce.rds" \
  --seed ${SEED} \
  --top_n ${TOP_N} \
  --overwrite
              