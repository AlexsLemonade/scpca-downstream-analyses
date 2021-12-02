#!/bin/bash
set -euo pipefail
 
###############################################################################
# Run the `01-filter-sce`, `02-normalize-sce`, and `03-dimension-reduction` 
# marker genes analysis scripts on the already pre-processed SingleCellExperiment
# RDS files associated with the Gawad project analysis.

#### Note that R 4.1 is required, along with bioconductor 3.14

# Usage
# Run truncated marker gene analysis for RDS files using miQC filtering
# run-truncated-marker-genes-analysis.sh --output_dir "data/results" \
# --sample_metadata "path/to/input-sample-file-metadata.tsv" \
# --mito_file "path/to/mitogenes-file.txt"

# Where `sample_metadata` should be a TSV file with the sample_id, library_id,
# and desired filtering_method for each sample whose RDS file should be read in.
###############################################################################

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

mito_file=${mito_file:-data/Homo_sapiens.GRCh38.103.mitogenes.txt}
SEED=${SEED:-2021}
TOP_N=${TOP_N:-2000}
output_dir="data/Gawad_processed_data/results"
sample_metadata="gawad-library-metadata.tsv"

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
while read -r sample_id library_id filtering_method; do 
  Rscript --vanilla 01-filter-sce.R \
    --sample_sce_filepath "data/Gawad_processed_data/${sample_id}/${library_id}_filtered.rds" \
    --sample_name ${library_id} \
    --mito_file ${mito_file} \
    --output_plots_directory "${output_dir}/${sample_id}/plots" \
    --output_filepath "${output_dir}/${sample_id}/${library_id}_filtered_${filtering_method}_sce.rds" \
    --seed ${SEED} \
    --gene_detected_row_cutoff 5 \
    --gene_means_cutoff 0.1 \
    --prob_compromised_cutoff 0.75 \
    --filtering_method ${filtering_method}
  
  # Run the normalization script on filtered SingleCellExperiment object
  Rscript --vanilla 02-normalize-sce.R \
    --sce "${output_dir}/${sample_id}/${library_id}_filtered_${filtering_method}_sce.rds" \
    --output_filepath "${output_dir}/${sample_id}/${library_id}_normalized_sce.rds" \
    --seed ${SEED}
  
  # Run the dimension reduction script on the normalized SingleCellExperiment
  # object
  Rscript --vanilla 03-dimension-reduction.R \
    --sce "${output_dir}/${sample_id}/${library_id}_normalized_sce.rds" \
    --seed ${SEED} \
    --top_n ${TOP_N} \
    --overwrite
done < "$sample_metadata"
