#!/bin/bash
set -euo pipefail
 
###############################################################################
# Run the marker genes analysis scripts including the optional marker genes
# mapping script in `utils`, and the prepare SingleCellExperiment object,
# filter, and normalize scripts (which are also optional based on the format of
# the data).

# Usage
# Run marker gene analysis for h5 files using miQC filtering
# run-marker-genes-analysis.sh --output_dir "data/results" \
# --sample_name "GSM4186961" \
# --sample_matrix "path/to/sample-file-matrix.h5" \
# --sample_metadata "path/to/sample-file-metadata.csv" \
# --marker_genes "data/marker-genes/nb_marker_genes.tsv" \
# --input_file_type "h5" \
# --ensembl_id_column_name "ensembl" \
# --gene_set_column_name "gene_set" \
# --plotting_identifier_type "symbol"

# Run marker gene analysis for raw cellranger files using manual filtering
# run-marker-genes-analysis.sh --output_dir "data/results" \
# --sample_name "GSM4186961" \
# --sample_matrix "path/to/sample-file-matrix.h5" \
# --marker_genes "data/marker-genes/nb_marker_genes.tsv" \
# --input_file_type "cellranger" \
# --filtering_method "miQC"
###############################################################################

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

RUN_MAPPING=${RUN_MAPPING:-FALSE}

mito_file=${mito_file:-data/Homo_sapiens.GRCh38.103.mitogenes.txt}
input_identifiers=${input_identifiers:-SYMBOL}
output_identifiers=${output_identifiers:-ENSEMBL}
identifier_column_name=${identifier_column_name:-gene_symbol}
organism=${organism:-Homo sapiens}
filtering_method=${filtering_method:-miQC}
SEED=${SEED:-2021}
TOP_N=${TOP_N:-2000}
plotting_identifier_type=${plotting_identifier_type:-symbol}
ensembl_id_column_name=${ensembl_id_column_name:-ensembl}

# grab variables from command line
while [ $# -gt 0 ]; do
    if [[ $1 == *'--'* ]]; then
        v="${1/--/}"
        declare $v="$2"
    fi
    shift
done

# Run the gene identifier mapping script if `RUN_MAPPING = TRUE`
if [[ ${RUN_MAPPING} == "TRUE" ]]
then
  Rscript --vanilla utils/map-marker-genes.R \
    --input_marker_gene_list ${marker_genes} \
    --input_identifiers ${input_identifiers} \
    --output_identifiers ${output_identifiers} \
    --identifier_column_name ${identifier_column_name} \
    --organism "${organism}" \
    --multi_mappings "list" \
    --output_file "${output_dir}/mapped_marker_genes.tsv"
else
  echo "Skipping gene identifier mapping."
fi

# Run the 00 prepare SingleCellExperiment object script
Rscript --vanilla 00-prepare-sce.R \
 --sample_matrix_filepath ${sample_matrix} \
 --sample_metadata_filepath ${sample_metadata} \
 --input_file_type ${input_file_type} \
 --output_dir "${output_dir}/${sample_name}" \
 --output_filename "${sample_name}_pre-filtered_sce.rds"

# Run the filtering script on pre-filtered SingleCellExperiment object (the
# implementation below incorporates `miQC` filtering)
Rscript --vanilla 01-filter-sce.R \
  --sample_sce_filepath "${output_dir}/${sample_name}/${sample_name}_pre-filtered_sce.rds" \
  --sample_name ${sample_name} \
  --mito_file ${mito_file} \
  --output_data_directory "${output_dir}/${sample_name}" \
  --output_plots_directory "${output_dir}/${sample_name}/plots" \
  --seed ${SEED} \
  --gene_detected_row_cutoff 5 \
  --gene_means_cutoff 0.1 \
  --filtering_method ${filtering_method}

# Run the normalization script on filtered SingleCellExperiment object
Rscript --vanilla 02-normalize-sce.R \
  --sce "${output_dir}/${sample_name}/${sample_name}_filtered_${filtering_method}_sce.rds" \
  --output_filepath "${output_dir}/${sample_name}/${sample_name}_normalized_sce.rds" \
  --seed ${SEED}

# Run the dimension reduction script on the normalized SingleCellExperiment
# object
Rscript --vanilla 03-dimension-reduction.R \
  --sce "${output_dir}/${sample_name}/${sample_name}_normalized_sce.rds" \
  --seed ${SEED} \
  --top_n ${TOP_N} \
  --overwrite

# Generate a html report displaying a hierarchical clustering heatmap and marker
# gene expression plots
if [[ ${plotting_identifier_type} == "symbol" ]]
then
Rscript -e "rmarkdown::render('marker-genes-report-template.Rmd', clean = TRUE,
              output_file = '${output_dir}/${sample_name}/${sample_name}_provided_markers_report.html',
              params = list(sample = '${sample_name}',
                            normalized_sce = '${output_dir}/${sample_name}/${sample_name}_normalized_sce.rds',
                            marker_genes = '${output_dir}/mapped_marker_genes.tsv',
                            ensembl_id_column = '${ensembl_id_column_name}',
                            gene_symbol_column = '${identifier_column_name}',
                            gene_set_column = '${gene_set_column_name}'),
              envir = new.env())"
elif [[ ${plotting_identifier_type} == "ensembl" ]]
then
Rscript -e "rmarkdown::render('marker-genes-report-template.Rmd', clean = TRUE,
              output_file = '${output_dir}/${sample_name}/${sample_name}_provided_markers_report.html',
              params = list(sample = '${sample_name}',
                            normalized_sce = '${output_dir}/${sample_name}/${sample_name}_normalized_sce.rds',
                            marker_genes = '${output_dir}/mapped_marker_genes.tsv',
                            ensembl_id_column = '${ensembl_id_column_name}',
                            gene_set_column = '${gene_set_column_name}'),
              envir = new.env())"
else
  echo "Plotting gene identifier specified is not supported. Specify 'symbol' or
  'ensembl' using the --plotting_identifier_type flag."
fi
