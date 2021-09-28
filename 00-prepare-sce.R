## Read in expression data and prepare SingleCellExperiment object for downstream
## analysis

# Command line usage:
# Rscript --vanilla 00-prepare-sce.R \
# --main_sample_dir data/GSM4186961 \
# --sample_matrix_file GSM4186961_raw_gene_bc_matrices_h5.h5 \
# --sample_metadata_file GSM4186961_metadata_HTAPP-312-SMP-901_fresh-T1_channel1.csv.gz \
# --input_file_type "h5" \
# --output_dir data/pre-filtered-sce \
# --output_filename GSM4186961_pre-filtered_sce.RDS

#### Set up --------------------------------------------------------------------

## Load libraries
library(magrittr)
library(scpcaTools)
library(optparse)

#### Command line arguments/options --------------------------------------------

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-r", "--main_sample_dir"),
    type = "character",
    default = NULL,
    help = "path to main sample directory",
  ),
  optparse::make_option(
    c("-d", "--sample_matrix_file"),
    type = "character",
    default = NULL,
    help = "path to sample matrix file(s)",
  ),
  optparse::make_option(
    c("-m", "--sample_metadata_file"),
    type = "character",
    default = NULL,
    help = "path to sample metadata file, if file type is h5",
  ),
  optparse::make_option(
    c("-y", "--input_file_type"),
    type = "character",
    default = NULL,
    help = "character string specifying the input file type",
  ),
  optparse::make_option(
    c("-o", "--output_dir"),
    type = "character",
    default = NULL,
    help = "output directory",
  ),
  optparse::make_option(
    c("-", "--output_filename"),
    type = "character",
    default = NULL,
    help = "filename for the filtered SCE RDS object",
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

opt$main_sample_dir <- "data/anderson-single-cell/GSE140819"
opt$sample_matrix_file <- "GSM4186961_HTAPP-312-SMP-901_fresh-T1_channel1_raw_gene_bc_matrices_h5.h5"
opt$sample_metadata_file <- "GSM4186961_metadata_HTAPP-312-SMP-901_fresh-T1_channel1.csv.gz"
opt$input_file_type <- "h5"
opt$output_dir <- "data/anderson-single-cell/pre-filtered-sce"
opt$output_filename <- "GSM4186961_pre-filtered_sce.rds"

# Directory to save output
output_dir <- file.path(opt$output_dir)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
#### Read in data --------------------------------------------------------------

matrix_file <- file.path(opt$main_sample_dir, opt$sample_matrix_file)

if(opt$input_file_type == "cellranger") {
  sce <- read_cellranger(matrix_file)
} else if (opt$input_file_type == "h5") {
  sce <- DropletUtils::read10xCounts(matrix_file)
  metadata <-
    data.table::fread(file.path(opt$main_sample_dir, opt$sample_metadata_file))
}

#### Filter based on metadata --------------------------------------------------

# filter out any instances where emptydrop or doublet is `TRUE`
metadata_filtered <- metadata %>%
  # create a `Barcode` column to filter the sce on
  dplyr::mutate(Barcode = paste0(gsub("^.*-", "", V1), "-1")) %>%
  dplyr::filter(emptydrop == FALSE & doublet == FALSE)

# filter using cell barcodes
filtered_sce <- sce[,colData(sce)$Barcode %in% metadata_filtered$Barcode]

#### Save filtered SCE ---------------------------------------------------------

readr::write_rds(filtered_sce, file.path(output_dir, opt$output_filename))
