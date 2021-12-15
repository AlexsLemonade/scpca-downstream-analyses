## Read in expression data and prepare SingleCellExperiment object for downstream
## analysis

# Command line usage:
# Rscript --vanilla 00-prepare-sce.R \
# --sample_matrix_filepath GSM4186961_raw_gene_bc_matrices_h5.h5 \
# --sample_metadata_filepath GSM4186961_metadata_HTAPP-312-SMP-901_fresh-T1_channel1.csv.gz \
# --input_file_type "h5" \
# --output_filepath data/pre-filtered-sce/GSM4186961_pre-filtered_sce.RDS

#### Set up --------------------------------------------------------------------

## Load project

# The `here::here()` function return the file path to the main top-level
# directory, hence we are able to provide this to `renv::load()` to find
# the project file
renv::load(here::here())

# Check that R version us at least 4.1
if (! (R.version$major == 4 && R.version$minor >= 1)){
  stop("R version must be at least 4.1")
}

# Check that Bioconductor version is 3.14
if (packageVersion("BiocVersion") < 3.14){
  stop("Bioconductor version is less than 3.14")
}

## Load libraries
library(magrittr)
library(scpcaTools)
library(optparse)
library(SingleCellExperiment)

#### Command line arguments/options --------------------------------------------

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-d", "--sample_matrix_filepath"),
    type = "character",
    default = NULL,
    help = "path to sample matrix file(s)",
  ),
  optparse::make_option(
    c("-m", "--sample_metadata_filepath"),
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
    c("-f", "--output_filename"),
    type = "character",
    default = NULL,
    help = "filename for the filtered SCE RDS object",
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Directory to save output
output_dir <- file.path(opt$output_dir)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
#### Read in data --------------------------------------------------------------

if(opt$input_file_type == "cellranger") {
  sce <- read_cellranger(opt$sample_matrix_filepath)
} else if (opt$input_file_type == "h5") {
  sce <-
    DropletUtils::read10xCounts(opt$sample_matrix_filepath)
  metadata <- data.table::fread(opt$sample_metadata_filepath)
  
  if ("emptydrop" %in% colnames(metadata)) {
    # filter out any instances where emptydrop or doublet is `TRUE`
    metadata_filtered <- metadata %>%
      # create a `Barcode` column to filter the sce on
      dplyr::mutate(Barcode = paste0(gsub("^.*-", "", V1), "-1")) %>%
      dplyr::filter(emptydrop == FALSE & doublet == FALSE)
  } else if (!("emptydrop") %in% colnames(metadata)) {
    # filter out any instances where emptydrop or doublet is `TRUE`
    metadata_filtered <- metadata %>%
      # create a `Barcode` column to filter the sce on
      dplyr::mutate(Barcode = paste0(gsub("^.*-", "", V1), "-1")) %>%
      dplyr::filter(doublet == FALSE)
  }
  # filter using cell barcodes
  sce <- sce[, colData(sce)$Barcode %in% metadata_filtered$Barcode]
  
}

#### Save output SCE -----------------------------------------------------------

readr::write_rds(sce, file.path(output_dir, opt$output_filename))
