## Read in SingleCellExperiment RDS object that has been filtered to remove
## cells on the basis of their mitochondrial read percentage and detected genes
## and to remove lowly expressed or undetected genes. This script will then
## normalize the filtered object using functions from the `scran` and `scater` 
## packages.

# Command line usage:

# Rscript --vanilla 02-normalize-sce.R \
#   --sce "data/anderson-single-cell/filtered/filtered_GSM4186961_miQC_sce.rds" \
#   --output_filepath data/anderson-single-cell/normalized/normalized_GSM4186961_sce.rds \
#   --seed 2021 \

## Set up -------------------------------------------------------------

## Load project

# `here::here()` looks at a number of criteria to identify the root 
# directory, including whether or not there is a .Rproj file present,
# so we can pass this to `renv::load()` to load the project file
renv::load(here::here())

# Check that R version us at least 4.1
if (! (R.version$major == 4 && R.version$minor >= 1)){
  stop("R version must be at least 4.1")
}

# Check that Bioconductor version is 3.14
if (packageVersion("BiocVersion") < 3.14){
  stop("Bioconductor version is less than 3.14")
}

## Command line arguments/options

# Library needed to declare command line options
library(optparse)

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-i", "--sce"),
    type = "character",
    default = NULL,
    help = "path to filtered sample SCE RDS file",
  ),
  optparse::make_option(
    c("-o", "--output_filepath"),
    type = "character",
    default = NULL,
    help = "path to output normalized RDS file"
  ),
  optparse::make_option(
    c("-s", "--seed"),
    type = "integer",
    default = 2021,
    help = "seed integer",
    metavar = "integer"
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## Load libraries
library(scater)
library(scran)
library(magrittr)

## Set the seed
set.seed(opt$seed)

## Define file paths

# Directory and file to save output
output_file <- opt$output_filepath
output_dir <- dirname(output_file)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#### Read in data --------------------------------------------------------------

# Read in filtered sce object
filtered_sce <- readr::read_rds(opt$sce)

#### Normalize the data --------------------------------------------------------

tryCatch(
  expr = {
    # Cluster similar cells
    qclust <- scran::quickCluster(filtered_sce)
  },
  error = function(e){ 
    print("Clustering similar cells failed. Skipping this sample.")
  },
  warning = function(w){
    print(paste0("Clustering similar cells failed. Skipping normalization and dimension reduction for filtered sample ", opt$sce))
  }
  
)

if (exists("qclust")) {
  # Compute sum factors for each cell cluster grouping
  filtered_sce <-
    scran::computeSumFactors(filtered_sce, clusters = qclust)
  
  # Include note in metadata re: clustering
  metadata(filtered_sce)$normalization <- "scater::logNormCounts clustered"
  
} else if (!exists("qclust")) {
  
  # Include note in metadata re: failed clustering
  metadata(filtered_sce)$normalization <- "scater::logNormCounts unclustered" 
  
  # Keep positive counts for `logNormCounts()`
  filtered_sce <- filtered_sce[, colSums(counts(filtered_sce)) > 0]
}
  
# Normalize and log transform
normalized_sce <- scater::logNormCounts(filtered_sce)

# Save output normalized file
readr::write_rds(normalized_sce, output_file)
