## Read in SingleCellExperiment RDS object that has been normalized using
## functions from the `scran` and `scater` packages. This script will add PCA
## and UMAP dimensionality results to the normalized SCE object, also using
## functions from the `scran` and `scater` packages.

# Command line usage:

# Rscript --vanilla 03-dimension-reduction.R \
#   --sce "data/anderson-single-cell/normalized/normalized_GSM4186961_sce.rds" \
#   --seed 2021 \
#   --top_n 2000 \
#   --overwrite "yes"

## Set up -------------------------------------------------------------

## Command line arguments/options

# Library needed to declare command line options
library(optparse)

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-i", "--sce"),
    type = "character",
    default = NULL,
    help = "path to normalized sample SCE RDS file",
  ),
  optparse::make_option(
    c("-n", "--top_n"),
    type = "double",
    default = 2000,
    help = "top number of high variance genes to use for dimension reduction;
            the default is top_n = 2000",
  ),
  optparse::make_option(
    c("-s", "--seed"),
    type = "integer",
    default = 2021,
    help = "seed integer",
    metavar = "integer"
  ),
  optparse::make_option(
    c("-o", "--overwrite"),
    action = "store_true",
    help = "specifies whether or not to overwrite any existing dimension reduction
            results"
  ),
  optparse::make_option(
    c("--output_filepath"),
    type = "character",
    default = NULL,
    help = "path to output dimensionality reduced RDS file"
  ),
  optparse::make_option(
    c("--project_root"),
    type = "character",
    help = "the path to the root directory for the R project and where the `utils` folder lives."
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# if project root is not provided use here::here()
if(is.null(opt$project_root)){
  project_root <- here::here()
} else {
  project_root <- opt$project_root
}

# Check that the input file exists
if (!file.exists(opt$sce)){
  stop(paste(opt$sce, "does not exist."))
}

# Check that `top_n` is an integer
if (opt$top_n %% 1 != 0){
  stop("The --top_n (-n) argument value must be an integer.")
}

# Source in set up function
source(file.path(project_root, "utils", "setup-functions.R"))

# Check R version
check_r_version()

# Set up renv
setup_renv(project_filepath = project_root)

## Load libraries
suppressPackageStartupMessages({
  library(scater)
  library(scran)
  library(magrittr)
})

## Set the seed
set.seed(opt$seed)

#### Read in data --------------------------------------------------------------

# Read in normalized sce object
normalized_sce <- readr::read_rds(opt$sce)

#### Add PCA and UMAP results --------------------------------------------------

# model gene variance using `scran:modelGeneVar()`
gene_variance <- scran::modelGeneVar(normalized_sce)

# select the most variable genes
subset_genes <- scran::getTopHVGs(gene_variance, n = opt$top_n)

# save the most variable genes to the metatadata
metadata(normalized_sce)$highly_variable_genes <- subset_genes

# make function to add dimensionality reduction to sce
dim_reduction <- function(normalized_sce, subset_genes) {

  # add PCA to normalized sce
  normalized_sce <- runPCA(normalized_sce, subset_row = subset_genes)

  # calculate a UMAP matrix using the PCA results
  normalized_sce <- runUMAP(normalized_sce, dimred = "PCA")
}

# if there are dimension reductions results in the sce and the user has not
# specified whether or not to overwrite the results, stop script with an error
if (!is.null(reducedDims(normalized_sce))) {
  if (!is.null(opt$overwrite)) {
    # perform dimensionality reduction
    message("Overwriting dimension reduction PCA and UMAP results.")
    normalized_sce <- dim_reduction(normalized_sce, subset_genes)
  } else {
    stop(
      "Dimensionality reduction results exist. Skipping dimensionality reduction
      steps. If you want to overwrite the existing results, use the --overwrite
      flag."
    )
  }
} else {
  # perform dimensionality reduction
  normalized_sce <- dim_reduction(normalized_sce, subset_genes)
}

#### Save normalized file with dimensionality results --------------------------

readr::write_rds(normalized_sce, opt$output_filepath)
