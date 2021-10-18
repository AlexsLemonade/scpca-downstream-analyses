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

## Load libraries
library(scater)
library(scran)
library(magrittr)
library(optparse)

#### Command line arguments/options --------------------------------------------

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
    type = "character",
    default = "yes",
    help = "specify whether or not to overwrite any existing dimension reduction
            results; the default is yes, overwrite the results",
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Set the seed
set.seed(opt$seed)

#### Read in data --------------------------------------------------------------

# Check that the file exists
if (!file.exists(opt$sce)){
  stop(paste(opt$sce, "does not exist."))
}

# Read in normalized sce object
normalized_sce <- readr::read_rds(opt$sce)

#### Add PCA and UMAP results --------------------------------------------------

# Check that `top_n` is an integer
if (!(opt$top_n %% 1 == 0)){
  stop(paste(opt$top_n, "must be an integer."))
}

# model gene variance using `scran:modelGeneVar()`
gene_variance <- scran::modelGeneVar(normalized_sce)

# select the most variable genes
subset_genes <- scran::getTopHVGs(gene_variance, n = opt$top_n)

# if there are dimension reductions results in the sce and the user has not
# specified whether or not to overwrite the results, stop script with an error
if (!is.null(reducedDims(normalized_sce)) & is.null(opt$overwrite)) {
  stop(
    "Do you want to overwrite the existing dimension reduction
               results? Specify 'yes' or 'no' using the --overwrite flag."
  )
} else if (opt$overwrite == "yes" |
           is.null(reducedDims(normalized_sce))) {
  # if `opt$overwrite == yes` or if there are no dimension reduction results in
  # the sce object, add PCA and UMAP results
  message("Overwriting dimension reduction PCA and UMAP results.")
  
  # calculate a PCA matrix using those genes
  normalized_sce <-
    runPCA(normalized_sce, subset_row = subset_genes)
  
  # calculate a UMAP matrix using the PCA results
  normalized_sce <- runUMAP(normalized_sce, dimred = "PCA")
  
} else if (opt$overwrite == "no") {
  message("Skipping dimension reduction steps.")
}

#### Save normalized file with dimensionality results --------------------------

readr::write_rds(normalized_sce, opt$sce)
