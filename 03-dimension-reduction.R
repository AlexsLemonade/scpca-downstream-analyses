## Read in SingleCellExperiment RDS object that has been normalized using
## functions from the `scran` and `scater` packages. This script will add PCA
## and UMAP dimensionality results to the normalized SCE object, also using
## functions from the `scran` and `scater` packages.

# Command line usage:

# Rscript --vanilla 03-dimension-reduction.R \
#   --sce "data/anderson-single-cell/normalized/normalized_GSM4186961_sce.rds" \
#   --seed 2021 \

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

# Set the seed
set.seed(opt$seed)

#### Read in data --------------------------------------------------------------

# Read in normalized sce object
normalized_sce <- readr::read_rds(opt$sce)

#### Add PCA and UMAP results --------------------------------------------------

# model gene variance using `scran:modelGeneVar()`
gene_variance <- scran::modelGeneVar(normalized_sce)

# select the most variable genes
subset_genes <- scran::getTopHVGs(gene_variance, n = 2000)

# calculate a PCA matrix using those genes
normalized_sce <-
  runPCA(normalized_sce, subset_row = subset_genes)

# calculate a UMAP matrix using the PCA results
normalized_sce <- runUMAP(normalized_sce,
                          dimred = "PCA")

#### Save normalized file with dimensionality results --------------------------

readr::write_rds(normalized_sce, opt$sce)
