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
