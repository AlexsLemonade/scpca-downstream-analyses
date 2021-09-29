## Read in and normalize filtered SingleCellExperiment data

# Command line usage:

# Rscript --vanilla 02-normalize-sce.R \
#   --sample_sce_filepath "data/anderson-single-cell/filtered/filtered_GSM4186961_miQC_sce.rds" \
#   --sample_name "GSM4186961" \
#   --output_data_directory data/anderson-single-cell/normalized \
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
    c("-a", "--filtered_sample_sce_filepath"),
    type = "character",
    default = NULL,
    help = "path to filtered sample SCE RDS file",
  ),
  optparse::make_option(
    c("-r", "--sample_name"),
    type = "character",
    default = NULL,
    help = "name of sample being analyzed",
  ),
  optparse::make_option(
    c("-o", "--output_data_directory"),
    type = "character",
    default = NULL,
    help = "output normalized data directory"
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

## Define file paths

# Directory and file to save output
normalized_dir <- file.path(opt$output_data_directory)
if (!dir.exists(normalized_dir)) {
  dir.create(normalized_dir, recursive = TRUE)
}

output_sce_file <- file.path(
  normalized_dir, paste0("normalized_", opt$sample_name, "_sce.rds"))

#### Read in data --------------------------------------------------------------

# Read in filtered sce object
filtered_sce <- readr::read_rds(opt$filtered_sample_sce_filepath)

#### Normalize the data --------------------------------------------------------

# Cluster similar cells
qclust <- scran::quickCluster(filtered_sce)

# Compute sum factors for each cell cluster grouping
filtered_sce <- scran::computeSumFactors(filtered_sce, clusters = qclust)

# Normalize and log transform
normalized_sce <- scater::logNormCounts(filtered_sce)

#### Add PCA and UMAP results --------------------------------------------------

# Model gene variance using `scran:modelGeneVar()`
gene_variance <- scran::modelGeneVar(normalized_sce)

# Select the most variable genes
subset_genes <- scran::getTopHVGs(gene_variance, n = 2000)

# Calculate a PCA matrix using those genes
normalized_sce <-
  runPCA(normalized_sce, subset_row = subset_genes)

# Calculate a UMAP matrix using the PCA results
normalized_sce <- runUMAP(normalized_sce,
                          dimred = "PCA")

#### Save output normalized file with dimensionality results -------------------

readr::write_rds(normalized_sce, output_sce_file)
