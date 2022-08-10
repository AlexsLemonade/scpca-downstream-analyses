## Read in SingleCellExperiment RDS object that has been normalized and has both PCA and UMAP embeddings.
## This script performs graph based clustering, by default Louvain,
## and outputs a SCE with the cluster assignments stored in the colData.

# Command line usage:

# Rscript --vanilla 04-clustering.R \
#   --sce "results/Gawad_processed_data/SCPCS000245/SCPCL000343_miQC_downstream_processed_normalized_reduced_sce.rds" \
#   --seed 2021 \
#   --cluster_type "louvain" \
#   --nearest_neighbors 10

## Set up -------------------------------------------------------------

# Check R and Bioconductor versions
check_r_bioc_versions()

library(optparse)

## Command line arguments/options
option_list <- list(
  optparse::make_option(
    c("-i", "--sce"),
    type = "character",
    default = NULL,
    help = "path to normalized SCE",
  ),
  optparse::make_option(
    c("-s", "--seed"),
    type = "integer",
    default = 2021,
    help = "seed integer",
    metavar = "integer"
  ),
  optparse::make_option(
    c("-c", "--cluster_type"),
    type = "character",
    default = "louvain",
    help = "Method used for clustering. Can be either louvain or walktrap.",
  ),
  optparse::make_option(
    c("-n", "--nearest_neighbors"),
    type = "integer",
    default = 10,
    help = "Number of nearest neighbors to include during graph construction."
  ),
  optparse::make_option(
    c("-o", "--output_filepath"),
    type = "character",
    default = NULL,
    help = "path to output RDS file containing SCE with cluster assignments added"
  ),
  optparse::make_option(
    c("--project_root"),
    type = "character",
    help = "the path to the root directory for the R project and where the `utils` folder lives."
  )
)

# Read the arguments passed
opt <- parse_args(OptionParser(option_list = option_list))

# if project root is not provided use here::here()
if(is.null(opt$project_root)){
  project_root <- here::here()
} else {
  project_root <- opt$project_root
}

# Source in set up function
source(file.path(project_root, "utils", "setup-functions.R"))

# source in clustering functions
source(file.path(project_root, "utils", "clustering-functions.R"))

# Load project
setup_renv(project_filepath = project_root)

## Load libraries
suppressPackageStartupMessages({
  library(magrittr)
  library(bluster)
  library(SingleCellExperiment)
})

## Set the seed
set.seed(opt$seed)

# Check that the input file exists
if (!file.exists(opt$sce)){
  stop(paste(opt$sce, "does not exist."))
}

# Check that clustering type is valid
if (!opt$cluster_type %in% c("louvain", "walktrap")) {
  stop("--cluster_type (-c) must be either louvain or walktrap.")
}

# Check that `nearest_neighbors` is an integer
if (opt$nearest_neighbors %% 1 != 0){
  stop("The --nearest_neighbors (-n) argument value must be an integer.")
}

# make sure that output file is provided and ends in rds
if (!is.null(opt$output_filepath)){
  if(!(stringr::str_ends(opt$output_filepath, ".rds"))){
    stop("output file name must end in .rds")
  }
} else {
  stop("--output_filepath (-o) must be provided.")
}

#### Read in data and check formatting -----------------------------------------

# Read in normalized sce object
sce <- readr::read_rds(opt$sce)

# check that data contains an SCE with PCA results
if(is(sce,"SingleCellExperiment")){
  if(!"PCA" %in% reducedDimNames(sce)) {
    stop("PCA results are not found in the provided SingleCellExperiment object.")
  }
} else {
  stop("file specified by --sce (-s) must contain a SingleCellExperiment object.")
}

#### Perform clustering --------------------------------------------------------

# perform clustering
sce <- graph_clustering(normalized_sce = sce,
                        params_range = opt$nearest_neighbors,
                        cluster_type = opt$cluster_type)

# write output file
readr::write_rds(sce, opt$output_filepath)
