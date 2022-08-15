## Read in SingleCellExperiment RDS object that has been normalized and has PCA
## embeddings.
## This script performs graph based clustering across a range of specified 
## nearest neighbors values and outputs a SCE with the cluster assignments 
## stored in the colData, as well as a separate data frame with cluster validity
## and stability statistics.

# Command line usage:

# Rscript --vanilla clustering-calculations.R \
#   --sce "example-results/sample01/library01_miQC_processed_sce.rds" \
#   --seed 2021 \
#   --cluster_type "louvain" \
#   --nearest_neighbors_range 5:25 \
#   --nearest_neighbors_increment 5 \
#   --output_sce_filepath "example-results/sample01/library01_miQC_processed_sce.rds" \
#   --output_stats_filepath "example-results/sample01/library01_cluster_stats.tsv"

## Set up -------------------------------------------------------------

# Check R and Bioconductor versions
check_r_bioc_versions()

## Command line arguments/options

# Library needed to declare command line options
library(optparse)

# Declare command line options
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
    c("-n", "--nearest_neighbors_range"),
    type = "integer",
    default = 5:25,
    help = "Range with number of nearest neighbors to include during graph construction."
  ),
  optparse::make_option(
    c("-i", "--nearest_neighbors_increment"),
    type = "integer",
    default = 5,
    help = "Increment to use when implementing the range number of nearest 
    neighbors for calculating the cluster stats."
  ),
  optparse::make_option(
    c("-o", "--output_sce_filepath"),
    type = "character",
    default = NULL,
    help = "path to output RDS file containing SCE with cluster assignments added"
  ),
  optparse::make_option(
    c("--output_stats_filepath"),
    type = "character",
    default = NULL,
    help = "path to output RDS file containing a data frame with cluster stats"
  ),
  optparse::make_option(
    c("--project_root"),
    type = "character",
    default = NULL,
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

# Make sure the output sce file is provided and ends in rds 
if (!is.null(opt$output_sce_filepath)) {
  if (!(stringr::str_ends(
    opt$output_sce_filepath,
    stringr::regex("\\.rds", ignore_case = TRUE)
  ))) {
    stop("output file name must end in .rds")
  }
} else {
  stop("--output_sce_filepath (-o) must be provided.")
}

# Make sure the output directory exists 
output_sce_dir <- dirname(opt$output_sce_filepath)
if(!dir.exists(output_sce_dir)){
  dir.create(output_sce_dir, recursive = TRUE)
}

#### Read in data and check formatting -----------------------------------------

# Read in normalized sce object
sce <- readr::read_rds(opt$sce)

# Check that the input data contains an SCE with PCA results 
if(is(sce,"SingleCellExperiment")){
  if(!"PCA" %in% reducedDimNames(sce)) {
    stop("PCA results are not found in the provided SingleCellExperiment object.")
  }
} else {
  stop("file specified by --sce (-s) must contain a SingleCellExperiment object.")
}

#### Perform clustering --------------------------------------------------------

# Perform graph-based clustering
sce <- graph_clustering(normalized_sce = sce,
                        params_range = opt$nearest_neighbors_range,
                        step_size = opt$nearest_neighbors_increment,
                        cluster_type = opt$cluster_type)

# Write output SCE file
readr::write_rds(sce, opt$output_sce_filepath)

### Calculate cluster validity stats -------------------------------------------

# Check the cluster validity stats for each of the clusters in the SCE object
# and return stats in a data frame
validity_stats_df <- create_metadata_stats_df(sce, opt$nearest_neighbors_range, 
                                              opt$nearest_neighbors_increment, opt$cluster_type)

# Summarize the stats and return in a data frame
summary_validity_stats_df <- summarize_clustering_stats(validity_stats_df) %>%
  dplyr::mutate(param_value = as.numeric(param_value))

### Calculate cluster stability stats ------------------------------------------

# Check the cluster stability stats and return the summary ari in a data frame
summary_stability_stats_df <- get_cluster_stability_summary(sce, opt$nearest_neighbors_range, 
                                                            opt$nearest_neighbors_increment, opt$cluster_type) %>%
  dplyr::mutate(param_value = as.numeric(param_value))

### Combine and save stats data frame ------------------------------------------

# Combined the summary stats data frames
combined_stats_df <- summary_validity_stats_df %>%
  dplyr::inner_join(summary_stability_stats_df, 
                    by = c("cluster_names_column", "param_value", "cluster_type"))

# Write output file
readr::write_tsv(combined_stats_df, opt$output_stats_filepath)
