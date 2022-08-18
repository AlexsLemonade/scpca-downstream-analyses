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

## Command line arguments/options

# Library needed to declare command line options
library(optparse)

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-i", "--sce"),
    type = "character",
    default = NULL,
    help = "path to RDS file containing SCE object with log-normalized counts 
    and PCA embeddings",
  ),
  optparse::make_option(
    c("-r", "--sample_id"),
    type = "character",
    default = NULL,
    help = "name of sample being analyzed",
  ),
  optparse::make_option(
    c("-l", "--library_id"),
    type = "character",
    default = NULL,
    help = "name of library being analyzed",
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
    help = "Range with number of nearest neighbors to include during graph 
    construction. Can be a range of values, e.g. 5:25, or a single value."
  ),
  optparse::make_option(
    c("-i", "--nearest_neighbors_increment"),
    type = "integer",
    default = 5,
    help = "Increment to use when implementing the range number of nearest 
    neighbors for calculating the cluster stats. If using a single value for the
    nearest neighbors range, set to NULL"
  ),
  optparse::make_option(
    c("-o", "--output_directory"),
    type = "character",
    default = NULL,
    help = "path to output directory that would hold the RDS file containing SCE
    with cluster assignments added, as well as the TSV files with the clustering stats"
  ),
  optparse::make_option(
    c("--project_root"),
    type = "character",
    default = NULL,
    help = "the path to the root directory for the R project and where the `utils` folder lives."
  ),
  optparse::make_option(
    c( "--overwrite"),
    action = "store_true",
    help = "specifies whether or not to overwrite any existing clustering results"
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

# Check R and Bioconductor versions
check_r_bioc_versions()

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

# Make sure the output directory exists 
if(!dir.exists(opt$output_directory)){
  dir.create(opt$output_directory, recursive = TRUE)
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

# Check for existing clustering results
if (!is.null(opt$overwrite)) {
  # Perform graph-based clustering
  message("Overwriting clustering results.")
  sce <- graph_clustering(
    normalized_sce = sce,
    params_range = opt$nearest_neighbors_range,
    step_size = opt$nearest_neighbors_increment,
    cluster_type = opt$cluster_type
  )
} else {
  stop(
    "Clustering results exist. Skipping clustering steps. If you want to
    overwrite the existing results, use the --overwrite flag."
  )
}

# Write output SCE file
readr::write_rds(sce, file.path(
  opt$output_directory,
  opt$sample_id,
  paste0(opt$library_id, "_processed_sce_with_clustering.rds")
))

### Calculate cluster validity stats -------------------------------------------

# Check the cluster validity stats for each of the clusters in the SCE object
# and return stats in a data frame
validity_stats_df <- create_metadata_stats_df(sce, 
                                              opt$nearest_neighbors_range, 
                                              opt$nearest_neighbors_increment, 
                                              opt$cluster_type)

# Write output file with all cluster validity stats
readr::write_tsv(validity_stats_df,
                 file.path(
                   opt$output_directory,
                   opt$sample_id,
                   paste0(opt$library_id, "_clustering_all_stats.tsv")
                 ))

# Summarize the stats and return in a data frame
summary_validity_stats_df <- summarize_clustering_stats(validity_stats_df) %>%
  dplyr::mutate(param_value = as.numeric(param_value))

# Write output file
readr::write_tsv(summary_validity_stats_df,
                 file.path(
                   opt$output_directory,
                   opt$sample_id,
                   paste0(opt$library_id, "_clustering_summary_validity_stats.tsv")
                 ))

### Calculate cluster stability stats ------------------------------------------

# Check the cluster stability stats and return the summary ari in a data frame
summary_stability_stats_df <-
  get_cluster_stability_summary(
    sce,
    opt$nearest_neighbors_range,
    opt$nearest_neighbors_increment,
    opt$cluster_type
  ) %>%
  dplyr::mutate(param_value = as.numeric(param_value))

# Write output file
readr::write_tsv(combined_stats_df,
                 file.path(
                   opt$output_directory,
                   opt$sample_id,
                   paste0(opt$library_id, "_clustering_summary_stability_stats.tsv")
                 ))
