## This script performs graph based clustering across a range of specified
## nearest neighbors values and outputs a SCE with the cluster assignments
## stored in the colData, as well as a separate data frame with cluster validity
## and stability statistics.
## Read in SingleCellExperiment RDS object that has been normalized and has PCA
## embeddings.

# Command line usage:

# Rscript --vanilla clustering-calculations.R \
#   --sce "example-results/sample01/library01_miQC_processed_sce.rds" \
#   --library_id "library01" \
#   --seed 2021 \
#   --cluster_types "louvain,walktrap" \
#   --nearest_neighbors_min 5
#   --nearest_neighbors_max 25 \
#   --nearest_neighbors_increment 5 \
#   --output_directory "example-results/sample01" \
#   --output_sce "example-results/sample01/library01_clustered_sce.rds"

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
    c("-c", "--cluster_types"),
    type = "character",
    default = "louvain,walktrap",
    help = "Method used for clustering. Can be either louvain or walktrap.",
  ),
  optparse::make_option(
    c("-n", "--nearest_neighbors_min"),
    type = "integer",
    default = 5,
    help = "The minimum number of a range of nearest neighbors values to include
    during graph construction. Default is 5."
  ),
  optparse::make_option(
    c("--nearest_neighbors_max"),
    type = "integer",
    default = 25,
    help = "The maximum number of a range of nearest neighbors values to include
    during graph construction. Default is 25."
  ),
  optparse::make_option(
    c("--nearest_neighbors_increment"),
    type = "integer",
    default = 1,
    help = "Increment to use when implementing the range number of nearest
    neighbors for calculating the cluster stats. If using a single value for the
    nearest neighbors range, set to 1"
  ),
  optparse::make_option(
    c("-o", "--output_directory"),
    type = "character",
    default = NULL,
    help = "path to output directory that would hold the TSV files with the
    clustering stats"
  ),
  optparse::make_option(
    c("--output_sce"),
    type = "character",
    default = NULL,
    help = "path to the RDS file containing SCE with cluster assignments added"
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

# Split up string of cluster types
cluster_types <-  stringr::str_split(opt$cluster_types, ",") %>%
  unlist() %>%
  stringr::str_trim()

# Check that clustering type is valid
for (cluster_type in cluster_types) {
  if (!cluster_type %in% c("louvain", "walktrap")) {
    stop("--cluster_types (-c) must be either louvain, walktrap, or both (comma-separated)")
  }
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
  stop("file specified by --sce (-i) must contain a SingleCellExperiment object.")
}

#### Perform clustering --------------------------------------------------------

# Define nearest neighbors range of values
nn_range <- define_nn_range(opt$nearest_neighbors_min,
                            opt$nearest_neighbors_max,
                            opt$nearest_neighbors_increment)

perform_clustering <- function(sce, nn_range, cluster_type, ...) {
  # Check for existing clustering results
  cluster_column_names <- sprintf("%s_%02d", cluster_type, nn_range)
  existing_columns <-
    intersect(cluster_column_names, colnames(colData(sce)))

  if (length(existing_columns) != 0) {
    if (!is.null(opt$overwrite)) {
      # Perform graph-based clustering
      message("Overwriting clustering results.")
      sce <- graph_clustering(
        normalized_sce = sce,
        nearest_neighbors_min = as.integer(opt$nearest_neighbors_min),
        nearest_neighbors_max = as.integer(opt$nearest_neighbors_max),
        step_size = opt$nearest_neighbors_increment,
        cluster_type = cluster_type
      )
    } else {
      message(glue::glue("
      Clustering results exist for {cluster_column_names}."))
      stop(
        "Skipping clustering steps. If you want to overwrite the existing results, use the --overwrite flag."
      )
    }
  } else {
    sce <- graph_clustering(
      normalized_sce = sce,
      nearest_neighbors_min = as.integer(opt$nearest_neighbors_min),
      nearest_neighbors_max = as.integer(opt$nearest_neighbors_max),
      step_size = opt$nearest_neighbors_increment,
      cluster_type = cluster_type
    )
  }

}

# Run the clustering wrapper function
for (cluster_type in cluster_types){
  sce <- perform_clustering(sce, nn_range, cluster_type,
                                  opt$nearest_neighbors_min,
                                  opt$nearest_neighbors_max,
                                  opt$nearest_neighbors_increment)
}


# Write output SCE file
readr::write_rds(sce, opt$output_sce)

### Calculate cluster validity stats -------------------------------------------

# Check the cluster validity stats and return the results in a data frame
validity_stats_df_list <- purrr::map(cluster_types,
                                     ~ create_metadata_stats_df(sce,
                                                                nn_range,
                                                                .x))
validity_stats_df <-
  dplyr::bind_rows(validity_stats_df_list) %>%
  dplyr::mutate(param_value = as.numeric(param_value))

# Write output file with all cluster validity stats
readr::write_tsv(validity_stats_df,
                 file.path(
                   opt$output_directory,
                   paste0(opt$library_id, "_clustering_all_validity_stats.tsv")
                 ))

# Summarize the validity stats and return in a data frame
summary_validity_stats_df <-
  summarize_clustering_stats(validity_stats_df) %>%
  dplyr::mutate(param_value = as.numeric(param_value))

# Write output file with summary cluster validity stats
readr::write_tsv(summary_validity_stats_df,
                 file.path(
                   opt$output_directory,
                   paste0(opt$library_id, "_clustering_summary_validity_stats.tsv")
                 ))

### Calculate cluster stability stats ------------------------------------------

# Check the cluster stability stats and return the summary ari in a data frame
summary_stability_stats_df_list <- purrr::map(cluster_types,
                                              ~ get_cluster_stability_summary(sce,
                                                                              nn_range,
                                                                              .x))
summary_stability_stats_df <-
  dplyr::bind_rows(summary_stability_stats_df_list) %>%
  dplyr::mutate(param_value = as.numeric(param_value))

# Write output file with summary cluster stability stats
readr::write_tsv(summary_stability_stats_df,
                 file.path(
                   opt$output_directory,
                   paste0(opt$library_id, "_clustering_summary_stability_stats.tsv")
                 ))
