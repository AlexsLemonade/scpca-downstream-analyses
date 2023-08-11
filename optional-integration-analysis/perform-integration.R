## This script performs integration on processed SCE objects that have been
## merged in the `merge-sce.R` script

## Set up ----------------------------------------------------------------------

if (!("harmony" %in% installed.packages())) {
  remotes::install_version("harmony", version = "0.1.1")
}

# Load libraries
library(optparse)
library(dplyr)
library(SingleCellExperiment)
library(scpcaTools)

# Declare command line options
option_list <- list(
  make_option(
    opt_str = c("--merged_sce_file"),
    type = "character",
    help = "Path to the RDS file containing the merged SCE object"
  ),
  make_option(
    opt_str = c("--integration_method"),
    type = "character",
    default = "fastMNN,harmony",
    help = "The integration method(s) to use when performing integration; 
    default is 'fastMNN,harmony'"
  ),
  make_option(
    opt_str = c("--fastmnn_auto_merge"),
    action = "store_true",
    default = FALSE,
    help = "Indicates whether or not to use the auto.merge option for `fastMNN` integration;
    to perform auto.merge, use `--fastmnn_auto_merge` and this will override any input order"
  ),
  make_option(
    opt_str = c("--fastmnn_merge_order"),
    type = "character",
    default = NULL,
    help = "Optional vector of library ids in the order of which they should be merged;
    will only be used if --fastmnn_auto_merge is FALSE"
  ),
  make_option(
    opt_str = c("-o", "--output_sce_file"),
    type = "character",
    help = "Path to output RDS file containing integrated object, must end in .rds"
  ),
  make_option(
    c("--project_root"),
    type = "character",
    default = NULL,
    help = "Path to the root directory for the R project and where the `utils` folder lives."
  )
)

# Read the arguments passed
opt <- parse_args(OptionParser(option_list = option_list))

# If project root is not provided use here::here()
if(is.null(opt$project_root)){
  project_root <- here::here()
} else {
  project_root <- opt$project_root
}

# Source in set up functions
source(file.path(project_root, "utils", "setup-functions.R"))
source(file.path(project_root, "utils", "integration-functions.R"))

# Check R version
check_r_version()

# Set up renv
setup_renv(project_filepath = project_root)

# Check that file extension for output file is correct
if(!(stringr::str_ends(opt$output_sce_file, ".rds"))){
  stop("output file name must end in .rds")
}

# Check that merged SCE file exists
if(!file.exists(opt$merged_sce_file)){
  stop(paste(opt$merged_sce_file, "does not exist."))
}

# Read in merged SCE -----------------------------------------------------------

merged_sce <- readr::read_rds(opt$merged_sce_file)

# Integrate SCEs ---------------------------------------------------------------

# Split up string of integration methods
integration_methods <-  stringr::str_split(opt$integration_method, ",") |>
  unlist() |>
  stringr::str_trim()

# Check that we support the provided integration method(s)
if (!all(integration_methods %in% c("fastMNN", "harmony"))) {
  stop("--integration_method must be either fastMNN, harmony, or both (comma-separated)")
}

# Check that `library_id` is in merged SCE object
if(!("library_id" %in% colnames(colData(merged_sce)))){
  stop("The 'library_id' column is missing from the `colData` of the merged SCE 
       object and is needed for integration.")
}

# Perform integration with specified method
if ("fastMNN" %in% integration_methods) {
  # Format `fastmnn_merge_order` if provided; can only be used when auto.merge is FALSE
  if(!opt$fastmnn_auto_merge) {
    if(!is.null(opt$fastmnn_merge_order)){
      fastmnn_merge_order <- stringr::str_split(opt$fastmnn_merge_order, ",") |>
        unlist() |>
        stringr::str_trim()
      if(!all(colData(merged_sce)$library_id %in% fastmnn_merge_order)){
        stop("The provided library ids do not match those in the SCE object. 
             Please re-run and provide --fastmnn_merge_order with all library ids 
             found in the `colData` of the SCE object.")
      }
    }
  } else {
    fastmnn_merge_order <- NULL
  }
  fastMNN_integrated_sce <- integrate_sces(merged_sce,
                                   integration_method = "fastMNN",
                                   batch_column = "library_id",
                                   auto.merge = opt$fastmnn_auto_merge,
                                   merge.order = fastmnn_merge_order)
  
  merged_sce <- add_integrated_pcs(merged_sce,
                                   integrated_pcs = reducedDim(fastMNN_integrated_sce, "fastMNN_PCA"),
                                   integration_method = "fastMNN")
  
}

if ("harmony" %in% integration_methods) {
  harmony_integrated_sce <- integrate_sces(merged_sce,
                                   integration_method = "harmony",
                                   batch_column = "library_id")
  
  merged_sce <- add_integrated_pcs(merged_sce,
                                   integrated_pcs = reducedDim(harmony_integrated_sce, "harmony_PCA"),
                                   integration_method = "harmony")
}


# Write integrated object to file ----------------------------------------------

readr::write_rds(merged_sce, opt$output_sce_file)
