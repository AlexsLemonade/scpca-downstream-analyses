## This script performs integration on processed SCE objects that have been
## merged in the `merge-sce.R` script

## Set up ----------------------------------------------------------------------

# Load libraries
library(optparse)
library(dplyr)
library(SingleCellExperiment)
library(scpcaTools)

# Declare command line options
option_list <- list(
  make_option(
    opt_str = c("--merged_sce_rds"),
    type = "character",
    help = "Path to the RDS file containing the merged SCE object"
  ),
  make_option(
    opt_str = c("--integration_method"),
    type = "character",
    default = c("fastMNN", "harmony"),
    help = "The integration method(s) to use when performing integration; 
    default is c('fastMNN', 'harmony')"
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

# Check R version
check_r_version()

# Set up renv
setup_renv(project_filepath = project_root)

# Check that file extension for output file is correct
if(!(stringr::str_ends(opt$output_sce_file, ".rds"))){
  stop("output file name must end in .rds")
}

# Read in merged SCE -----------------------------------------------------------

merged_sce <- readr::read_rds(opt$merged_sce_rds)

# Integrate SCEs ---------------------------------------------------------------

# Perform integration with specified method
if ("fastMNN" %in% opt$integration_method) {
  integrated_sce <- integrate_sces(merged_sce,
                                   integration_method = "fastMNN",
                                   batch_column = "library_id")
}

if ("harmony" %in% opt$integration_method) {
  integrated_sce <- integrate_sces(merged_sce,
                                   integration_method = "harmony",
                                   batch_column = "library_id")
}

# Write integrated object to file ----------------------------------------------

write_rds(integrated_sce, opt$output_sce_file)

