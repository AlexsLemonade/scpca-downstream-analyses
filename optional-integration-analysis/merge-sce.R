## This script takes in individual SCE objects and merges them into one SCE object
## for a data integration analysis.

## Set up ----------------------------------------------------------------------

# Load libraries
library(optparse)
library(SingleCellExperiment)

# Declare command line options
option_list <- list(
  make_option(
    opt_str = c("--input_metadata_tsv"),
    type = "character",
    help = "Path to the TSV file containing the list of library IDs corresponding to the libraries being integrated, must have a `library_id` column."
  ),
  make_option(
    opt_str = c("--input_directory"),
    type = "character",
    help = "File path to the directory storing the input SCE object RDS files."
  ),
  make_option(
    opt_str = c("-o", "--output_sce_file"),
    type = "character",
    help = "Path to output RDS file, must end in .rds"
  ),
  optparse::make_option(
    c("--project_root"),
    type = "character",
    default = NULL,
    help = "Path to the root directory for the R project and where the `utils` folder lives."
  ),
  make_option(
    opt_str = c("-t", "--threads"),
    type = "integer",
    default = 1,
    help = "Number of multiprocessing threads to use"
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

# Check that input metadata file is provided
if(is.null(opt$input_metadata_tsv)){
  stop("TSV file with the list of library IDs associated with the individual SCE objects to merge is missing.")
} else {
  # List of library ids
  input_metadata <- readr::read_tsv(opt$input_metadata_tsv)
  library_ids <- input_metadata$library_id
}

# Check that there is more than one library id associated with an SCE to be merged
if(length(library_ids) == 1){
  stop("Only 1 input file provided, no merging or integration will be performed for this group")
}

# List of SCE filepaths
sce_files <- library_ids |> purrr::map(~ file.path(opt$input_directory, paste0(.x, "_processed.rds")))

# Check that input files exist
missing_sce_files <- c()
for (sce in sce_files) {
  if(!file.exists(sce)){
    missing_sce_files <- missing_sce_files + sce
  }
}
if(length(missing_sce_files) > 0){
  stop(
    glue::glue(
      "\nCannot find input file: {missing_sce_files}."
    )
  )
}

# Set up multiprocessing params
if(opt$threads > 1){
  bp_param = BiocParallel::MulticoreParam(opt$threads)
} else {
  bp_param = BiocParallel::SerialParam()
}

# Merge SCEs -------------------------------------------------------------------

# Read in list of SCEs
sce_list <- purrr::map(sce_files, readr::read_rds)

# Set names of SCE list
names(sce_list) <- library_ids

# Check that all input RDS files contain SCE objects
sce_checks <- purrr::map(sce_list,
                         \(x) is(x, "SingleCellExperiment"))
if(!all(sce_checks)){
  stop(
    "All input files must contain a `SingleCellExperiment` object."
  )
}

# Create combined SCE object
merged_sce <- scpcaTools::merge_sce_list(sce_list,
                                         batch_column = "library_id",
                                         preserve_rowdata_cols = "gene_symbol",
                                         cell_id_column = "cell_id")

# Save combined SCE object
readr::write_rds(merged_sce, opt$output_sce_file)
