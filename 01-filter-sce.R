## Read in and filter SingleCellExperiment data

# Command line usage:

# Rscript --vanilla 01-filter-sce.R \
#   --sample_sce_filepath "data/anderson-single-cell/pre-filtered-sce/GSM4186961_pre-filtered_sce.rds" \
#   --sample_id "GSM4186961" \
#   --library_id "GSL4186960" \
#   --mito_file data/Homo_sapiens.GRCh38.103.mitogenes.txt \
#   --output_filepath data/anderson-single-cell/results/sample_filtered_sce.rds \
#   --seed 2021 \
#   --gene_detected_row_cutoff 5 \
#   --gene_means_cutoff 0.1 \
#   --prob_compromised_cutoff 0.75 \
#   --filtering_method "miQC"

## Set up -------------------------------------------------------------

# Check R and Bioconductor versions
check_r_bioc_versions()

## Command line arguments/options

# Library needed to declare command line options
library(optparse)

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-a", "--sample_sce_filepath"),
    type = "character",
    default = NULL,
    help = "path to sample SCE RDS file",
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
    c("-i", "--mito_file"),
    type = "character",
    default = NULL,
    help = "the path to the mito file"
  ),
  optparse::make_option(
    c("--output_filepath"),
    type = "character",
    default = NULL,
    help = "filepath to output filtered file"
  ),
  optparse::make_option(
    c("-s", "--seed"),
    type = "integer",
    default = 2021,
    help = "seed integer",
    metavar = "integer"
  ),
  optparse::make_option(
    c("-g", "--gene_detected_row_cutoff"),
    type = "integer",
    default = 5,
    help = "cutoff for percent of cells genes must be detected in",
    metavar = "integer"
  ),
  optparse::make_option(
    c("-m", "--gene_means_cutoff"),
    type = "double",
    default = 0.1,
    help = "gene mean expression cutoff",
    metavar = "double"
  ),
  optparse::make_option(
    c("-c", "--mito_percent_cutoff"),
    type = "double",
    default = 20,
    help = "cell mitochondrial percent cutoff -- not needed for miQC filtering",
    metavar = "double"
  ),
  optparse::make_option(
    c("-p", "--detected_gene_cutoff"),
    type = "integer",
    default = 500,
    help = "cell detected genes cutoff -- not needed for miQC filtering",
    metavar = "integer"
  ),
  optparse::make_option(
    c("-u", "--umi_count_cutofff"),
    type = "integer",
    default = 500,
    help = "cell UMI counts cutoff -- not needed for miQC filtering",
    metavar = "integer"
  ),
  optparse::make_option(
    c("--prob_compromised_cutoff"),
    type = "double",
    default = 0.75,
    help = "probability compromised cutoff used for filtering cells if using miQC",
    metavar = "double"
  ),
  optparse::make_option(
    c("-f", "--filtering_method"),
    type = "character",
    default = "miQC",
    help = "the selected filtered method -- can be miQC or manual; will be miQC
            by default"
  ),
  optparse::make_option(
    c("--project_root"),
    type = "character",
    help = "the path to the root directory for the R project and where the `utils` folder lives."
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# if project root is not provided use here::here()
if(is.null(opt$project_root)){
  project_root <- here::here()
} else {
  project_root <- opt$project_root
}

# Source in set up function
source(file.path(project_root, "utils", "setup-functions.R"))

# Load project
setup_renv(project_filepath = project_root)

# Check that input arguments are valid
if (!opt$filtering_method %in% c("manual", "miQC")) {
  stop("Incorrect filtering method. Specify `manual` or `miQC` filtering.")
}

## Load libraries
library(scater)
library(scran)
library(ggplot2)
library(magrittr)
library(scpcaTools)
library(SingleCellExperiment)
library(cowplot)
library(scuttle)
library(tryCatchLog)

# source filtering functions
source(file.path(project_root, "utils", "filtering-functions.R"))

## Set the seed
set.seed(opt$seed)

## Define file paths

# Directory and file to save output
output_file <- opt$output_filepath
output_dir <- dirname(output_file)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#### Filter data ---------------------------------------------------------------

# Read in sce object
sce_qc <- readr::read_rds(opt$sample_sce_filepath)

# Read in mito genes
mito_genes <- unique(readLines(opt$mito_file))

if (is.null(sce_qc$detected)) {
  # We will filter by cells using `scpcaTools::add_cell_mito_qc()`
  sce_qc <- add_cell_mito_qc(sce_qc, mito_genes)
}

# Perform filtering based on specified method
if (opt$filtering_method == "manual") {

  # manually filter the cells
  filtered_sce <- manual_cell_filtering(sce = sce_qc,
                                        mito_percent_cutoff = opt$mito_percent_cutoff,
                                        detected_gene_cutoff = opt$detected_gene_cutoff,
                                        umi_count_cutoff = opt$umi_count_cutoff)


} else if (opt$filtering_method == "miQC") {

  model <- NULL
  filtered_sce <- NULL
  model_attempt <- 0

  # This can fail in a few ways, so we will wrap the next steps in a while/try loop
  while(model_attempt < 3 &&
        (!is(model, "flexmix") || length(model@components) < 2 ) &&
        is.null(filtered_sce)){
    model_attempt <- model_attempt + 1
    try({
      model <- miQC::mixtureModel(sce_qc)
      # filter step can fail too
      filtered_sce <-
        miQC::filterCells(sce_qc,
                          model = model,
                          posterior_cutoff = opt$prob_compromised_cutoff,
                          verbose = FALSE)

      # Save model in metadata for plotting later
      metadata(filtered_sce)$miQC_model <- model

      # Include note in metadata re: filtering
      metadata(filtered_sce)$filtering <- "miQC filtered"
      metadata(filtered_sce)$probability_compromised_cutoff <- opt$prob_compromised_cutoff

    }, silent = TRUE)
  }

  if(is.null(filtered_sce) || is.null(model)) {
    # if miQC failed after 3 attempts then do manual filtering instead
    warning(glue::glue("
                       miQC filtering failed for {opt$sample_sce_filepath}.
                       Using manual filtering instead.
                       "))

    # manually filter the cells
    filtered_sce <- manual_cell_filtering(sce = sce_qc,
                                          mito_percent_cutoff = opt$mito_percent_cutoff,
                                          detected_gene_cutoff = opt$detected_gene_cutoff,
                                          umi_count_cutoff = opt$umi_count_cutoff)

  }
}

# Remove old gene-level rowData statistics and recalculate
drop_cols = colnames(rowData(filtered_sce)) %in% c('mean', 'detected')
rowData(filtered_sce) <- rowData(filtered_sce)[!drop_cols]
filtered_sce <- scater::addPerFeatureQC(filtered_sce)

# Filter the genes (rows)
detected <-
  rowData(filtered_sce)$detected > opt$gene_detected_row_cutoff
expressed <- rowData(filtered_sce)$mean > opt$gene_means_cutoff
filtered_sce <- filtered_sce[detected & expressed, ]

# Save sample and library id in metadata of filtered object
metadata(filtered_sce)$sample <- opt$sample_id
metadata(filtered_sce)$library <- opt$library_id

# Save output filtered sce
readr::write_rds(filtered_sce, output_file)
