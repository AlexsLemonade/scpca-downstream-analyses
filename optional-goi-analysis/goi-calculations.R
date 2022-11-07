## This script performs optional gene identifier mapping, as well as genes of 
#  interest calculations in preparation for plotting.
## Read in SingleCellExperiment RDS object that has been normalized and has PCA
## embeddings.

# Command line usage:

# Rscript --vanilla goi-calculations.R \
#   --sce "example-results/sample01/library01_miQC_processed_sce.rds" \
#   --library_id "library01" \
#   --seed 2021 \
#   --input_goi_list "example-data/goi-lists/sample01_goi_list.tsv" \
#   --perform_mapping \
#   --sce_rownames_identifier "ENSEMBL" \
#   --provided_identifier "SYMBOL" \
#   --organism "Homo sapiens" \
#   --multi_mappings "first" \
#   --output_directory "example-results/sample01"

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
    help = "Unique ID corresponding to library being analyzed",
  ),
  optparse::make_option(
    c("-s", "--seed"),
    type = "integer",
    default = 2021,
    help = "Random seed to set",
    metavar = "integer"
  ),
  optparse::make_option(
    c("-g", "--input_goi_list"),
    type = "character",
    default = NULL,
    help = "file path to the input genes of interest list, must contain gene 
    identifiers to plot in a column named `gene_id`."
  ),
  optparse::make_option(
    c( "--perform_mapping"),
    action = "store_true",
    default = FALSE,
    help = "specifies whether or not to perform gene identifier mapping; the 
    default is FALSE as mapping is not required when ENSEMBL ids are provided"
  ),
  optparse::make_option(
    c("--sce_rownames_identifier"),
    type = "character",
    default ="ENSEMBL",
    help = "gene identifier type associated with rownames of provided SCE object"
  ),
  optparse::make_option(
    c("--provided_identifier"),
    type = "character",
    default = NULL,
    help = "type of gene identifiers that are provided in the GOI list"
  ),
  optparse::make_option(
    c("-o", "--organism"),
    type = "character",
    default = NULL,
    help = "the scientific name of the organism relevant to the genes to be
            mapped; this workflow currently supports 'Homo sapiens', 
            'Mus musculus', and 'Danio rerio' as options here"
  ),
  optparse::make_option(
    c("-m", "--multi_mappings"),
    type = "character",
    default = "first",
    help = "how should instances of multiple gene identifier mappings be handled
            - may specify 'list' to return all the mappings or 'first' to keep
            only the first mapping"
  ),
  optparse::make_option(
    c("--output_directory"),
    type = "character",
    default = NULL,
    help = "path to output directory"
  ),
  optparse::make_option(
    c("--project_root"),
    type = "character",
    default = NULL,
    help = "the path to the root directory for the R project and where the 
    `utils` folder lives."
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

# Source in goi functions 
source(file.path(project_root, "utils", "provided-goi-analysis-functions.R"))

# Check R and Bioconductor versions
check_r_bioc_versions()

# Load project
setup_renv(project_filepath = project_root)

## Set the seed
set.seed(opt$seed)

# Check that the input file exists
if (!file.exists(opt$sce)){
  stop(paste(opt$sce, "does not exist."))
}

# Check that the outout directory exists
if (!dir.exists(opt$output_directory)) {
  dir.create(opt$output_directory, recursive = TRUE)
}

#### Load libraries ------------------------------------------------------------

suppressPackageStartupMessages({
  library(magrittr)
  library(readr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ComplexHeatmap)
  library(SingleCellExperiment)
})

#### Read in data --------------------------------------------------------------

goi_list <- data.table::fread(opt$input_goi_list, stringsAsFactors = FALSE)
sce <- read_rds(file.path(opt$sce))

# Check for normalized data in SingleCellExperiment object
if(is.null(logcounts(sce))){
  stop("There is no normalized data in the provided SingleCellExperiment object.
       You may want to run the core workflow analysis before re-running this GOI analysis.")
}

# Check for dimensionality reduction results
if(!(any(c("PCA", "UMAP") %in% reducedDimNames(sce)))){
  stop("There are no PCA or UMAP results in the provided SingleCellExperiment object.
       You may want to run the core workflow analysis before re-running this GOI analysis.")
}

# Check that `gene_id` column is in GOI list
if(is.null(goi_list$gene_id)){
  stop("There is no column named `gene_id` in the provided genes of interest list.
       Please rename the column holding your gene identifiers `gene_id` before 
       re-running this GOI analysis.")
}
#### Perform mapping -----------------------------------------------------------

if (opt$perform_mapping == TRUE) {
  # turn column name into symbol for pulling the column info out of data frame
  ids_for_mapping <- goi_list %>%
    dplyr::pull(gene_id)
  
  # Define the annotation packages based on the specified organism
  annotation_list <- list(
    'Homo sapiens' = "org.Hs.eg.db",
    'Mus musculus' = "org.Mm.eg.db",
    'Danio rerio' = "org.Dr.eg.db"  #,
    # new organisms would go here
  )
  # Error handling
  if (!(opt$organism %in% names(annotation_list))) {
    stop(paste(opt$organism, "is not supported"))
  }
  # Load library and assign annotation package
  library(annotation_list[[opt$organism]], character.only = TRUE)
  annotation_package <-
    eval(parse(text = annotation_list[[opt$organism]]))
  
  # Check that provided keytypes are present in the annotation package
  if (!(opt$provided_identifier %in% keytypes(annotation_package))) {
    stop(paste(opt$provided_identifier, "is not a supported type of gene identifiers."))
  }
  
  # Perform mapping
  goi_list <- mapIds(
    # organism annotation package
    annotation_package,
    # the provided gene identifiers
    keys = ids_for_mapping,
    # the type of provided gene identifiers
    keytype = opt$provided_identifier,
    # the type of gene identifiers to map to
    column = opt$sce_rownames_identifier,
    multiVals = opt$multi_mappings
  ) %>%
    tibble::enframe(name = 'gene_id',
                    value = tolower(opt$sce_rownames_identifier)) %>%
    # enframe() makes a `list` column; we will simplify it with unnest()
    # This will result in one row of our data frame per list item
    tidyr::unnest(cols = tolower(opt$sce_rownames_identifier)) %>%
    # grab only the unique rows
    dplyr::distinct() %>%
    dplyr::left_join(goi_list, by = "gene_id")
  
  # define mapped goi associated with rownames of SCE object
  goi_rownames_column <- tolower(opt$sce_rownames_identifier)
  goi_rownames <- goi_list %>%
    dplyr::pull(goi_rownames_column)
  
} else {
  # if no mapping is performed then set goi to input id column
  goi_rownames <- goi_list %>%
    dplyr::pull(gene_id)
}

# Save goi list to file to be used as input for template 
write_tsv(goi_list, file.path(
  opt$output_directory,
  paste0(opt$library_id, "_mapped_genes.tsv")
))

#### Prepare data for plotting -------------------------------------------------

# get logcounts from normalized sce object
normalized_logcounts_matrix <- as.matrix(t(logcounts(sce)))

if(any(colnames(normalized_logcounts_matrix) %in% goi_rownames)) {
  # filter counts matrix to only data associated with the provided genes of
  # interest
  normalized_logcounts_matrix <-
    normalized_logcounts_matrix[, colnames(normalized_logcounts_matrix) %in% goi_rownames]
} else {
  stop(
    "Provided gene identifiers cannot be found in the column names of the 
    logcounts matrix. You may need to re-run script with `--perform_mapping` to 
    map to Ensembl gene identifiers."
  )
}

# transform counts into z-scores
normalized_zscores_matrix <- scale(normalized_logcounts_matrix,
                                       center = TRUE, scale = TRUE)

# prepare the matrix and column annotation for heatmap plotting
if(!is.null(goi_list$gene_set)) {
  if (!is.null(opt$provided_identifier)) {
    normalized_zscores_matrix <- colnames_to_gene_symbols(
      normalized_zscores_matrix,
      goi_list,
      goi_rownames_column,
      "gene_id")
  }
  column_annotation <-
    prepare_heatmap_annotation(normalized_zscores_matrix,
                               goi_list,
                               "gene_id",
                               "gene_set")
  
} else {
  if (!is.null(opt$provided_identifier)) {
    normalized_zscores_matrix <-
      colnames_to_gene_symbols(
        normalized_zscores_matrix,
        goi_list,
        goi_rownames_column,
        "gene_id"
      )
  } 
  column_annotation <- NULL
}

# Convert heatmap matrix to a sparse matrix and save to file
normalized_zscores_matrix <- as(normalized_zscores_matrix, "sparseMatrix")
Matrix::writeMM(normalized_zscores_matrix,
          file.path(
            opt$output_directory,
            paste0(opt$library_id, "_normalized_zscores.mtx")
          ))

# Save heatmap column annotation to file
write_rds(column_annotation, file.path(
  opt$output_directory,
  paste0(opt$library_id, "_heatmap_annotation.rds")
))
