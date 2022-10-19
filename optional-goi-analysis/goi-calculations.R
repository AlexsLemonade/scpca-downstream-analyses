## This script performs optional gene identifier mapping, as well as genes of 
#  interest calculations in preparation for plotting.
## Read in SingleCellExperiment RDS object that has been normalized and has PCA
## embeddings.

# Command line usage:

# Rscript --vanilla goi-calculations.R \
#   --sce "example-results/sample01/library01_miQC_processed_sce.rds" \
#   --library_id "library01" \
#   --seed 2021 \
#   --input_goi_list "data/anderson-single-cell/goi-lists/nb_goi_list.tsv" \
#   --perform_mapping TRUE \
#   --input_identifiers "SYMBOL" \
#   --output_identifiers "ENSEMBL" \
#   --identifier_column_name "gene_symbol" \
#   --organism "Homo sapiens" \
#   --multi_mappings "list" \
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
    c("-g", "--input_goi_list"),
    type = "character",
    default = NULL,
    help = "file path to the unmapped input genes of interest list"
  ),
  optparse::make_option(
    c( "--perform_mapping"),
    action = "store_false",
    help = "specifies whether or not to perform gene identifier mapping"
  ),
  optparse::make_option(
    c("--input_identifiers"),
    type = "character",
    default = NULL,
    help = "type of gene identifiers provided in the list -- 'SYMBOL' or
           'ENSEMBL' for example, see `keytypes()` in the organism's annotation
            package for more options"
  ),
  optparse::make_option(
    c("--output_identifiers"),
    type = "character",
    default = NULL,
    help = "type of gene identifiers to be mapped to and returned"
  ),
  optparse::make_option(
    c("-c", "--identifier_column_name"),
    type = "character",
    default = NULL,
    help = "the name of the column in the input file that contains the gene
            identifiers to be converted from"
  ),
  optparse::make_option(
    c("--gene_set_column_name"),
    type = "character",
    default = NULL,
    help = "the name of the column in the input file that contains the gene set
            information"
  ),
  optparse::make_option(
    c("-o", "--organism"),
    type = "character",
    default = NULL,
    help = "the scientific name of the organism relevant to the genes to be
            mapped, 'Homo sapiens' or 'Mus musculus' for examples"
  ),
  optparse::make_option(
    c("-m", "--multi_mappings"),
    type = "character",
    default = NULL,
    help = "how should instances of multiple gene identifier mappings be handled
            - may specify 'list' to return all the mappings or 'first' to keep
            only the first mapping"
  ),
  optparse::make_option(
    c("--output_directory"),
    type = "character",
    default = NULL,
    help = "path to output directory that would hold the TSV files with the 
    clustering stats"
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
})

#### Read in data --------------------------------------------------------------

goi_list <- data.table::fread(opt$input_goi_list, stringsAsFactors = FALSE)
sce <- read_rds(file.path(opt$sce))

#### Perform mapping -----------------------------------------------------------

if (opt$perform_mapping == TRUE) {
  # turn column name into symbol for pulling the column info out of data frame
  identifier_column_name <- rlang::sym(opt$identifier_column_name)
  ids_for_mapping <- goi_list %>%
    dplyr::pull(identifier_column_name)
  
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
  
  # Perform mapping
  goi_list <- mapIds(
    # organism annotation package
    annotation_package,
    # the provided gene identifiers
    keys = ids_for_mapping,
    # the type of provided gene identifiers
    keytype = opt$input_identifiers,
    # the type of gene identifiers to map to
    column = opt$output_identifiers,
    multiVals = opt$multi_mappings
  ) %>%
    tibble::enframe(
      name = opt$identifier_column_name,
      value = tolower(opt$output_identifiers)
    ) %>%
    # enframe() makes a `list` column; we will simplify it with unnest()
    # This will result in one row of our data frame per list item
    tidyr::unnest(cols = tolower(opt$output_identifiers)) %>%
    # grab only the unique rows
    dplyr::distinct() %>%
    # join the remaining columns
    dplyr::left_join(goi_list, by = opt$identifier_column_name)
  
  # Save mapped object to file
  write_tsv(goi_list, file.path(
    opt$output_directory,
    paste0(opt$library_id, "_mapped_genes.tsv")
  ))
}

#### Prepare data for plotting -------------------------------------------------

identifier_column_name_sym <- rlang::sym(opt$identifier_column_name)

# get logcounts from normalized sce object
normalized_sce_logcounts_matrix <- as.matrix(t(logcounts(sce)))

# filter counts matrix to only data associated with the provided genes of
# interest
normalized_sce_logcounts_matrix <- normalized_sce_logcounts_matrix[,colnames(normalized_sce_logcounts_matrix) %in% goi_list[[tolower(opt$output_identifiers)]]]

# transform counts into z-scores
normalized_sce_zscores_matrix <- scale(normalized_sce_logcounts_matrix,
                                       center = TRUE, scale = TRUE)

# prepare the matrix and column annotation for heatmap plotting
if(!is.null(opt$gene_set_column_name)) {
  if (!is.null(opt$output_identifiers)) {
    normalized_sce_zscores_matrix <- colnames_to_gene_symbols(
      normalized_sce_zscores_matrix,
      goi_list,
      tolower(opt$output_identifiers),
      opt$identifier_column_name
    )
    gene_id_column <- tolower(opt$output_identifiers)
    column_annotation <-
      prepare_heatmap_annotation(
        normalized_sce_zscores_matrix,
        goi_list,
        gene_id_column,
        opt$gene_set_column_name
      )
  } else {
    gene_id_column <- opt$identifier_column_name
    column_annotation <-
      prepare_heatmap_annotation(
        normalized_sce_zscores_matrix,
        goi_list,
        gene_id_column,
        opt$gene_set_column_name
      )
  }
} else {
  if (!is.null(opt$output_identifiers)) {
    normalized_sce_zscores_matrix <-
      colnames_to_gene_symbols(
        normalized_sce_zscores_matrix,
        goi_list,
        tolower(opt$output_identifiers),
        opt$identifier_column_name
      )
    column_annotation <- NULL
  } else {
    column_annotation <- NULL
  }
}

# Save matrix object to file
write_rds(normalized_sce_zscores_matrix,
          file.path(
            opt$output_directory,
            paste0(opt$library_id, "_normalized_sce_zscores.rds")
          ))

# Save heatmap column annotation to file
write_rds(column_annotation, file.path(
  opt$output_directory,
  paste0(opt$library_id, "_heatmap_annotation.rds")
))
