# Map provided gene identifiers for downstream use in the marker genes analysis.
# The type of gene identifiers you are converting *from*, the column in which
# they can be found in the input file, and the type of gene identifiers you are
# converting *to* should all be specified using the command line flags outlined
# below. 
# 
# Run `Rscript utils/map-marker-genes.R --help` for more on what values should
# be specified for each of the command line flags.

# Command line usage:
# Rscript --vanilla utils/map-marker-genes.R \
# --input_marker_gene_list "data/anderson-single-cell/marker-genes/nb_marker_genes.tsv" \
# --input_identifiers "SYMBOL" \
# --output_identifiers "ENSEMBL" \
# --identifier_column_name "gene_symbol" \
# --organism "Homo sapiens" \
# --multi_mappings "list" \
# --output_file "data/anderson-single-cell/marker-genes/mapped_nb_marker_genes.tsv"

#### Command line arguments/options --------------------------------------------

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-l", "--input_marker_gene_list"),
    type = "character",
    default = NULL,
    help = "file path to the unmapped input marker gene list"
  ),
  optparse::make_option(
    c("-i", "--input_identifiers"),
    type = "character",
    default = NULL,
    help = "type of gene identifiers provided in the list -- 'SYMBOL' or
           'ENSEMBL' for example, see `keytypes()` in the organism's annotation
            package for more options"
  ),
  optparse::make_option(
    c("-e", "--output_identifiers"),
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
    c("-p", "--output_file"),
    type = "character",
    default = NULL,
    help = "file path to the mapped output marker gene file"
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Define file paths
output_file <- opt$output_file
output_dir <- dirname(output_file)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#### Load libraries ------------------------------------------------------------

library(magrittr)
library(readr)
library(AnnotationDbi)

#### Read in data --------------------------------------------------------------

marker_genes <- data.table::fread(opt$input_marker_gene_list, stringsAsFactors = FALSE)

#### Perform mapping -----------------------------------------------------------

# turn column name into symbol for pulling the column info out of data frame
identifier_column_name <- rlang::sym(opt$identifier_column_name)
ids_for_mapping <- marker_genes %>%
  dplyr::pull(identifier_column_name)

# Define the annotation packages based on the specified organism
annotation_list <- list(
  'Homo sapiens' = "org.Hs.eg.db",
  'Mus musculus' = "org.Mm.eg.db",
  'Danio rerio' = "org.Dr.eg.db"  #,
  # new organisms would go here
)
# Error handling
if (!(opt$organism %in% names(annotation_list))){
  stop(paste(opt$organism, "is not supported"))
}
# Load library and assign annotation package
library(annotation_list[[opt$organism]], character.only = TRUE)
annotation_package <- eval(parse(text = annotation_list[[opt$organism]]))

# Perform mapping
marker_genes_mapped <- mapIds(
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
  dplyr::left_join(marker_genes, by = opt$identifier_column_name)

#### Save mapped object to file ------------------------------------------------

write_tsv(marker_genes_mapped, opt$output_file)
