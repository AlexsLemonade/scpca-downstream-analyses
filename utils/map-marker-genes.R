# Map provided gene identifiers for downstream use in the marker genes analysis.
# Provided gene identifiers may be in the format of gene symbols or Ensembl IDs
# and should be specified using the --type_of_identifier flag.

# Command line usage:
# Rscript --vanilla 00-map-marker-genes.R \
# --input_marker_gene_list "data/anderson-single-cell/marker-genes/nb_marker_genes.tsv" \
# --input_identifiers "SYMBOL" \
# --output_identifiers "ENSEMBL" \
# --identifier_column_name "gene_symbol" \
# --organism "human" \
# --multi_mappings "list" \
# --output_file "data/anderson-single-cell/marker-genes/mapped_nb_marker_genes.tsv"

#### Load libraries ------------------------------------------------------------

library(magrittr)
library(readr)
library(AnnotationDbi)

#### Command line arguments/options --------------------------------------------

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-l", "--input_marker_gene_list"),
    type = "character",
    default = NULL,
    help = "file path to the unmapped input marker gene list(s)"
  ),
  optparse::make_option(
    c("-i", "--input_identifiers"),
    type = "character",
    default = NULL,
    help = "type of gene identifiers provided in the list(s) -- 'SYMBOL' or
           'ENSEMBL' for example, see `keytypes(org.Xx.eg.db)` for more options"
  ),
  optparse::make_option(
    c("-e", "--output_identifiers"),
    type = "character",
    default = NULL,
    help = "type of gene identifiers to be mapped to and returned-- see
            `keytypes(org.Xx.eg.db)` for more options"
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
    help = "the name of the organism annotation Bioconductor package relevant to
            the genes to be mapped"
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
    default = "all_marker_genes_mapped.tsv",
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

#### Read in data --------------------------------------------------------------

marker_genes <- data.table::fread(opt$input_marker_gene_list, stringsAsFactors = FALSE)

#### Perform mapping -----------------------------------------------------------

# turn column name into symbol for pulling the column info out of data frame
identifier_column_name <- rlang::sym(opt$identifier_column_name)
ids_for_mapping <- marker_genes %>%
  dplyr::pull(identifier_column_name)

# define the annotation packages based on the specified organism
if(opt$organism == "human") {
  library(org.Hs.eg.db)
  annotation_package <- org.Hs.eg.db
} else if (opt$organism == "mouse") {
  library(org.Mm.eg.db)
  annotation_package <- org.Mm.eg.db
} else if (opt$organism == "zebrafish") {
  library(org.Dr.eg.db)
  annotation_package <- org.Dr.eg.db
}

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
