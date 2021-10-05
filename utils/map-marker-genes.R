# Map provided gene identifiers for downstream use in the marker genes analysis.
# Provided gene identifiers may be in the format of gene symbols or Ensembl IDs
# and should be specified using the --type_of_identifier flag.

# Command line usage:
# Rscript --vanilla 00-map-marker-genes.R \
#         --input_marker_gene_list "all_marker_genes.tsv" \
#         --type_of_identifier "gene symbol" \
#         --output_file "mapped_all_marker_genes.tsv"

#### Load libraries ------------------------------------------------------------

library(magrittr)
library(readr)
library(AnnotationDbi)
library(org.Hs.eg.db)

#### Command line arguments/options --------------------------------------------

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-i", "--input_marker_gene_list"),
    type = "character",
    default = NULL,
    help = "file path to the unmapped input marker gene list(s)"
  ),
  optparse::make_option(
    c("-t", "--type_of_identifier"),
    type = "character",
    default = NULL,
    help = "type of gene identifiers provided in the list(s) -- can be gene
            symbol or ensembl gene id"
  ),
  optparse::make_option(
    c("-o", "--output_file"),
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

marker_genes <- data.table::fread(opt$input_marker_gene_list)

#### Perform mapping -----------------------------------------------------------

if (opt$type_of_identifier == "gene symbol") {
  marker_genes_mapped <- mapIds(
    # organism annotation package
    org.Hs.eg.db,
    # the provided gene identifiers
    keys = marker_genes$gene_symbol,
    # the type of provided gene identifiers
    keytype = "SYMBOL",
    # the type of gene identifiers to map to
    column = "ENSEMBL",
    multiVals = "list"
  ) %>%
    tibble::enframe(name = "gene_symbol", value = "gene_id") %>%
    # enframe() makes a `list` column; we will simplify it with unnest()
    # This will result in one row of our data frame per list item
    tidyr::unnest(cols = gene_id) %>%
    # grab only the unique rows
    dplyr::distinct() %>%
    # join the remaining columns
    dplyr::left_join(marker_genes, by = "gene_symbol")
} else if (opt$type_of_identifier == "ensembl id") {
  marker_genes_mapped <- mapIds(
    # organism annotation package
    org.Hs.eg.db,
    # the provided gene identifiers
    keys = marker_genes$gene_id,
    # the type of provided gene identifiers
    keytype = "ENSEMBL",
    # the type of gene identifiers to map to
    column = "SYMBOL",
    multiVals = "list"
  ) %>%
    tibble::enframe(name = "gene_id", value = "gene_symbol") %>%
    # enframe() makes a `list` column; we will simplify it with unnest()
    # This will result in one row of our data frame per list item
    tidyr::unnest(cols = gene_symbol) %>%
    # grab only the unique rows
    dplyr::distinct() %>%
    # join the remaining columns
    dplyr::left_join(marker_genes, by = "gene_id")
} else {
  stop(
    "Incorrect input gene identifier columns. Columns containing gene symbols
       should be named `gene_symbol` while those containing Ensembl gene IDs
       should be named `gene_id`"
  )
}

#### Save mapped object to file ------------------------------------------------

write_tsv(marker_genes_mapped, opt$output_file)

