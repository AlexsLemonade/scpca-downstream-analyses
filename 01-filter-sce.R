## Read in and filter SingleCellExperiment data

# Command line usage:

# Rscript --vanilla 01-filter-sce.R \
#   --sample_sce_filepath "data/anderson-single-cell/pre-filtered-sce/GSM4186961_pre-filtered_sce.rds" \
#   --sample_name "GSM4186961" \
#   --mito_file data/Homo_sapiens.GRCh38.103.mitogenes.txt \
#   --output_filepath data/anderson-single-cell/results/sample_filtered_sce.rds \
#   --output_plots_directory plots \
#   --seed 2021 \
#   --gene_detected_row_cutoff 5 \
#   --gene_means_cutoff 0.1 \
#   --prob_compromised_cutoff 0.75 \
#   --filtering_method "miQC"

## Set up -------------------------------------------------------------

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
    c("-r", "--sample_name"),
    type = "character",
    default = NULL,
    help = "name of sample being analyzed",
  ),
  optparse::make_option(
    c("-i", "--mito_file"),
    type = "character",
    default = NULL,
    help = "the path to the mito file"
  ),
  optparse::make_option(
    c("--output_plots_directory"),
    type = "character",
    default = NULL,
    help = "output plots directory"
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
    type = "integer",
    default = 0.1,
    help = "gene mean expression cutoff",
    metavar = "integer"
  ),
  optparse::make_option(
    c("-c", "--mito_percent_cutoff"),
    type = "integer",
    default = 20,
    help = "cell mitochondrial percent cutoff -- not needed for miQC filtering",
    metavar = "integer"
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
    type = "integer",
    default = 0.75,
    help = "cell probability compromised cutoff -- this should be provided if
            `prob_compromised` already exists in the object",
    metavar = "integer"
  ),
  optparse::make_option(
    c("-f", "--filtering_method"),
    type = "character",
    default = "miQC",
    help = "the selected filtered method -- can be miQC or manual; will be miQC
            by default"
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## Load project

# `here::here()` looks at a number of criteria to identify the root 
# directory, including whether or not there is a .Rproj file present,
# so we can pass this to `renv::load()` to load the project file
renv::load(here::here())

# Check that R version us at least 4.1
if (! (R.version$major == 4 && R.version$minor >= 1)){
  stop("R version must be at least 4.1")
}

# Check that Bioconductor version is 3.14
if (packageVersion("BiocVersion") < 3.14){
  stop("Bioconductor version is less than 3.14")
}

## Load libraries
library(scater)
library(scran)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(scpcaTools)
library(SingleCellExperiment)
library(cowplot)
library(scuttle)
library(tryCatchLog)

if (!("miQC" %in% installed.packages())) {
  if (!requireNamespace("BiocManager")) {
    install.packages("BiocManager")
  }
  BiocManager::install("miQC")
}

## Set the seed
set.seed(opt$seed)

## Define file paths

# Directory and file to save output
output_file <- opt$output_filepath
output_dir <- dirname(output_file)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

plots_dir <- file.path(opt$output_plots_directory)
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

output_filtered_cell_plot <- file.path(
  plots_dir,
  paste0(opt$sample_name, "_", opt$filtering_method, "_cell_filtering.png"))

#### Filter data ---------------------------------------------------------------

# Read in sce object
sce_qc <- readr::read_rds(opt$sample_sce_filepath)

# Read in mito genes
mito_genes <- unique(readLines(opt$mito_file))

if (is.null(sce_qc$detected)) {
  # We will filter by features (genes) using `scater::addPerFeatureQC()` and by
  # cells using `scpcaTools::add_cell_mito_qc()`
  sce_qc <- add_cell_mito_qc(sce_qc, mito_genes)
}

if (!is.null(sce_qc$prob_compromised)) {
  
  # Remove `prob_compromised` if it exists, as this will cause errors with
  # plotModel
  sce_qc$prob_compromised <- NULL
}

# Perform filtering based on specified method
if (!opt$filtering_method %in% c("manual", "miQC")) {
  
  stop("Incorrect filtering method. Specify `manual` or `miQC` filtering.")
  
} else if (opt$filtering_method == "manual") {
  
  # Filter the cells
  mito_filter <-
    colData(sce_qc)$mito_percent < opt$mito_percent_cutoff
  gene_filter <-
    colData(sce_qc)$detected > opt$detected_gene_cutoff
  sum_filter <- colData(sce_qc)$sum > opt$umi_count_cutoff
  filtered_sce <-
    sce_qc[, mito_filter & gene_filter & sum_filter]
  
  coldata_qc <- data.frame(colData(sce_qc))
  filtered_cell_plot <-
    ggplot(coldata_qc, aes(x = sum, y = detected, color = mito_percent)) +
    geom_point(alpha = 0.5) +
    scale_color_viridis_c() +
    labs(x = "Total Count", y = "Number of Genes Expressed", color = "Mitochondrial \nFraction") +
    theme_classic() +
    geom_hline(yintercept = opt$detected_gene_cutoff) +
    geom_vline(xintercept = opt$umi_count_cutoff)
  
  # Include note in metadata re: filtering
  metadata(filtered_sce)$filtering <- "manually filtered"
  
} else if (opt$filtering_method == "miQC") {
  # Rename columns for `mixtureModel()`
  names(colData(sce_qc)) <-
    stringr::str_replace(names(colData(sce_qc)),
                         "^mito_",
                         "subsets_mito_")
  tryCatch(
    expr = {
      # Generate miQC model
      sce_model <- miQC::mixtureModel(sce_qc)
      
    },
    error = function(e) {
      print(
        paste0(
          "miQC filtering failed. Skipping filtering for sample ",
          opt$sample_sce_filepath,
          ". Try `manual` filtering instead."
        )
      )
    }
    
  )
  
  if (exists("sce_model")) {
    # Filter cells
    filtered_sce <- miQC::filterCells(sce_qc, sce_model)
    
    # Plot model
    filtered_model_plot <- miQC::plotModel(sce_qc, sce_model)
    
    # Plot filtering
    filtered_cell_plot <- miQC::plotFiltering(sce_qc, sce_model)
    
    # Combine plots
    filtered_cell_plot <-
      ggarrange(filtered_model_plot,
                filtered_cell_plot,
                ncol = 1,
                nrow = 2)
    # Include note in metadata re: filtering
    metadata(filtered_sce)$filtering <- "miQC filtered"
    
  } else {
    filtered_sce <- sce_qc
    # Include note in metadata re: failed filtering
    metadata(filtered_sce)$filtering <- "miQC filtering failed"
    # Implement `plotMetrics` when a model cannot be generated
    filtered_cell_plot <- miQC::plotMetrics(filtered_sce)
  }
  
  # Save plot
  ggsave(output_filtered_cell_plot, filtered_cell_plot)
  
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

# Save output filtered sce
readr::write_rds(filtered_sce, output_file)
