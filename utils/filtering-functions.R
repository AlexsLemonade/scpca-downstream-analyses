## Custom functions for filtering.

# #### USAGE
# This script is intended to be sourced in the script as follows:
#
# source(file.path("utils", "filtering-functions.R"))

manual_cell_filtering <- function(sce,
                                  mito_percent_cutoff,
                                  detected_gene_cutoff,
                                  umi_count_cutoff){
  # Purpose: Remove low quality cells from a SingleCellExperiment object using 
  # provided manual thresholds. 
  
  # Args:
  #   sce: Unfiltered SingleCellExperiment object
  #   mito_percent_cutoff: Maximum percentage of mitochondrial reads per cell after filtering. 
  #   detected_gene_cutoff: Minimum number of genes detected per cell after filtering. 
  #   umi_count_cutoff: Minimum number of UMI per cell after filtering. 
  
  # filter cells based on provided cutoffs
  mito_filter <-
    colData(sce)$subsets_mito_percent < mito_percent_cutoff
  gene_filter <-
    colData(sce)$detected > detected_gene_cutoff
  sum_filter <- colData(sce)$sum > umi_count_cutoff
  filtered_sce <-
    sce[, mito_filter & gene_filter & sum_filter]

  # Include note in metadata re: filtering
  metadata(filtered_sce)$filtering <- "manually filtered"
  metadata(filtered_sce)$mito_percent_cutoff <- mito_percent_cutoff
  metadata(filtered_sce)$umi_count_cutoff <- umi_count_cutoff
  
  return(filtered_sce)
  
}

plot_manual_filtering <- function(sce,
                                  detected_gene_cutoff,
                                  umi_count_cutoff){
  
  # Purpose: Create summary plot of manual filtering. 
  
  # Args: 
  #   sce: Unfiltered SingleCellExperiment object.
  #   detected_gene_cutoff: Minimum number of genes detected per cell after filtering. 
  #   umi_count_cutoff: Minimum number of UMI per cell after filtering. 
  
  # grab coldata as dataframe from sce
  coldata_qc <- data.frame(colData(sce))
  
  filtered_cell_plot <-
    ggplot(coldata_qc, aes(x = sum, y = detected, color = subsets_mito_percent)) +
    geom_point(alpha = 0.5) +
    scale_color_viridis_c() +
    labs(x = "Total Count",
         y = "Number of Genes Expressed",
         color = "Mitochondrial\nFraction") +
    theme_classic() +
    geom_hline(yintercept = detected_gene_cutoff) +
    geom_vline(xintercept = umi_count_cutoff)
  
  return(filtered_cell_plot)
  
}
