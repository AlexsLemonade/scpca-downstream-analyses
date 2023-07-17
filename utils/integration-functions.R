## Custom functions for data integration.

# #### USAGE
# This script is intended to be sourced in the script as follows:
#
# source(file.path("utils", "integration-functions.R"))

plot_integration_umap <- function(sce,
                                  integration_method,
                                  cell_label_column,
                                  legend_title = NULL, 
                                  legend_labels = NULL,
                                  plot_colors = NULL,
                                  plot_title = NULL,
                                  seed = NULL) {
  
  # Purpose: Generate UMAP plots for a merged or integrated SingleCellExperiment object containing multiple batches
  
  # Args:
  #   sce: SingleCellExperiment object containing integration results
  #   integration_method: the method used to perform integration. i.e. "fastMNN", "harmony"
  #   cell_label_column: the name of the colData column to label cells by
  #   legend_title: the title of the plot legend
  #   legend_labels: a vector of names to use to label the legend
  #   plot_colors: a vector of colors to use when plotting; default is NULL
  #   plot_title: the title of the plot
  #   seed: an integer to set the seed as for reproducibility
  
  # set seed for reproducibility 
  set.seed(seed)
  
  # check that column to label cells by is present in colData
  if(!cell_label_column %in% colnames(colData(sce))){
    stop("Provided cell_label_column should be present in the SCE object.")
  }
  
  
  # Define colors if not provided, or if provided check the size
  if (is.null(plot_colors)) {
    num_colors <- length(unique(sce[[cell_label_column]]))
    plot_colors <- rainbow(num_colors)    
  } else {
    if (!(length(plot_colors)) == length(unique(sce[[cell_label_column]]))) {
      stop("The number of provided colors does not match the number of labels.")
    }
  }
  
  if(integration_method == "unintegrated"){
    umap_name <- "UMAP"
  } else {
    # grab dim reduction name to use for plotting
    umap_name <- paste0(integration_method, "_UMAP")
  }
  
  # randomly shuffle cells prior to plotting
  col_order <- sample(ncol(sce))
  shuffled_sce <- sce[,col_order]
  
  # create umap and label with provided cell label column
  umap <- scater::plotReducedDim(shuffled_sce,
                                 dimred = umap_name,
                                 colour_by = cell_label_column,
                                 point_size = 0.1,
                                 point_alpha = 0.4) +
    scale_color_manual(values = plot_colors, 
                       name = legend_title, 
                       labels = legend_labels) +
    # relabel legend and resize dots
    guides(color = guide_legend(override.aes = list(size = 3),
                                label.theme = element_text(size = 16))) +
    theme(legend.position = "none",
          text = element_text(size = 14),
          legend.title = element_text(size = 16)) +
    labs(
      title = plot_title
    )
  
  return(umap)
}

add_integrated_pcs <- function(merged_sce,
                               integrated_pcs,
                               integration_method){  
  # Purpose: Add integrated PCs as calculated by the integration method to the 
  #          merged SingleCellExperiment object 
  
  # Args:
  #   merged_sce: a merged SingleCellExperiment object
  #   integrated_pcs: the full set of integrated pcs
  #   integration_method: the method used to perform integration. i.e. "fastMNN", "harmony"
  
  # check that integration method is provided
  if(is.null(integration_method)){
    stop("Integration method is missing.")
  } else if (!integration_method %in% c("fastMNN", "harmony")) {
    stop("Integration method can be 'fastMNN' or 'harmony'.")
  }
  
  # create pca and umap names
  pca_name <- glue::glue("{integration_method}_PCA")
  umap_name <- glue::glue("{integration_method}_UMAP")
  
  # add UMAP
  reducedDim(merged_sce, pca_name) <- integrated_pcs
  merged_sce <- scater::runUMAP(merged_sce, dimred = pca_name, name = umap_name)
  
  return(merged_sce)
  
}
