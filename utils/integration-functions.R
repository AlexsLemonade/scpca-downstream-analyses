## Custom functions for data integration.

# #### USAGE
# This script is intended to be sourced in the script as follows:
#
# source(file.path("utils", "integration-functions.R"))

plot_umap_panel <- function(sce,
                            cell_label_column,
                            umap_name,
                            plot_colors,
                            plot_title = NULL, 
                            legend_title = NULL,
                            legend_labels = NULL,
                            seed = NULL){
  
  set.seed(seed)
  
  # check that plot_colors is equal to the number of categories present in the cell label column
  color_categories <- unique(colData(sce)[,cell_label_column])
  if(length(plot_colors) != length(color_categories)){
    stop("Number of colors provided must be equal to the number of categories used to classify cells in
         the specified cell_label_column.")
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

plot_integration_umap <- function(sce,
                                  integration_method,
                                  cell_label_column,
                                  include_legend_counts = TRUE,
                                  legend_title = NULL, 
                                  legend_labels = NULL,
                                  plot_colors = NULL,
                                  seed = NULL) {
  
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
  
  umap <- plot_umap_panel(sce = sce,
                          cell_label_column,
                          umap_name = umap_name, 
                          plot_colors = plot_colors,
                          plot_title = integration_method,
                          legend_title = legend_title,
                          legend_labels = legend_labels,
                          seed = seed)
  
  return(umap)
}

perform_dim_reduction <- function(merged_sce,
                                  prefix = NULL){
  
  # create pca and umap names
  pca_name <- "PCA"
  umap_name <- "UMAP"
  if(!is.null(prefix)){
    pca_name <- paste(prefix, pca_name, sep = "_")
    umap_name <- paste(prefix, umap_name, sep = "_")
  }
  
  # add UMAP
  merged_sce <- scater::runUMAP(merged_sce, dimred = pca_name, name = umap_name)
  
  return(merged_sce)
  
}
