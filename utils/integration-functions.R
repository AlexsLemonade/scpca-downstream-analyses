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
                                label.theme = element_text(size = 12))) +
    theme(legend.position = "none",
          text = element_text(size = 14),
          legend.title = element_text(size = 16)) +
    theme_bw() +
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

plot_asw <- function(asw_df,
                     seed = NULL, 
                     batch_label) {
  # Purpose: Function to plot average silhouette width (ASW) metric
  #
  # Args: 
  #   asw_df: data frame containing silhouette width values calculated on both
  #   integrated and unintegrated SCEs. Expected columns are at least
  #   `rep`, `silhouette_width`, `silhouette_cluster`, and `integration_method`.
  #   seed: for sina plot reproducibility
  #   batch_label: label to include in plot for batch, if by_batch is TRUE
  
  # Set seed if given
  set.seed(seed)
  
  # Check that all expected columns are present in dataframe
  expected_columns <- c("pc_name", "rep", "silhouette_width", "silhouette_cluster")
  if(!all(expected_columns%in% colnames(asw_df))){
    stop("Required columns are missing from input dataframe, make sure that `calculate_silhouette_width` has been run successfully.")
  }
  
  # Prepare and plot data by batch
  asw_plot <- asw_df |>
    dplyr::mutate(
      integration_method = dplyr::if_else(
        pc_name == "PCA",
        "Pre-Integration",
        stringr::str_remove(pc_name, "_PCA")
      ),
      integration_method_factor = forcats::fct_relevel(integration_method, "Pre-Integration")
    ) |>
    # the `silhouette_cluster` column contains the true identity; rename for ease
    dplyr::group_by(rep, integration_method_factor, silhouette_cluster) |>
    dplyr::summarize(asw = mean(abs(silhouette_width))) |>
    dplyr::ungroup() |>
    ggplot() +
    aes(x = integration_method_factor,
        y = asw,
        color = silhouette_cluster) +
    ggforce::geom_sina(size = 1,
                       alpha = 0.5,
                       position = position_dodge(width = 0.5)) +
    # add median/IQR pointrange to plot
    stat_summary(
      aes(group = silhouette_cluster),
      color = "black",
      fun = "median",
      fun.min = function(x) {
        quantile(x, 0.25)
      },
      fun.max = function(x) {
        quantile(x, 0.75)
      },
      geom = "pointrange",
      position = position_dodge(width = 0.5),
      size = 0.2
    ) +
    guides(color = guide_legend(override.aes = list(size = 3),
                                label.theme = element_text(size = 12)))
  
  # Add shared labeling
  asw_plot <- asw_plot + 
    labs(
      x = "Integration method",
      y = "Average silhouette width",
      color = "Batch"
    )
  
  # return the plot
  return(asw_plot)
  
}
