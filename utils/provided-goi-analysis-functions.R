## Custom functions to be sourced in the provided goi analysis reports template
## notebook.

# #### USAGE
# This script is intended to be sourced in the script as follows:
#
# source(file.path("utils", "provided-goi-analysis-functions.R"))

colnames_to_plotting_symbols <- function(normalized_sce_matrix,
                                    goi_list,
                                    sce_rownames_column,
                                    plotting_column) {
  # Given a normalized SingleCellExperiment matrix, the name of the column with
  # the SCE rowname identifiers, the name of the associated plotting identifiers,
  # convert the matrix column names into the plotting ids for heatmap plotting.
  #
  # Args:
  #   normalized_sce_matrix: matrix retrieved from a normalized
  #                          SingleCellExperiment object
  #   goi_list: data frame with provided genes of interest relevant to the data
  #             in the SingleCellExperiment object
  #   sce_rownames_column: name of the column with the SCE row identifiers
  #   plotting_column: name of the column with the identifiers to be used for 
  #                    plotting
  
  
  # turn the gene symbol column name into a symbol for use when subsetting
  plotting_column_sym <- rlang::sym(plotting_column)
  
  # we will want to replace the column names our matrix with the relevant gene
  # symbols for plotting
  map_ensembl_symbols <-
    data.frame(ensembl = colnames(normalized_sce_matrix)) %>%
    dplyr::left_join(goi_list, by = c("ensembl" = sce_rownames_column))
  
  colnames(normalized_sce_matrix) <-
    map_ensembl_symbols[[plotting_column_sym]]
  
  return(normalized_sce_matrix)
  
}

prepare_heatmap_annotation <- function(normalized_sce_matrix,
                                       goi_list,
                                       gene_id_column,
                                       gene_set) {
  # Given a normalized SingleCellExperiment matrix, the name of the column with
  # the gene identifiers (can be ensembl or gene symbols if present), and
  # optionally (at command line) the name of column with the gene set info,
  # prepare the column annotation object for the heatmap.
  #
  # Args:
  #   normalized_sce_matrix: matrix retrieved from a normalized
  #                          SingleCellExperiment object
  #   goi_list: data frame with genes of interest relevant to the data in the
  #                 SingleCellExperiment object
  #   gene_id_column: name of the column with the gene identifiers (can be
  #                   ensembl or gene symbols if provided)
  #   gene_set: name of the column with the gene set information
  
  gene_id_column_sym <- rlang::sym(gene_id_column)
  gene_set_sym <- rlang::sym(gene_set)
  
  # create column annotation data frame
  annotation_df <- goi_list %>%
    dplyr::select(gene_id_column_sym, gene_set_sym) %>%
    dplyr::distinct()
  
  rownames(annotation_df) <- NULL
  
  # filter to ensure only those genes of the goi list that had data stored
  # in the sce object are kept in the annotation object
  annotation_df <- annotation_df %>%
    dplyr::filter(!!gene_id_column_sym %in% colnames(normalized_sce_matrix)) %>%
    tibble::column_to_rownames(gene_id_column)
  
  gene_set_names <- annotation_df %>%
    dplyr::arrange(gene_set) %>%
    dplyr::pull(gene_set) %>%
    unique()
  
  if (length(gene_set_names) > 8) {
    stop(
      "There are more than eight unique gene set values. Please choose a different gene set column with eight or less unique values or re-run without providing one."
    )
  } else {
    # get colors for the annotation object
    annotation_colors <-
      palette.colors(palette = "Okabe-Ito")[2:(length(gene_set_names) + 1)]

    gene_set_colors = annotation_colors[1:length(gene_set_names)]
    names(gene_set_colors) = gene_set_names

    # create the column annotation for the ComplexHeatmap
    column_annotation <- HeatmapAnnotation(
      df = annotation_df,
      col = list(gene_set_names = gene_set_colors),
      annotation_label = "Gene Set"
    )
  }
}

prepare_expression_df <- function(normalized_sce,
                                  goi_list,
                                  sce_rownames_column) {
  # Given a normalized SingleCellExperiment object and a goi list,
  # prepare a gene expression data frame for plotting.
  #
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object
  #   goi_list: data frame with genes of interest relevant to the data in the
  #                 SingleCellExperiment object
  #   sce_rownames_column: name of the column with the SCE rowname identifiers
  #                              used for mapping
  
  # Turn the SCE rownames into symbols for use when subsetting
  sce_rownames_column_sym <- rlang::sym(sce_rownames_column)
  
  # Prepare a data frame containing the logcounts expression values associated
  # with each marker gene symbol
  expression_means_df <- logcounts(normalized_sce[rownames(normalized_sce) %in% goi_list[[sce_rownames_column_sym]],]) %>%
    t() %>%
    as.matrix() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell_barcode") %>%
    tidyr::pivot_longer(
      cols = -c("cell_barcode"),
      names_to = sce_rownames_column,
      values_to = "gene_expression"
    ) %>%
    dplyr::left_join(goi_list, by = sce_rownames_column)
  
  return(expression_means_df)
  
}

plot_goi_expression_sina <- function(normalized_sce,
                                     goi_list,
                                     sce_rownames_column,
                                     optional_plotting_column = NULL,
                                     use_rownames = TRUE) {
  
  # Given a normalized SingleCellExperiment object and a goi list,
  # plot the gene expression of each gene of interest on a sina plot against the
  # average mean expression of all the genes of interest.
  #
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object
  #   goi_list: data frame with genes of interest relevant to the data in the
  #                 SingleCellExperiment object
  #   sce_rownames_column: name of the column with the SCE rowname identifiers
  #                              used for mapping
  #   optional_plotting_column: name of the column with the identifiers that
  #                             would be used for plotting if provided
  #   use_rownames: indicates whether or not the SCE rowname identifiers should
  #                 be used when plotting; default is TRUE
  #
  # Returns:
  #   sina_expression_plot: sina plots displaying the logcounts expression
  #                         values for each gene of interest symbol relative to
  #                         the mean logcounts expression of all the genes of
  #                         interest
  
  # Turn the SCE rownames into symbols for use when subsetting
  sce_rownames_column_sym <- rlang::sym(sce_rownames_column)
  
  # Run the `prepare_expression_df` function
  expression_means_df <- prepare_expression_df(normalized_sce,
                                               goi_list,
                                               sce_rownames_column)
  
  # Check for optional plotting identifiers if `use_rownames` is FALSE
  if(!use_rownames && is.null(optional_plotting_column)){
    stop("When `use_rownames = FALSE`, a column name that contains optional 
    plotting identifiers should be provided. Please provide this column name or
    implement `use_rownames = TRUE`.")
  }
  
  if(!is.null(optional_plotting_column)) {
    optional_plotting_column_sym <- rlang::sym(optional_plotting_column)
    
    expression_means_df <- expression_means_df %>%
      dplyr::select(-cell_barcode) %>%
      dplyr::distinct()

    # calculate average gene expression
    avg_gene_exp = mean(colMeans(logcounts(normalized_sce)))
    
    sina_expression_plot <- ggplot(expression_means_df, aes(
      x = !!optional_plotting_column_sym,
      y = gene_expression,
      color = gene_expression)) +
      ggforce::geom_sina(size = 0.2) +
      labs(x = "gene id", y = "gene expression", color = "gene expression") +
      stat_summary(
        aes(group = !!optional_plotting_column_sym),
        color = "red",
        # median and quartiles for point range
        fun = "mean",
        fun.min = function(x) {
          quantile(x, 0.25)
        },
        fun.max = function(x) {
          quantile(x, 0.75)
        },
        geom = "pointrange",
        position = position_dodge(width = 0.7),
        size = 0.2,
        shape = 21
      ) +
      geom_hline(yintercept = avg_gene_exp, linetype = "dashed")
    
  } else {
    
    sina_expression_plot <- ggplot(expression_means_df, aes(
      x = !!sce_rownames_column_sym,
      y = gene_expression,
      color = gene_expression)) +
      ggforce::geom_sina(size = 0.2) +
      labs(x = gsub('_', ' ', sce_rownames_column_sym), y = "gene expression", color = "gene expression") +
      stat_summary(
        aes(group = !!sce_rownames_column_sym),
        color = "red",
        # median and quartiles for point range
        fun = "mean",
        fun.min = function(x) {
          quantile(x, 0.25)
        },
        fun.max = function(x) {
          quantile(x, 0.75)
        },
        geom = "pointrange",
        position = position_dodge(width = 0.7),
        size = 0.2,
        shape = 21
      )
  }
  
  return(sina_expression_plot)
}

plot_goi_expression_umap <- function(normalized_sce,
                                     goi_list,
                                     sce_rownames_column,
                                     optional_plotting_column = NULL,
                                     use_rownames = TRUE){
  # Given a normalized SingleCellExperiment object and a vector of marker genes,
  # plot the gene expression on the UMAP results from the SingleCellExperiment
  # object.
  #
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object containing UMAP
  #                   results
  #   goi_list: data frame with marker genes relevant to the data in the
  #                 SingleCellExperiment object
  #   sce_rownames_column: name of the column with the SCE rowname identifiers
  #                              used for mapping
  #   optional_plotting_column: name of the column with the gene symbols that would be
  #                              used for plotting if provided
  #   use_rownames: indicates whether or not the SCE rowname identifiers should
  #                 be used when plotting; default is TRUE
  #
  # Returns:
  #   umap_plot: UMAP plots colored by gene expression for each of the
  #              individual marker gene symbols
  
  # Run the `prepare_expression_df` function
  expression_means_df <- prepare_expression_df(normalized_sce,
                                               goi_list,
                                               sce_rownames_column)
  
  # Get the UMAP matrix and then join the UMAP results with the expression data using the cell barcodes
  expression_umap_df <- data.frame(reducedDim(normalized_sce, "UMAP")) %>%
    tibble::rownames_to_column("cell_barcode") %>%
    dplyr::left_join(expression_means_df, by = "cell_barcode")
  
  # Check for optional plotting identifiers if `use_rownames` is FALSE
  if(!use_rownames && is.null(optional_plotting_column)){
    stop("When `use_rownames = FALSE`, a column name that contains optional 
    plotting identifiers should be provided. Please provide this column name or
    implement `use_rownames = TRUE`.")
  }
  
  # Plot UMAP, color by marker gene expression
  if(!is.null(optional_plotting_column)) {
    # Turn the gene symbol column name into a symbol for use when filtering
    optional_plotting_column_sym <- rlang::sym(optional_plotting_column)
    expression_umap_df <- expression_umap_df %>%
      dplyr::filter(!is.na(!!optional_plotting_column_sym))
    
    umap_plot <- ggplot(expression_umap_df,
                        aes(x = X1, y = X2, color = gene_expression)) +
      geom_point(size = 0.01) +
      facet_wrap(as.formula(paste("~", optional_plotting_column))) +
      scale_color_viridis_c() +
      labs(x = "UMAP1", y = "UMAP2", color = "gene expression")
  } else {
    # Turn the SCE rownames column into a symbol for use when filtering
    sce_rownames_column_sym <- rlang::sym(sce_rownames_column)
    expression_umap_df <- expression_umap_df %>%
      dplyr::filter(!is.na(!!sce_rownames_column_sym))
    
    umap_plot <- ggplot(expression_umap_df,
                        aes(x = X1, y = X2, color = gene_expression)) +
      geom_point(size = 0.01) +
      facet_wrap(as.formula(paste("~", sce_rownames_column))) +
      scale_color_viridis_c() +
      labs(x = "UMAP1", y = "UMAP2", color = "gene expression")
  }
  
  return(umap_plot)
}

plot_goi_expression_pca <- function(normalized_sce,
                                    goi_list,
                                    sce_rownames_column,
                                    optional_plotting_column = NULL,
                                    use_rownames = TRUE){
  
  # Given a normalized SingleCellExperiment object and a vector of marker genes,
  # plot the gene expression on the PCA results from the SingleCellExperiment
  # object.
  #
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object containing PCA
  #                   results
  #   goi_list: data frame with marker genes relevant to the data in the
  #                 SingleCellExperiment object
  #   sce_rownames_column: name of the column with the SCE rowname identifiers
  #                              used for mapping
  #   optional_plotting_column: name of the column with the gene symbols that would be
  #                              used for plotting if provided
  #   use_rownames: indicates whether or not the SCE rowname identifiers should
  #                 be used when plotting; default is TRUE
  #
  # Returns:
  #   pca_plot: PCA plots colored by gene expression for each of the
  #              individual marker gene symbols
  
  # Run the `prepare_expression_df` function
  expression_means_df <- prepare_expression_df(normalized_sce,
                                               goi_list,
                                               sce_rownames_column)
  
  # Get the PCA matrix and join the PCA results with the expression data using the cell barcodes
  expression_pca_df <- data.frame(reducedDim(normalized_sce, "PCA")) %>%
    tibble::rownames_to_column("cell_barcode") %>%
    dplyr::left_join(expression_means_df, by = "cell_barcode")
  
  # Check for optional plotting identifiers if `use_rownames` is FALSE
  if(use_rownames == FALSE && is.null(optional_plotting_column)){
    stop("When `use_rownames` = FALSE, a column name that contains optional 
    plotting identifiers should be provided. Please provide this column name or
    implement `use_rownames` = TRUE.")
  }
  
  # Plot first two principal components, color by the expression values
  # associated with the provided genes of interest
  if (!is.null(optional_plotting_column)) {
    # Turn the gene symbol column name into a symbol for use when filtering
    optional_plotting_column_sym <- rlang::sym(optional_plotting_column)
    expression_pca_df <- expression_pca_df %>%
      dplyr::filter(!is.na(!!optional_plotting_column_sym))
    
    pca_plot <- ggplot(expression_pca_df,
                       aes(x = PC1, y = PC2, color = gene_expression)) +
      geom_point(size = 0.3) +
      facet_wrap(as.formula(paste("~", optional_plotting_column))) +
      scale_color_viridis_c() +
      labs(color = "gene expression")
  } else {
    # Turn the SCE rownames column into a symbol for use when filtering
    sce_rownames_column_sym <- rlang::sym(sce_rownames_column)
    expression_pca_df <- expression_pca_df %>%
      dplyr::filter(!is.na(!!sce_rownames_column_sym))
    
    pca_plot <- ggplot(expression_pca_df,
                       aes(x = PC1, y = PC2, color = gene_expression)) +
      geom_point(size = 0.3) +
      facet_wrap(as.formula(paste("~", sce_rownames_column))) +
      scale_color_viridis_c() +
      labs(color = "gene expression")
  }
  
  return(pca_plot)
}
