## Custom functions to be sourced in the marker genes analysis reports template
## notebook.

# #### USAGE
# This script is intended to be sourced in the script as follows:
#
# source(file.path("utils", "marker-genes-analysis-functions.R"))

prepare_matrix_colnames <- function(normalized_sce_matrix,
                                    marker_genes,
                                    ensembl_id_column,
                                    gene_symbol_column) {
  # Given a normalized SingleCellExperiment matrix, the name of the column with
  # the Ensembl gene identifiers, and optionally (at command line) the name of
  # the associated gene symbols, prepare the matrix column names for heatmap
  # plotting.
  #
  # Args:
  #   normalized_sce_matrix: matrix retrieved from a normalized
  #                          SingleCellExperiment object
  #   marker_genes: data frame with marker genes relevant to the data in the
  #                 SingleCellExperiment object
  #   ensembl_id_column: name of the column with the ensembl gene identifiers
  #   gene_symbol_column: name of the column with the gene symbols that would be
  #                       used for plotting if provided at the command line
  
  # turn the ensembl column names into a symbol for use when subsetting
  ensembl_id_column_sym <- rlang::sym(ensembl_id_column)
  
  if (!is.null(gene_symbol_column)) {
    gene_symbol_column_sym <- rlang::sym(gene_symbol_column)
    
    # we will want to replace the column names our matrix with the relevant gene
    # symbols for plotting
    map_ensembl_symbols <-
      data.frame(ensembl = colnames(normalized_sce_matrix)) %>%
      dplyr::left_join(marker_genes, by = c("ensembl" = ensembl_id_column))
    
    colnames(normalized_sce_matrix) <-
      map_ensembl_symbols[[gene_symbol_column_sym]]
  } else {
    return(normalized_sce_matrix)
  }
  
  return(normalized_sce_matrix)
  
}

prepare_heatmap_annotation <- function(normalized_sce_matrix,
                                       marker_genes,
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
  #   marker_genes: data frame with marker genes relevant to the data in the
  #                 SingleCellExperiment object
  #   gene_id_column: name of the column with the gene identifiers (can be
  #                   ensembl or gene symbols if provided)
  #   gene_set: name of the column with the gene set information
  
  if (!is.null(gene_set)) {
    gene_id_column_sym <- rlang::sym(gene_id_column)
    gene_set_sym <- rlang::sym(gene_set)
    
    # create column annotation data frame
    annotation_df <- marker_genes %>%
      dplyr::select(gene_id_column_sym, gene_set_sym) %>%
      dplyr::distinct()
    
    rownames(annotation_df) <- NULL
    
    # filter to ensure only those genes of that marker genes that had data stored
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
      gene_set_colors <- as.list(gene_set_colors)
      
      # create the column annotation for the ComplexHeatmap
      column_annotation <- HeatmapAnnotation(
        df = annotation_df,
        col = gene_set_colors,
        annotation_label = c("Gene Set")
      )
    }
  } else if (is.null(gene_set)) {
    column_annotation <- NULL
  }
}

prepare_expression_df <- function(normalized_sce,
                                  marker_genes,
                                  ensembl_id_column) {
  # Given a normalized SingleCellExperiment object and a vector of marker genes,
  # prepare a gene expression data frame for plotting.
  #
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object
  #   marker_genes: data frame with marker genes relevant to the data in the
  #                 SingleCellExperiment object
  #   ensembl_id_column: name of the column with the ensembl gene identifiers
  #                              used for mapping
  
  # Turn the ensembl and gene symbol column names into symbols for use when
  # subsetting
  ensembl_id_column_sym <- rlang::sym(ensembl_id_column)
  
  # Prepare a data frame containing the logcounts expression values associated
  # with each marker gene symbol
  expression_means_df <- logcounts(normalized_sce[rownames(normalized_sce) %in% marker_genes[[ensembl_id_column_sym]],]) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell_barcode") %>%
    # add column containing average gene expression of each cell for all genes
    dplyr::mutate(A_all_mean_exp = colMeans(logcounts(normalized_sce))) %>%
    tidyr::pivot_longer(
      cols = -c("cell_barcode"),
      names_to = ensembl_id_column,
      values_to = "gene_expression"
    ) %>%
    dplyr::left_join(marker_genes, by = ensembl_id_column)
  
  return(expression_means_df)
  
}

plot_markers_expression_sina <- function(normalized_sce,
                                         marker_genes,
                                         ensembl_id_column,
                                         gene_symbol_column = NULL) {
  
  # Given a normalized SingleCellExperiment object and a vector of marker genes,
  # plot the gene expression of each marker gene on a sina plot against the
  # average mean expression of all the marker genes.
  #
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object
  #   marker_genes: data frame with marker genes relevant to the data in the
  #                 SingleCellExperiment object
  #   ensembl_id_column: name of the column with the ensembl gene identifiers
  #                              used for mapping
  #   gene_symbol_column: name of the column with the gene symbols that would be
  #                              used for plotting if provided
  #
  # Returns:
  #   sina_expression_plot: sina plots displaying the logcounts expression
  #                         values for each marker gene symbol relative to the
  #                         mean logcounts expression of all the marker genes
  
  # Turn the ensembl and gene symbol column names into symbols for use when
  # subsetting
  ensembl_id_column_sym <- rlang::sym(ensembl_id_column)
  
  # Run the `prepare_expression_df` function
  expression_means_df <- prepare_expression_df(normalized_sce,
                                               marker_genes,
                                               ensembl_id_column)
  
  if(!is.null(gene_symbol_column)) {
    gene_symbol_column_sym <- rlang::sym(gene_symbol_column)
    
    expression_means_df <- expression_means_df %>%
      dplyr::select(-cell_barcode) %>%
      dplyr::distinct() %>%
      # ensure that the symbol column contains the all_mean_expressed name
      # rather than NA
      dplyr::mutate(
        !!gene_symbol_column := ifelse(
          !!ensembl_id_column_sym == "A_all_mean_exp",
          "A_all_mean_exp",
          !!gene_symbol_column_sym
        )
      )
    
    sina_expression_plot <- ggplot(expression_means_df, aes(
      x = !!gene_symbol_column_sym,
      y = gene_expression,
      color = gene_expression)) +
      ggforce::geom_sina() +
      theme(axis.text.x = element_text(angle = 90))
    
  } else {
    
    sina_expression_plot <- ggplot(expression_means_df, aes(
      x = !!ensembl_id_column_sym,
      y = gene_expression,
      color = gene_expression)) +
      ggforce::geom_sina() +
      theme(axis.text.x = element_text(angle = 90))
  }
  
  return(sina_expression_plot)
}

plot_markers_expression_umap <- function(normalized_sce,
                                         marker_genes,
                                         ensembl_id_column,
                                         gene_symbol_column = NULL){
  # Given a normalized SingleCellExperiment object and a vector of marker genes,
  # plot the gene expression on the UMAP results from the SingleCellExperiment
  # object.
  #
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object containing UMAP
  #                   results
  #   marker_genes: data frame with marker genes relevant to the data in the
  #                 SingleCellExperiment object
  #   ensembl_id_column: name of the column with the ensembl gene identifiers
  #                              used for mapping
  #   gene_symbol_column: name of the column with the gene symbols that would be
  #                              used for plotting if provided
  #
  # Returns:
  #   umap_plot: UMAP plots colored by gene expression for each of the
  #              individual marker gene symbols
  
  # Run the `prepare_expression_df` function
  expression_means_df <- prepare_expression_df(normalized_sce,
                                               marker_genes,
                                               ensembl_id_column)
  
  # Get the UMAP matrix and then join the UMAP results with the expression data using the cell barcodes
  expression_umap_df <- data.frame(reducedDim(normalized_sce, "UMAP")) %>%
    tibble::rownames_to_column("cell_barcode") %>%
    dplyr::left_join(expression_means_df, by = "cell_barcode")
  
  
  # Plot UMAP, color by marker gene expression
  if(!is.null(gene_symbol_column)) {
    umap_plot <- ggplot(expression_umap_df,
                        aes(x = X1, y = X2, color = gene_expression)) +
      geom_point(size = 0.5) +
      facet_wrap(as.formula(paste("~", gene_symbol_column))) +
      scale_color_viridis_c()
  } else {
    umap_plot <- ggplot(expression_umap_df,
                        aes(x = X1, y = X2, color = gene_expression)) +
      geom_point(size = 0.5) +
      facet_wrap(as.formula(paste("~", ensembl_id_column))) +
      scale_color_viridis_c()
  }
  
  return(umap_plot)
}

plot_markers_expression_pca <- function(normalized_sce,
                                        marker_genes,
                                        ensembl_id_column,
                                        gene_symbol_column = NULL){
  
  # Given a normalized SingleCellExperiment object and a vector of marker genes,
  # plot the gene expression on the PCA results from the SingleCellExperiment
  # object.
  #
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object containing PCA
  #                   results
  #   marker_genes: data frame with marker genes relevant to the data in the
  #                 SingleCellExperiment object
  #   ensembl_id_column: name of the column with the ensembl gene identifiers
  #                              used for mapping
  #   gene_symbol_column: name of the column with the gene symbols that would be
  #                              used for plotting if provided
  #
  # Returns:
  #   pca_plot: PCA plots colored by gene expression for each of the
  #              individual marker gene symbols
  
  # Run the `prepare_expression_df` function
  expression_means_df <- prepare_expression_df(normalized_sce,
                                               marker_genes,
                                               ensembl_id_column)
  
  # Get the PCA matrix and join the PCA results with the expression data using the cell barcodes
  expression_pca_df <- data.frame(reducedDim(normalized_sce, "PCA")) %>%
    tibble::rownames_to_column("cell_barcode") %>%
    dplyr::left_join(expression_means_df, by = "cell_barcode")
  
  # Plot first two principal components, color by marker gene expression
  if (!is.null(gene_symbol_column)) {
    pca_plot <- ggplot(expression_pca_df,
                       aes(x = PC1, y = PC2, color = gene_expression)) +
      geom_point(size = 0.5) +
      facet_wrap(as.formula(paste("~", gene_symbol_column))) +
      scale_color_viridis_c()
  } else {
    pca_plot <- ggplot(expression_pca_df,
                       aes(x = PC1, y = PC2, color = gene_expression)) +
      geom_point(size = 0.5) +
      facet_wrap(as.formula(paste("~", ensembl_id_column))) +
      scale_color_viridis_c()
  }
  
  return(pca_plot)
}
