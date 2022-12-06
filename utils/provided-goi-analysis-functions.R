## Custom functions to be sourced in the provided goi analysis reports template
## notebook.

# #### USAGE
# This script is intended to be sourced in the script as follows:
#
# source(file.path("utils", "provided-goi-analysis-functions.R"))

colnames_to_gene_symbols <- function(normalized_sce_matrix,
                                    goi_list,
                                    ensembl_id_column,
                                    gene_symbol_column) {
  # Given a normalized SingleCellExperiment matrix, the name of the column with
  # the Ensembl gene identifiers, the name of the associated gene symbols,
  # convert the matrix column names into symbols for heatmap plotting.
  #
  # Args:
  #   normalized_sce_matrix: matrix retrieved from a normalized
  #                          SingleCellExperiment object
  #   goi_list: data frame with provided genes of interest relevant to the data
  #             in the SingleCellExperiment object
  #   ensembl_id_column: name of the column with the ensembl gene identifiers
  #   gene_symbol_column: name of the column with the gene symbols that would be
  #                       used for plotting if provided at the command line
  
  # turn the gene symbol column name into a symbol for use when subsetting
  gene_symbol_column_sym <- rlang::sym(gene_symbol_column)
  
  # we will want to replace the column names our matrix with the relevant gene
  # symbols for plotting
  map_ensembl_symbols <-
    data.frame(ensembl = colnames(normalized_sce_matrix)) %>%
    dplyr::left_join(goi_list, by = c("ensembl" = ensembl_id_column))
  
  colnames(normalized_sce_matrix) <-
    map_ensembl_symbols[[gene_symbol_column_sym]]
  
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
                                  ensembl_id_column) {
  # Given a normalized SingleCellExperiment object and a goi list,
  # prepare a gene expression data frame for plotting.
  #
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object
  #   goi_list: data frame with genes of interest relevant to the data in the
  #                 SingleCellExperiment object
  #   ensembl_id_column: name of the column with the ensembl gene identifiers
  #                              used for mapping
  
  # Turn the ensembl and gene symbol column names into symbols for use when
  # subsetting
  ensembl_id_column_sym <- rlang::sym(ensembl_id_column)
  
  # Prepare a data frame containing the logcounts expression values associated
  # with each marker gene symbol
  expression_means_df <- logcounts(normalized_sce[rownames(normalized_sce) %in% goi_list[[ensembl_id_column_sym]],]) %>%
    t() %>%
    as.matrix() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell_barcode") %>%
    tidyr::pivot_longer(
      cols = -c("cell_barcode"),
      names_to = ensembl_id_column,
      values_to = "gene_expression"
    ) %>%
    dplyr::left_join(goi_list, by = ensembl_id_column)
  
  return(expression_means_df)
  
}

plot_goi_expression_sina <- function(normalized_sce,
                                     goi_list,
                                     ensembl_id_column,
                                     gene_symbol_column = NULL) {
  
  # Given a normalized SingleCellExperiment object and a goi list,
  # plot the gene expression of each gene of interest on a sina plot against the
  # average mean expression of all the genes of interest.
  #
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object
  #   goi_list: data frame with genes of interest relevant to the data in the
  #                 SingleCellExperiment object
  #   ensembl_id_column: name of the column with the ensembl gene identifiers
  #                              used for mapping
  #   gene_symbol_column: name of the column with the gene symbols that would be
  #                              used for plotting if provided
  #
  # Returns:
  #   sina_expression_plot: sina plots displaying the logcounts expression
  #                         values for each gene of interest symbol relative to
  #                         the mean logcounts expression of all the genes of
  #                         interest
  
  # Turn the ensembl and gene symbol column names into symbols for use when
  # subsetting
  ensembl_id_column_sym <- rlang::sym(ensembl_id_column)
  
  # Run the `prepare_expression_df` function
  expression_means_df <- prepare_expression_df(normalized_sce,
                                               goi_list,
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
    
    # calculate average gene expression
    avg_gene_exp = mean(colMeans(logcounts(normalized_sce)))
    
    sina_expression_plot <- ggplot(expression_means_df, aes(
      x = !!gene_symbol_column_sym,
      y = gene_expression,
      color = gene_expression)) +
      ggforce::geom_sina(size = 0.2) +
      theme(axis.text.x = element_text(angle = 90)) +
      stat_summary(
        aes(group = !!gene_symbol_column_sym),
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
      geom_hline(yintercept = avg_gene_exp, type = "dashed")
    
  } else {
    
    sina_expression_plot <- ggplot(expression_means_df, aes(
      x = !!ensembl_id_column_sym,
      y = gene_expression,
      color = gene_expression)) +
      ggforce::geom_sina(size = 0.2) +
      theme(axis.text.x = element_text(angle = 90)) +
      stat_summary(
        aes(group = !!ensembl_id_column_sym),
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
                                         ensembl_id_column,
                                         gene_symbol_column = NULL){
  # Given a normalized SingleCellExperiment object and a vector of marker genes,
  # plot the gene expression on the UMAP results from the SingleCellExperiment
  # object.
  #
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object containing UMAP
  #                   results
  #   goi_list: data frame with marker genes relevant to the data in the
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
                                               goi_list,
                                               ensembl_id_column)
  
  # Get the UMAP matrix and then join the UMAP results with the expression data using the cell barcodes
  expression_umap_df <- data.frame(reducedDim(normalized_sce, "UMAP")) %>%
    tibble::rownames_to_column("cell_barcode") %>%
    dplyr::left_join(expression_means_df, by = "cell_barcode")
  
  
  # Plot UMAP, color by marker gene expression
  if(!is.null(gene_symbol_column)) {
    # Turn the gene symbol column name into a symbol for use when filtering
    gene_symbol_column_sym <- rlang::sym(gene_symbol_column)
    expression_umap_df <- expression_umap_df %>%
      dplyr::filter(!is.na(!!gene_symbol_column_sym))
    
    umap_plot <- ggplot(expression_umap_df,
                        aes(x = X1, y = X2, color = gene_expression)) +
      geom_point(size = 0.01) +
      facet_wrap(as.formula(paste("~", gene_symbol_column))) +
      scale_color_viridis_c()
  } else {
    # Turn the ensembl column name into a symbol for use when filtering
    ensembl_id_column_sym <- rlang::sym(ensembl_id_column)
    expression_umap_df <- expression_umap_df %>%
      dplyr::filter(!is.na(!!ensembl_id_column_sym))
    
    umap_plot <- ggplot(expression_umap_df,
                        aes(x = X1, y = X2, color = gene_expression)) +
      geom_point(size = 0.01) +
      facet_wrap(as.formula(paste("~", ensembl_id_column))) +
      scale_color_viridis_c()
  }
  
  return(umap_plot)
}

plot_goi_expression_pca <- function(normalized_sce,
                                        goi_list,
                                        ensembl_id_column,
                                        gene_symbol_column = NULL){
  
  # Given a normalized SingleCellExperiment object and a vector of marker genes,
  # plot the gene expression on the PCA results from the SingleCellExperiment
  # object.
  #
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object containing PCA
  #                   results
  #   goi_list: data frame with marker genes relevant to the data in the
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
                                               goi_list,
                                               ensembl_id_column)
  
  # Get the PCA matrix and join the PCA results with the expression data using the cell barcodes
  expression_pca_df <- data.frame(reducedDim(normalized_sce, "PCA")) %>%
    tibble::rownames_to_column("cell_barcode") %>%
    dplyr::left_join(expression_means_df, by = "cell_barcode")
  
  # Plot first two principal components, color by the expression values
  # associated with the provided genes of interest
  if (!is.null(gene_symbol_column)) {
    # Turn the gene symbol column name into a symbol for use when filtering
    gene_symbol_column_sym <- rlang::sym(gene_symbol_column)
    expression_pca_df <- expression_pca_df %>%
      dplyr::filter(!is.na(!!gene_symbol_column_sym))
    
    pca_plot <- ggplot(expression_pca_df,
                       aes(x = PC1, y = PC2, color = gene_expression)) +
      geom_point(size = 0.3) +
      facet_wrap(as.formula(paste("~", gene_symbol_column))) +
      scale_color_viridis_c()
  } else {
    # Turn the gene symbol column name into a symbol for use when filtering
    ensembl_id_column_sym <- rlang::sym(ensembl_id_column)
    expression_pca_df <- expression_pca_df %>%
      dplyr::filter(!is.na(!!ensembl_id_column_sym))
    
    pca_plot <- ggplot(expression_pca_df,
                       aes(x = PC1, y = PC2, color = gene_expression)) +
      geom_point(size = 0.3) +
      facet_wrap(as.formula(paste("~", ensembl_id_column))) +
      scale_color_viridis_c()
  }
  
  return(pca_plot)
}
