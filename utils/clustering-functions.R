## Custom functions to be sourced in the clustering R notebook.

# #### USAGE
# This script is intended to be sourced in the script as follows:
#
# source(file.path("utils", "clustering-functions.R"))

kmeans_clustering <- function(normalized_sce, k, cluster_name,
                              check_stability = FALSE, seed = 2021) {
  # Purpose: Perform the k-means clustering on a normalized SingleCellExperiment
  # object
  
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object
  #   k: the desired number of centers
  #   cluster_name: name for storing the cluster results in the
  #                 SingleCellExperiment object
  #   check_stability: if 'TRUE', the stability of the clusters will be
  #                    calculated; the default is 'FALSE'
  #   seed: an integer to set the seed as for reproducibility
  
  # set the seed for reproducible results
  set.seed(seed)
  
  # extract the principal components matrix
  pca_matrix <- reducedDim(normalized_sce, "PCA")
  
  # perform k-means clustering
  clusters <- clusterRows(pca_matrix, KmeansParam(centers = k))
  
  # store cluster results in the SCE object
  normalized_sce[[cluster_name]] <- factor(clusters)
  
  # check cluster stability if the `check_stability == TRUE`
  if (check_stability) {
    # create name for the bootstrapping results to be stored in
    bootstrapping_name <-
      paste("bootstrapping_results", cluster_name, sep = "_")
    
    metadata(normalized_sce)[[bootstrapping_name]] <-
      check_cluster_stability(normalized_sce,
                              pca_matrix,
                              bootstrapping_name,
                              cluster_name,
                              k)
  }
  
  return(normalized_sce)
  
}

graph_clustering <- function(normalized_sce, k, cluster_name, type,
                             cluster_function, check_stability = FALSE, seed = 2021) {
  # Purpose: Perform the k-means clustering on a normalized SingleCellExperiment
  # object
  
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object
  #   k: the desired number of centers
  #   cluster_name: name for storing the cluster results in the
  #                 SingleCellExperiment object
  #   type: the type of weighting scheme -- can by "jaccard" or "rank"
  #   cluster_function: the name of the community detection algorithm that is
  #                     being tested -- can be "walktrap" or "louvain"
  #   check_stability: if 'TRUE', the stability of the clusters will be
  #                    calculated; the default is 'FALSE'
  #   seed: an integer to set the seed as for reproducibility
  
  
  # set the seed for reproducible results
  set.seed(seed)
  
  # turn arguments into symbols
  type_sym <- rlang::sym(type)
  cluster_function_sym <- rlang::sym(cluster_function)
  
  # extract the principal components matrix
  pca_matrix <- reducedDim(normalized_sce, "PCA")
  
  # perform graph-based clustering
  clusters <- clusterRows(pca_matrix,
                          NNGraphParam(
                            k = k,
                            type = type,
                            cluster.fun = cluster_function
                          ))
  
  # store cluster results in the SCE object
  normalized_sce[[cluster_name]] <- factor(clusters)
  
  # check cluster stability if the `check_stability == TRUE`
  if (check_stability) {
    # create name for the bootstrapping results to be stored in
    bootstrapping_name <-
      paste("bootstrapping_results", cluster_name, sep = "_")
    
    metadata(normalized_sce)[[bootstrapping_name]] <-
      check_cluster_stability(normalized_sce,
                              pca_matrix,
                              bootstrapping_name,
                              cluster_name,
                              k)
  }
  
  return(normalized_sce)
  
}

check_cluster_stability <- function(normalized_sce, pca_matrix,
                                    bootstrapping_name, cluster_name, k) {
  # Purpose: To use bootstrapping to check the stability of given clusters
  
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object
  #   pca_matrix: the matrix object that contains the PCA results that were
  #               extracted from the noramlized SingleCellExperiment object
  #   bootstrapping_name: name for storing the bootstrapping results in the
  #                       SingleCellExperiment object
  #   cluster_name: name associated with where the cluster results in the
  #                 SingleCellExperiment object are stored
  #   k: the desired number of centers
  
  # check cluster stability if the `check_stability == TRUE`
  if (startsWith(cluster_name, "kcluster")) {
    # Create subfunctions for `bootstrapStability()`
    cluster_stability_subfunction <- function(x) {
      clusterRows(x, KmeansParam(centers = k))
    }
  } else if (startsWith(cluster_name, "walktrap_cluster")) {
    cluster_stability_subfunction <- function(x) {
      clusterRows(x, NNGraphParam(k = k, cluster.fun = "walktrap"))
    }
  } else if (startsWith(cluster_name, "louvain_cluster")) {
    cluster_stability_subfunction <- function(x) {
      clusterRows(x,
                  NNGraphParam(
                    k = k,
                    type = "jaccard",
                    cluster.fun = "louvain"
                  ))
    }
  }
  
  # run the `bootstapStability()` function
  metadata(normalized_sce)[[bootstrapping_name]] <-
    bootstrapStability(pca_matrix,
                       FUN = cluster_stability_subfunction,
                       clusters = normalized_sce[[cluster_name]])
  
}

cluster_validity_stats <- function(clustered_sce, cluster_name) {
  # Purpose: Check the cluster purity of the clusters in the clustered
  # SingleCellExperiment object
  
  # Args:
  #   clustered_sce: SingleCellExperiment object with clustered results
  #   cluster_name: name associated with where the cluster results in the
  #                 SingleCellExperiment object are stored
  
  # isolate the clusters stored in the `cluster_name` slot of the SCE object
  clusters <- clustered_sce[[(cluster_name)]]
  
  # use `neighborPurity` to check the purity of clusters
  purity_df <-
    neighborPurity(reducedDim(clustered_sce, "PCA"), clusters = clusters) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("cell_barcode")
  purity_df$maximum <- factor(purity_df$maximum)
  
  # use `approxSilhouette()` to approximate silhouettes
  sil_df <- approxSilhouette(reducedDim(clustered_sce, "PCA"),
                             clusters = clusters) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell_barcode")
  
  # account for negative widths
  sil_df$closest <-
    factor(ifelse(sil_df$width > 0, clusters, sil_df$other))
  
  cluster_stats_df <- purity_df %>%
    dplyr::left_join(sil_df) %>%
    dplyr::mutate(cluster = factor(clusters))
  
  return(cluster_stats_df)
}

plot_clustering_validity <- function(cluster_validity_df, measure, colour_var,
                                     facet_var) {
  # Purpose: Plot the provided clustering data frame
  
  # Args:
  #   clustered_validity_df: data frame with cluster validity stats
  #   measure: string associated with the column whose values should be on the 
  #            y-axis
  #   colour_var: string associated with the column whose values should be used
  #               color the points on the plot
  #   facet_var: string associated with the column whose values should be used
  #              to facet on, to create individual plots
  
  # convert into symbols for plotting
  measure <- rlang::sym(measure)
  colour_var <- rlang::sym(colour_var)
  facet_var <- rlang::sym(facet_var)
  
  # plot the cluster validity data frames
  ggplot(cluster_validity_df, aes(x = cluster, y = !!measure, colour = !!colour_var)) +
    ggbeeswarm::geom_quasirandom(method = "smiley") +
    facet_wrap(facet_var, ncol = 1)
}

plot_cluster_stability <- function(normalized_sce, cluster_name) {
  # Purpose: Plot the cluster bootstrapping stability values of the clusters
  # stored in the SingleCellExperiment object
  
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object
  #   cluster_name: name associated with where the cluster results in the
  #                 SingleCellExperiment object are stored
  
  # set the name from which the bootstrapping results can be retrieved from the SCE object
  bootstrapping_name <- paste("bootstrapping_results", cluster_name, sep = "_")
  
  # plot cluster stability
  stability_pheatmap <- pheatmap(
    metadata(normalized_sce)[[bootstrapping_name]],
    cluster_row = FALSE,
    cluster_col = FALSE,
    color = viridis::magma(100),
    breaks = seq(-1, 1, length.out = 101),
    main = cluster_name
  )
}
