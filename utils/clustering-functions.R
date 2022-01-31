## Custom functions to be sourced in the clustering R notebook.

# #### USAGE
# This script is intended to be sourced in the script as follows:
#
# source(file.path("utils", "clustering-functions.R"))

kmeans_clustering <- function(normalized_sce,
                              k_range,
                              check_stability = FALSE,
                              seed = 2021) {
  # Purpose: Perform the k-means clustering on a normalized SingleCellExperiment
  # object
  
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object
  #   k_range: the range of the desired number of centers
  #   check_stability: if 'TRUE', the stability of the clusters will be
  #                    calculated; the default is 'FALSE'
  #   seed: an integer to set the seed as for reproducibility
  
  # first check that the normalized object is a SingleCellExperiment object
  if(!is(normalized_sce,"SingleCellExperiment")){
    stop("normalized_sce must be a SingleCellExperiment object.")
  }
  
  # Perform k-means clustering
  for (k in k_range) {
    cluster_name <- paste("kcluster", k, sep = "")
    
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
      normalized_sce <- check_cluster_stability(normalized_sce,
                                                pca_matrix,
                                                cluster_type = "kmeans",
                                                cluster_name,
                                                k)
    }
    
    
  }
  
  return(normalized_sce)
  
}

graph_clustering <- function(normalized_sce,
                             nn_range,
                             weighting_type = "rank",
                             cluster_function = "walktrap",
                             check_stability = FALSE,
                             seed = 2021) {
  # Purpose: Perform the graph based clustering on a normalized SingleCellExperiment
  # object
  
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object
  #   nn_range: the range of the number of nearest neighbors to consider during
  #             graph construction
  #   weighting_type: the type of weighting scheme -- can be "rank", "number", or "jaccard"
  #   cluster_function: the name of the community detection algorithm that is
  #                     being tested -- can be "walktrap" or "louvain"
  #   check_stability: if 'TRUE', the stability of the clusters will be
  #                    calculated; the default is 'FALSE'
  #   seed: an integer to set the seed as for reproducibility
  
  # first check that the normalized object is a SingleCellExperiment object
  if(!is(normalized_sce,"SingleCellExperiment")){
    stop("normalized_sce must be a SingleCellExperiment object.")
  }
  
  # perform the graph-based clustering
  for (nearest_neighbors in nn_range) {
    # set cluster name
    cluster_name <- paste(cluster_function,"_cluster", nearest_neighbors, sep = "")
    
    # set the seed for reproducible results
    set.seed(seed)
    
    # extract the principal components matrix
    pca_matrix <- reducedDim(normalized_sce, "PCA")
    
    # perform graph-based clustering
    clusters <- clusterRows(
      pca_matrix,
      NNGraphParam(
        k = nearest_neighbors,
        type = weighting_type,
        cluster.fun = cluster_function
      )
    )
    
    # store cluster results in the SCE object
    normalized_sce[[cluster_name]] <- factor(clusters)
    
    # check cluster stability if the `check_stability == TRUE`
    if (check_stability) {
      normalized_sce <- check_cluster_stability(
        normalized_sce,
        pca_matrix,
        cluster_type = "graph",
        cluster_name,
        nearest_neighbors,
        weighting_type,
        cluster_function
      )
    }
    
  }
  
  return(normalized_sce)
  
}

check_cluster_stability <- function(normalized_sce,
                                    pca_matrix,
                                    cluster_type,
                                    cluster_name,
                                    k,
                                    weighting_type = "rank",
                                    cluster_function = "walktrap") {
  # Purpose: To use bootstrapping to check the stability of given clusters
  
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object
  #   pca_matrix: the matrix object that contains the PCA results that were
  #               extracted from the noramlized SingleCellExperiment object
  #   cluster_type: the type of clustering method performed - can be "kmeans or graph"
  #   cluster_name: name associated with where the cluster results in the
  #                 SingleCellExperiment object are stored
  #   k: the desired number of centers
  #   weighting_type: the type of weighting scheme -- can be "jaccard" or "rank";
  #                   "rank" being the default
  #   cluster_function: the name of the community detection algorithm that is
  #                     being tested -- can be "walktrap" or "louvain";
  #                     "walktrap" being the default
  
  # create name for the bootstrapping results to be stored in
  bootstrapping_name <- paste("bootstrapping_results", cluster_name, sep = "_")
  
  # check cluster stability if the `check_stability == TRUE`
  if (cluster_type == "kmeans") {
    # Create subfunctions for `bootstrapStability()`
    cluster_stability_subfunction <- function(x) {
      clusterRows(x, KmeansParam(centers = k))
    }
  } else if (cluster_type == "graph") {
    cluster_stability_subfunction <- function(x) {
      clusterRows(x,
                  NNGraphParam(
                    k = k,
                    type = weighting_type,
                    cluster.fun = cluster_function
                  ))
    }
  } else {
    stop("Clustering type is not valid. Please use --cluster_type to specify 'kmeans' or 'graph' clustering.")
  }
  
  # run the `bootstapStability()` function
  metadata(normalized_sce)[[bootstrapping_name]] <-
    bootstrapStability(pca_matrix,
                       FUN = cluster_stability_subfunction,
                       clusters = normalized_sce[[cluster_name]])
  
  return(normalized_sce)
  
}

add_metadata_clustering_stats <- function(clustered_sce, cluster_names, cluster_type) {
  # Purpose: Check the cluster purity of the clusters in the clustered
  # SingleCellExperiment object and store results within the object
  
  # Args:
  #   clustered_sce: SingleCellExperiment object with clustered results
  #   cluster_names: vector of names associated with where the cluster results in the
  #                 SingleCellExperiment object are stored
  #   cluster_type: the type of clustering method performed - can be "kmeans or graph"
  
  get_cluster_stats <- function(clustered_sce, cluster_names) {
    
    for (cluster_name in cluster_names) {
      # isolate the clusters stored in the `cluster_name` slot of the SCE object
      clusters <- clustered_sce[[(cluster_name)]]
    }
    
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
    
    # join purity and silhoette info in one data.frame
    cluster_stats_df <- purity_df %>%
      dplyr::left_join(sil_df) %>%
      dplyr::mutate(cluster = factor(clusters))
  }
  
  # save data.frame to the cluster validity list of data.frames
  cluster_validity_df_list <- cluster_names %>% 
    purrr::map(~ get_cluster_stats(clustered_sce = clustered_sce, .x))
  
  # now bind the rows of all the cluster validity data.frames in the list
  cluster_validity_df <- dplyr::bind_rows(cluster_validity_df_list,
                                          .id = "cluster_names")
  
  # create a summary data.frame of the results across the individual clusters
  validity_summary_df <- cluster_validity_df %>%
    dplyr::group_by(cluster_names) %>%
    dplyr::summarize(avg_purity = median(purity),
                     avg_maximum = median(as.numeric(maximum)),
                     avg_width = median(width),
                     avg_closest = median(as.numeric(closest))) %>%
    dplyr::select(cluster_names, avg_purity, avg_maximum, avg_width, avg_closest)
  
  metadata(clustered_sce)$all_stats[[cluster_type]] <- cluster_validity_df
  metadata(clustered_sce)$summary_stats[[cluster_type]]  <- validity_summary_df
  
  return(clustered_sce)
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
