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
    cluster_name <- paste("kmeans", k, sep = "_")
    
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
    cluster_name <- paste(cluster_function, nearest_neighbors, sep = "_")
    
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

get_cluster_stats <- function(clustered_sce, cluster_column_name) {
  # Purpose: Check the validity stats of the clusters in the clustered
  # SingleCellExperiment object
  
  # Args:
  #   clustered_sce: SingleCellExperiment object with clustered results
  #   cluster_column_name: name of the column with the associated cluster names
  
  # isolate the clusters stored in the `cluster_name` slot of the SCE object
  clusters <- clustered_sce[[(cluster_column_name)]]
  
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

create_metadata_stats_df <- function(clustered_sce, params_range, increments, cluster_type) {
  # Purpose: Calculate and return a data frame with the validity stats of the
  # clusters in the SingleCellExperiment object
  
  # Args:
  #   clustered_sce: SingleCellExperiment object with clustered results
  #   params_range: the range of numeric parameters to test for clustering
  #   increments: a numeric value representing the increments by which to explore
  #               the params range of values
  #   cluster_type: the type of clustering method performed - can be "kmeans or graph"
  
  # define cluster names
  param_values <- seq(min(params_range), max(params_range),increments)
  cluster_names_column <- paste(cluster_type, param_values, sep = "_")
  
  # save data.frame to the cluster validity list of data.frames
  cluster_validity_df_list <- cluster_names_column %>% 
    purrr::map(~ get_cluster_stats(clustered_sce = clustered_sce, .x))
  names(cluster_validity_df_list) <- cluster_names_column
  
  # now bind the rows of all the cluster validity data.frames in the list
  cluster_validity_df <- dplyr::bind_rows(cluster_validity_df_list,
                                          .id = "cluster_names_column") %>%
    tidyr::separate(cluster_names_column, 
                    c("cluster_type", "param_value"),
                    remove = FALSE) # keep original column for grouping later
    
  
  return(cluster_validity_df)
}


summarize_clustering_stats <- function(cluster_validity_df) {
  # Purpose: Calculate and return a summary data frame of the provided cluster
  # validity stats
  
  # Args:
  #   clustered_validity_df: data.frame with cluster validity stats associated
  #                          with their relevant cluster names
  
  # create a summary data.frame of the results across the individual clusters
  validity_summary_df <- cluster_validity_df %>%
    dplyr::group_by(cluster_names_column) %>%
    dplyr::summarize(avg_purity = median(purity),
                     avg_maximum = median(as.numeric(maximum)),
                     avg_width = median(width),
                     avg_closest = median(as.numeric(closest))) %>%
    dplyr::select(cluster_names_column, avg_purity, avg_maximum, avg_width, avg_closest)
  
  return(validity_summary_df)
}

plot_cluster_purity <- function(cluster_validity_df) {
  
  # Purpose: Generate a plot displaying the cluster purity stats of the clusters
  # in the SingleCellExperiment object
  
  # Args:
  #   clustered_validity_df: data.frame with cluster validity stats associated
  #                          with their relevant cluster names
  
  # prepare data frame for plotting
  metadata <- cluster_validity_df %>%
    # create a column for the color scale to make things easier to control the color later
    # for cluster purity, do the majority of neighboring cells come from the assigned cluster (yes) or a different cluster (no)
    # for silhouette width if the cluster matches, then the silhouette width is positive, if not it's negative
    dplyr::mutate(
      color_scale = ifelse(maximum == cluster,  "yes",  "no"),
      param_value = as.numeric(param_value),
      cluster_param_assignment = paste(cluster_type, param_value, cluster, sep = "_")
    )
  
  # set colors and title for plotting
  colors = c("gray", "red")
  names(colors) = levels(metadata$color_scale)
  legend_title = "Neighboring cells \nbelong to assigned cluster"
  
  # plot the cluster validity data frames
  plot <-
    ggplot(metadata, aes(x = cluster, y = purity, colour = color_scale)) +
    ggbeeswarm::geom_quasirandom(method = "smiley", size = 0.2) +
    scale_color_manual(values = c("yes" = "gray",
                                  "no" = "red")) +
    stat_summary(
      aes(group = cluster_param_assignment),
      color = "black",
      # median and quartiles for point range
      fun = "median",
      fun.min = function(x) {
        quantile(x, 0.25)
      },
      fun.max = function(x) {
        quantile(x, 0.75)
      },
      geom = "pointrange",
      position = position_dodge(width = 0.9),
      size = 0.2
    ) +
    theme(text = element_text(size = 18)) +
    labs(x = paste0(unique(metadata$cluster_type), "Parameters"),
         color = legend_title) +
    facet_wrap( ~ param_value, scale="free") + 
    theme_bw()
  
  return(plot)
}

plot_cluster_silhouette_width <- function(cluster_validity_df) {
  # Purpose: Calculate and return a data frame with the validity stats of the
  # clusters in the SingleCellExperiment object
  
  # Args:
  #   clustered_validity_df: data.frame with cluster validity stats associated
  #                          with their relevant cluster names
  
  # prepare data frame for plotting
  metadata <- cluster_validity_df %>%
    # create a column for the color scale to make things easier to control the color later
    # for cluster purity, do the majority of neighboring cells come from the assigned cluster (yes) or a different cluster (no)
    # for silhouette width if the cluster matches, then the silhouette width is positive, if not it's negative
    dplyr::mutate(
      param_value = as.numeric(param_value),
      cluster_param_assignment = paste(cluster_type, param_value, cluster, sep = "_")
    )
  
  # set title for plotting
  legend_title = "Positive Silhouette Width"
  
  # plot the cluster validity data frames
  plot <-
    ggplot(metadata, aes(x = cluster, y = width)) +
    ggbeeswarm::geom_quasirandom(method = "pseudorandom", size = 0.2) +
    geom_hline(yintercept = 0) +
    theme(text = element_text(size = 18)) +
    labs(x = paste0(unique(metadata$cluster_type), "Parameters"),
         color = legend_title) +
    facet_wrap( ~ param_value, scale="free_x") +
    stat_summary(
      aes(group = cluster_param_assignment),
      color = "red",
      # median and quartiles for point range
      fun = "median",
      fun.min = function(x) {
        quantile(x, 0.25)
      },
      fun.max = function(x) {
        quantile(x, 0.75)
      },
      geom = "pointrange",
      position = position_dodge(width = 0.9),
      size = 0.5
    ) +
    theme_bw()
  
  return(plot)
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
