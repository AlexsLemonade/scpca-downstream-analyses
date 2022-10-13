## Custom functions to be sourced when performing clustering.

# #### USAGE
# This script is intended to be sourced in the script as follows:
#
# source(file.path("utils", "clustering-functions.R"))

kmeans_clustering <- function(normalized_sce,
                              params_range,
                              step_size = NULL,
                              seed = 2021) {
  # Purpose: Perform the k-means clustering on a normalized SingleCellExperiment
  # object
  
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object
  #   params_range: the range of numeric parameters to test for clustering, 
  #     can be a range of values, e.g. c(1:10), or a single value
  #   step_size: if a range of values is provided to `params_range`, 
  #     a numeric value representing the step size to use within the params range of values
  #   seed: an integer to set the seed as for reproducibility
  
  # first check that the normalized object is a SingleCellExperiment object
  if(!is(normalized_sce,"SingleCellExperiment")){
    stop("`normalized_sce` must be a SingleCellExperiment object.")
  }
  
  # check that params_range is an integer
  if(!is.integer(params_range)){
    stop("`params_range` must be an integer.")
  }
  
  # if a step size exists then create a sequence of params for clustering 
  if(!is.null(step_size)){
    k_range <- seq(min(params_range), max(params_range), step_size) 
    # if no step size has been input then the params range is directly used for clustering 
  } else {
    k_range <- params_range
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
    
    
  }
  
  return(normalized_sce)
  
}

define_nn_range <- function(nearest_neighbors_min,
                            nearest_neighbors_max,
                            nearest_neighbors_increment = NULL){
  # Purpose: Define a nearest neighbors range sequence with the provided
  # range and increment values
  
  # Args:
  #   nearest_neighbors_min: minimum number of a range of nearest neighbors 
  #                          values to include when calculating/plotting the 
  #                          clustering results.
  #   nearest_neighbors_max: maximum number of a range of nearest neighbors 
  #                          values to include when calculating/plotting the 
  #                          clustering results.
  #   nearest_neighbors_increment: increment to use when implementing the range 
  #                                number of nearest neighbors for cluster stats.
  
  # If no nearest neighbors increment has been input then the provided range is 
  # directly used for clustering and nearest neighbors increment is set to 1
  if (is.null(nearest_neighbors_increment)) {
    nearest_neighbors_increment <- 1
  }
  nn_range <- seq(
    nearest_neighbors_min,
    nearest_neighbors_max,
    nearest_neighbors_increment
  )
  
  return(nn_range)
}

graph_clustering <- function(normalized_sce,
                             nearest_neighbors_min,
                             nearest_neighbors_max,
                             step_size = NULL,
                             cluster_type = "walktrap",
                             seed = 2021) {
  # Purpose: Perform the graph based clustering on a normalized SingleCellExperiment
  # object
  
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object
  #   nearest_neighbors_min: minimum number of a range of nearest neighbors 
  #                          values to include when calculating/plotting the 
  #                          clustering results.
  #   nearest_neighbors_max: maximum number of a range of nearest neighbors 
  #                          values to include when calculating/plotting the 
  #                          clustering results.
  #   step_size: if a range of values is provided to `params_range`, 
  #     a numeric value representing the step size to use within the params range of values
  #   cluster_type: the type of graph-based clustering method that is being 
  #     tested -- can be "walktrap" or "louvain"; the default is "walktrap"
  #   seed: an integer to set the seed as for reproducibility
  
  # first check that the normalized object is a SingleCellExperiment object
  if(!is(normalized_sce,"SingleCellExperiment")){
    stop("normalized_sce must be a SingleCellExperiment object.")
  }
  
  # check that nearest_neighbors_min and nearest_neighbors_max are integers
  if(!is.integer(nearest_neighbors_min)){
    stop("`nearest_neighbors_min` must be an integer.")
  }
  if(!is.integer(nearest_neighbors_max)){
    stop("`nearest_neighbors_max` must be an integer.")
  }
  
  # define nearest neighbors range
  nn_range <- define_nn_range(nearest_neighbors_min,
                              nearest_neighbors_max,
                              step_size)
  
  # determine weighting type to use based on graph detection algorithm specified 
  # if louvain is used, use jaccard 
  # if walktrap is used, use rank 
  weighting_type <- ifelse(cluster_type == "louvain", "jaccard", "rank")
  
  # perform the graph-based clustering
  for (nearest_neighbors in nn_range) {
    # set cluster name
    cluster_name <- sprintf("%s_%02d", cluster_type, nearest_neighbors)
    
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
        cluster.fun = cluster_type
      )
    )
    
    # store cluster results in the SCE object
    normalized_sce[[cluster_name]] <- factor(clusters)
    
  }
  
  return(normalized_sce)
  
}

check_cluster_stability <- function(pca_matrix,
                                    cluster_assignments,
                                    cluster_type,
                                    k,
                                    iterations = 20,
                                    seed = 2021) {
  # Purpose: To use bootstrapping to check the stability of given clusters
  
  # Args:
  #   pca_matrix: the matrix object that contains the PCA results that were
  #               extracted from a normalized SingleCellExperiment object
  #   cluster_assignments: original cluster assignments from the SingleCellExperiment
  #   cluster_type: the type of clustering method performed - can be "kmeans", 
  #                 "walktrap", or "louvain"
  #   k: the desired number of centers
  #   iterations: The number of iterations to perform for bootstrapping, default is 20.
  #   seed: an integer to set the seed as for reproducibility
  
  # set the seed for reproducible results
  set.seed(seed)
  
  # check cluster stability if the `check_stability == TRUE`
  if (cluster_type == "kmeans") {
    # Create subfunctions needed for bootstrapping
    cluster_stability_subfunction <- function(x) {
      clusterRows(x, KmeansParam(centers = k))
    }
  } else if (cluster_type %in% c("walktrap", "louvain")) {
    # grab weighting type based on cluster type
    weighting_type <- ifelse(cluster_type == "louvain", "jaccard", "rank")
    
    # create cluster stability subfunction
    cluster_stability_subfunction <- function(x) {
      clusterRows(x,
                  NNGraphParam(
                    k = k,
                    type = weighting_type,
                    cluster.fun = cluster_type,
                  ))
    }
  } else {
    stop("Clustering type is not valid. Please use --cluster_type to specify 'kmeans', 'walktrap' or 'louvain' clustering.")
  }
  
  # perform bootstrapping across a given number of iterations
  ari <- c()
  for (iter in 1:iterations){
      # sample cells with replacement 
      sample_cells <- sample(nrow(pca_matrix), nrow(pca_matrix), replace=TRUE)
      resampled_pca <- pca_matrix[sample_cells,,drop=FALSE]
      
      # perform clustering on sampled cells 
      resampled_clusters <- cluster_stability_subfunction(resampled_pca)
      
      # calculate ARI between new clustering and original clustering 
      ari[iter] <- pdfCluster::adj.rand.index(resampled_clusters, cluster_assignments[sample_cells])
    }
  return(ari)
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

create_metadata_stats_df <- function(clustered_sce, params_range, cluster_type) {
  # Purpose: Calculate and return a data frame with the validity stats of the
  # clusters in the SingleCellExperiment object
  
  # Args:
  #   clustered_sce: SingleCellExperiment object with clustered results
  #   params_range: the range of numeric parameters to test for clustering
  #   cluster_type: the type of clustering method performed - can be "kmeans", 
  #                 "walktrap", or "louvain"
  
  # define cluster names
  cluster_names_column <- sprintf("%s_%02d", cluster_type, params_range)
  
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
    dplyr::group_by(cluster_names_column, cluster_type, param_value) %>%
    # here we calculate and store the median values of the cluster stats in
    # columns beginning with `avg_`, while we calculate and store the median 
    # absolute deviation (MAD) values in columns beginning with `mad_`
    dplyr::summarize(avg_purity = median(purity),
                     mad_purity = mad(purity),
                     avg_width = median(width),
                     mad_width = mad(width))
  
  return(validity_summary_df)
}

plot_cluster_purity <- function(cluster_validity_df, num_col, point_size = 0.7) {
  
  # Purpose: Generate a plot displaying the cluster purity stats of the clusters
  # in the SingleCellExperiment object
  
  # Args:
  #   clustered_validity_df: data.frame with cluster validity stats associated
  #                          with their relevant cluster names
  #   num_col: number of columns to use when facetting
  #   point_size: desired size of each plotted point; default is 0.7
  
  # prepare data frame for plotting
  metadata <- cluster_validity_df %>%
    # create a column for the color scale to make things easier to control the color later
    # for cluster purity, do the majority of neighboring cells come from the assigned cluster (yes) or a different cluster (no)
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
    ggbeeswarm::geom_quasirandom(method = "smiley", size = point_size) +
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
    labs(x = "Cluster Assignment",
         color = legend_title) +
    facet_wrap( ~ cluster_names_column, scale="free_x", ncol = num_col, dir = "v")
  
  return(plot)
}

plot_cluster_silhouette_width <- function(cluster_validity_df, num_col, point_size = 0.7) {
  # Purpose: Calculate and return a data frame with the validity stats of the
  # clusters in the SingleCellExperiment object
  
  # Args:
  #   clustered_validity_df: data.frame with cluster validity stats associated
  #                          with their relevant cluster names
  #   num_col: number of columns to use when facetting
  #   point_size: desired size of each plotted point; default is 0.7
  
  # prepare data frame for plotting
  metadata <- cluster_validity_df %>%
    # for silhouette width if the cluster matches, then the silhouette width is positive, if not it's negative
    dplyr::mutate(
      param_value = as.numeric(param_value),
      cluster_param_assignment = paste(cluster_type, param_value, cluster, sep = "_")
    )
  
  # plot the cluster validity data frames
  plot <-
    ggplot(metadata, aes(x = cluster, y = width)) +
    ggbeeswarm::geom_quasirandom(method = "pseudorandom", size = point_size) +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    labs(x = "Cluster Assignment") +
    facet_wrap( ~ cluster_names_column, scale="free_x", ncol = num_col, dir = "v") +
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
    )
  
  return(plot)
}                                

plot_avg_validity_stats <- function(cluster_validity_summary_df_list,
                                    measure){
  # Purpose: Plot the summary stats for each clustering type using the list
  #          of provided cluster validity summary data frames
  
  # Args:
  #   cluster_validity_summary_df_list: list of data frames with cluster 
  #                                     validity stats associated with their 
  #                                     relevant cluster names
  #   measure: the cluster validity measure to use for plotting, can be 
  #            "avg_purity" or "avg_width"
  
  # prepare a data frame for plotting
  cluster_validity_summary_df <- dplyr::bind_rows(cluster_validity_summary_df_list) %>%
    dplyr::mutate(param_value = as.numeric(param_value))
  
  # convert summary symbols for plotting
  summary_y <- rlang::sym(measure)
  
  # grab column with median absolute deviation info
  if (measure == "avg_purity") {
    mad_column <- "mad_purity"
    # define y-axis range by finding the largest/smallest possible values of the median +/- mad 
    # if those values are outside the possible range of the metric then use them to define the limits of the y-axis
    y_lower <- min(0, min(cluster_validity_summary_df$avg_purity - cluster_validity_summary_df$mad_purity))
    y_upper <-  max(1, max(cluster_validity_summary_df$avg_purity + cluster_validity_summary_df$mad_purity))
    y_range <- c(y_lower, y_upper)
  } else if (measure == "avg_width") {
    mad_column <- "mad_width"
    y_lower <- min(-1, min(cluster_validity_summary_df$avg_width - cluster_validity_summary_df$mad_width))
    y_upper <- max(1, max(cluster_validity_summary_df$avg_width + cluster_validity_summary_df$mad_width))
    y_range <- c(y_lower, y_upper)
  } else {
    stop("Please specify 'avg_purity' or 'avg_width' to the `measure` argument.")
  }
  
  mad_column <- rlang::sym(mad_column)
  
  # plot the summary stats
  summary_plot <- ggplot(
    cluster_validity_summary_df,
    aes(
      x = param_value,
      y = !!summary_y,
      color = cluster_type)) +
    geom_pointrange(aes(x = param_value, y = !!summary_y, 
                        ymin = !!summary_y - !!mad_column,
                        ymax = !!summary_y + !!mad_column),
                    color = "black",
                    position = position_dodge2(width = 0.6)) +
    geom_line() +
    ylim(y_range) + 
    labs(x = "Parameter value",
         y = gsub("_", " ", measure),
         color = "Cluster type") +
    guides(color = guide_legend(override.aes = list(size = 5)))
  
  return(summary_plot)
}

get_cluster_stability_summary <- function(normalized_sce,
                                   params_range,
                                   cluster_type) {
  # Purpose: Calculate and return a data frame of ARI values of the bootstrapping
  # replicates and the associated original clusters stored in the 
  # SingleCellExperiment object for the specified clustering parameters
  
  # Args:
  #   normalized_sce: normalized SingleCellExperiment object
  #   params_range: the range of numeric parameters to test for clustering
  #   cluster_type: the type of clustering method performed - can be "kmeans",
  #                 "walktrap", or "louvain"
  
  # extract the principal components matrix
  pca_matrix <- reducedDim(normalized_sce, "PCA")
  
  # define cluster names
  if (cluster_type == "kmeans") {
    cluster_names <- sprintf("%s_%02d", cluster_type, params_range)
  } else if (cluster_type %in% c("walktrap", "louvain")) {
    cluster_names <- sprintf("%s_%02d", cluster_type, params_range)
  }
  
  ari_results <- list()
  
  for (name in cluster_names) {
    # grab k from cluster name
    k <- as.numeric(sub('.+_(.+)', '\\1', name))
    
    # run custom `check_cluster_stability` function and save results to a list
    ari_results[[name]] <- check_cluster_stability(
      pca_matrix,
      cluster_assignments = normalized_sce[[name]],
      cluster_type,
      k
    )
  }
  
  # combine ARI results from each of the cluster assignments and get data frame
  # ready for plotting
  plot_ari_df <- data.frame(ari_results) %>%
    tibble::rownames_to_column("bootstrap_iteration") %>%
    tidyr::pivot_longer(-c("bootstrap_iteration"), values_to = "ARI", names_to = "cluster_names_column") %>%
    tidyr::separate(cluster_names_column, 
                    c("cluster_type", "param_value"),
                    remove = FALSE) %>% # keep original column for grouping later
    dplyr::mutate(param_value = as.numeric(param_value))
  
  return(plot_ari_df)
}

plot_cluster_stability_ari <- function(ari_plotting_df, point_size = 0.7) {
  # Purpose: Plot the ARI values of the bootstrapping replicates and the 
  # associated original clusters stored in the SingleCellExperiment object for 
  # the specified clustering parameters
  
  # Args:
  #   ari_plotting_df: data frame with ARI values for plotting
  #   point_size: desired size of each plotted point; default is 0.7
  
  # set params as factors
  ari_plotting_df <- ari_plotting_df %>%
    dplyr::mutate(param_value = as.factor(param_value))
  
  # set params range
  params_range <- sort(unique(ari_plotting_df$param_value))
  
  plot <-
    ggplot(ari_plotting_df, aes(x = param_value, y = ARI, group = param_value)) +
    geom_violin() +
    ggforce::geom_sina(size = point_size) +
    stat_summary(
      aes(group = param_value),
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
      size = 0.2
    ) +
    scale_x_discrete(name = "Parameter value",
                     limits = params_range) +
    facet_wrap( ~ cluster_type)
  
  return(plot)
}

plot_summary_cluster_stability_ari <- function(ari_df_list) {
  # Purpose: Calculate a summary data frame of the provided cluster
  # stability ARI values and return a plot displaying these summary values
  
  # Args:
  #   ari_df_list: list of data frames with cluster stability ARI values 
  #                associated with their relevant clustering type and param values
  
  # prepare a data frame for plotting
  ari_combined_df <- dplyr::bind_rows(ari_df_list)
  
  # create a summary data.frame of the results across the individual clusters
  ari_summary_df <- ari_combined_df %>%
    dplyr::group_by(cluster_names_column, cluster_type, param_value) %>%
    # here we calculate and store the median of the ARI values along with the
    # median absolute deviation (MAD) values
    dplyr::summarize(median_ARI = median(ARI),
                     mad_ARI = mad(ARI),
                     .groups = 'drop')

  
  # plot the summary ARI values
  ari_summary_plot <- ggplot(
    ari_summary_df,
    aes(
      x = param_value,
      y = median_ARI,
      color = cluster_type)) +
    geom_pointrange(aes(x = param_value, y = median_ARI, 
                        ymin = median_ARI - mad_ARI,
                        ymax = median_ARI + mad_ARI),
                    color = "black",
                    position = position_dodge2(width = 0.6)) +
    geom_line() +
    labs(x = "Parameter value",
         y = "Median ARI",
         color = "Cluster type")
  
  return(ari_summary_plot)
}
