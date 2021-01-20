
#' Title
#'
#' @param variable
#' @param n_quantiles
#' @param mark_below
#'
#' @return
#' @export
#'
#' @examples
mark_with_na <- function(variable, n_qntls = 5, keep_qntls = 5){

  var_quantiles <-
    dplyr::ntile(x = variable, n = n_qntls) %>%
    base::as.character()

  set_to_na <- !var_quantiles %in% base::as.character(keep_qntls)

  variable[set_to_na] <- NA

  base::return(variable)

}

#' @rdname mark_with_na
#' @export

mark_all_with_na <- function(coords_df, n_qntls = 5, keep_qntls = 5){

  check_coords_df(coords_df)

  numeric_vars <-
    dplyr::select(coords_df, -dplyr::all_of(x = coords_df_vars)) %>%
    dplyr::select_if(.predicate = base::is.numeric) %>%
    base::colnames()

  dplyr::mutate(
    .data = coords_df,
    dplyr::across(.cols = numeric_vars, .fns = mark_with_na,
                  n_qntls = n_qntls, keep_qntls = keep_qntls)
  )

}



# use for visualization




#' Title
#'
#' @param x
#' @param l
#'
#' @return
#' @export
#'
#' @examples
generate_sample <- function(x, l = 250){

  base::seq(x[1], x[2], length.out = l) %>%
    base::sample(size = l, replace = TRUE)

}

#
#' Title
#'
#' @param dropped_df dropped df has variables c("x", "y", value)
#' and is the result of tidyr::drop_na()
#'
#' @return
#' @export
#'
#' @examples
generate_reference_df <- function(n, x_coords, y_coords){

  n_barcode_spots <- n

  n_corner <- base::round(n_barcode_spots/4, digits = 0)

  s_min_x <-
    x_coords %>% stats::quantile(probs = base::seq(0,1, length.out = 10)) %>%
    utils::head(2) %>% generate_sample(l = n_corner)

  s_min_y <-
    y_coords %>% stats::quantile(probs = base::seq(0,1, length.out = 10)) %>%
    utils::head(2) %>% generate_sample(l = n_corner)

  s_max_x <-
    x_coords %>% stats::quantile(probs = base::seq(0,1, length.out = 10)) %>%
    utils::tail(2) %>% generate_sample(l = n_corner)

  s_max_y <-
    y_coords %>% stats::quantile(probs = base::seq(0,1, length.out = 10)) %>%
    utils::tail(2) %>% generate_sample(l = n_corner)

  reference_df <-
    data.frame(
      x = c(s_min_x, s_min_x, s_max_x, s_max_x),
      y = c(s_max_y, s_min_y, s_min_y, s_max_y)
    )

  base::return(reference_df)

}

#' Title
#'
#' @param dropped_df
#' @param max_cluster
#' @param n_start
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
iterate_over_kmeans <- function(dropped_df,
                                vars = c("x", "y"),
                                max_cluster = 9,
                                n_start = 10,
                                ...){

  purrr::map_df(.x = 1:max_cluster,
                df = dropped_df[, vars],
                n_start = n_start,
                .f = function(k, df, n_start){

                  model <- stats::kmeans(x = df, centers = k, nstart = n_start)

                  tibble::tibble(
                    k = k,
                    tot_wss = model$tot.withinss,
                    centers = list(model$centers)
                  )

                })

}


iterate_over_pam <- function(df,
                             vars = c("x", "y"),
                             max_cluster = 9,
                             n_start = 10,
                             ...){

  purrr::map(.x = 2:max_cluster,
             df = df[, vars],
             .f = function(k, df){

               cluster::pam(
                 x = df,
                 k = k,
                 metric = "euclidean"
               )


             })

}


#' Title
#'
#' @param variance
#'
#' @return
#' @export
#'
#' @examples
find_elbow <- function(variance) {
  if (base::is.unsorted(-variance)) {

    base::return(0)

  } else {

    # Finding distance from each point on the curve to the diagonal.
    dy <- -base::diff(base::range(variance))
    dx <- base::length(variance) - 1
    l2 <- base::sqrt(dx^2 + dy^2)
    dx <- dx/l2
    dy <- dy/l2

    dy0 <- variance - variance[1]
    dx0 <- base::seq_along(variance) - 1

    parallel.l2 <- base::sqrt((dx0 * dx)^2 + (dy0 * dy)^2)
    normal.x <- dx0 - dx * parallel.l2
    normal.y <- dy0 - dy * parallel.l2
    normal.l2 <- base::sqrt(normal.x^2 + normal.y^2)

    # Picking the maximum normal that lies below the line.
    # If the entire curve is above the line, we just pick the last point.
    below.line <- normal.x < 0 & normal.y < 0
    if (!base::any(below.line)) {

      base::length(variance) %>%
        base::return()

    } else {

      base::which(below.line)[base::which.max(normal.l2[below.line])] %>%
        base::return()

    }
  }
}

#' Title
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
compute_k_optimum <- function(df){

  row_index <-
    find_elbow(variance = df$tot_wss)

  base::return(df$k[row_index])

}


compute_pam_with_k_optimum <- function(dropped_df,
                                       k,
                                       ...){

  cluster::pam(
    x = dropped_df[, c("x", "y")],
    k = k,
    ...
  )

}


#' Title
#'
#' @param marked_df
#' @param verbose
#' @param pb
#'
#' @return
#' @export
#'
#' @examples
evaluate_gene_cluster_tendency <- function(marked_df,
                                           verbose = TRUE,
                                           pb = NULL,
                                           ref_totwss,
                                           ...){

  if(base::isTRUE(verbose)){ pb$tick() }

  dropped_df <-
    tidyr::drop_na(marked_df[, c("x", "y", "values", "barcodes")])

  var_twss <-
    iterate_over_kmeans(dropped_df)

  var_twss_n <-
    tibble::add_row(
      .data = var_twss,
      k = 0,
      tot_wss = ref_totwss,
      .before = 1
    )

  optimal_k <-
    compute_k_optimum(df = var_twss_n)

  if(optimal_k > 1){

    dist_mtr <-
      stats::dist(x = dropped_df[, c("x", "y", "values")])

    hclust <-
      stats::hclust(d = dist_mtr)

    hclust_vec <-
      stats::cutree(tree = hclust, k = optimal_k)

    clustered_df <-
      dplyr::mutate(dropped_df, cluster = base::as.factor(x = hclust_vec))

  } else {

    clustered_df <-
      dplyr::mutate(dropped_df, cluster = base::as.factor(x = "1"))

  }

  res_df <-
    assess_intra_cluster_distance(clustered_df = clustered_df)

  base::return(res_df)

}

#' @rdname evaluate_gene_cluster_tendency
#' @export
evaluate_gene_cluster_tendency_dbscan <- function(marked_df,
                                                  verbose = TRUE,
                                                  pb = NULL){

  if(base::isTRUE(verbose)){ pb$tick() }

  {
    dropped_df <- marked_df

    # use density based clustering to filter out noisy points
    dbc_res <-
      dbscan::dbscan(
        x = base::as.matrix(dropped_df[,c("x", "y")]),
        eps = dbscan::kNNdist(x = base::as.matrix(dropped_df[,c("x", "y")]), k = 3) %>% base::mean(),
        minPts = 3
      )

    dropped_df <-
      dplyr::mutate(.data = dropped_df, cluster = base::as.character(dbc_res$cluster)) %>%
      dplyr::filter(cluster != "0")

    n_bcs <- base::nrow(dropped_df)
    n_clusters <- dplyr::n_distinct(dropped_df$cluster)

    if(n_clusters == 1){

      cluster_keep <- "1"

    } else {

      cluster_count <-
        dplyr::count(x = dropped_df, cluster) %>%
        dplyr::filter(n >= n_bcs*1.5/n_clusters)

      cluster_keep <- cluster_count$cluster

    }

    dropped_df <-
      dplyr::filter(dropped_df, cluster %in% {{cluster_keep}})

    center_df <-
      dplyr::group_by(dropped_df, cluster) %>%
      dplyr::mutate(remaining_bcs = )
      dplyr::summarise(
        size = dplyr::n(),
        center_x = base::mean(x),
        center_y = base::mean(y)
      ) %>%
      dplyr::mutate(cluster = base::as.character(dplyr::row_number())) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(n_cluster = dplyr::n()) %>%
      dplyr::select(n_cluster, cluster, size, center_x, center_y)

  }

  base::return(center_df)

}

#' @rdname evaluate_gene_cluster_tendency
#' @export
evaluate_gene_cluster_tendency_dbscan2 <- function(marked_df, #!!! changed center_df compuation
                                                   verbose = TRUE,
                                                   pb = NULL){

  if(base::isTRUE(verbose)){ pb$tick() }

  {
    dropped_df <- marked_df

    size_total <- base::nrow(dropped_df)

    # use density based clustering to filter out noisy points
    dbc_res <-
      dbscan::dbscan(
        x = base::as.matrix(dropped_df[,c("x", "y")]),
        eps = dbscan::kNNdist(x = base::as.matrix(dropped_df[,c("x", "y")]), k = 3) %>% base::mean(), # arbitry threshold
        minPts = 3
      )

    dropped_df <-
      dplyr::mutate(.data = dropped_df, cluster = base::as.character(dbc_res$cluster)) %>%
      dplyr::filter(cluster != "0")

    dropped_df_first_cluster <- dropped_df

    n_bcs <- base::nrow(dropped_df)
    n_clusters <- dplyr::n_distinct(dropped_df$cluster)

    if(n_clusters == 1){

      cluster_keep <- "1"

    } else {

      cluster_count <-
        dplyr::count(x = dropped_df, cluster) %>%
        dplyr::filter(n >= n_bcs*1.5/n_clusters) # arbitrary

      cluster_keep <- cluster_count$cluster

    }

    dropped_df <-
      dplyr::filter(dropped_df, cluster %in% {{cluster_keep}})

    dropped_df_second_cluster <- dropped_df

    size_noisless <- base::nrow(dropped_df)

    center_df <-
      dplyr::group_by(dropped_df, cluster) %>%
      tidyr::nest() %>%
      dplyr::mutate(
        remaining_barcodes = purrr::map(data, .f = ~ dplyr::select(.x, barcodes)),
        center_x = purrr::map_dbl(data, .f = ~ base::mean(.x[["x"]], na.rm = TRUE)),
        center_y = purrr::map_dbl(data, .f = ~ base::mean(.x[["y"]], na.rm = TRUE)),
        size_cluster = purrr::map_dbl(data, .f = ~ base::nrow(.x))
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        cluster = forcats::as_factor(dplyr::row_number()),
        n_cluster = dplyr::n(),
        size_total = {{size_total}},
        size_noisless = {{size_noisless}},
        noise_total = size_total - size_noisless
        ) %>%
      dplyr::select(cluster, n_cluster, size_cluster, size_noisless, size_total, center_x, center_y,remaining_barcodes, -data )

  }

  base::return(center_df)

}

#' @rdname evaluate_gene_cluster_tendency
#' @export
evaluate_gene_cluster_tendency_hdbscan <- function(marked_df, #!!! every barcode spot gets sum of all cluster belonging propabilities
                                                   verbose = TRUE,
                                                   pb = NULL){

  if(base::isTRUE(verbose)){ pb$tick() }

  {
    dropped_df <- marked_df

    # use density based clustering to filter out noisy points
    hdbscan_res <- dbscan::hdbscan(x = dropped_df[,c("x","y")], minPts = 17)

    dropped_df_res <-
      dplyr::mutate(
        .data = dropped_df,
        cluster = base::as.character(hdbscan_res$cluster),
        membership_prob = base::as.numeric(hdbscan_res$membership_prob)
        ) %>%
      dplyr::filter(cluster != "0")

  }

  base::return(dropped_df_res)

}




#' Title
#'
#' @description Assesses the gene's tendency to cluster together.
#' The lower the \emph{intra_cluster_dist} the higher the tendency
#' to cluter.
#'
#' @param clustered_df dropped_df that contains an additional
#' cluster variable indicating the cluster belonging of
#' every barcode spot.
#'
#' @return
#' @export
#'
#' @details The distances between all barcode spots of each cluster are
#' computed and normalized for the size of the cluster.
#'
#' @examples
assess_intra_cluster_distance <- function(clustered_df){

  assessment_res_df <-
    purrr::map_df(.x = base::levels(clustered_df$cluster),
                  df = clustered_df,
                  .f = function(cl, df){

                    size <- base::nrow(df)

                    dist_mtr <-
                      dplyr::filter(df, cluster == {{cl}}) %>%
                      tibble::column_to_rownames(var = "barcodes") %>%
                      dplyr::select(x, y) %>%
                      stats::dist() %>%
                      base::as.matrix()

                    dist_mtr[base::upper.tri(x = dist_mtr, diag = TRUE)] <- NA

                    n_barcodes <- dplyr::n_distinct(base::rownames(dist_mtr))

                    hlpr_dist_mtr_to_df(dist_mtr, varnames = c("barcodes", "barcodes2")) %>%
                      dplyr::group_by(barcodes) %>%
                      dplyr::summarise(
                        mean_dist = base::mean(distance),
                        .groups = "drop") %>%
                      dplyr::ungroup() %>%
                      dplyr::filter(mean_dist == base::min(mean_dist)) %>%
                      dplyr::mutate(
                        intra_cluster_dist = mean_dist / n_barcodes,
                        size = {{size}}
                        )

                  }) %>%
    dplyr::mutate(
      n_cluster = dplyr::n(),
      cluster = base::levels(clustered_df$cluster),
      mean_intra_cluster_dist = base::mean(intra_cluster_dist),
      median_intra_cluster_dist = stats::median(intra_cluster_dist)
    )

  original_renamed_df <-
    dplyr::select(.data = clustered_df, center_x = x, center_y = y, barcodes)

  joined_df <-
    dplyr::left_join(x = assessment_res_df, y = original_renamed_df, by = "barcodes") %>%
    dplyr::select(n_cluster, cluster, size, median_intra_cluster_dist,
                  mean_intra_cluster_dist, intra_cluster_dist,
                  center_x, center_y, center_barcodes = barcodes)

  base::return(joined_df)

}




