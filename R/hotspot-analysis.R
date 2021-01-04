

# Main function -----------------------------------------------------------

#' Title
#'
#' @param object
#' @param of_sample
#' @param threshold_stw
#' @param threshold_stpv
#' @param smooth_span
#' @param with_ties
#' @param normalize
#' @param max_cluster
#' @param n_start
#'
#' @return
#' @export
#'
#' @examples
runHotspotAnalysis <- function(object,
                               of_sample = "",
                               threshold_stw = 0.75,
                               threshold_stpv = 0.1,
                               smooth_span = 0.01,
                               max_cluster = 9,
                               n_start = 10,
                               verbose = TRUE){


  # 1. Control --------------------------------------------------------------

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  # 2. Filter genes for normally distributed ones ---------------------------

  genes <-
    getGeneMetaDf(object = object, of_sample = of_sample) %>%
    dplyr::filter(st_W >= threshold_stw & st_pv <= threshold_stpv) %>%
    dplyr::pull(var = "genes")


  # 3. Join gene variables and smooth spatially  ----------------------------

  gene_df <-
    joinWith(
      object = object,
      spata_df = getCoordsDf(object, of_sample = of_sample),
      genes = genes,
      smooth = TRUE,
      smooth_span = smooth_span,
      normalize = TRUE
    )


  # 4. Mark barcode spots with expression levels below threshold with NA ----

  marked_df <-
    dplyr::mutate(gene_df, dplyr::across(.cols = dplyr::all_of(x = genes), .fns = mark_with_na))


  # 5. Shift the focus to the genes and nest the data.frame -----------------

  nested_df <-
    tidyr::pivot_longer(
      data = marked_df,
      cols = dplyr::all_of(genes),
      names_to = "genes",
      values_to = "values"
    ) %>%
    dplyr::group_by(genes) %>%
    tidyr::nest()


  # 6. Iterate over all genes  ----------------------------------------------

  pb <-
    progress::progress_bar$new(
      format = "Progress: [:bar] :percent eta: :eta",
      total = base::nrow(nested_df), clear = FALSE, width = 100)

  evaluated_df <-
    dplyr::mutate(
      .data = nested_df,
      hotspot_evaluation = purrr::map(.x = data,
                                       pb = pb,
                                       verbose = verbose,
                                       .f = evaluate_gene_cluster_tendency)
    )

  unnested_df <-
    dplyr::transmute(
      .data = evaluated_df,
      genes = genes,
      n_cluster = purrr::map(.x = hotspot_evaluation, .f = "n_cluster") %>%
                  purrr::flatten_dbl(),
      avg_sil_width = purrr::map(.x = hotspot_evaluation, .f = ~ .x[["silinfo"]][["avg.width"]]) %>%
                      purrr::map(.x = ., .f = ~ base::ifelse(base::is.null(.x), NA, .x)) %>%
                      purrr::flatten_dbl(),
      cluster_sil_width = purrr::map(.x = hotspot_evaluation, .f = ~ .x[["silinfo"]][["clus.avg.widths"]]),
      centers = purrr::map(.x = hotspot_evaluation, .f = ~ .x[["medoids"]])
    ) %>%
    tidyr::unnest(dplyr::all_of(c("centers", "cluster_sil_width"))) %>%
    dplyr::group_by(genes) %>%
    dplyr::mutate(cluster = dplyr::row_number()) %>%
    dplyr::select(genes, n_cluster, avg_sil_width, cluster, cluster_sil_width, center_x = x, center_y = y)


  # last. Store results and return object -----------------------------------

  mtr_name <- getActiveMatrixName(object, of_sample = of_sample)

  add_to_gmdf <-
    dplyr::select(unnested_df, genes, n_cluster, avg_sil_width)

  gmd_list <-
    getGeneMetaData(object, of_sample = of_sample, mtr_name = mtr_name)

  gmd_list$df <-
    dplyr::left_join(x = gmd_list$df, y = add_to_gmdf, by = "genes") %>%
    dplyr::mutate(hotspot_evaluated = genes %in% add_to_gmdf$genes)

  object <- addGeneMetaData(object = object,
                            of_sample = of_sample,
                            meta_data_list = gmd_list)

  hotspot_list <-
    list(
      "df" = unnested_df,
      "mtr_name" = getActiveMatrixName(object, of_sample = of_sample),
      "threshold_stw" = threshold_stw,
      "threshold_stpv" = threshold_stpv,
      "smooth_span" = smooth_span,
      "max_cluster" = max_cluster,
      "n_start" = n_start
    )

  object <- setHotspotList(object = object,
                           of_sample = of_sample,
                           hotspot_list = hotspot_list)




  base::return(object)

}




# Subfunctions ------------------------------------------------------------

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
mark_with_na <- function(variable, n_quantiles = 5, mark_below = 5){

  var_quantiles <-
    stats::quantile(x = variable, probs = seq(0, 1, length.out = n_quantiles))

  threshold <- var_quantiles[mark_below - 1]

  variable[variable < threshold] <- NA

  base::return(variable)

}



# use for visualization

#' Title
#'
#' @param marked_df
#' @param var
#' @param pt_size
#'
#' @return
#' @export
#'
#' @examples
#'
plotSurfaceQuantile <- function(marked_df, var, pt_size = 2){

  var <- marked_df[,var] %>% base::names()

  plot_df <-
    dplyr::mutate(.data = df,
                  color = dplyr::if_else(base::is.na(!!rlang::sym(var)), "discard", "keep")
    ) %>%
    dplyr::select(x, y, color, dplyr::all_of(x = var))

  ggplot(plot_df, aes(x = x, y = y, color = color, alpha = color)) +
    geom_point(size = pt_size) +
    theme_bw() +
    theme(legend.position = "right", panel.grid = element_blank()) +
    scale_color_manual(values = c("keep" = "forestgreen", "discard" = "lightgrey")) +
    scale_alpha_manual(values = c("keep" = 1, "discard" = 0.1))

}


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
generate_reference_df <- function(dropped_df, marked_df){

  n_barcode_spots <- base::nrow(dropped_df)

  n_corner <- base::round(n_barcode_spots/4, digits = 0)

  s_min_x <-
    marked_df$x %>% stats::quantile(probs = base::seq(0,1, length.out = 10)) %>%
    utils::head(2) %>% generate_sample(l = n_corner)

  s_min_y <-
    marked_df$y %>% stats::quantile(probs = base::seq(0,1, length.out = 10)) %>%
    utils::head(2) %>% generate_sample(l = n_corner)

  s_max_x <-
    marked_df$x %>% stats::quantile(probs = base::seq(0,1, length.out = 10)) %>%
    utils::tail(2) %>% generate_sample(l = n_corner)

  s_max_y <-
    marked_df$y %>% stats::quantile(probs = base::seq(0,1, length.out = 10)) %>%
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
                                max_cluster = 9,
                                n_start = 10,
                                ...){

  purrr::map_df(.x = 1:max_cluster,
                df = dropped_df[, c("x", "y")],
                n_start = n_start,
                .f = function(k, df, n_start){

    model <- stats::kmeans(x = df, centers = k, nstart = n_start, ...)

    tibble::tibble(
      k = k,
      tot_wss = model$tot.withinss,
      centers = list(model$centers)
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
evaluate_gene_cluster_tendency <- function(marked_df, verbose = TRUE, pb = NULL){

  if(base::isTRUE(verbose)){ pb$tick() }

  dropped_df <-
    tidyr::drop_na(marked_df[, c("x", "y", "values")])

  reference_df <-
    generate_reference_df(dropped_df)

  var_twss <-
    iterate_over_kmeans(dropped_df)

  reference_twss <-
    iterate_over_kmeans(reference_df, max_cluster = 1)

  var_twss_n <-
    tibble::add_row(
      .data = var_twss,
      k = 0,
      tot_wss = reference_twss$tot_wss[1],
      .before = 1
    )

  optimal_k <-
    compute_k_optimum(df = var_twss_n)


  pam_obj <- compute_pam_with_k_optimum(dropped_df, k = optimal_k)

  res_list <- list(
    "n_cluster" = optimal_k,
    "silinfo" = pam_obj$silinfo[c("clus.avg.widths", "avg.width")],
    "cluster_info" = pam_obj$clusinfo,
    "medoids" = base::as.data.frame(pam_obj$medoids)
  )

}

