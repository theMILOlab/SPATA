

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
                               smooth_span = 0.02,
                               max_hotspots = 9,
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
                                      max_cluster = max_hotspots,
                                      n_start = n_start,
                                      .f = evaluate_gene_cluster_tendency)
    )

  unnested_df <-
    tidyr::unnest(evaluated_df, cols = "hotspot_evaluation") %>%
    dplyr::select(-data) %>%
    dplyr::ungroup()


  # 7. Evaluate sample hotspots ---------------------------------------------

  # compile hotspot df
  icld_quantiles <-
    stats::quantile(unnested_df$intra_cluster_dist)

  threshold_icld <- icld_quantiles[2]

  hotspot_df <-
    dplyr::filter(unnested_df, intra_cluster_dist < threshold_icld)

  hotspot_k_total_wss <-
    iterate_over_kmeans(hotspot_df, vars = c("center_x", "center_y"))

  hotspot_optimal_k <-
    compute_k_optimum(df = hotspot_k_total_wss)

  hotspot_pam <-
    dplyr::select(hotspot_df, center_x, center_y) %>%
    cluster::pam(k = hqc_optimal_k)

  hotspot_df$hotspots <-
    base::as.factor(stringr::str_c("htsp", hqc_pam$clustering, sep = "_"))

  # compile hotspot info df
  hotspot_info <-
    dplyr::group_by(.data = hotspot_df, hotspots) %>%
    dplyr::summarise(
      n_genes = dplyr::n_distinct(genes),
      n_barcode_spots = dplyr::n_distinct(center_barcodes),
      mean_intra_hotspot_dist = base::mean(intra_cluster_dist)
    ) %>%
    base::cbind(., hotspot_pam$medoids) %>%
    dplyr::mutate(
      center_barcodes = hotspot_df$center_barcodes[hotspot_pam$id.med],
      center_genes = hotspot_df$genes[hotspot_pam$id.med],
      average_sil_widths = hotspot_pam$silinfo$clus.avg.widths
    ) %>%
    base::cbind(., hotspot_pam$clusinfo[, c("max_diss", "av_diss", "diameter", "separation")])

  # compute hotspot enrichtment df
  summarised_df <-
    dplyr::left_join(x = getGeneSetDf(object), y = hotspot_df, by = c("gene" = "genes")) %>%
    dplyr::group_by(ont) %>%
    dplyr::mutate(exists = !base::is.na(hotspots), n_genes = dplyr::n()) %>%
    dplyr::group_by(ont, hotspots) %>%
    dplyr::summarise(
      n_genes_found = base::sum(exists),
      p_genes_found = n_genes_found/n_genes,
      n_genes = base::mean(n_genes)
    ) %>%
    tidyr::drop_na() %>%
    dplyr::ungroup()

  hotspot_enrichment_df <-
    dplyr::group_by(summarised_df, hotspots, ont) %>%
    dplyr::summarise(enrichment_score = base::mean(p_genes_found * n_genes)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(hotspots) %>%
    dplyr::arrange(dplyr::desc(enrichment_score), .by_group = TRUE)


  # last. Store results and return object -----------------------------------

  mtr_name <- getActiveMatrixName(object, of_sample = of_sample)

  add_to_gmdf <-
    dplyr::select(unnested_df, genes, n_cluster, mean_intra_cluster_dist, median_intra_cluster_dist)

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
      "max_hotspots" = max_hotspots,
      "mtr_name" = getActiveMatrixName(object, of_sample = of_sample),
      "n_start" = n_start,
      "sample" = of_sample,
      "smooth_span" = smooth_span,
      "suggestion" = list("df" = dplyr::arrange(hotspot_df, hotspots),
                          "gse_df" = hotspot_enrichment_df,
                          "info" = hotspot_summary),
      "threshold_icld" = threshold_icld,
      "threshold_stpv" = threshold_stpv,
      "threshold_stw" = threshold_stw,
      "tot_wss" = hotspot_k_total_wss
    )

  object <- setHotspotList(object = object,
                           of_sample = of_sample,
                           hotspot_list = hotspot_list)

  base::return(object)

}




