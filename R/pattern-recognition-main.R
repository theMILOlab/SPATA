

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
#' @param max_pattern
#' @param n_start
#'
#' @return
#' @export

runPatternRecognition <- function(object,
                                  genes = 1500,
                                  genes_additional = NULL,
                                  with_ties = TRUE,
                                  method_pr = "hotspot",
                                  threshold_qntl = 0.8,
                                  threshold_stw = 0.5,
                                  threshold_stpv = 0.1,
                                  smooth = TRUE,
                                  smooth_span = 0.02,
                                  verbose = TRUE,
                                  of_sample = NA){

  confuns::are_values(c("genes", "threshold_stw", "threshold_stpv"), mode = "numeric")

  # 1. Control --------------------------------------------------------------

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  # 2. Filter genes for normally distributed ones ---------------------------

  genes <-
    getGeneMetaDf(object = object, of_sample = of_sample) %>%
    dplyr::filter(stw >= threshold_stw & stpv <= threshold_stpv) %>%
    dplyr::slice_max(order_by = var, n = genes, with_ties = with_ties) %>%
    dplyr::pull(var = "genes") %>%
    base::unique()

  # 3. Join gene variables and smooth spatially  ----------------------------

  gene_df <-
    joinWith(
      object = object,
      spata_df = getCoordsDf(object, of_sample = of_sample),
      genes = genes,
      smooth = smooth,
      smooth_span = smooth_span,
      normalize = TRUE
    )

  # 4. Mark barcode spots with expression levels below threshold with NA ----
  #    for every joined gene

  confuns::give_feedback(
    msg = "Preparing data for cluster evaluation.",
    verbose = verbose
  )

  marked_df <-
    mark_all_with_na(gene_df, n_qntls = threshold_qntl*10, keep_qntls = threshold_qntl*10)

  # 5. Shift the focus to the genes and nest the data.frame -----------------

  nested_df <-
    tidyr::pivot_longer(
      data = marked_df,
      cols = dplyr::all_of(genes),
      names_to = "genes",
      values_to = "values"
    ) %>%
    tidyr::drop_na() %>%
    dplyr::group_by(genes) %>%
    tidyr::nest()

  # keep total number of barcode spots
  n_bcsp <- base::nrow(nested_df$data[[1]])


  # 6. Iterate over all genes  ----------------------------------------------

  confuns::give_feedback(
    msg = glue::glue("Evaluating cluster tendency of {base::nrow(nested_df)} genes."),
    verbose = verbose
  )

  pb <-
    progress::progress_bar$new(
      format = "Progress: [:bar] :percent eta: :eta",
      total = base::nrow(nested_df), clear = FALSE, width = 100)

  pattern_evaluation_list <-
    purrr::map(.x = nested_df$data,
               pb = pb, verbose = verbose,
               .f = purrr::safely(.f = evaluate_gene_cluster_tendency_dbscan, otherwise = NA))

  # lgl vector where evaluation failed
  failed_evaluation <-
    purrr::map(.x = pattern_evaluation_list, .f = ~ base::all(base::is.na(.x[["result"]]))) %>%
    purrr::flatten_lgl()

  ignored_genes <-
    nested_df[failed_evaluation, ] %>%
    dplyr::pull(genes)

  confuns::give_feedback(
    msg = glue::glue("Evaluation failed for {f_genes} {ref_genes}.",
                     f_genes = base::length(ignored_genes),
                     ref_genes = "genes"),
    verbose = verbose
  )

  # keep only successful evaluations
  evaluated_df <-
    tibble::as_tibble(x = nested_df[!failed_evaluation, ]) %>%
    dplyr::mutate(
      pattern_evaluation = purrr::map(.x = pattern_evaluation_list[!failed_evaluation],
                                      .f = "result")
      )

  cluster_eval_df <-
    tidyr::unnest(evaluated_df, cols = "pattern_evaluation") %>%
    dplyr::select(-data) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(n_bcsp = {{n_bcsp}}, new_cluster_exclusivity = size / n_bcsp)


  # 7. Evaluate sample hotspots ---------------------------------------------

  if(method_pr == "hotspot"){

    pr_list <-
      hlpr_assess_hotspot_results(
        object = object,
        cluster_eval_df = cluster_eval_df,
        ignored_genes = ignored_genes,
        of_sample = of_sample,
        smooth_span = smooth_span,
        threshold_qntl = threshold_qntl,
        threshold_stpv = threshold_stpv,
        threshold_stw = threshold_stw,
        verbose = TRUE)

    object <-
      setPrResults(object = object,
                   of_sample = of_sample,
                   pr_list = pr_list,
                   method = method_pr)

  }

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  base::return(object)

}




