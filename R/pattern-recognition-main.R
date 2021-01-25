

# Main function -----------------------------------------------------------

#' @title Find and reconstruct underlying spatial gene expression patterns
#'
#' @inherit check_sample params
#' @inherit check_smooth params
#' @param threshold_stw Numeric value. The minimal shapiro test W-value
#' a gene must have in order to be included.
#' @param threshold_stpv Numeric value. The maximal shapiro test p-value
#' a gene must have in order to be included.
#' @param with_ties Logical value. Indicates whether ties are kept. If set
#' to TRUE this might result in more genes than specified in \code{n_genes}.
#'
#' @return An updated spata-object containing the analysis results.
#' @export

runPatternRecognition <- function(object,
                                  n_genes = 5000,
                                  genes_additional = NULL,
                                  with_ties = TRUE,
                                  method_pr = "hotspot",
                                  threshold_qntl = 0.8,
                                  threshold_stpv = 0.1,
                                  threshold_stw = 0.5,
                                  smooth = FALSE,
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
    dplyr::slice_max(order_by = var, n = n_genes, with_ties = with_ties) %>%
    dplyr::pull(var = "genes") %>%
    base::unique()

  # 3. Join gene variables and smooth spatially  ----------------------------

  coords_df <- getCoordsDf(object, of_sample = of_sample)

  gene_df <-
    joinWith(
      object = object,
      spata_df = coords_df,
      genes = genes,
      smooth = smooth,
      smooth_span = smooth_span,
      normalize = TRUE
    )

  genes <- genes[genes %in% base::colnames(gene_df)]

  # 4. Mark barcode spots with expression levels below threshold with NA ----
  #    for every joined gene

  confuns::give_feedback(
    msg = "Preparing gene data for cluster evaluation.",
    verbose = verbose
  )

  marked_df <-
    mark_all_with_na2(gene_df, percentile = threshold_qntl) # arbitrary threshold

  # 5. Shift focus to the genes and nest the data.frame ---------------------

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
               .f = purrr::safely(.f = evaluate_gene_cluster_tendency_dbscan2, otherwise = NA))

  # lgl vector where evaluation failed
  failed_evaluation <-
    purrr::map(.x = pattern_evaluation_list, .f = ~ base::all(base::is.na(.x[["result"]]))) %>%
    purrr::flatten_lgl()

  failed_genes <-
    nested_df[failed_evaluation, ] %>%
    dplyr::pull(genes)

  if(base::length(failed_genes) >= 1){

    confuns::give_feedback(
      msg = glue::glue("Evaluation failed for {f_genes} {ref_genes}.",
                       f_genes = base::length(failed_genes),
                       ref_genes = confuns::adapt_reference(failed_genes, sg = "gene", pl = "genes")),
      verbose = verbose
    )

  }

  # keep only successful evaluations

  successful_evaluations <-
    purrr::map(.x = pattern_evaluation_list[!failed_evaluation], .f = "result")

  cluster_eval_df <-
    tibble::as_tibble(x = nested_df[!failed_evaluation, ]) %>%
    dplyr::mutate(pattern_evaluation = successful_evaluations) %>%
    tidyr::unnest(cols = "pattern_evaluation") %>%
    dplyr::select(-data) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(noise_total = size_total - size_noisless)


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




