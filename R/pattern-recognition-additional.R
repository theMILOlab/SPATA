

#' Title
#'
#' @param object
#' @param of_sample
#' @param n_patterns
#' @param plotPrSummary
#'
#' @return
#' @export
#'

findNPatterns <- function(object, n_patterns, plotPrSummary = TRUE,  of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  confuns::is_value(n_patterns, mode = "numeric")

  pr_list <-
    getPrResults(object, of_sample)

  pr_list_new <-
    hlpr_assess_pattern_res(
      object = object,
      cluster_eval_df = pr_list$df,
      max_patterns = pr_list$max_patterns,
      n_start = pr_list$n_start,
      of_sample = of_sample,
      smooth_span = pr_list$smooth_span,
      threshold_stpv = pr_list$threshold_stpv,
      threshold_stw = pr_list$threshold_stw,
      n_patterns = n_patterns
    )

  if(!base::isFALSE(plotPrSummary)){

    object <-
      setPrResults(object, of_sample, pr_list = pr_list_new)

    summary_plot <-
      confuns::call_flexibly(
        fn = "plotPrSummary",
        fn.ns = "SPATA",
        default = list(object = object)
      )

    plot(summary_plot)

  }

  base::return(pr_list)

}
