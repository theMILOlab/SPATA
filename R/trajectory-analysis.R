#' @title Trajectory expression analysis
#'
#' @description These functions analyze the expression dynamics along
#' a specified trajectory.
#'
#' @inherit check_sample params
#' @inherit check_trajectory params
#' @inherit check_gene_sets params
#' @inherit check_genes params
#' @inherit check_method params
#' @inherit verbose params
#'
#' @return A nested data.frame with information about the dynamics of each gene
#' or gene set.
#'
#' @export

rankTrajectoryGeneSets <- function(object,
                                   trajectory_name,
                                   of_sample,
                                   gene_sets,
                                   method_gs = "mean",
                                   method_padj = "fdr",
                                   verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  # all checks
  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  check_trajectory(object, trajectory_name = trajectory_name, of_sample = of_sample)
  check_method(method_gs = method_gs, method_padj = method_padj)

  gene_sets <- check_gene_sets(object = object, gene_sets = gene_sets)

  # -----


  # 2. Prepare source data --------------------------------------------------

  t_object <- getTrajectoryObject(object = object,
                                  trajectory_name = trajectory_name,
                                  of_sample = of_sample)

  # nest compiled trajectory data.frame
  nested_ctdf <-
    hlpr_summarize_trajectory_df(
      object = object,
      ctdf = t_object@compiled_trajectory_df,
      variables = gene_sets,
      method_gs = method_gs,
      accuracy = 5,
      verbose = verbose) %>%
    dplyr::group_by(!!rlang::sym("gene_sets")) %>%
    tidyr::nest()

  # -----


  # 3. Modeling pipelines ---------------------------------------------------

  # linear modeling

  if(base::isTRUE(verbose)){base::message("Computing linear models.")}

  lm_ctdf <-
    dplyr::mutate(.data = nested_ctdf,
                  lm = purrr::map(.x = data, .f = hlpr_lm_trajectory),
                  lm_info = purrr::map(.x = lm, .f = hlpr_lm_info)) %>%
    tidyr::unnest(lm_info) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(lm_adj_pvalue = stats::p.adjust(p = lm_pvalue,
                                                  method = method_padj)) %>%
    dplyr::arrange(dplyr::desc(x = lm_adj_rsquared))

  # -----

  base::return(lm_ctdf)

}


#' @rdname rankTrajectoryGeneSets
#' @export

rankTrajectoryGenes <- function(object,
                                trajectory_name,
                                of_sample,
                                genes,
                                method_padj = "fdr",
                                verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  # all checks
  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  check_trajectory(object, trajectory_name = trajectory_name, of_sample = of_sample)
  check_method(method_gs = method_gs, method_padj = method_padj)

  genes <- check_genes(object = object, genes = genes)

  # -----

  # 2. Prepare source data --------------------------------------------------

  t_object <- getTrajectoryObject(object = object,
                                  trajectory_name = trajectory_name,
                                  of_sample = of_sample)

  # nest compiled trajectory data.frame
  nested_ctdf <-
    hlpr_summarize_trajectory_df(
      object = object,
      ctdf = t_object@compiled_trajectory_df,
      variables = genes,
      method_gs = method_gs,
      accuracy = 5,
      verbose = verbose) %>%
    dplyr::group_by(!!rlang::sym("genes")) %>%
    tidyr::nest()

  # -----

  # 3. Modeling pipelines ---------------------------------------------------

  # linear modeling

  if(base::isTRUE(verbose)){base::message("Computing linear models.")}

  lm_ctdf <-
    dplyr::mutate(.data = nested_ctdf,
                  lm = purrr::map(.x = data, .f = hlpr_lm_trajectory),
                  lm_info = purrr::map(.x = lm, .f = hlpr_lm_info)) %>%
    tidyr::unnest(lm_info) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(lm_adj_pvalue = stats::p.adjust(p = lm_pvalue,
                                                  method = method_padj)) %>%
    dplyr::arrange(dplyr::desc(x = lm_adj_rsquared))

  # -----

  base::return(lm_ctdf)

}
