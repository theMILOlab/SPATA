
#' @title Filter gene-set data.frame
#'
#' @description Checks the objects gene-set data.frame for gene-sets that
#' are composed of genes that exist in the given expression matrix.
#'
#' @inherit check_object params
#' @param limit Numeric value between 1 and 100. The minimum percentage of gene-set genes
#' that have to exist in the given expression matrix in order for a gene set to stay in the
#' gene-set data.frame.
#'
#' @return An updated spata-object and an informative message about how many
#' gene-sets have been discarded and how many gene-sets remain.
#'
#' @details Example: Assuming that gene-set 'x' is composed of 30 genes. The expression matrix
#' however contains only 15 of them. If argument \code{limit} is set to 75 gene-set 'x'
#' is removed since the percentage of genes of which the given expression matrix
#' contains information about is only 50.
#'
#' @export

adjustGeneSetDf <- function(object, limit = 50){

  # 1. Control --------------------------------------------------------------

  check_object(object)
  confuns::is_value(limit, mode = "numeric", ref = "limit")
  if(!dplyr::between(limit, left = 1, right = 99)){

    base::stop("Argument 'limit' needs to be a numeric value between 1 and 99.")

  }

  limit <- limit/100

  # -----

  # 2. Cleaning -------------------------------------------------------------

  base::message(glue::glue("Calculating percentage of genes found in expression matrix for {dplyr::n_distinct(object@used_genesets$ont)} gene sets."))

  all_genes <- getGenes(object, simplify = TRUE, in_sample = "all")

  filtered_df <-
    dplyr::group_by(.data = object@used_genesets, ont) %>%
    dplyr::mutate(
      gene_count = dplyr::n(),
      gene_found = gene %in% all_genes,
      n_found = base::sum(gene_found),
      p_found = base::round(n_found/gene_count, digits = 2)
    ) %>%
    dplyr::filter(p_found > {{limit}}) %>%
    dplyr::ungroup()

  n_all_gs <-
    getGeneSets(object) %>%
    base::length()

  n_remaining_gs <-
    dplyr::pull(filtered_df, var = ont) %>%
    base::unique() %>%
    base::length()

  n_removed_gs <- n_all_gs - n_remaining_gs

  base::message(glue::glue("Removed {n_removed_gs} gene-sets. Number of remaining gene-sets: {n_remaining_gs} "))

  object@used_genesets <-
    dplyr::select(filtered_df, ont, gene)

  base::return(object)

}


