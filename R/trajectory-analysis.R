
# Objects -----------------------------------------------------------------

trajectory_patterns <- c("Linear descending", "Linear ascending", "Gradient descending", "Logarithmic descending",
                         "Logarithmic ascending", "Gradient ascending","Sinus",  "Sinus (reversed)", "One peak",
                         "One peak (reversed)", "Two peaks (reversed)", "Two peaks")

# -----


# Functions ---------------------------------------------------------------


#' @title Trajectory expression analysis
#'
#' @description These functions analyze the expression dynamics along
#' a specified trajectory by fitting a variety of models to the trajectory's
#' expression trend.
#'
#' @param stdf A summarized trajectory data.frame. (e.g. obtained by
#' \code{getSummarizedTrajectoryDf()}).
#'
#' @return A ranked trajectory data.frame. These nested data.frames contain
#' the following variables:
#'
#' @return A nested data.frame with information about the dynamics of each gene
#' or gene set.
#'
#' @export

rankTrajectoryGeneSets <- function(stdf){

  # 1. Control --------------------------------------------------------------

  check_summarized_trajectory_df(stdf)
  base::stopifnot("gene_sets" %in% base::colnames(stdf))

  # -----

  # 2. Ranking --------------------------------------------------------------

  ranked_df <-
    dplyr::group_by(.data = stdf, !!rlang::sym("gene_sets")) %>%
    dplyr::mutate(values = confuns::normalize(x = values)) %>%
    tidyr::nest() %>>%
    "Adding models." %>>%
    dplyr::mutate(models = purrr::map(.x = data, .f = hlpr_add_models)) %>>%
    "Fitting models." %>>%
    dplyr::mutate(residuals = purrr::map(.x = data, .f = hlpr_add_residuals)) %>>%
    "Calculating residuals." %>>%
    dplyr::mutate(auc = purrr::map(.x = residuals, .f = hlpr_summarize_residuals))

  # -----

  return(ranked_df)

}

#' @rdname rankTrajectoryGeneSets
#' @export
rankTrajectoryGenes <- function(stdf){

  # 1. Control --------------------------------------------------------------

  check_summarized_trajectory_df(stdf)
  base::stopifnot("genes" %in% base::colnames(stdf))

  # -----

  # 2. Ranking --------------------------------------------------------------

  ranked_df <-
    dplyr::group_by(.data = stdf, !!rlang::sym("genes")) %>%
    dplyr::mutate(values = confuns::normalize(x = values)) %>%
    tidyr::nest() %>>%
    "Adding models." %>>%
    dplyr::mutate(models = purrr::map(.x = data, .f = hlpr_add_models)) %>>%
    "Fitting models." %>>%
    dplyr::mutate(residuals = purrr::map(.x = data, .f = hlpr_add_residuals)) %>>%
    "Calculating residuals." %>>%
    dplyr::mutate(auc = purrr::map(.x = residuals, .f = hlpr_summarize_residuals))

  # -----

  return(ranked_df)

}



#' @title Assess trajectory ranking.
#'
#' @description Takes a ranked trajectory data.frame and returns a data.frame
#' that informs about how well the ranked gene- or gene set exprssion-trends
#' fitted certain patterns.
#'
#' @param rtdf A ranked trajectory data.frame.
#' @param pattern The pattern(s) you are interested in specified as a character
#' vector. If set to NULL all patterns are included.
#' @param max_auc Numeric value. The maximum area-under-the-curve-value allowed.
#' @param names_only Logical. If set to TRUE only the names of the assessed ranking
#' are returned as a character vector. (Convenient to use as input for functions
#' taking gene set- or gene vectors as input.)
#'
#' @return A data.frame arranged by the residuals area-under-the-curve-values describing
#' how well a model fitted the expression trend of a gene or gene set.
#'
#' @export
#'

assessTrajectoryTrends <- function(rtdf, pattern = NULL, max_auc = NULL, names_only = FALSE){

  # 1. Control --------------------------------------------------------------

  check_rtdf(rtdf = rtdf)

  if(!base::is.null(pattern)){confuns::is_vec(pattern, "character", "pattern")}
  if(!base::is.null(max_auc)){confuns::is_value(max_auc, "numeric", "max_auc")}

  confuns::is_value(names_only, mode = "logical", "names_only")

  # -----

  # 2. Data wrangling -------------------------------------------------------

  arranged_df <-
    dplyr::select(.data = rtdf, -data, - models, - residuals) %>%
    tidyr::unnest(cols = dplyr::all_of("auc")) %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with("p_"),
      names_to = "pattern",
      names_prefix = "p_",
      values_to = "auc"
    ) %>%
    dplyr::arrange(auc) %>%
    dplyr::filter(auc <= base::ifelse(base::is.numeric(max_auc),
                                     yes = max_auc,
                                     no = base::max(auc))) %>%
    dplyr::mutate(
      pattern = hlpr_name_models(pattern)
    )


  # filter
  if(!base::is.null(pattern)){

    confuns::is_vec(x = pattern, mode = "character", ref = "pattern")

    flt_df <-
      dplyr::filter(arranged_df, pattern %in% {{pattern}})

    if(base::nrow(flt_df) == 0){

      base::stop("Unknown input for 'pattern'.")

    }

  } else {

    flt_df <- arranged_df

  }

  # -----

  if(base::isTRUE(names_only)){

    dplyr::pull(.data = flt_df, var = 1) %>%
    base::return()

  } else {

    base::return(flt_df)

  }

}

# -----


