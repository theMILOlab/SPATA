
# Objects -----------------------------------------------------------------

trajectory_patterns <- c("Linear descending", "Linear ascending", "Gradient descending", "Logarithmic descending",
                         "Logarithmic ascending", "Gradient ascending","Sinus",  "Sinus (reversed)", "One peak",
                         "One peak (reversed)", "Two peaks (reversed)", "Two peaks")

# -----


# Functions ---------------------------------------------------------------


#' @title Trajectory expression analysis
#'
#' @description Analyze the expression dynamics along
#' a specified trajectory by fitting a variety of models to the genes or
#' gene sets expression trends.
#'
#' @param stdf A summarized trajectory data.frame. (e.g. obtained by
#' \code{getSummarizedTrajectoryDf()}).
#'
#' @return A nested data.frame with information about the dynamics of each gene
#' or gene set.
#'
#' @export


rankTrajectoryTrends <- function(stdf){

  # 1. Control --------------------------------------------------------------

  check_summarized_trajectory_df(stdf)

  var <- base::colnames(stdf)[base::colnames(stdf) %in% c("genes", "gene_sets")]

  if(base::length(var) != 1){

    base::stop("Data.frame 'stdf' must have one column of name 'genes' or 'gene_sets'.")

  }

  # -----

  # 2. Ranking --------------------------------------------------------------

  ranked_df <-
    dplyr::group_by(.data = stdf, !!rlang::sym(var)) %>%
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

assessTrajectoryTrends <- function(rtdf){

  # 1. Control --------------------------------------------------------------

  check_rtdf(rtdf = rtdf)

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
    dplyr::mutate(pattern = hlpr_name_models(pattern))

  # -----

  base::return(arranged_df)


}



#' @title Trajectory trends
#'
#' @description Extracts the trajectories of a desired trend.
#'
#' @param atdf An assessed trajectory data.frame.
#' @param limit Numeric value. The maximum area-under-the-curve value the
#' trajectory-trend-assessment might have.
#' @param pattern Character vector. The patterns of interest.
#'
#' @return A character vector of gene or gene-set names that follow the specified
#' patterns to the specified degree.
#' @export
#'

trajectoryTrend <- function(atdf,
                            limit = 2,
                            pattern){

  confuns::is_vec(x = pattern, mode = "character", "pattern")

  res <-
    hlpr_filter_trend(atdf = atdf,
                      limit = limit,
                      poi = pattern)

  base::return(res)

}


# -----


