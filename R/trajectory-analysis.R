
# Objects -----------------------------------------------------------------

#' @title Trajectory patterns
#'
#' @description Character vectors containing the names of valid trajectory patterns.
#'
#' @export
trajectory_patterns <- c("Linear descending", "Linear ascending", "Gradient descending", "Logarithmic descending",
                         "Logarithmic ascending", "Gradient ascending","Sinus",  "Sinus (reversed)", "One peak",
                         "One peak (reversed)", "Two peaks (reversed)", "Two peaks", "Early peak", "Late peak")
#' @export
linear_trends <- c("Linear descending", "Linear ascending")

#' @export
gradient_trends <- c("Gradient descending", "Gradient ascending")

#' @export
peak_trends <- c("One peak", "Late peak", "Early peak")

#' @export
logarithmic_trends <- c("Logarithmic descending", "Logarithmic ascending")

# -----



# Helper functions --------------------------------------------------------

#' @title Rank trajectory trends.
#'
#' @description Analyze the expression dynamics along
#' a specified trajectory by fitting a variety of models to the genes or
#' gene sets expression trends.
#'
#' @param stdf A summarized trajectory data.frame. (e.g. obtained by
#' \code{getTrajectoryDf()}).
#'
#' @return A nested data.frame with information about the dynamics of each gene
#' or gene set.
#'
#' @export

hlpr_rank_trajectory_trends <- function(stdf, verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  check_stdf(stdf)

  var <- "variables"

  if(base::length(var) != 1){

    base::stop("Data.frame 'stdf' must have one column either name 'genes' or 'gene_sets'.")

  }

  # -----

  # 2. Ranking --------------------------------------------------------------

  # nest data.frame
  nested_df <-
    dplyr::group_by(.data = stdf, !!rlang::sym(var)) %>%
      dplyr::mutate(values = confuns::normalize(x = values)) %>%
      tidyr::nest()


  # add residuals to data.frame
  if(base::isTRUE(verbose)){

    base::message("Fitting models.")

    pb_add <- progress::progress_bar$new(
      format = "Progress: [:bar] :percent eta: :eta",
      total = base::nrow(nested_df),
      clear = FALSE,
      width = 100
    )

  } else {

    pb_add <- NULL

  }

  w_residuals <-
    dplyr::mutate(.data = nested_df,
                  residuals = purrr::map(.x = data, .f = hlpr_add_residuals, pb = pb_add))

  # rank data.frame
  if(base::isTRUE(verbose)){

    base::message("Calculating residuals.")

    pb_calc <- progress::progress_bar$new(
      format = "Progress: [:bar] :percent eta: :eta",
      total = base::nrow(nested_df),
      clear = FALSE,
      width = 100
    )

  } else {

    pb_calc <- NULL

  }

  ranked_df <-
    dplyr::mutate(.data = w_residuals,
                  auc = purrr::map(.x = residuals, .f = hlpr_summarize_residuals, pb = pb_calc))

  # -----

  if(base::isTRUE(verbose)){base::message("Done.")}

  return(ranked_df)

}



#' @title Assess trajectory ranking.
#'
#' @description Takes a ranked trajectory data.frame and returns a data.frame
#' that informs about how well the ranked gene- or gene set expression-trends
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

hlpr_assess_trajectory_trends <- function(rtdf, verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  check_rtdf(rtdf = rtdf)

  # -----

  # 2. Data wrangling -------------------------------------------------------

  if(base::isTRUE(verbose)){base::message("Asessing trajectory trends.")}

  arranged_df <-
    dplyr::select(.data = rtdf, -data, -residuals) %>%
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

  if(base::isTRUE(verbose)){base::message("Done.")}

  base::return(arranged_df)

}


#' @title Filter variables of a certain trend
#'
#' @description Extracts the genes or gene sets that follow a desired trend.
#'
#' @param atdf An assessed trajectory data.frame (easily accessed via
#' \code{assessTrajectoryTrends()}).
#' @param limit Numeric value. The maximum area-under-the-curve value the
#' trajectory-trend-assessment might have.
#' @param trend Character vector. The patterns of interest.
#' @param variables_only Logical. If set to TRUE a character of variable-names is returned.
#' If set to FALSE the filtered data.frame is returned.
#'
#' @return A character vector of gene or gene-set names that follow the specified
#' patterns to the specified degree.
#' @export

filterTrends <- function(atdf, limit = 5, trends = "all", variables_only = TRUE){

  if(base::all(trends == "all")){

    trends <- trajectory_patterns

  }

  confuns::is_vec(x = trends, mode = "character", "trends")
  trends <- confuns::check_vector(input = trends,
                        against = trajectory_patterns,
                        verbose = TRUE,
                        ref.input = "argument 'trends'",
                        ref.against = "known trajectory trends")

  if(base::isTRUE(variables_only)){

    res <-
      hlpr_filter_trend(atdf = atdf,
                        limit = limit,
                        poi = trends) # poi = patterns of interest

  } else {
    res <-
      dplyr::filter(.data = atdf, pattern %in% trends) %>%
      dplyr::filter(auc <= limit) %>%
      dplyr::group_by(variables) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(pattern) %>%
      dplyr::arrange(auc, .by_group = TRUE)

  }

  base::return(res)

}

# -----



# Main functions ----------------------------------------------------------

#' @title Trajectory trend analysis
#'
#' @description Analyzes the trend of gene and gene-set-expressions along
#' trajectories by fitting a variety of mathematical models to them and
#' by assessing the quality of the fit.
#'
#' \itemize{
#'  \item{\code{assessTrajectoryTrends()}: Takes a valid spata-object and constructs
#'  the subsequent data.frame from scratch.}
#'  \item{\code{assessTrajectoryTrends2()}: Takes a summarized trajectory data.frame
#'  returned by \code{getTrajectoryDf()}.}
#'  }
#'
#' @inherit check_sample params
#' @inherit check_trajectory params
#' @inherit check_variables params
#' @inherit verbose params
#' @inherit hlpr_rank_trajectory_trends params
#'
#' @return A data.frame arranged by the residuals area-under-the-curve-values describing
#' how well a model fitted the expression trend of a gene or gene set.
#'
#' @export

assessTrajectoryTrends <- function(object,
                                   trajectory_name,
                                   of_sample = "",
                                   variables,
                                   accuracy = 5,
                                   verbose = TRUE){


  # 1. Control --------------------------------------------------------------

  check_object(object)

  of_sample <- check_sample(object, of_sample, desired_length = 1)

  check_trajectory(object, trajectory_name, of_sample)

  # -----


  # 2. Main part ------------------------------------------------------------

  # get trajectory data.frame

    stdf <- getTrajectoryDf(object = object,
                           trajectory_name = trajectory_name,
                           of_sample = of_sample,
                           variables = variables,
                           accuracy = accuracy,
                           verbose = verbose)

  rtdf <- hlpr_rank_trajectory_trends(stdf = stdf, verbose = verbose)

  atdf <- hlpr_assess_trajectory_trends(rtdf = rtdf, verbose = verbose)

  # -----

  base::return(atdf)

}

#' @rdname assessTrajectoryTrends
#' @export
assessTrajectoryTrends2 <- function(stdf, verbose = TRUE){

  # 2. Main part ------------------------------------------------------------

  check_stdf(stdf = stdf)

  rtdf <- hlpr_rank_trajectory_trends(stdf = stdf, verbose = verbose)

  atdf <- hlpr_assess_trajectory_trends(rtdf = rtdf, verbose = verbose)

  # -----

  base::return(atdf)

}

# -----

