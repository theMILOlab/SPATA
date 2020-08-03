#' @title Extract coordinates
#'
#' @description Extracts spatial information of whole samples, trajectories,
#' segments etc. in form of coordinates. (See 'Value').
#'
#' @param object A valid object of class \emph{spata}.
#' @inherit check_sample params
#' @inherit check_segment params
#' @inherit check_trajectory params
#' @inherit check_method params
#'
#' @return A data.frame that contains the unique identifiers
#' (keys): \emph{barcodes, sample} and the coordinates.
#'
#'  \itemize{
#'   \item{ \code{coordsSpatial()}:  \emph{x, y}}
#'   \item{ \code{coordsSegment()}: \emph{x, y}}
#'   \item{ \code{coordsTrajectory()}: \emph{x, y}}
#'   \item{ \code{coordsTSNE()}: \emph{tsne1, tsne2}}
#'   \item{ \code{coordsUMAP()}: \emph{umap1, umap2}}
#'   }.
#'
#' @export

coordsSpatial <- function(object,
                          of_sample){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)

  # adjusting check
  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  # -----

  # 2. Data extraction ------------------------------------------------------

  coords_df <-
    object@coordinates %>%
    dplyr::filter(sample %in% of_sample)

  # -----

  base::return(coords_df)


}

#' @rdname coordsSpatial
#' @export
coordsSegment <- function(object,
                          of_sample,
                          segment_name){

  # 1. Control --------------------------------------------------------------

  # adjusting check
  bc_segm <- check_segment(object = object,
                           segment_name = segment_name)

  # -----

  # 2. Data extraction ------------------------------------------------------

  coords <-
    coordsSpatial(object = object,of_sample = of_sample) %>%
    dplyr::filter(barcodes %in% bc_segm) %>%
    dplyr::mutate(segment_name = {{segment_name}})

  # -----

  base::return(coords)


}

#' @rdname coordsSpatial
#' @export
coordsTrajectory <- function(object,
                             of_sample,
                             trajectory_name){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_trajectory(object, trajectory_name, of_sample)

  # adjusting check
  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  # -----

  # 2. Data extraction ------------------------------------------------------

  t_object <-
    trajectory(object, trajectory_name = trajectory_name, of_sample = of_sample)

  # -----

  return(t_object@compiled_trajectory_df)

}

#' @rdname coordsSpatial
#' @export
coordsDimRed <- function(object,
                         of_sample,
                         method_dr = c("UMAP", "TSNE")
                         ){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_method(method_dr = method_dr)

  # adjusting check
  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  # -----

  # 2. Data extraction ------------------------------------------------------

  dr_strings <- stringr::str_c(base::tolower(x = method_dr), 1:2, sep = "")

  dim_red_df <-
    methods::slot(object = object@dim_red, name = method_dr) %>%
    dplyr::filter(sample %in% of_sample) %>%
    dplyr::select(dplyr::all_of(x = c("barcodes", "sample", dr_strings))) %>%
    tibble::remove_rownames()

  # -----

  if(base::nrow(dim_red_df) == 0){

    base::stop("There seems to be no data for method: ", method_dr)

  }

  base::return(dim_red_df)

}

#' @rdname coordsSpatial
#' @export
coordsUMAP <- function(object,
                       of_sample){

  coordsDimRed(object = object,
               of_sample = of_sample,
               method_dr = "UMAP")

}

#' @rdname coordsSpatial
#' @export
coordsTSNE <- function(object,
                       of_sample){

  coordsDimRed(object = object,
               of_sample = of_sample,
               method_dr = "TSNE")

}


