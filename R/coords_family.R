#' Extract coordinates or positions
#'
#' @param object A valid object of class \emph{spata}.
#' @param of_sample The sample(s) name specified as a character.
#' @param of_segment The segment name specified as a character of length one.
#' @param of_trajectory The trajectory name specified as a character of length one.
#' @param method_dr The dimensional reduction method from which to extract the
#' positions specified as a character of length one.

#' @return A data.frame that contains in all cases the unique identifiers
#' (keys): \emph{barcodes, sample}.
#'
#' Depending on the function called additional variables such as \emph{x, y,
#' umap1, umap2, tsne1 or tsne2}.
#'
#' @export

coordsSpatial <- function(object,
                          of_sample,
                          of_segment = NULL){


  # 1. Control valid input --------------------------------------------------

  # object
  validation(x = object)

  # sample
  of_sample <- check_sample(object, sample_input = of_sample, desired_length = 1)

  # segment
  if(!is.null(of_segment) && !is.character(of_segment) && length(of_segment) != 1){

    base::stop("Argument 'of_segment' needs to be either NULL or a character vector of length one.")

  } else if(!is.null(of_segment)){

    segm_bc <-
      featureData(object, of_sample = of_sample) %>%
      dplyr::filter(segment == of_segment) %>%
      dplyr::pull(barcodes)

    if(base::length(segm_bc) == 0){

      base::stop(stringr::str_c("There is no segment of name' ", of_segment,
                                "' in sample '", of_sample, "'.", sep = ""))

    }

  }


  # 2. Data extraction ------------------------------------------------------

  coords_df <-
    object@coordinates %>%
    dplyr::filter(sample %in% of_sample)

  if(!is.null(of_segment)){

    coords_df %>%
      dplyr::filter(barcodes %in% segm_bc) %>%
      base::return()

  } else{

    base::return(coords_df)

  }

}

#' @rdname coordsSpatial
#' @export
coordsTrajectory <- function(object,
                             of_sample,
                             of_trajectory){


  # 1. Control valid input --------------------------------------------------

  # object
  validation(x = object)

  # sample
  of_sample <- check_sample(object, sample_input = of_sample, desired_length = 1)

  # trajectory
  if(!is.character(of_trajectory) && !base::length(of_trajectory) == 1){

    base::stop("Argument 'of_trajectory' needs to be a character vector of length one.")

  } else if(!of_trajectory %in% getTrajectoryNames(object, of_sample = of_sample)){

    base::stop(stringr::str_c("There is no trajectory of name '", of_trajectory,
                              "' in sample '", of_sample, "'.", sep = ""))

  }


  # 2. Data extraction ------------------------------------------------------

  t_object <-
    getTrajectoryObject(object,
                        trajectory_name = of_trajectory,
                        of_sample = of_sample)

  return(t_object@compiled_trajectory_df)

}

#' @rdname coordsSpatial
#' @export
coordsDimRed <- function(object,
                         of_sample,
                         method_dr = c("UMAP", "TSNE")){



  # 1. Control valid input --------------------------------------------------

  # object
  validation(x = object)

  # method
  if(base::length(method_dr) != 1){

    stop("Argument 'method_dr' needs to be of length one.")

  } else if(!method_dr %in% c("UMAP", "TSNE")) {

    stop("Argument 'method_dr' needs to be  'UMAP' or 'TSNE'.")

  } else {

    dr_strings <- stringr::str_c(base::tolower(x = method_dr), 1:2, sep = "")

  }

  # sample
  of_sample <- check_sample(object = object, sample_input = of_sample)

  # 2. Data extraction ------------------------------------------------------

  dim_red_df <-
    methods::slot(object = object@dim_red, name = method_dr) %>%
    dplyr::filter(sample %in% of_sample) %>%
    dplyr::select(dplyr::all_of(x = c("barcodes", "sample", dr_strings))) %>%
    tibble::remove_rownames()

  base::return(dim_red_df)

}


