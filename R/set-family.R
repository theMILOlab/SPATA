
# Slot: Information -------------------------------------------------------

#' @title
#'
#' @inherit check_object param
#' @param name Character value. The name of the matrix that is to be set as
#' the active expression matrix.
#'
#' @return
#' @export

setActiveMatrix <- function(object, name){

  check_object(object)
  confuns::is_value(x = name, mode = "character")

  # check if 'name' is slot in @data
  mtr_names <- base::names(object@data)
  confuns::check_one_of(input = name,
                        against = mtr_names[mtr_names != "counts"],
                        ref.input = "name")

  base::message(glue::glue("Active expression matrix set to '{name}'."))

  # set name
  object@information$active_mtr <- name

  base::return(object)

}


#' Title
#'
#' @inherit check_object params
#' @param info_list A named list with slots \code{$activation, $assesment, $bottleneck, $dropout, $epochs, $layers}.
#'
#' @return A spata-object.
#' @export
#'

setAutoencoderInfo <- function(object, info_list){

  check_object(object)
  confuns::is_list(info_list)

  object@information$autoencoder <- info_list

  base::return(object)

}

#' Title
#'
#' @inherit check_object params
#' @param assessment_list Named list with slots \code{$df} and \code{$set_up}.
#'
#' @return A spata-object.
#' @export

setAutoencoderAssessment <- function(object, assessment_list){

  check_object(object)
  confuns::check_data_frame(
    df = assessment_list$df,
    var.class = list("activation" = c("character", "factor"),
                     "bottleneck" = c("character", "factor"),
                     "total_var" = c("numeric", "integer", "double")),
    ref = "assessment_list$df"
  )

  object@information$autoencoder$assessment <- assessment_list

  base::return(object)

}


# Slot: Data --------------------------------------------------------------

setDenoisedMatrix <- function(object, denoised_mtr){

  check_object(object)

  object@data$denoised <- denoised_mtr

  base::return(object)

}

setScaledMatrix <- function(object, scaled_mtr){

  check_object(object)

  object@data$scaled <- scaled_mtr

  base::return(object)

}

addExpressionMatrix <- function(object, expr_mtr, name){

  check_object(object)

  object@data[[name]] <- expr_mtr

  base::return(object)

}




