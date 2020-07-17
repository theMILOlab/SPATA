#' @title Spatial Segmentation
#'
#' @description The function \code{createSegmentation()} provides access to an
#' interactive mini-shiny application that allows to separate a sample into
#' several segments.
#'
#' @param object A valid spata-object.
#'
#' @return An updated version of the spata-object specified as \code{object}
#' now containing the information about all drawn segments.
#'

createSegmentation <- function(object = NULL){

  validation(x = object)

  ##----- launch application
  new_object <-
    shiny::runApp(shiny::shinyApp(ui = shiny_ui_spatial_segmentation,
                                  server = shiny_server_spatial_segmentation)
    )

  return(new_object)

}



#' @title Spatial Trajectories
#'
#' @description The function \code{createTrajectories()} provides access to an
#' interactive mini-shiny application that allows to draw trajectories.
#'
#' @param object A valid spata-object.
#'
#' @return An updated version of the spata-object specified as \code{object}
#' now containing the information about all drawn trajectories.
#' @export
#'
createTrajectories <- function(object){

  validation(x = object)

  new_object <-
    shiny::runApp(
      shiny::shinyApp(ui = shiny_ui_spatial_trajectories,
                      server = shiny_server_spatial_trajectories)
    )

  return(new_object)

}
