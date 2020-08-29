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
#' @export

createSegmentation <- function(object = NULL){

  validation(x = object)

  ##----- launch application
  new_object <-
    shiny::runApp(
      shiny::shinyApp(ui = {shiny::fluidPage(

        ##----- title
        shiny::titlePanel(title = "Spatial Segmentation"),

        shinybusy::add_busy_spinner(spin = "cube-grid", margins = c(0, 10), color = "red"),

        ##----- main panel
        shiny::mainPanel(
          shiny::fluidRow(
            shiny::column(width = 4,
                          shiny::wellPanel(
                            shiny::tags$h3(shiny::strong("Instructions")),
                            shiny::HTML("<br>"),
                            shiny::helpText("1. Click on 'Plot & Update' to display the sample according to the adjustments you set up or changed."),
                            shiny::HTML("<br>"),
                            shiny::helpText("2. Determine the vertices of the segment by 'double - clicking' the position on the plot. Be aware of the order of the vertices."),
                            shiny::HTML("<br>"),
                            shiny::helpText("3. Highlight or reset the segment by clicking the respective button below."),
                            shiny::splitLayout(
                              shiny::actionButton("highlight_segment", label = "Highlight", width = "100%"),
                              shiny::actionButton("reset_segment", label = "Reset ", width = "100%"),
                              cellWidths = c("50%", "50%")
                            ),
                            shiny::HTML("<br><br>"),
                            shiny::helpText("4. Enter the name you want to give the highlighted segment and click the 'Save'-button."),
                            shiny::splitLayout(
                              shiny::actionButton("save_segment", "Save Segment", width = "100%"),
                              shiny::textInput("name_segment", label = NULL, placeholder = "Name segment", value = ""),
                              cellWidths = c("50%", "50%")
                            ),
                            shiny::HTML("<br>"),
                            shiny::helpText("5. If you want to remove certain segments type in the respective name and click the 'Remove'-button.
                                          Choose 'Feature' as color code and then 'segment' to see the current segmentation."),
                            shiny::splitLayout(
                              shiny::actionButton("remove_segment", "Remove Segment", width = "100%"),
                              shiny::textInput("name_segment_rmv", label = NULL, placeholder = "Name segment", value = "")
                            ),
                            shiny::HTML("<br>"),
                            shiny::helpText("6. If you are done click on 'Close application'."),
                            shiny::HTML("<br>"),
                            shiny::fluidRow(
                              shiny::column(width = 12, align = "center",
                                            shiny::actionButton("close_app", label = "Close application", width = "50%")
                              )
                            )
                          )),
            shiny::column(width = 8,
                          moduleSurfacePlotUI(id = "segmentation")
            )

          )
        )

      )},
      server = function(input, output, session){

        # Reactive values -----------------------------------------------------------

        # a reactive spata object
        spata_obj <- shiny::reactiveVal(value = object)

        # df and ggplot layer of the currently drawn segment
        vertices_df <-
          shiny::reactiveVal(value = data.frame(x = numeric(0),
                                                y = numeric(0)))

        vertices_layer <- shiny::reactiveVal(value = list())

        # a list about the parameters of the currently displayed surface plot
        current <- reactiveVal(value = list())

        #
        segmentation_df <- reactive({

          segm_df <-
            featureData(object = spata_obj()) %>%
            dplyr::filter(sample == current()$sample) %>%
            dplyr::filter(segment != "") %>%
            dplyr::select(barcodes, segment)

          return(segm_df)

        })


        # Modularized plot surface part -------------------------------------------

        module_return <- moduleSurfacePlotServer(id = "segmentation",
                                                 object = object,
                                                 final_plot = shiny::reactive(final_plot()),
                                                 reactive_object = shiny::reactive(spata_obj()))

        # update current()
        oe <- shiny::observeEvent(module_return()$current_setting(), {

          current(module_return()$current_setting())

        })

        # final plot
        final_plot <- shiny::reactive({

          module_return()$assembled_plot() +
            vertices_layer()

        })

        # Observe events ----------------------------------------------------------

        ##--- 1. grow vertices data and update vertices layer frame with every click
        oe <- shiny::observeEvent(module_return()$dblclick(), {

          ## 1. computation
          vrtcs_list <- module_return()$dblclick()
          new_df <- dplyr::add_row(.data = vertices_df(),
                                   x = vrtcs_list$x,
                                   y = vrtcs_list$y)

          ## 2.1 update vertices df
          vertices_df(new_df)

          ## 2.2 update vertices geom layer
          if(nrow(vertices_df()) != 0){

            new_layer <- list(ggplot2::geom_point(data = vertices_df(), mapping = ggplot2::aes(x = x, y = y), size = 3.5, color = "black"),
                              ggplot2::geom_path(data = vertices_df(), mapping = ggplot2::aes(x = x, y = y), size = 1, color = "black"))

            vertices_layer(new_layer)

          } else {

            new_layer <- NULL
            vertices_layer(list())

          }

        })

        ##--- 2.1 convert vertices layer to geom_polygon to highlight the segment
        oe <- shiny::observeEvent(input$highlight_segment, {

          checkpoint(evaluate = base::nrow(vertices_df()) > 2, case_false = "insufficient_n_vertices")

          new_layer <- list(ggplot2::geom_polygon(data = vertices_df(),
                                                  mapping = ggplot2::aes(x = x, y = y),
                                                  alpha = 0.75, colour = "orange", fill = "orange",
                                                  size = 1))
          vertices_layer(new_layer)

        })

        ##--- 2.2 reset current() vertices
        oe <- shiny::observeEvent(input$reset_segment, {

          vertices_df(data.frame(x = numeric(0), y = numeric(0)))
          vertices_layer(list())

        })

        ##--- 3. save the highlighted segment
        oe <- shiny::observeEvent(input$save_segment, {

          checkpoint(evaluate = input$name_segment != "", case_false = "invalid_segment_name")
          checkpoint(evaluate = !input$name_segment %in% segmentation_df()$segment, case_false = "occupied_segment_name")
          checkpoint(evaluate = nrow(vertices_df()) > 2, case_false = "insufficient_n_vertices")

          sample_coords <- coordinates(objec = spata_obj(), of_sample = current()$sample)

          ## 1. determine positions of each point with respect to the defined segment
          positions <-  sp::point.in.polygon(point.x = sample_coords$x, # x coordinates of all spatial positions
                                             point.y = sample_coords$y, # y coordaintes of all spatial positions
                                             pol.x = vertices_df()$x, # x coordinates of the segments vertices
                                             pol.y = vertices_df()$y) # y coordinates of the segments vertices

          ## 2. update spata obj

          # 2.1 extract object
          spata_obj <- spata_obj()

          # 2.2 update fdata

          # extract feature data
          fdata <- featureData(spata_obj)

          # update sample subset
          fdata_subset <-
            fdata %>%
            dplyr::filter(sample == current()$sample) %>%
            dplyr::mutate(
              positions = positions,
              segment = dplyr::if_else(condition = positions %in% c(1,2,3), true = input$name_segment, false = segment)
            ) %>%
            dplyr::select(-positions)

          # exchange sample subset
          fdata[fdata$sample == current()$sample, ] <- fdata_subset

          featureData(spata_obj) <- fdata


          # 2.4 update and check
          spata_obj(spata_obj)

          if(input$name_segment %in% base::unique(featureData(spata_obj(), of_sample = current()$sample)$segment)){

            shiny::showNotification(ui = stringr::str_c(input$name_segment, "has been saved.", sep = " "), type = "message")

          }


          ## 3. reset vertices values
          vertices_df(data.frame(x = numeric(0), y = numeric(0)))
          vertices_layer(list())

        })

        ##--- 4. remove segments
        oe <- shiny::observeEvent(input$remove_segment, {

          spata_obj <- spata_obj()
          fdata <- featureData(spata_obj)

          checkpoint(evaluate = input$name_segment_rmv %in% base::unique(fdata$segment), case_false = "segment_name_not_found")

          fdata_new <-
            fdata %>%
            dplyr::filter(sample == current()$sample) %>%
            dplyr::mutate(segment = dplyr::if_else(segment == input$name_segment_rmv, "", segment))

          fdata[fdata$sample == current()$sample, ] <- fdata_new
          featureData(spata_obj) <- fdata

          spata_obj(spata_obj)

          if(!input$name_segment_rmv %in% featureData(spata_obj(), of_sample = current()$sample)$segment){

            shiny::showNotification(ui = stringr::str_c("Segment '", input$name_segment_rmv, "' has been successfully removed.", sep = ""), type = "message")

          }

        })

        ##--- 5. close application and return spata object
        oe <- shiny::observeEvent(input$close_app, {

          shiny::stopApp(returnValue = spata_obj())

        })

      })
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
      shiny::shinyApp(
        ui = {shiny::fluidPage(

          ##----- title
          shiny::titlePanel(title = "Spatial Trajectories"),

          shinybusy::add_busy_spinner(spin = "cube-grid", margins = c(0, 10), color = "red"),

          ##----- main panel
          shiny::mainPanel(
            shiny::fluidRow(
              shiny::column(width = 4,
                            shiny::wellPanel(
                              shiny::tags$h3(shiny::strong("Instructions")),
                              shiny::HTML("<br>"),
                              shiny::helpText("1. Click on 'Plot & Update' to display the sample according to the adjustments you set up or changed."),
                              shiny::HTML("<br>"),
                              shiny::helpText("2. Determine the vertices of the trajectory by 'double - clicking' the position on the plot."),
                              shiny::HTML("<br>"),
                              shiny::helpText("3. Highlight or reset the trajectory by clicking the respective button below."),
                              shiny::sliderInput("trajectory_width", label = "Determine width of trajectory", value = 20, min = 5, max = 100, step = 1),
                              shiny::HTML("<br>"),
                              shiny::splitLayout(
                                shiny::actionButton("highlight_trajectory", label = "Highlight", width = "100%"),
                                shiny::actionButton("reset_trajectory", label = "Reset ", width = "100%"),
                                cellWidths = c("50%", "50%")
                              ),
                              shiny::HTML("<br>"),
                              shiny::helpText("4. Enter the name you want to give the trajectory as well as a 'guiding comment' and click the 'Save'-button."),
                              shiny::splitLayout(
                                shiny::actionButton("save_trajectory", "Save Trajectory", width = "100%"),
                                shiny::textInput("name_trajectory", label = NULL, placeholder = "Name trajectory", value = ""),
                                cellWidths = c("50%", "50%")
                              ),
                              shiny::textInput("comment_trajectory", label = NULL, placeholder = "A guiding comment.", value = ""),
                              shiny::HTML("<br>"),
                              shiny::helpText("5. If you are done click on 'Close application'."),
                              shiny::HTML("<br>"),
                              shiny::fluidRow(
                                shiny::column(width = 12, align = "center",
                                              shiny::actionButton("close_app", label = "Close application", width = "50%")
                                )
                              )
                            )),
              shiny::column(width = 8,
                            moduleSurfacePlotUI(id = "trajectories")
              )

            )
          )

        )},
        server = function(input, output, session){


          # Reactive values ---------------------------------------------------------
          spata_obj <- shiny::reactiveVal(value = object)

          vertices_df <-
            shiny::reactiveVal(value = data.frame(x = numeric(0),
                                                  y = numeric(0)))

          segment_trajectory_df <- shiny::reactiveVal(
            value = data.frame(x = numeric(0),
                               y = numeric(0),
                               xend = numeric(0),
                               yend = numeric(0),
                               part = character(0))
          )

          compiled_trajectory_df <- shiny::reactiveVal(
            value = data.frame(barcodes = character(0),
                               sample = character(0),
                               x = numeric(0),
                               y = numeric(0),
                               projection_length = numeric(0),
                               trajectory_part = character(0),
                               stringsAsFactors = F)
          )

          current <- shiny::reactiveVal(value = list())

          # -----

          # Modularized plot surface part -------------------------------------------

          module_return <- moduleSurfacePlotServer(id = "trajectories",
                                                   object = object,
                                                   final_plot = shiny::reactive(final_plot()),
                                                   reactive_object = shiny::reactive(spata_obj()))

          # update current()
          oe <- shiny::observeEvent(module_return()$current_setting(), {

            current(module_return()$current_setting())

          })

          # final plot
          final_plot <- shiny::reactive({

            module_return()$assembled_plot() +
              trajectory_point_add_on() +
              trajectory_segment_add_on()

          })



          # trjectory add ons
          trajectory_segment_add_on <- shiny::reactive({

            new_layer <- list()

            # update geom_point layer
            if(nrow(vertices_df()) >= 1){

              new_layer[[1]] <-
                ggplot2::geom_point(data = vertices_df(),
                                    mapping = ggplot2::aes(x = x, y = y),
                                    size = 3.5, color = "black")

            }

            # update geom_segment layer
            if(base::nrow(segment_trajectory_df()) >= 1){

              new_layer[[2]] <-
                ggplot2::geom_segment(data = segment_trajectory_df(),
                                      mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                                      size = 1.25, color = "black",
                                      arrow = ggplot2::arrow(length = ggplot2::unit(0.125, "inches"))
                )

            }

            return(new_layer)

          })

          # highlight points of trajectory
          trajectory_point_add_on <- shiny::reactive({

            if(!base::nrow(compiled_trajectory_df()) == 0){

              joined_traj_df <-
                dplyr::left_join(x = compiled_trajectory_df(),
                                 y = dplyr::select(module_return()$smoothed_df(), -x, -y),
                                 by = "barcodes")

              color_var <- dplyr::pull(.data = joined_traj_df, module_return()$variable())
              size <- module_return()$current_setting()$size

              add_on_layer <-
                list(
                  #ggplot2::geom_point(data = module_return()$smoothed_df(), size = size,
                  #                    mapping = ggplot2::aes(x = x, y = y), alpha = current()$pt_alpha, color = "lightgrey"),
                  ggplot2::geom_point(data = joined_traj_df, size = size,
                                      mapping = ggplot2::aes(x = x, y = y, color = color_var))
                )

            } else {

              add_on_layer <- list()

            }

            return(add_on_layer)

          })

          # -----


          # Observe events and reactive events --------------------------------------

          # 1. add trajectory vertice consecutively
          oe <- shiny::observeEvent(module_return()$dblclick(), {

            # 1. prolong and update data.frame
            vrtcs_list <- module_return()$dblclick()
            new_df <- dplyr::add_row(.data = vertices_df(),
                                     x = vrtcs_list$x,
                                     y = vrtcs_list$y)

            vertices_df(new_df)

            # 2. update trajectory df
            n_vrt <- nrow(vertices_df())

            if(n_vrt >= 2){

              stdf <-
                segment_trajectory_df() %>%
                dplyr::add_row(
                  x = base::as.numeric(vertices_df()[(n_vrt-1), 1]),
                  y = base::as.numeric(vertices_df()[(n_vrt-1), 2]),
                  xend = base::as.numeric(vertices_df()[(n_vrt), 1]),
                  yend = base::as.numeric(vertices_df()[(n_vrt), 2]),
                  part = stringr::str_c("part", n_vrt-1 , sep = "_")
                )

              print(stdf)
              segment_trajectory_df(stats::na.omit(stdf))

            } else {

              segment_trajectory_df(data.frame(
                x = numeric(0),
                y = numeric(0),
                xend = numeric(0),
                yend = numeric(0),
                part = character(0),
                stringsAsFactors = F))

            }

          })

          # 2.1
          oe <- shiny::observeEvent(input$highlight_trajectory, {

            checkpoint(evaluate = nrow(segment_trajectory_df()) >= 1, case_false = "insufficient_n_vertices2")

            compiled_trajectory_df <-
              hlpr_compile_trajectory(segment_trajectory_df = segment_trajectory_df(),
                                      trajectory_width = input$trajectory_width,
                                      object = spata_obj(),
                                      sample = current()$sample)

            compiled_trajectory_df(compiled_trajectory_df)

          })


          # 2.2 reset current() vertices
          oe <- shiny::observeEvent(input$reset_trajectory, {

            vertices_df(data.frame(x = numeric(0),
                                   y = numeric(0)))

            segment_trajectory_df(data.frame(x = numeric(0),
                                             y = numeric(0),
                                             xend = numeric(0),
                                             yend = numeric(0),
                                             part = character(0),
                                             stringsAsFactors = F))

            compiled_trajectory_df(data.frame(barcodes = character(0),
                                              sample = character(0),
                                              x = numeric(0),
                                              y = numeric(0),
                                              projection_length = numeric(0),
                                              trajectory_part = character(0),
                                              stringsAsFactors = F))

          })

          ##--- 3. save the highlighted trajectory
          oe <- shiny::observeEvent(input$save_trajectory, {

            ## control
            checkpoint(evaluate = base::nrow(compiled_trajectory_df()) > 0, case_false = "insufficient_n_vertices2")
            checkpoint(evaluate = shiny::isTruthy(x = input$name_trajectory), case_false = "invalid_trajectory_name")
            checkpoint(evaluate = !input$name_trajectory %in% getTrajectoryNames(spata_obj(), current()$sample),
                       case_false = "occupied_trajectory_name")

            ## save trajectory
            spata_obj <- spata_obj()

            spata_obj@trajectories[[current()$sample]][[input$name_trajectory]] <-
              methods::new("spatialTrajectory",
                           compiled_trajectory_df = compiled_trajectory_df(),
                           segment_trajectory_df = segment_trajectory_df(),
                           comment = input$comment_trajectory,
                           name = input$name_trajectory,
                           sample = current()$sample)

            spata_obj(spata_obj)


            ## control
            check <- trajectory(spata_obj(), trajectory_name = input$name_trajectory, of_sample = current()$sample)

            ## feedback and reset
            if(base::identical(check@compiled_trajectory_df, compiled_trajectory_df())){

              shiny::showNotification(ui = "Trajectory has been stored.", type = "message", duration = 7)


              vertices_df(data.frame(x = numeric(0),
                                     y = numeric(0)))

              segment_trajectory_df(data.frame(x = numeric(0),
                                               y = numeric(0),
                                               xend = numeric(0),
                                               yend = numeric(0),
                                               part = character(0),
                                               stringsAsFactors = F))

              compiled_trajectory_df(data.frame(barcodes = character(0),
                                                sample = character(0),
                                                x = numeric(0),
                                                y = numeric(0),
                                                projection_length = numeric(0),
                                                trajectory_part = character(0),
                                                stringsAsFactors = F))

            } else {

              shiny::showNotification(ui = "Could not save trajectory.")

            }

          })

          ##--- 5. close application and return spata object
          oe <- shiny::observeEvent(input$close_app, {

            stopApp(returnValue = spata_obj())

          })

        }))

  return(new_object)

}



#' @title Monocle3 Pseudotime
#'
#' @description Calculates the pseudotime values for every barcode in a given sample
#' and adds it under the specified name to the feature data in the provided
#' spata-object.
#'
#' @param object A valid spata object.
#' @param use_cds_file A directory leading to a .rds file containing a valid
#' cell_data_set-object previously calculated for the specified object. Specified
#' as a character value. If set to FALSE the cell_data_set object will be created
#' from scratch (which requires interactive node-choosing).
#' @param save_cds_file A filename/directory (that does not already exists) under which the used or created cell_data_set-object
#' is stored specified as a character value. Should end with \emph{'.rds'} or  \emph{'.RDS'}.
#' @param preprocess_method Given to \code{monocle3::preprocess_cds()} if \code{use_cds_file} isn't a character string.
#' @param cluster_method Given to \code{monocle3::cluster_cells()} if \code{use_cds_file} isn't a character string.
#' @param feature_name The name under which the created pseudotime-variable is stored in the provided object. Will overwrite
#' already existing features of the same name!
#' @param verbose Logical value. If set to TRUE informative messages with respect
#' to the computational progress made will be printed.
#'
#' (Warning messages will always be printed.)
#'
#' @return An updated spata-object.
#' @export
#'

createPseudotime <- function(object,
                             use_cds_file = FALSE,
                             save_cds_file = FALSE,
                             preprocess_method = "PCA",
                             cluster_method = "leiden",
                             feature_name = "pseudotime",
                             verbose = TRUE){

  check_object(object)

  cds <-
    hlpr_compile_cds(object = object,
                     use_cds_file = use_cds_file,
                     save_cds_file = save_cds_file,
                     preprocess_method = preprocess_method,
                     cluster_method = cluster_method,
                     verbose = verbose)

  ps_time <-
    as.data.frame(monocle3::pseudotime(cds)) %>%
    magrittr::set_colnames(value = feature_name) %>%
    tibble::rownames_to_column(var = "barcodes")

  ps_time[base::is.infinite(ps_time[,feature_name]), feature_name] <- NA

  object <- addFeature(object, feature_df = ps_time, feature_name = feature_name, overwrite = TRUE)

  base::return(object)

}





















