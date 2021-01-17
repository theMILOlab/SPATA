
#' @title Plot the sample
#'
#' @description Displays the spatial dimension of the sample and colors the
#' surface according to the expression of genes, gene sets or features.
#'
#' \itemize{
#'
#'  \item{ \code{plotSurface()} Takes the spata-object as the starting point and creates the
#'  necessary data.frame from scratch according to additional parameters.}
#'  \item{ \code{plotSurface2()} Takes a data.frame as the starting point.}
#'  \item{ \code{plotSurfaceInteractive()} Takes only the spata-object and opens a shiny
#'  application which allows for interactive plotting.}
#'
#' }
#'
#' @inherit check_coords_df params
#' @inherit check_sample params
#' @inherit check_color_to params
#' @inherit check_smooth params
#' @inherit check_method params
#' @inherit check_pt params
#' @inherit check_display params
#' @inherit check_assign params
#' @inherit verbose params
#' @inherit image_dummy params
#' @param complete Logical. If the provided spata-object has been subsetted by
#'  \code{subsetBySegment()} the original sample is completed with grey barcode
#'  spots.
#'
#' @inherit plot_family params
#'
#' @export

plotSurface <- function(object,
                        color_to = NULL,
                        method_gs = NULL,
                        normalize = NULL,
                        smooth = NULL,
                        smooth_span = NULL,
                        pt_alpha = NULL,
                        pt_clr = NULL,
                        pt_clrp = NULL,
                        pt_clrsp = NULL,
                        pt_size = NULL,
                        display_image = NULL,
                        display_title = NULL,
                        complete = NULL,
                        verbose = NULL,
                        of_sample = NA,
                        ...){

  # 1. Control --------------------------------------------------------------

  # lazy check
  hlpr_assign_arguments(object)


  check_pt(pt_size, pt_alpha, pt_clrsp)
  check_display(display_title, display_image)

  # adjusting check
  of_sample <- check_sample(object = object,
                            of_sample = of_sample,
                            desired_length = 1)


  if(!base::is.null(color_to)){

    color_to <- check_color_to(
      color_to = color_to,
      all_genes = getGenes(object, of_sample = of_sample),
      all_gene_sets = getGeneSets(object),
      all_features = getFeatureNames(object, of_sample = of_sample)
      )

  } else {

    color_to <- list("color" = pt_clr)

  }


  # -----

  # 2. Data extraction and plot preparation ---------------------------------

  coords_df <- getCoordsDf(object, of_sample = of_sample)

  plot_list <-
    hlpr_scatterplot(object = object,
                     spata_df = coords_df,
                     color_to = color_to,
                     pt_size = pt_size,
                     pt_alpha = pt_alpha,
                     pt_clrp = pt_clrp,
                     pt_clrsp = pt_clrsp,
                     method_gs = method_gs,
                     normalize = normalize,
                     smooth = smooth,
                     smooth_span = smooth_span,
                     verbose = verbose,
                     complete = complete,
                     ...)

  # -----


  # 5. Plotting --------------------------------------------------------------

  ggplot2::ggplot(data = plot_list$data, mapping = ggplot2::aes(x = x, y = y)) +
    hlpr_image_add_on(object, display_image, of_sample) +
    plot_list$add_on +
    ggplot2::coord_equal() +
    ggplot2::theme_void()


  # -----

}


#' @rdname plotSurface
#' @export
plotSurface2 <- function(coords_df,
                         color_to,
                         pt_size = 2,
                         pt_alpha = 0.9,
                         pt_clrsp = "inferno",
                         pt_clrp = "milo",
                         image = NULL,
                         ...){

   # 1. Control --------------------------------------------------------------
  confuns::is_value(color_to, "character", "color_to")

  # check lazy
  if(!color_to %in% base::colnames(coords_df)){

    base::stop("Argument 'color_to' needs be a variable of 'coords_df' specified as a character value.")

  }

  check_pt(pt_size, pt_alpha, pt_clrsp)
  check_coords_df(coords_df)

  # -----

  # 2. Plotting -------------------------------------------------------------


  ggplot2::ggplot(data = coords_df) +
    hlpr_image_add_on2(image) +
    ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = y,
                                               color = .data[[color_to]]),
                        size = pt_size, alpha = pt_alpha) +
    confuns::scale_color_add_on(
      clrp = pt_clrp,
      clrsp = pt_clrsp,
      variable = dplyr::pull(coords_df, {{color_to}}),
      ...) +
    ggplot2::theme_void() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL, y = NULL)

  # -----

}


#' @rdname plotSurface
#' @export
plotSurfaceInteractive <- function(object){

  check_object(object)

  surface_plots <-
    shiny::runApp(
      shiny::shinyApp(
        ui = function(){

          shinydashboard::dashboardPage(

            shinydashboard::dashboardHeader(title = "Surface Plots"),

            shinydashboard::dashboardSidebar(
              collapsed = TRUE,
              shinydashboard::sidebarMenu(
                shinydashboard::menuItem(
                  text = "Surface Plots",
                  tabName = "surface_plots",
                  selected = TRUE
                )
              )
            ),

            shinydashboard::dashboardBody(

              #----- busy indicator
              shinybusy::add_busy_spinner(spin = "cube-grid", margins = c(0,10), color = "red"),

              #----- tab items
              shinydashboard::tabItems(
                tab_surface_plots_return()
              )

            )

          )

        },
        server = function(input, output, session){

          # render uis

          output$saved_plots <- shiny::renderUI({

            saved_plots <- base::names(plot_list())

            shiny::validate(
              shiny::need(base::length(saved_plots) != 0, message = "No plots have been saved yet.")
            )

            shinyWidgets::checkboxGroupButtons(
              inputId = "saved_plots",
              label = "Choose plots to export",
              choices = saved_plots,
              selected = saved_plots,
              checkIcon = list(yes = icon("ok", lib = "glyphicon")))

          })

          #  reactive
          plot_list <- shiny::reactiveVal(value = list())
          plot_df <- shiny::reactiveVal(value = data.frame())

          # module return list
          module_return <-
            moduleSurfacePlotServer(id = "isp",
                                    object = object,
                                    final_plot = shiny::reactive(module_return()$assembled_plot()),
                                    reactive_object = shiny::reactive(object))

          # final plot
          final_plot <- shiny::reactive({

            module_return()$assembled_plot()

          })

          # store plot in list
          oe <- shiny::observeEvent(input$save_plot, {

            plot_list <- plot_list()

            if(input$plot_name %in% base::names(plot_list) | input$plot_name == ""){

              shiny::showNotification(ui = "Plot name is already taken or invalid.", type = "error")

            } else {

              plot_list[[input$plot_name]] <- final_plot()
              plot_list(plot_list)
              shiny::showNotification(ui = "Plot has been saved.", type = "message")

            }

          })


          # return last plot
          oe <- shiny::observeEvent(input$return_plot, {

            plot_list <- plot_list()

            shiny::stopApp(returnValue = plot_list[base::names(plot_list) %in% input$saved_plots])

          })



          # Distribution plotting ---------------------------------------------------

          output$surface_variable <- shiny::renderPlot({

            plot_df <- module_return()$smoothed_df()
            var_name <- base::colnames(plot_df)[5]

            if(base::is.numeric(dplyr::pull(plot_df, var_name))){

              plot_type <- input$surface_variable_plot_type

              if(plot_type == "violin"){

                add_on <- ggplot2::theme(
                  axis.text.x = ggplot2::element_blank(),
                  axis.ticks.x = ggplot2::element_blank()
                )

              } else {

                add_on <- list()
              }

              plotDistribution2(df = plot_df,
                                plot_type = plot_type,
                                binwidth = 0.05,
                                verbose = FALSE) + add_on

            } else {

              ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(x = .data[[var_name]])) +
                ggplot2::geom_bar(mapping = ggplot2::aes(fill = .data[[var_name]]), color = "black") +
                ggplot2::theme_classic() +
                ggplot2::theme(legend.position = "none") +
                confuns::scale_color_add_on(aes = "fill",
                                            variable = "discrete",
                                            clrp = module_return()$current_setting()$pt_clrp) +
                ggplot2::labs(y = "Count")

            }

          })


        }
      )
    )

  # return surface plot
  return(surface_plots)

}





#' @title Plot several surface plots at the same time
#'
#' @description Plots a surface plot for every valid element in argument
#' \code{variables} and colors it according to each element.
#'
#' \itemize{
#'  \item{ \code{plotSurfaceComparison()} Takes the spata-object as the starting point and creates the
#'  necessary data.frame from scratch according to additional parameters.}
#'  \item{ \code{plotSurfaceComparison2()} Takes a data.frame as the starting point. }
#'  }
#'
#' @param data A data.frame containing the numeric variables of interest.
#' @inherit check_variables params
#' @inherit check_pt params
#' @inherit check_method params
#' @inherit check_smooth params
#' @inherit check_display params
#' @inherit check_assign params
#' @inherit verbose params
#' @inherit image_dummy params
#' @param ... Additional arguments given to \code{ggplot2::facet_wrap()}.
#'
#' @export


plotSurfaceComparison <- function(object,
                                  variables,
                                  method_gs = NULL,
                                  normalize = NULL,
                                  smooth = NULL,
                                  smooth_span = NULL,
                                  pt_size = NULL,
                                  pt_alpha = NULL,
                                  pt_clrsp = NULL,
                                  display_image = NULL,
                                  verbose = NULL,
                                  of_sample = NA,
                                  ...){

  # 1. Control --------------------------------------------------------------

  # lazy check
  hlpr_assign_arguments(object)

  # adjusting check
  of_sample <- check_sample(object, of_sample, 1)

  all_gene_sets <- getGeneSets(object, "all")
  all_genes <- getGenes(object, of_sample = of_sample)
  all_features <-
    base::suppressWarnings(
      check_features(object, features = getFeatureNames(object, of_sample = of_sample),
                     valid_classes = "numeric")
    )

  variables <- check_variables(variables = variables,
                               all_gene_sets = all_gene_sets,
                               all_genes = all_genes,
                               all_features = all_features,
                               simplify = FALSE)

  # -----

  # 2. Extract and join data ------------------------------------------------

  joined_df <-
    joinWithVariables(object = object,
                      spata_df = getCoordsDf(object, of_sample = of_sample),
                      variables = variables,
                      average_genes = FALSE,
                      smooth = smooth,
                      smooth_span = smooth_span,
                      normalize = normalize,
                      verbose = verbose)
  # -----

  # adjust data.frame for use of ggplot2::facets

  variables <- base::unname(base::unlist(variables))
  n_variables <- base::length(variables)

  plot_df <-
    tidyr::pivot_longer(
      data = joined_df,
      cols = dplyr::all_of(variables),
      names_to = "variables",
      values_to = "values"
    )

  # plotting

  if(base::isTRUE(verbose)){base::message(glue::glue("Plotting {n_variables} different variables. (This can take a few seconds.)"))}

  ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(x = x, y = y)) +
    hlpr_image_add_on(object, display_image, of_sample) +
    ggplot2::geom_point(mapping = ggplot2::aes(color = values),
                        size = pt_size, alpha = pt_alpha) +
    confuns::scale_color_add_on(variable = plot_df$values, clrsp = pt_clrsp) +
    ggplot2::theme_void() +
    ggplot2::facet_wrap(facets = ~ variables, ...) +
    ggplot2::labs(color = "Expr.\nscore")

}

#' @rdname plotSurfaceComparison
#' @export

plotSurfaceComparison2 <- function(coords_df,
                                   variables = "all",
                                   pt_size = 2,
                                   pt_alpha = 0.9,
                                   pt_clrsp = "inferno",
                                   image = NULL,
                                   verbose = TRUE,
                                   ...){


    # 1. Control --------------------------------------------------------------

    stopifnot(base::is.data.frame(coords_df))
    confuns::is_vec(variables, "character", "variables")

    check_pt(pt_size, pt_alpha, pt_clrsp)

    if(base::all(variables == "all")){

      if(base::isTRUE(verbose)){base::message("Argument 'variables' set to 'all'. Extracting all valid, numeric variables.")}

      cnames <- base::colnames(dplyr::select_if(.tbl = coords_df, .predicate = base::is.numeric))

      valid_variables <- cnames[!cnames %in% c("x", "y", "umap1", "umap2", "tsne1", "tsne2")]

    } else {

      check_list <-
        purrr::map(c("x", "y", variables), function(i){c("numeric", "integer")}) %>%
        magrittr::set_names(value = c("x", "y", variables))

      confuns::check_data_frame(
        df = coords_df,
        var.class = check_list,
        ref = "coords_df"
      )

      valid_variables <- variables

      if(base::isTRUE(verbose)){"All specified variables found."}

    }

    # -----


    n_valid_variables <- base::length(valid_variables)
    ref <- base::ifelse(n_valid_variables > 1,
                        yes = "different variables. (This can take a few seconds.)",
                        no = "variable.")

    if(base::isTRUE(verbose)){base::message(glue::glue("Plotting {n_valid_variables} {ref}"))}


    # adjust data.frame for use of ggplot2::facets
    shifted_data <-
      tidyr::pivot_longer(
        data = coords_df,
        cols = dplyr::all_of(valid_variables),
        names_to = "aspects",
        values_to = "values"
      )

    # plotting

    ggplot2::ggplot(data = shifted_data, mapping = ggplot2::aes(x = x, y = y)) +
      hlpr_image_add_on2(image) +
      ggplot2::geom_point(mapping = ggplot2::aes(color = values),
                          size = pt_size, alpha = pt_alpha) +
      confuns::scale_color_add_on(variable = shifted_data$values, clrsp = pt_clrsp) +
      ggplot2::theme_void() +
      ggplot2::facet_wrap(facets = ~ aspects, ...) +
      ggplot2::labs(color = "Values")

  }



#' Title
#'
#' @param object
#' @param of_sample
#' @param filled
#' @param clrsp
#' @param pt_alpha
#' @param pt_color
#' @param pt_size
#' @param display_image
#' @param display_labels
#' @param ...
#'
#' @return
#' @export
#'
plotSurfaceHotspots <- function(object,
                                plot_type = "expression",
                                clrsp = NULL,
                                pt_alpha = NULL,
                                pt_clr = NULL,
                                pt_clrp = NULL,
                                pt_size = NULL,
                                pt_clrsp = NULL,
                                smooth = NULL,
                                smooth_span = NULL,
                                display_points = NULL,
                                display_image = NULL,
                                display_labels = NULL,
                                verbose = NULL,
                                of_sample = NA,
                                ...){


  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  confuns::check_one_of(
    input = plot_type,
    against = hotspot_plot_types,
    ref.input = "input for argument 'plot_type'"
  )

  # -----


  # 2. Data extraction -----------------------------------------------------

  pr_sugg <-
    getPrSuggestion(object = object, of_sample = of_sample)

  # -----

  # 3. Set up add ons -------------------------------------------------------

  # density vs. expression

  if(plot_type == "expression"){

    plot_df <-
      purrr::map_df(.x = check_pattern(object, of_sample = of_sample),
                    .f = function(htsp){

                      genes <-
                        dplyr::filter(.data = pr_sugg$df, hotspot == {{htsp}}) %>%
                        dplyr::pull(genes)

                      df <-
                        joinWith(object = object,
                                 spata_df = getCoordsDf(object, of_sample = of_sample),
                                 genes = genes,
                                 smooth = smooth,
                                 smooth_span = smooth_span,
                                 average_genes = TRUE,
                                 verbose = FALSE,
                                 normalize = TRUE) %>%
                        dplyr::mutate(hotspot = {{htsp}})

                    }) %>%
      dplyr::mutate(hotspot = forcats::as_factor(hotspot))

    hotspot_add_on <-
      list(
        ggplot2::geom_point(data = plot_df, mapping = ggplot2::aes(x = x, y = y, color = mean_genes),
                            size = pt_size, alpha = pt_alpha),
        ggplot2::facet_wrap(facets = . ~ hotspot, ...),
        scale_color_add_on(clrsp = pt_clrsp),
        ggplot2::theme_void(),
        ggplot2::theme(legend.position = "none")
      )

    panel_grid_major <- ggplot2::element_blank()
    display_points <- FALSE
    display_labels <- TRUE


  } else {

    coords_df <-
      getCoordsDf(object = object, of_sample = of_sample) %>%
      dplyr::select(x, y) %>%
      dplyr::mutate(include = TRUE)

    plot_df <-
      dplyr::select(.data = pr_sugg$df, x = center_x, y = center_y) %>%
      dplyr::mutate(include = TRUE) %>%
      base::rbind(., coords_df)

    if(plot_type == "density_filled"){

      hotspot_add_on <-
        list(
          ggplot2::geom_density2d_filled(
            data = plot_df, mapping = ggplot2::aes(x = x, y = y), ...
          ),
          ggplot2::scale_fill_viridis_d(option = clrsp)
        )

      panel_grid_major <- ggplot2::element_blank()
      display_points <- FALSE

    } else if(plot_type == "density_2d"){

      hotspot_add_on <-
        list(
          ggplot2::geom_density2d(
            data = plot_df, mapping = ggplot2::aes(x = x, y = y), ...
          ),
          ggplot2::scale_fill_viridis_d(option = clrsp)
        )

      if(base::isTRUE(display_points)){

        panel_grid_major <- ggplot2::element_blank()

      } else {

        panel_grid_major <- ggplot2::element_line()

      }

    }

  }

  # background points
  if(base::isTRUE(display_points)){

    point_add_on <-
      ggplot2::geom_point(data = coords_df,
                          mapping = ggplot2::aes(x = x, y = y),
                          size = pt_size,
                          alpha = pt_alpha,
                          color = pt_clr)

  } else {

    point_add_on <- NULL

  }

  # labels
  if(base::isTRUE(display_labels)){

    label_add_on <-
      ggrepel::geom_label_repel(
        data = pr_sugg$info_df,
        mapping = ggplot2::aes(x = center_x, y = center_y, label = hotspot)
      )

  } else {

    label_add_on <- NULL

  }



  # 4. Plotting -------------------------------------------------------------

  ggplot2::ggplot() +
    hlpr_image_add_on(object = object,
                      display_image = display_image,
                      of_sample = of_sample) +
    point_add_on +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = panel_grid_major,
      panel.grid.minor = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(),
      legend.position = "none"
    ) +
    hotspot_add_on +
    label_add_on +
    ggplot2::labs(x = NULL, y = NULL)

}


#' Title
#'
#' @param object
#' @param of_sample
#' @param color_to
#' @param n_qntls
#' @param keep_qntls
#' @param pt_alpha
#' @param pt_clrp
#' @param pt_size
#' @param smooth
#' @param smooth_span
#'
#' @return
#' @export

plotSurfaceQuantiles <- function(object,
                                 color_to,
                                 n_qntls = 5,
                                 keep_qntls = 1:n_qntls,
                                 pt_alpha = NULL,
                                 pt_clrp = NULL,
                                 pt_size = NULL,
                                 smooth = NULL,
                                 smooth_span = NULL,
                                 verbose = NULL,
                                 of_sample = NA,
                                 ...){

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  confuns::is_value(x = n_qntls, mode = "numeric")
  confuns::is_vec(x = keep_qntls, mode = "numeric")

  color_to_list <-
    check_color_to(
      color_to = color_to,
      all_features = getFeatureNames(object, of_sample = of_sample, of_class = numeric_classes),
      all_genes = getGenes(object, of_sample = of_sample),
      all_gene_sets = getGeneSets(object)
    )

  plot_df <-
    joinWithVariables(
      object = object,
      spata_df = getCoordsDf(object, of_sample = of_sample),
      variables = color_to_list,
      smooth = smooth,
      smooth_span = smooth_span,
      verbose = verbose
      ) %>%
    confuns::bin_numeric_variable(
      df = .,
      num_variable = color_to,
      discr_variable = stringr::str_c(color_to, " "),
      n_bins = n_qntls
    ) %>%
    tidyr::pivot_longer(
      cols = stringr::str_c(color_to, " "),
      names_to = "variables",
      values_to = "values"
    )

  values <-
    dplyr::pull(plot_df, var = "values") %>%
    base::levels()

  keep_values <- values[keep_qntls]

  discard <-
    dplyr::filter(plot_df, !(values %in% base::as.character(keep_values))) %>%
    dplyr::pull(var = "values") %>%
    base::unique() %>%
    base::as.character()

  adjust_clrs <- base::rep("lightgrey", base::length(discard))

  base::names(adjust_clrs) <- discard

  plotSurface2(coords_df = plot_df,
               color_to = "values",
               pt_alpha = pt_alpha,
               pt_clrp = pt_clrp,
               pt_size = pt_size,
               adjust = adjust_clrs) +
    ggplot2::facet_wrap(facets = . ~ variables) +
    ggplot2::labs(color = "Quantiles")

}


