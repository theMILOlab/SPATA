
#' @title Plot the sample
#'
#' @description Displays the sample and colors the surface according to the
#' expression of genes and gene sets or other featured characteristics.
#'
#' \itemize{
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
#'
#' @inherit plot_family params
#'
#' @export

plotSurface <- function(object,
                        of_sample = "",
                        color_to = NULL,
                        method_gs = "mean",
                        normalize = TRUE,
                        smooth = FALSE,
                        smooth_span = 0.02,
                        pt_size = 2,
                        pt_alpha = 1,
                        pt_clrsp = "inferno",
                        pt_clrp = "milo",
                        pt_clr = "lightgrey",
                        display_image = FALSE,
                        display_title = FALSE,
                        assign = FALSE,
                        assign_name,
                        verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_pt(pt_size, pt_alpha, pt_clrsp)
  check_display(display_title, display_image)
  check_assign(assign, assign_name)

  # adjusting check
  of_sample <- check_sample(object = object,
                            of_sample = of_sample,
                            desired_length = 1)

  all_genes <- getGenes(object, in_sample = of_sample)
  all_gene_sets <- getGeneSets(object)
  all_features <- getFeatureNames(object)

  if(!base::all(color_to %in% all_genes)){
    color_to <- check_color_to(color_to = color_to,
                               all_features = all_features,
                               all_gene_sets = all_gene_sets,
                               all_genes = all_genes)
  } else if(base::is.null(color_to)){

    base::message("No specification for argument 'color_to'. Color according to sample belonging.")

  } else {

    color_to <- list("genes" = color_to)

  }

  # -----

  spata_df <- getCoordinates(object, of_sample = of_sample)

  # 2. Join data and prepare ggplot add-ons ---------------------------------

  # if of length one and feature
  if("features" %in% base::names(color_to)){

    spata_df <- joinWithFeatures(object = object,
                                  spata_df = spata_df,
                                  features = color_to$features,
                                  smooth = smooth,
                                  smooth_span = smooth_span,
                                  verbose = verbose)

    labs_add_on <- hlpr_labs_add_on(input = color_to, input_str = "Feature:",
                                    color_str = color_to,
                                    display_title = display_title)

    # assemble ggplot add on
    ggplot_add_on <- list(
      ggplot2::geom_point(data = spata_df, size = pt_size, alpha = pt_alpha,
                          mapping = ggplot2::aes(x = x, y = y, color = .data[[color_to$features]])),
      confuns::scale_color_add_on(clrsp = pt_clrsp, clrp = pt_clrp, variable = dplyr::pull(spata_df, color_to$feature)),
      labs_add_on
    )

    # if of length one and gene set
  } else if("gene_sets" %in% base::names(color_to)){

    spata_df <- joinWithGeneSets(object = object,
                                  spata_df = spata_df,
                                  gene_sets = color_to$gene_sets,
                                  method_gs = method_gs,
                                  smooth = smooth,
                                  smooth_span = smooth_span,
                                  normalize = normalize,
                                  verbose = verbose)

    labs_add_on <- hlpr_labs_add_on(input = color_to$gene_sets, input_str = "Gene set:",
                                    color_str = hlpr_gene_set_name(color_to$gene_sets),
                                    display_title = display_title)

    # assemble ggplot add-on
    ggplot_add_on <- list(
      ggplot2::geom_point(data = spata_df, size = pt_size, alpha = pt_alpha,
                          mapping = ggplot2::aes(x = x, y = y, color = .data[[color_to$gene_sets]])),
      confuns::scale_color_add_on(clrsp = pt_clrsp),
      labs_add_on
    )

  } else if("genes" %in% base::names(color_to)){

    spata_df <- joinWithGenes(object = object,
                               spata_df = spata_df,
                               genes = color_to$genes,
                               average_genes = TRUE,
                               smooth = smooth,
                               smooth_span = smooth_span,
                               normalize = normalize,
                               verbose = verbose)

    color_str <- base::ifelse(test = base::length(color_to$genes) == 1,
                              yes = color_to$genes,
                              no = "Mean expr.\nscore")

    labs_add_on <- hlpr_labs_add_on(input = color_to$genes, input_str = "Genes:",
                                    color_str = color_str,
                                    display_title = display_title)

    # assemble ggplot add-on
    ggplot_add_on <- list(
      ggplot2::geom_point(data = spata_df, size = pt_size, alpha = pt_alpha,
                          mapping = ggplot2::aes_string(x = "x",
                                                        y = "y",
                                                        color = "mean_genes")),
      confuns::scale_color_add_on(clrsp = pt_clrsp),
      labs_add_on
    )


  } else if(base::is.null(color_to)){

    ggplot_add_on <- list(ggplot2::geom_point(data = spata_df, size = pt_size, alpha = pt_alpha, color = pt_clr,
                                              mapping = ggplot2::aes(x = x, y = y)))

  }

  # -----


  # 5. Plotting --------------------------------------------------------------

  hlpr_assign(assign = assign,
              object = list("point" = spata_df),
              name = assign_name)

  ggplot2::ggplot() +
    hlpr_image_add_on2(object, display_image, of_sample) +
    ggplot_add_on +
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
                         image = NULL){

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
    hlpr_image_add_on(image) +
    ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = y,
                                               color = .data[[color_to]]),
                        size = pt_size, alpha = pt_alpha) +
    confuns::scale_color_add_on(clrp = pt_clrp, clrsp = pt_clrsp, variable = dplyr::pull(coords_df, {{color_to}})) +
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
                                  of_sample = "",
                                  variables,
                                  method_gs = "mean",
                                  normalize = TRUE,
                                  smooth = FALSE,
                                  smooth_span = 0.02,
                                  pt_size = 2,
                                  pt_alpha = 1,
                                  pt_clrsp = "inferno",
                                  display_image = FALSE,
                                  assign = FALSE,
                                  assign_name,
                                  verbose = TRUE,
                                  ...){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_assign(assign, assign_name)

  # adjusting check
  of_sample <- check_sample(object, of_sample, 1)

  all_gene_sets <- getGeneSets(object, "all")
  all_genes <- getGenes(object, in_sample = of_sample)
  all_features <-
    base::suppressWarnings(
      check_features(object, features = getFeatureNames(object), valid_classes = "numeric")
    )

  variables <- check_variables(variables = variables,
                               all_gene_sets = all_gene_sets,
                               all_genes = all_genes,
                               all_features = all_features,
                               simplify = FALSE)

  # -----

  # 2. Extract and join data ------------------------------------------------

  data <- getCoordinates(object, of_sample = of_sample)


  # join data.frame with variables to compare
  if("gene_sets" %in% base::names(variables)){

    data <- joinWithGeneSets(object,
                             spata_df = data,
                             gene_sets = variables$gene_sets,
                             method_gs = method_gs,
                             smooth = smooth,
                             smooth_span = smooth_span,
                             normalize = normalize,
                             verbose = verbose)

  }

  if("genes" %in% base::names(variables)){

    data <- joinWithGenes(object,
                          spata_df = data,
                          genes = variables$genes,
                          average_genes = FALSE,
                          smooth = smooth,
                          smooth_span = smooth_span,
                          normalize = normalize,
                          verbose = verbose)

  }

  if("features" %in% base::names(variables)){

    data <- joinWithFeatures(object,
                             spata_df = data,
                             features = variables$features,
                             smooth = smooth,
                             smooth_span = smooth_span,
                             verbose = verbose)

  }

  # -----

  # adjust data.frame for use of ggplot2::facets

  variables <- base::unname(base::unlist(variables))
  n_variables <- base::length(variables)

  shifted_data <-
    tidyr::pivot_longer(
      data = data,
      cols = dplyr::all_of(variables),
      names_to = "aspects",
      values_to = "values"
    )

  # plotting

  hlpr_assign(assign = assign,
              object = list("point" = shifted_data),
              name = assign_name)

  if(base::isTRUE(verbose)){base::message(glue::glue("Plotting {n_variables} different variables. (This can take a few seconds.)"))}

  ggplot2::ggplot(data = shifted_data, mapping = ggplot2::aes(x = x, y = y)) +
    hlpr_image_add_on2(object, display_image, of_sample) +
    ggplot2::geom_point(mapping = ggplot2::aes(color = values),
                        size = pt_size, alpha = pt_alpha) +
    confuns::scale_color_add_on(variable = shifted_data$values, clrsp = pt_clrsp) +
    ggplot2::theme_void() +
    ggplot2::facet_wrap(facets = ~ aspects, ...) +
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
      hlpr_image_add_on(image) +
      ggplot2::geom_point(mapping = ggplot2::aes(color = values),
                          size = pt_size, alpha = pt_alpha) +
      confuns::scale_color_add_on(variable = shifted_data$values, clrsp = pt_clrsp) +
      ggplot2::theme_void() +
      ggplot2::facet_wrap(facets = ~ aspects, ...) +
      ggplot2::labs(color = "Values")

  }






