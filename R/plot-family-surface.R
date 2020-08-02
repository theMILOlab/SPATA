
# Surface related ---------------------------------------------------------

#' @title Plot the sample
#'
#' @description Displays the sample and colors the surface according to the
#' expression of genes and gene sets or other featured characteristics.
#'
#' \itemize{
#'  \item{ \code{plotSurface()} Takes a data.frame as the starting point.}
#'  \item{ \code{plotSurface2()} Takes the spata-object as the starting point and creates the
#'  necessary data.frame from scratch according to additional paramters.}
#'  \item{ \code{plotSurfaceInteractive()} Takes only the spata-object and opens a mini-shiny
#'  application which allows for interactive plotting.}
#'
#' }
#'
#' @param data A data.frame containing at least the variables \emph{'x'} and \emph{'y'}.
#'
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

plotSurface <- function(data,
                        color_to,
                        pt_size = 2,
                        pt_alpha = 0.9,
                        pt_clrsp = "inferno",
                        pt_clrsp_dir = 1,
                        image = NULL){

  # 1. Control --------------------------------------------------------------

  # check lazy
  if(!base::is.character(color_to) |
     !base::length(color_to) == 1 |
     !color_to %in% base::colnames(data)){

    base::stopp("Argument 'color_to' needs be a variable of 'data' specified as a character value.")

  }

  check_pt(pt_size, pt_alpha, pt_clrsp)
  check_coordinate_variables(data = data, x = "x", y = "y")

  # -----

  # 2. Plotting -------------------------------------------------------------

  ggplot2::ggplot(data = data) +
    hlpr_image_add_on(image) +
    ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = y,
                                               color = .data[[color_to]]),
                        size = pt_size, alpha = pt_alpha) +
    ggplot2::scale_color_viridis_c(option = pt_clrsp) +
    ggplot2::theme_void() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL, y = NULL)

  # -----

}

#' @rdname plotSurface
#' @export
plotSurface2 <- function(object,
                         of_sample,
                         color_to,
                         method_gs = "mean",
                         smooth = FALSE,
                         smooth_span = 0.02,
                         pt_size = 2,
                         pt_alpha = 1,
                         pt_clrsp = "inferno",
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

  all_genes <- getGenes(object)
  all_gene_sets <- getGeneSets(object)
  all_features <- getFeatureNames(object)

  if(!base::all(color_to %in% all_genes)){
    color_to <- check_color_to(color_to = color_to,
                               all_features = all_features,
                               all_gene_sets = all_gene_sets,
                               all_genes = all_genes,
                               max_length = 1)
  } else {

    color_to <- list("genes" = color_to)

  }

  # adjusting check
  of_sample <- check_sample(object = object,
                            of_sample = of_sample,
                            desired_length = 1)

  # -----

  coords_df <- coordsSpatial(object, of_sample = of_sample)

  # 2. Join data and prepare ggplot add-ons ---------------------------------

  # if of length one and feature
  if("features" %in% base::names(color_to)){

    coords_df <- joinWithFeatures(object = object,
                                  coords_df = coords_df,
                                  features = color_to$features,
                                  smooth = smooth,
                                  smooth_span = smooth_span,
                                  verbose = verbose)

    labs_add_on <- hlpr_labs_add_on(input = color_to, input_str = "Feature:",
                                    color_str = color_to,
                                    display_title = display_title)

    # colour spectrum
    if(base::is.numeric(coords_df$feature)){

      scale_color_add_on <- ggplot2::scale_color_viridis_c(option = pt_clrsp)

    } else {

      scale_color_add_on <- NULL

    }

    # assemble ggplot add on
    ggplot_add_on <- list(
      ggplot2::geom_point(data = coords_df, size = pt_size, alpha = pt_alpha,
                          mapping = ggplot2::aes(x = x, y = y, color = .data[[color_to$features]])),
      scale_color_add_on,
      labs_add_on
    )

    # if of length one and gene set
  } else if("gene_sets" %in% base::names(color_to)){

    coords_df <- joinWithGeneSets(object = object,
                                  coords_df = coords_df,
                                  gene_sets = color_to$gene_sets,
                                  method_gs = method_gs,
                                  smooth = smooth,
                                  smooth_span = smooth_span,
                                  verbose = verbose)

    labs_add_on <- hlpr_labs_add_on(input = color_to$gene_sets, input_str = "Gene set:",
                                    color_str = hlpr_gene_set_name(color_to$gene_sets),
                                    display_title = display_title)

    # assemble ggplot add-on
    ggplot_add_on <- list(
      ggplot2::geom_point(data = coords_df, size = pt_size, alpha = pt_alpha,
                          mapping = ggplot2::aes(x = x, y = y, color = .data[[color_to$gene_sets]])),
      ggplot2::scale_color_viridis_c(option = pt_clrsp),
      labs_add_on
    )

  } else if("genes" %in% base::names(color_to)){

    coords_df <- joinWithGenes(object = object,
                               coords_df = coords_df,
                               genes = color_to$genes,
                               average_genes = TRUE,
                               smooth = smooth,
                               smooth_span = smooth_span,
                               verbose = verbose)

    color_str <- base::ifelse(test = base::length(color_to$genes) == 1,
                              yes = color_to$genes,
                              no = "Mean expr.\nscore")

    labs_add_on <- hlpr_labs_add_on(input = color_to, input_str = "Genes:",
                                    color_str = color_str,
                                    display_title = display_title)

    # assemble ggplot add-on
    ggplot_add_on <- list(
      ggplot2::geom_point(data = coords_df, size = pt_size, alpha = pt_alpha,
                          mapping = ggplot2::aes_string(x = "x",
                                                        y = "y",
                                                        color = "mean_genes")),
      ggplot2::scale_color_viridis_c(option = pt_clrsp),
      labs_add_on
    )


  } else if(base::is.null(color_to)){

    ggplot_add_on <- list(ggplot2::geom_point(data = coords_df, size = pt_size, alpha = 0,
                                              mapping = ggplot2::aes(x = x, y = y)))

  }

  # -----


  # 5. Plotting --------------------------------------------------------------

  hlpr_assign(assign = assign,
              object = list("point" = coords_df),
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
plotSurfaceInteractive <- function(object){

  validation(object)

  surface_plots <-
    shiny::runApp(
      shiny::shinyApp(
        ui = function(){

          shiny::fluidPage(

            #----- title
            shiny::titlePanel(title = "Surface Plot"),

            #----- shiny add-ons
            shinybusy::add_busy_spinner(spin = "cube-grid", margins = c(0,10), color = "red"),

            #----- main panel
            shiny::mainPanel(
              shiny::column(width = 12, align = "center",
                            moduleSurfacePlotUI(id = "isp"),
                            shiny::textInput("plot_name", label = NULL, value = "", placeholder = "Plot name"),
                            shiny::actionButton("save_plot", label = "Save Plot"),
                            shiny::actionButton("return_plot", label = "Return Plots")
              )
            )

          )},
        server = function(input, output, session){

          # plot list
          plot_list <- shiny::reactiveVal(value = list())

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

            if(base::length(plot_list) == 1){

              shiny::stopApp(returnValue = plot_list[[1]])

            } else {

              shiny::stopApp(returnValue = plot_list)

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
#' \code{color_to} and colors it according to the element.
#'
#' \itemize{
#'  \item{ \code{plotSurfaceComparison()} Takes a data.frame as the starting point. }
#'  \item{ \code{plotSurfaceComparison2()} Takes the spata-object as the starting point and creates the
#'  necessary data.frame from scratch according to additional parameters.}
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

plotSurfaceComparison <- function(data,
                                  variables = NULL,
                                  pt_size = 2,
                                  pt_alpha = 0.9,
                                  pt_clrsp = "inferno",
                                  image = NULL,
                                  ...){


  # 1. Control --------------------------------------------------------------

  stopifnot(base::is.data.frame(data))
  stopifnot(base::is.null(variables) | base::is.character(variables))

  check_pt(pt_size, pt_alpha, pt_clrsp)

  num_data <- data[,base::sapply(data, base::is.numeric)]
  num_variables <- base::colnames(num_data)

  if(base::is.null(variables)){

    valid_variables <- num_variables[!num_variables %in% c("x", "y", "umap1", "umap2", "tsne1", "tsne2")]

  } else {

    valid_variables <- variables[variables %in% num_variables]
    valid_variables <- valid_variables[!valid_variables %in% c("x", "y", "umap1", "umap2", "tsne1", "tsne2")]

    invalid_variables <- variables[!variables %in% valid_variables]

    if(base::length(invalid_variables) > 0){

      invalid_variables <- stringr::str_c(invalid_variables, collapse = "', '")

      base::warning(stringr::str_c("Ignoring non-numeric, invalid or not found variables: '",
                                   invalid_variables, "'", sep = "" ))

    }

  }

  # -----

  # adjust data.frame for use of ggplot2::facets
  shifted_data <-
    tidyr::pivot_longer(
      data = data,
      cols = dplyr::all_of(valid_variables),
      names_to = "aspects",
      values_to = "values"
    )

  # plotting

  ggplot2::ggplot(data = shifted_data, mapping = ggplot2::aes(x = x, y = y)) +
    hlpr_image_add_on(image) +
    ggplot2::geom_point(mapping = ggplot2::aes(color = values),
                        size = pt_size, alpha = pt_alpha) +
    ggplot2::scale_color_viridis_c(option = pt_clrsp) +
    ggplot2::theme_void() +
    ggplot2::facet_wrap(facets = ~ aspects, ...) +
    ggplot2::labs(color = "Values")

}

#' @rdname plotSurfaceComparison
#' @export
plotSurfaceComparison2 <- function(object,
                                  of_sample,
                                  variables,
                                  method_gs = "mean",
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
  all_genes <- getGenes(object)
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

  data <- coordsSpatial(object = object,
                        of_sample = of_sample)


  # join data.frame with variables to compare
  if("gene_sets" %in% base::names(variables)){

    data <- joinWithGeneSets(object,
                             coords_df = data,
                             gene_sets = variables$gene_sets,
                             method_gs = method_gs,
                             smooth = smooth,
                             smooth_span = smooth_span,
                             verbose = verbose)

  }

  if("genes" %in% base::names(variables)){

    data <- joinWithGenes(object,
                          coords_df = data,
                          genes = variables$genes,
                          average_genes = FALSE,
                          smooth = smooth,
                          smooth_span = smooth_span,
                          verbose = verbose)

  }

  if("features" %in% base::names(variables)){

    data <- joinWithFeatures(object,
                             coords_df = data,
                             features = variables$features,
                             smooth = smooth,
                             smooth_span = smooth_span,
                             verbose = verbose)

  }

  # -----

  # adjust data.frame for use of ggplot2::facets
  shifted_data <-
    tidyr::pivot_longer(
      data = data,
      cols = dplyr::all_of(base::unname(base::unlist(variables))),
      names_to = "aspects",
      values_to = "values"
    )

  # plotting

  hlpr_assign(assign = assign,
              object = list("point" = shifted_data),
              name = assign_name)

  ggplot2::ggplot(data = shifted_data, mapping = ggplot2::aes(x = x, y = y)) +
    hlpr_image_add_on2(object, display_image, of_sample) +
    ggplot2::geom_point(mapping = ggplot2::aes(color = values),
                        size = pt_size, alpha = pt_alpha) +
    ggplot2::scale_color_viridis_c(option = pt_clrsp) +
    ggplot2::theme_void() +
    ggplot2::facet_wrap(facets = ~ aspects, ...) +
    ggplot2::labs(color = "Expr.\nscore")

}




