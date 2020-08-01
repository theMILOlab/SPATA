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
#' @param color_to The variable to be displayed by color. Specified
#' as a character value. (Must be one of \code{base::colnames(data)}).
#' @param image An image of class \emph{Image} to be displayed in the background.
#' Easily accessible via \code{SPATA::image()}.
#'
#' @param object A valid object of class \emph{spata}.
#' @param of_sample The samples name specified as a character of length one.
#' @param color_to The information to be displayed by color specified as a
#' character vector. If you specify a feature or a gene set this vector needs
#' to be of length one. If you specify more than one gene the average
#' expression of these genes will be calculated.
#' @param smooth Logical value. If set to TRUE \code{plotSurface()} will smooth
#'  the values displayed by color deploying \code{stats::loess()}.
#' @param smooth_span Numeric value, given to \code{stats::loess()} if
#'  \code{smooth} is set to TRUE.
#' @param method_gs The method according to which gene sets will be handled
#' specified as a character of length one. This can be either \emph{mean} or one
#' of \emph{gsva, ssgsea, zscore, or plage}. The latter four will be given to
#' \code{gsva::GSVA()}. Ignored if \code{color_to} isn't a gene set.
#' @param pt_size The size of the points specified as a numeric value.
#' @param pt_alpha The transparency of the points specified as a numeric value.
#' @param pt_clrsp The colour spectrum used to display \code{color_to} if the
#' specified variable is continuous. Needs to be one of \emph{'inferno', 'magma',
#' 'plasma', 'cividis' or 'viridis'}.
#' @param pt_clrsp_dir The direction of the color spectrum specified as either \emph{1}
#' or \emph{-1}.
#' @param display_image Logical value to specify whether the histology image
#' of the sample is supposed to be displayed underneath the plot.
#'
#' (Can be made visible with low \code{pt_alpha} or \code{color_to} set to NULL.)
#' @param display_title Logical value to specify whether a title is supposed to
#' be displayed.
#' @param assign Logical value. If set to TRUE a list of the data needed to plot the
#' returned plot will be assigned to the global environment.
#' @param assign_name The name of the list that is assigned to the global environment
#' specified as a character value that does not already exist.
#' @param verbose Logical value to specify whether informative messages are
#' supposed to be printed or not.
#'
#' @return Returns a ggplot-object that can be additionally customized according
#' to the rules of the ggplot2-framework.
#'
#' @export
#'

plotSurface <- function(data,
                        color_to,
                        pt_size = 2,
                        pt_alpha = 0.9,
                        pt_clrsp = "inferno",
                        pt_clrsp_dir = 1,
                        image = NULL){

  # Control -----------------------------------------------------------------

  if(!base::is.character(color_to) |
     !base::length(color_to) == 1 |
     !color_to %in% base::colnames(data)){

    base::stopp("Argument 'color_to' needs be a variable of 'data' specified as a character value.")

  }

  check_pt_input(pt_size, pt_alpha, pt_clrsp)
  check_coordinates(data = data, x = "x", y = "y")


  # Plotting ----------------------------------------------------------------

  ggplot2::ggplot(data = data) +
    hlpr_image_add_on1(image) +
    ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = y,
                                               color = .data[[color_to]]),
                        size = pt_size, alpha = pt_alpha) +
    ggplot2::scale_color_viridis_c(option = pt_clrsp, limits = limits_clrsp, direction = pt_clrsp_dir) +
    ggplot2::theme_void() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL, y = NULL)

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
                         display_image = F,
                         display_title = F,
                         assign = FALSE,
                         assign_name,
                         verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  validation(x = object)

  check_pt_input(pt_size, pt_alpha, pt_clrsp)
  check_assign(assign, assign_name)

  of_sample <- check_sample(object = object,
                            sample_input = of_sample,
                            desired_length = 1)

  if(!is.character(color_to) & !is.null(color_to)){

    stop("Argument 'color_to' needs to be either NULL or a character vector.")

  }

  if(!is.logical(display_image)){

    stop("Argument 'img_bckgrd' needs to be logical.")

  }

  if(!is.logical(display_title)){

    stop("Argument 'display_title' needs to be logical.")

  }


  # 2. Extract coordinates --------------------------------------------------

  coords_df <- coordsSpatial(object, of_sample = of_sample)


  # 3. Join data and prepare ggplot add-ons ---------------------------------

  # if of length one and feature
  if(base::length(color_to) == 1 && color_to %in% getFeatureNames(object = object)){

    color_to <- check_features(object, features = color_to)

    coords_df <- joinWithFeatures(object = object,
                                  coords_df = coords_df,
                                  features = color_to,
                                  smooth = smooth,
                                  smooth_span = smooth_span,
                                  verbose = verbose)

    # ensure bugless ggplot2::aes_string
    base::colnames(coords_df)[base::colnames(coords_df) == color_to] <- "feature"

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
                          mapping = ggplot2::aes_string(x ="x",
                                                        y = "y",
                                                        color = "feature")),
      scale_color_add_on,
      labs_add_on
    )

    # if of length one and gene set
  } else if(length(color_to) == 1 && color_to %in% getGeneSets(object = object)){

    color_to <- check_gene_sets(object, gene_sets = color_to)

    coords_df <- joinWithGeneSets(object = object,
                                  coords_df = coords_df,
                                  gene_sets = color_to,
                                  method_gs = method_gs,
                                  smooth = smooth,
                                  smooth_span = smooth_span,
                                  verbose = verbose)

    labs_add_on <- hlpr_labs_add_on(input = color_to, input_str = "Gene set:",
                                    color_str = "Expr.\nscore",
                                    display_title = display_title)

    # ensure bugless ggplot2::aes_string
    base::colnames(coords_df)[base::colnames(coords_df) == color_to] <- "gene_set"

    # assemble ggplot add-on
    ggplot_add_on <- list(
      ggplot2::geom_point(data = coords_df, size = pt_size, alpha = pt_alpha,
                          mapping = ggplot2::aes_string(x = "x",
                                                        y = "y",
                                                        color = "gene_set")),
      ggplot2::scale_color_viridis_c(option = pt_clrsp),
      labs_add_on
    )



  } else if(base::any(color_to %in% getGenes(object = object))){

    rna_assay <- exprMtr(object, of_sample = of_sample)
    color_to <- check_genes(object, genes = color_to, rna_assay = rna_assay)

    coords_df <- joinWithGenes(object = object,
                               coords_df = coords_df,
                               genes = color_to,
                               average_genes = TRUE,
                               smooth = smooth,
                               smooth_span = smooth_span,
                               verbose = verbose)

    labs_add_on <- hlpr_labs_add_on(input = color_to, input_str = "Genes:",
                                    color_str = "Mean expr.\nscore",
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

  } else {

    base::warning("Could not map color to specified argument 'color_to'! (Hint: Features and Gene sets need to be of length one.)")

    ggplot_add_on <- list(ggplot2::geom_point(data = coords_df, size = pt_size, alpha = pt_alpha,
                                              mapping = ggplot2::aes(x = x, y = y)))

  }


  # 5. Plotting --------------------------------------------------------------

  hlpr_assign(assign = assign,
              object = list("point" = coords_df),
              name = assign_name)

  ggplot2::ggplot() +
    hlpr_image_add_on2(object, display_image, of_sample) +
    ggplot_add_on +
    ggplot2::coord_equal() +
    ggplot2::theme_void()

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
#'  necessary data.frame from scratch according to additional paramters.}
#'  }
#'
#' @param color_to The genes and gene sets whose expression you want to compare. If
#' specified as NULL in \code{plotSurfaceComparison()} the function will take all
#' appropriate, numeric variables of \code{data}.
#' @param ... Additional arguments given to \code{ggplot2::facet_wrap()}.
#'
#' @inherit plotSurface params return
#'
#' @export
#'

plotSurfaceComparison <- function(data,
                                  color_to = NULL,
                                  pt_size = 2,
                                  pt_alpha = 0.9,
                                  pt_clrsp = "inferno",
                                  image = NULL,
                                  ...){


  # Control -----------------------------------------------------------------

  stopifnot(base::is.data.frame(data))
  stopifnot(base::is.null(color_to) | base::is.character(color_to))

  check_pt_input(pt_size, pt_alpha, pt_clrsp)

  num_data <- data[,base::sapply(data, base::is.numeric)]
  num_color_to <- base::colnames(num_data)

  if(base::is.null(color_to)){

    valid_color_to <- num_color_to[!num_color_to %in% c("x", "y", "umap1", "umap2", "tsne1", "tsne2")]

  } else {

    valid_color_to <- color_to[color_to %in% num_color_to]
    valid_color_to <- valid_color_to[!valid_color_to %in% c("x", "y", "umap1", "umap2", "tsne1", "tsne2")]

    invalid_color_to <- color_to[!color_to %in% valid_color_to]

    if(base::length(invalid_color_to) > 0){

      invalid_color_to <- stringr::str_c(invalid_color_to, collapse = "', '")

      base::warning(stringr::str_c("Ignoring non-numeric, invalid or not found color_to: '",
                                   invalid_color_to, "'", sep = "" ))

    }

  }

  # Shift data --------------------------------------------------------------

  # adjust data.frame for use of ggplot2::facets
  shifted_data <-
    tidyr::pivot_longer(
      data = data,
      cols = dplyr::all_of(valid_color_to),
      names_to = "aspects",
      values_to = "values"
    )


  # Plotting ----------------------------------------------------------------

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
                                  color_to,
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

  # Control -----------------------------------------------------------------

  validation(object)
  check_assign(assign, assign_name)
  of_sample <- check_sample(object, of_sample, 1)

  all_gene_sets <- getGeneSets(object, "all")
  all_genes <- getGenes(object)
  all_features <- check_features(object, features = getFeatureNames(object), valid_classes = "numeric")

  color_to <- check_color_to(color_to = color_to,
                             all_gene_sets = all_gene_sets,
                             all_genes = all_genes,
                             all_features = all_features)



  # Extract and join data ---------------------------------------------------

  data <- coordsSpatial(object = object,
                        of_sample = of_sample)


  # join data.frame with variables to compare
  if("gene_sets" %in% base::names(color_to)){

    data <- joinWithGeneSets(object,
                             coords_df = data,
                             gene_sets = color_to$gene_sets,
                             method_gs = method_gs,
                             smooth = smooth,
                             smooth_span = smooth_span,
                             verbose = verbose)

  }

  if("genes" %in% base::names(color_to)){

    data <- joinWithGenes(object,
                          coords_df = data,
                          genes = color_to$genes,
                          average_genes = FALSE,
                          smooth = smooth,
                          smooth_span = smooth_span,
                          verbose = verbose)

  }

  if("features" %in% base::names(color_to)){

    data <- joinWithFeatures(object,
                             coords_df = data,
                             features = color_to$features,
                             smooth = smooth,
                             smooth_span = smooth_span,
                             verbose = verbose)

  }

  # adjust data.frame for use of ggplot2::facets
  shifted_data <-
    tidyr::pivot_longer(
      data = data,
      cols = dplyr::all_of(base::unname(base::unlist(color_to))),
      names_to = "aspects",
      values_to = "values"
    )


  # Plotting ----------------------------------------------------------------

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




