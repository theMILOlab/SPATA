
# Dimensional reduction related -------------------------------------------
plotDimRed <- function(object,
                       method_dr,
                       of_sample,
                       color_to = NULL,
                       method_gs = "mean",
                       pt_size = 2,
                       pt_alpha = 1,
                       pt_clrsp = "inferno",
                       verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  validation(x = object)

  check_pt_input(pt_size, pt_alpha, pt_clrsp)

  of_sample <- check_sample(object = object, sample_input = of_sample)

  if(!is.character(color_to) & !is.null(color_to)){

    stop("Argument 'color_to' needs to be either NULL or a character vector.")

  }


  # 2. Extract dimensional reduction ----------------------------------------

  dimRed_df <- coordsDimRed(object, method_dr = method_dr, of_sample = of_sample)

  if(base::nrow(dimRed_df) == 0){

    base::stop("There seems to be no data for method: ", method_dr)

  }

  # 3. Join data and prepare ggplot add-ons ---------------------------------

  # if of length one and feature
  if(base::length(color_to) == 1 && color_to %in% getFeatureNames(object = object)){

    color_to <- check_features(object, features = color_to)

    dimRed_df <- joinWithFeatures(object = object,
                                  coords_df = dimRed_df,
                                  features = color_to,
                                  smooth = FALSE,
                                  verbose = verbose)

    # ensure bugless ggplot2::aes_string
    base::colnames(dimRed_df)[base::colnames(dimRed_df) == color_to] <- "feature"

    # colour spectrum
    if(base::is.numeric(dimRed_df$feature)){

      scale_color_add_on <- ggplot2::scale_color_viridis_c(option = pt_clrsp)

    } else {

      scale_color_add_on <- NULL

    }

    # assemble ggplot add on
    ggplot_add_on <- list(
      ggplot2::geom_point(data = dimRed_df, size = pt_size, alpha = pt_alpha,
                          mapping = ggplot2::aes_string(x = stringr::str_c(base::tolower(method_dr), 1, sep = ""),
                                                        y = stringr::str_c(base::tolower(method_dr), 2, sep = ""),
                                                        color = "feature")),
      scale_color_add_on,
      ggplot2::labs(color = color_to)
    )

    # if of length one and gene set
  } else if(length(color_to) == 1 && color_to %in% getGeneSets(object = object)){

    color_to <- check_gene_sets(object, gene_sets = color_to)

    dimRed_df <- joinWithGeneSets(object = object,
                                  coords_df = dimRed_df,
                                  gene_sets = color_to,
                                  method_gs = method_gs,
                                  smooth = FALSE,
                                  verbose = verbose)

    # ensure bugless ggplot2::aes_string
    base::colnames(dimRed_df)[base::colnames(dimRed_df) == color_to] <- "gene_set"

    # assemble ggplot add-on
    ggplot_add_on <- list(
      ggplot2::geom_point(data = dimRed_df, size = pt_size, alpha = pt_alpha,
                          mapping = ggplot2::aes_string(x = stringr::str_c(base::tolower(method_dr), 1, sep = ""),
                                                        y = stringr::str_c(base::tolower(method_dr), 2, sep = ""),
                                                        color = "gene_set")),
      ggplot2::scale_color_viridis_c(option = pt_clrsp),
      ggplot2::labs(color = "Expr.\nscore", title = stringr::str_c("Gene set: ", color_to, " (", method_gs, ")", sep = ""))
    )



  } else if(base::any(color_to %in% getGenes(object = object))){

    rna_assay <- exprMtr(object, of_sample = of_sample)
    color_to <- check_genes(object, genes = color_to, rna_assay = rna_assay)

    dimRed_df <- joinWithGenes(object = object,
                               coords_df = dimRed_df,
                               genes = color_to,
                               average_genes = TRUE,
                               smooth = FALSE,
                               verbose = verbose)

    # assemble ggplot add-on
    ggplot_add_on <- list(
      ggplot2::geom_point(data = dimRed_df, size = pt_size, alpha = pt_alpha,
                          mapping = ggplot2::aes_string(x = stringr::str_c(base::tolower(method_dr), 1, sep = ""),
                                                        y = stringr::str_c(base::tolower(method_dr), 2, sep = ""),
                                                        color = "mean_genes")),
      ggplot2::scale_color_viridis_c(option = pt_clrsp),
      ggplot2::labs(color = "Mean expr.\nscore")
    )


  } else {

    base::warning("Could not map color to specified argument 'color_to'! (Hint: Features and Gene sets need to be of length one.)")

    ggplot_add_on <- NULL

  }



  # 4. Plotting -------------------------------------------------------------

  if(base::is.list(ggplot_add_on)){

    ggplot2::ggplot(data = dimRed_df) +
      ggplot_add_on +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )

  } else if(base::is.null(ggplot_add_on)) {

    ggplot2::ggplot(data = dimRed_df) +
      ggplot2::geom_point(data = dimRed_df, size = pt_size, alpha = pt_alpha,
                          mapping = ggplot2::aes_string(x = stringr::str_c(base::tolower(method_dr), 1, sep = ""),
                                                        y = stringr::str_c(base::tolower(method_dr), 2, sep = ""))
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )

  }

}


#' @title Plot dimensional reduction
#'
#' @param object A valid spata-object.
#' @param of_sample The sample(s) name specified as a character.
#' @param color_to The information you want to display by color specified as a
#' character vector. If you specify a feature or a gene set this vector needs
#' to be of length one. If you specify more than one gene the average
#' expression of these genes will be calculated and displayed by color.
#' @param method_gs The method according to which gene sets will be handled
#' specified as a character of length one. This can be either \emph{mean} or one
#' of \emph{gsva, ssgsea, zscore, or plage}. The latter four will be given to
#' \code{gsva::GSVA()}.
#' @param pt_size The size of the points specified as a numeric value.
#' @param pt_alpha The transparency of the points specified as a numeric value.
#' @param pt_clrsp The colour spectrum used to display \code{color_to} if the
#' specified variable is continuous. Needs to be one of \emph{inferno, magma,
#' plasma, cividis or viridis}.
#' @param verbose Logical value. If set to TRUE informative messages with respect
#' to the computational progress made will be printed.
#'
#' (Warning messages will always be printed.)
#'
#' @return Returns a ggplot-object that can be additionally customized according
#' to the rules of the ggplot2-framework.
#'
#' @export
#'

plotUMAP <- function(object,
                     of_sample,
                     color_to = NULL,
                     method_gs = "mean",
                     pt_size = 2,
                     pt_alpha = 1,
                     pt_clrsp = "inferno",
                     verbose = TRUE){

  plotDimRed(object = object,
             method_dr = "UMAP",
             of_sample = of_sample,
             color_to = color_to,
             method_gs = method_gs,
             pt_size = pt_size,
             pt_alpha = pt_alpha,
             pt_clrsp = pt_clrsp,
             verbose = verbose)

}

#' @rdname plotUMAP
#'
#' @export
plotTSNE <- function(object,
                     of_sample,
                     color_to = NULL,
                     method_gs = "mean",
                     pt_size = 2,
                     pt_alpha = 1,
                     pt_clrsp = "inferno",
                     verbose = TRUE){

  plotDimRed(object = object,
             method_dr = "TSNE",
             of_sample = of_sample,
             color_to = color_to,
             method_gs = method_gs,
             pt_size = pt_size,
             pt_alpha = pt_alpha,
             pt_clrsp = pt_clrsp,
             verbose = verbose)

}



# Miscellaneous -----------------------------------------------------------


#' @title Gene set state plot
#'
#' @description This function takes four gene sets and visualizes the relative
#' expression of these four gene sets for every barcode by computing it's respective
#' x- and y- coordinates in the state plot. (See details.)
#'
#' (\code{plotFourStates2()} generates the data.frame that needs to be specified
#' in \code{data} from scratch. It then calls \code{plotFourStates()}).
#'
#' @param object A valid spata-object.
#' @param of_sample The sample from which to extract the data specified as a
#' character value.
#' @param data A data.frame containing at least the variables \emph{barcodes, \code{states.}}.
#' Whereby the states-variables contain the respective expression values of the specified
#' gene sets. See 'See also' for how to easily obtain these data.frames.
#' @param states The gene sets defining the four states specified as a character vector
#' of length 4.
#' @param color_to The variable in the data frame that is supposed to be displayed by color
#' specified as a character value.
#' @param pt_size The size of the points specified as a numeric value.
#' @param pt_alpha The transparency of the points specified as a numeric value.
#' @param pt_clrsp The colour spectrum used to display \code{color_to} if the
#' specified variable is continuous. Needs to be one of \emph{inferno, magma,
#' plasma, cividis or viridis}.
#' @param display_labels Logical value.
#' @param verbose Logical value. If set to TRUE informative messages with respect
#' to the computational progress made will be printed.
#'
#' (Warning messages will always be printed.)
#'
#' @seealso Combine \code{coordsSpatial()} and \code{joinWithGeneSets()} to obtain
#' a valid input data.frame for \code{data}.
#'
#'
#' @return Returns a ggplot-object that can be additionally customized according
#' to the rules of the ggplot2-framework.
#'
#' @export
#'

plotFourStates <- function(data,
                           states,
                           color_to = NULL,
                           pt_size = 1.5,
                           pt_alpha = 0.9,
                           pt_clrsp = "inferno",
                           display_labels = TRUE){


  # Control -----------------------------------------------------------------

  if(!base::is.data.frame(data)){

    base::stop("Argument 'data' needs to be of type data.frame.")

  } else if(!"barcodes" %in% base::colnames(data)){

    base::stop("Data.frame 'data' needs to have a variable named 'barcodes'.")

  }

  if(!base::is.null(color_to)){

    if(!base::length(color_to) == 1 |
       !color_to %in% base::colnames(data)){

      base::stop("Argument 'color_to' is not a variable of 'data'.")

    }

  }

  check_pt_clrsp(pt_clrsp)

  # Data wrangling ----------------------------------------------------------

  sym <- rlang::sym
  max <- base::max
  abs <- base::abs
  log2 <- base::log2

  plot_df <-
    tidyr::pivot_longer(
      data = data,
      cols = dplyr::all_of(states),
      names_to = "gene_set",
      values_to = "gene_set_expr"
    ) %>%
    dplyr::group_by(barcodes) %>%
    dplyr::filter(gene_set_expr == max(gene_set_expr)) %>%
    dplyr::ungroup() %>%
    dplyr::select(barcodes, max_gene_set = gene_set, max_expr = gene_set_expr) %>%
    dplyr::mutate(max_loc = dplyr::if_else(max_gene_set %in% states[1:2], true = "top", false = "bottom")) %>%
    dplyr::left_join(x = dplyr::select(data, -x, -y), y = ., by = "barcodes") %>%
    dplyr::mutate(
      pos_x = dplyr::case_when(
        max_loc == "top" & !!sym(states[1]) > !!sym(states[2]) ~ (log2(abs((!!sym(states[1]) - !!sym(states[2])) + 1)) * -1),
        max_loc == "top" & !!sym(states[2]) > !!sym(states[1]) ~ log2(abs((!!sym(states[2]) - !!sym(states[1])) + 1)),
        max_loc == "bottom" & !!sym(states[3]) > !!sym(states[4]) ~ (log2(abs((!!sym(states[3]) - !!sym(states[4])) + 1)) * -1),
        max_loc == "bottom" & !!sym(states[4]) > !!sym(states[3]) ~ log2(abs((!!sym(states[4]) - !!sym(states[3])) + 1)))
    ) %>%
    dplyr::group_by(barcodes) %>%
    dplyr::mutate(
      pos_y = dplyr::case_when(
        max_loc == "bottom" ~ (log2(abs(max(c(!!sym(states[3]), !!sym(states[4]))) - max(!!sym(states[1]), !!sym(states[2])) + 1)) * -1),
        max_loc == "top" ~ log2(abs(max(c(!!sym(states[1]), !!sym(states[2]))) - max(!!sym(states[3]), !!sym(states[4])) + 1))
      )
    ) %>%
    dplyr::filter(!base::is.na(pos_x) & !is.na(pos_y))



  # Additional add ons ------------------------------------------------------

  states <- hlpr_gene_set_name(states)
  color_to_lab <- hlpr_gene_set_name(color_to)

  #xlab <- base::bquote(paste("log2(GSV-Score "[.(states[3])]*" - GSV-Score "[.(states[4])]*")"))
  #ylab <- base::bquote(paste("log2(GSV-Score "[.(states[2])]*" - GSV-Score "[.(states[1])]*")"))


  # scale color add-on
  if(!is.null(color_to) && base::is.numeric(dplyr::pull(plot_df, var = {{color_to}}))){

    scale_color_add_on <- ggplot2::scale_colour_viridis_c(option = pt_clrsp)

  } else {

    scale_color_add_on <- NULL

  }

  max <- base::max(plot_df$pos_x, plot_df$pos_y)

  # geom text add-on
  if(base::isTRUE(display_labels)){

    tpx <- max * 0.7
    tpy <- max * 1.05

    # assemble text data.frame
    text_df <- data.frame(
      "x" = as.numeric(c(-tpx, tpx, -tpx, tpx)),
      "y" = as.numeric(c(tpy , tpy, -tpy, -tpy)),
      "states" = states
    )

    geom_text_add_on <-
      ggplot2::geom_text(mapping = ggplot2::aes(x = x, y = y, label = states),
                         data = text_df)

    print(text_df)

  } else {

    geom_text_add_on <- NULL

  }


  # plotting
  ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "lightgrey") +
    ggplot2::geom_hline(yintercept = 0,  linetype = "dashed", color = "lightgrey") +
    geom_text_add_on +
    ggplot2::geom_point(mapping = ggplot2::aes_string(x = "pos_x", y = "pos_y", color = color_to),
                        size = pt_size, alpha = pt_alpha, data = plot_df) +
    ggplot2::scale_x_continuous(limits = c(-max*1.1, max*1.1), expand = c(0,0)) +
    ggplot2::scale_y_continuous(limits = c(-max*1.1, max*1.1), expand = c(0,0)) +
    scale_color_add_on +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    ggplot2::labs(x = NULL, y = NULL, color = color_to_lab)

}


#' @rdname plotFourStates
#' @export
#'

plotFourStates2 <- function(object,
                            of_sample,
                            states,
                            color_to = NULL,
                            method_gs = "gsva",
                            pt_size = 1.5,
                            pt_alpha = 0.9,
                            pt_clrsp = "inferno",
                            display_labels = TRUE,
                            assign = FALSE,
                            assign_name,
                            verbose = TRUE){


  # Control -----------------------------------------------------------------

  validation(object)

  check_pt_input(pt_size, pt_alpha, pt_clrsp)
  check_assign(assign, assign_name)

  of_sample <- check_sample(object, sample_input = of_sample, desired_length = 1)
  states <- check_gene_sets(object, gene_sets = states, max_length = 4)

  if(base::length(states) != 4){

    base::stop(stringr::str_c(base::length(states), "valid gene sets provided.",
                              "Need four.",sep = " "))

  }

  all_genes <- getGenes(object)
  all_gene_sets <- getGeneSets(object)
  all_features <- getFeatureNames(object)


  if(!base::is.null(color_to)){
    color_to <- check_color_to(color_to = color_to,
                               all_features = all_features,
                               all_gene_sets = all_gene_sets,
                               all_genes = all_genes,
                               max_length = 1)
  }



  # Data extraction ---------------------------------------------------------

  data <-
    coordsSpatial(object = object,
                  of_sample = of_sample) %>%
    joinWithGeneSets(object,
                     coords_df = .,
                     gene_sets = states,
                     normalize = TRUE,
                     method_gs = method_gs,
                     verbose = verbose)

  if(!base::is.null(color_to)){

    if(color_to %in% all_genes){

      data <-
        joinWithGenes(object,
                      coords_df = data,
                      genes = color_to,
                      average_genes = FALSE,
                      normalize = TRUE,
                      verbose = verbose)

    } else if(color_to %in% all_gene_sets){

      data <-
        joinWithGeneSets(object,
                         coords_df = data,
                         gene_sets = color_to,
                         method_gs = method_gs,
                         normalize = TRUE,
                         verbose = verbose)

    } else if(color_to %in% all_features){

      data <-
        joinWithFeatures(object,
                         coords_df = data,
                         features = color_to,
                         normalize = TRUE,
                         verbose = verbose)

    }

  }


  # Plotting ----------------------------------------------------------------

  hlpr_assign(assign = assign,
              object = list("point" = data),
              name = assign_name)

  plotFourStates(data = data,
                 states = states,
                 color_to = color_to,
                 pt_size = pt_size,
                 pt_alpha = pt_alpha,
                 pt_clrsp = pt_clrsp,
                 display_labels = display_labels)


}



#' Plot segmentation
#'
#' @description Displays the segmentation of a specified sample that was drawn with
#' \code{SPATA::createSegmentation()}.
#'
#' @param object A valid spata-object.
#' @param of_sample The sample to be displayed specified as a character vector of
#' length one.
#' @param pt_size The size of the points specified as a numeric value.
#' @param verbose Logical value. If set to TRUE informative messages with respect
#' to the computational progress made will be printed.
#'
#' @return Returns a ggplot-object that can be additionally customized according
#' to the rules of the ggplot2-framework.
#'
#' @export
#'
plotSegmentation <- function(object,
                             of_sample,
                             pt_size = 2,
                             verbose = FALSE){

  validation(x = object)
  of_sample <- check_sample(object, of_sample, desired_length = 1)
  check_pt_input(pt_size = pt_size)


  plot_df <-
    coordsSpatial(object, of_sample = of_sample) %>%
    joinWithFeatures(object, coords_df = ., features = "segment", verbose = verbose)

  segment_df <- dplyr::filter(plot_df, segment != "")

  ggplot2::ggplot() +
    ggplot2::geom_point(data = plot_df, mapping = ggplot2::aes(x = x, y = y), size = pt_size, color = "lightgrey") +
    ggplot2::geom_point(data = segment_df, mapping = ggplot2::aes(x = x, y = y, color = segment)) +
    ggplot2::theme_void() +
    ggplot2::labs(color = "Segments") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5)))


}



# Surface related ---------------------------------------------------------

#' @title Plot surface
#'
#' @description Displays the sample and color the surface according to the
#' expression of genes and gene sets or other featured characteristics.
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
#' specified variable is continuous. Needs to be one of \emph{inferno, magma,
#' plasma, cividis or viridis}.
#' @param display_image Logical value to specify whether the histological image
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

plotSurface <- function(object,
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


  # 4. Set up additional add-ons --------------------------------------------

  # set up background
  if(base::isTRUE(display_image)){

    image_raster <-
      image(object, of_sample) %>%
      grDevices::as.raster()

    img_info <-
      image_raster %>%
      magick::image_read() %>%
      magick::image_info()

    st_image <-
      image_raster %>%
      magick::image_read() %>%
      magick::image_flip()

    image_add_on <-
      ggplot2::annotation_raster(raster = st_image,
                                 xmin = 0, ymin = 0,
                                 xmax = img_info$width,
                                 ymax = img_info$height)

  } else {

    image_add_on <- NULL

  }



  # 5. Plotting --------------------------------------------------------------

  hlpr_assign(assign = assign,
              object = list("point" = coords_df),
              name = assign_name)

  ggplot2::ggplot() +
    image_add_on +
    ggplot_add_on +
    ggplot2::coord_equal() +
    ggplot2::theme_void()

}

#' @rdname plotSurface
#'
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
#' \code{color_to} and colors the surface according to the expression of
#' genes and gene sets.
#'
#' @param object A valid object of class \emph{spata}.
#' @param of_sample The samples name specified as a character of length one.
#' @param color_to The genes and gene sets whoose expression you want to
#' compare specified as a character vector or a list.
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
#' specified variable is continuous. Needs to be one of \emph{inferno, magma,
#' plasma, cividis or viridis}.
#' @param display_image Logical value to specify whether the histological image
#' of the sample is supposed to be displayed underneath the plot.
#' (Can be made visible with low \code{pt_alpha} or \code{color_to} set to NULL.)
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

plotSurfaceComparison <- function(object,
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
                                  verbose = TRUE){

  # control
  validation(object)
  check_assign(assign, assign_name)
  of_sample <- check_sample(object, of_sample, 1)

  all_gene_sets <- getGeneSets(object, "all")
  all_genes <- getGenes(object, "all", "all")

  color_to <- check_color_to(color_to = color_to,
                             all_gene_sets = all_gene_sets,
                             all_genes = all_genes)

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

  # adjust data.frame for use of ggplot2::facets
  shifted_data <-
    tidyr::pivot_longer(
      data = data,
      cols = dplyr::all_of(base::unname(base::unlist(color_to))),
      names_to = "aspects",
      values_to = "values"
    )

  # display sample image
  image_add_on <- hlpr_image_add_on(object = object,
                                    display_image = display_image)


  # plot
  hlpr_assign(assign = assign,
              object = list("point" = shifted_data),
              name = assign_name)

  ggplot2::ggplot(data = shifted_data, mapping = ggplot2::aes(x = x, y = y)) +
    image_add_on +
    ggplot2::geom_point(mapping = ggplot2::aes(color = values),
                        size = pt_size, alpha = pt_alpha) +
    ggplot2::scale_color_viridis_c(option = pt_clrsp) +
    ggplot2::theme_void() +
    ggplot2::facet_wrap(facets = ~ aspects, scales = "fixed") +
    ggplot2::labs(color = "Expr.\nscore")

}




# Trajectory related ------------------------------------------------------


#' @title Plot trajectory
#'
#' @description Displays a trajectory of a specified sample that was
#' drawn with \code{SPATA::createTrajectories()}.
#'
#' @param object A valid object of class \emph{spata}.
#' @param trajectory_name The trajectory to plot specified as a character vector
#' of length one.
#' @param of_sample The samples name specified as a character of length one.
#' @param color_to The information to be displayed by color specified as a
#' character vector. If you specify a feature or a gene set this vector needs
#' to be of length one. If you specify more than one gene the average
#' expression of these genes will be calculated..
#' @param method_gs The method according to which gene sets will be handled
#' specified as a character of length one. This can be either \emph{mean} or one
#' of \emph{gsva, ssgsea, zscore, or plage}. The latter four will be given to
#' \code{gsva::GSVA()}. Ignored if \code{color_to} isn't a gene set.
#' @param smooth Logical value. If set to TRUE \code{plotSurface()} will smooth
#'  the values displayed by color deploying \code{stats::loess()}.
#' @param smooth_span Numeric value, given to \code{stats::loess()} if
#'  \code{smooth} is set to TRUE.
#' @param pt_size The size of the points specified as a numeric value.
#' @param pt_alpha The transparency of the points specified as a numeric value.
#' @param pt_clr The color of the points if \code{color_to} is set to NULL.
#' @param pt_clrsp The color spectrum used to display \code{color_to} if the
#' specified variable is continuous. Needs to be one of \emph{inferno, magma,
#' plasma, cividis or viridis}.
#' @param sgmt_size The size of the segment arrrow specified as a numeric value.
#' @param display_image Logical value. If set to TRUE the image will be displayed
#' as the background. If set to FALSE barcodes that do not fall into the trajectory
#' will be displayed in grey.
#' @param display_title Logical value.
#' @param verbose Logical value. If set to TRUE informative messages with respect
#' to the computational progress made will be printed.
#'
#' @return Returns a ggplot-object that can be additionally customized according
#' to the rules of the ggplot2-framework.
#' @export
#'

plotTrajectory <- function(object,
                           trajectory_name,
                           of_sample,
                           color_to = NULL,
                           method_gs = "mean",
                           smooth = FALSE,
                           smooth_span = 0.02,
                           pt_size = 2.5,
                           pt_alpha = 1,
                           pt_clr = "steelblue",
                           pt_clrsp = "inferno",
                           sgmt_size = 1,
                           display_image = F,
                           display_title = F,
                           verbose = T){


  # 1. Control --------------------------------------------------------------

  validation(x = object)

  of_sample <- check_sample(object = object, sample_input = of_sample, desired_length = 1)

  if(!is.character(color_to) & !is.null(color_to)){

    stop("Argument 'color_to' needs to be either NULL or a character vector.")

  }
  if(!pt_clrsp %in% c("inferno", "viridis", "cividis", "plasma", "magma")){

    stop("Please enter a valid color spectrum. (inferno, viridis, cividis, plasma, magma)")

  }
  if(length(method_gs) != 1 | !method_gs %in% c("mean", "gsva", "ssgsea", "zscore", "plage")){

    stop("Invalid input for argument 'method_gs'.")

  }
  if(isTRUE(smooth) && !is.numeric(smooth_span)){

    message("Argument 'span' needs to be numeric.")
    smooth <- FALSE

  } else if(isTRUE(smooth) && smooth_span < 0.01){

    stop("Argument 'smooth_span' needs to be higher than '0.01'.")

  }

  if(!is.logical(smooth)){

    message("Argument 'smooth' needs to be logical.")
    smooth <- FALSE

  }

  if(!is.logical(display_image)){

    stop("Argument 'display_image' needs to be logical.")

  }
  if(!is.logical(display_title)){

    stop("Argument 'display_title' needs to be logical.")

  }



  # 2. Extract data ---------------------------------------------------------

  t_object <-
    trajectory(object = object, trajectory = trajectory_name, of_sample = of_sample)

  trajectory_bc <- dplyr::pull(t_object@compiled_trajectory_df, barcodes)
  trajectory_sgmt_df <- t_object@segment_trajectory_df

  sample_coords <-
    coordsSpatial(object, of_sample = of_sample)


  # 3. Determine trajectory geom_point layer --------------------------------

  if(base::length(color_to) == 1 && color_to %in% getFeatureNames(object = object)){

    # check feature
    feature <- check_features(object, features = color_to)

    # join coordinates with feature
    plot_df <-
      joinWithFeatures(object = object,
                       coords_df = sample_coords,
                       features = feature,
                       smooth = smooth,
                       smooth_span = smooth_span,
                       verbose = verbose) %>%
      dplyr::filter(barcodes %in% trajectory_bc)

    # prepare ggplot add on
    base::colnames(plot_df)[base::colnames(plot_df) == feature] <- "feature"

    if(base::is.numeric(plot_df$feature)){

      scale_color_add_on <- ggplot2::scale_color_viridis_c(option = pt_clrsp)

    } else {

      scale_color_add_on <- NULL

    }

    if(base::isTRUE(display_title)){

      plot_title <- stringr::str_c("Feature:", feature, sep = " ")

      labs_add_on <- ggplot2::labs(title = plot_title, color = feature)

    } else {

      labs_add_on <- ggplot2::labs(color = feature)

    }

    ggplot_add_on <-
      list(
        ggplot2::geom_point(data = plot_df, size = pt_size, alpha = pt_alpha,
                            mapping = ggplot2::aes_string(x = "x", y = "y", color = "feature")),
        labs_add_on,
        scale_color_add_on
      )


  } else if(base::length(color_to) == 1 && color_to %in% getGeneSets(object = object)){

    # check gene set
    gene_set <- check_gene_sets(object, gene_sets = color_to)

    # join coordinates with gene set
    plot_df <-
      joinWithGeneSets(object = object,
                       coords_df = sample_coords,
                       gene_sets = gene_set,
                       method_gs = method_gs,
                       smooth = smooth,
                       smooth_span = smooth_span,
                       verbose = verbose) %>%
      dplyr::filter(barcodes %in% trajectory_bc)

    # prepare ggplot add on
    base::colnames(plot_df)[base::colnames(plot_df) == gene_set] <- "gene_set"

    if(isTRUE(display_title)){

      gene_set_string <- stringr::str_c(gene_set, " (", method_gs, ")", sep = "")

      plot_title <- stringr::str_c("Gene set:", gene_set_string, sep = " ")

      labs_add_on <- ggplot2::labs(title = plot_title, color = "Expr.\nscore")

    } else {

      labs_add_on <- ggplot2::labs(color = "Expr.\nscore")

    }

    ggplot_add_on <-
      list(
        ggplot2::geom_point(data = plot_df, size = pt_size, alpha = pt_alpha,
                            mapping = ggplot2::aes_string(x = "x", y = "y", color = "gene_set")),
        ggplot2::scale_color_viridis_c(option = pt_clrsp),
        labs_add_on
      )



  } else if(base::any(color_to %in% getGenes(object = object))){

    # check genes
    rna_assay <- exprMtr(object, of_sample = of_sample)
    genes <- check_genes(object, genes = color_to, rna_assay = rna_assay)

    # join coordinates with genes
    plot_df <-
      joinWithGenes(object = object,
                    coords_df = sample_coords,
                    genes = genes,
                    average_genes = TRUE,
                    smooth = smooth,
                    smooth_span = smooth_span,
                    verbose = verbose) %>%
      dplyr::filter(barcodes %in% trajectory_bc)

    # prepare ggplot add on
    if(base::isTRUE(display_title)){

      if(base::length(genes) > 5){

        genes <- c(genes[1:5], stringr::str_c("... +", (length(genes)-5), sep = " "))

      }

      genes_string <- stringr::str_c(genes, collapse = ", ")

      plot_title <- stringr::str_c("Genes:", genes_string, sep = " ")

      labs_add_on <- ggplot2::labs(title = plot_title, color = "Mean expr.\n score")

    } else {

      labs_add_on <- ggplot2::labs(color = "Mean expr.\n score")
    }

    ggplot_add_on <-
      list(
        ggplot2::geom_point(data = plot_df, size = pt_size, alpha = pt_alpha,
                            mapping = ggplot2::aes_string(x = "x", y = "y", color = "mean_genes")),
        ggplot2::scale_color_viridis_c(option = pt_clrsp),
        labs_add_on
      )


  } else {

    if(!base::is.null(color_to)){

      warning("Could not map argument 'color_to'. Does input match the requirements? (Hint: Features and gene sets need to be of length one.)")

    }

    ggplot_add_on <- ggplot2::geom_point(data = t_object@compiled_trajectory_df,
                                         mapping = ggplot2::aes(x = x, y = y), color = pt_clr, size = pt_size)

  }


  # 4. Set up additional add-ons --------------------------------------------

  # set up background
  if(isTRUE(display_image)){

    image_raster <-
      image(object, of_sample) %>%
      grDevices::as.raster()

    img_info <-
      image_raster %>%
      magick::image_read() %>%
      magick::image_info()

    st_image <-
      image_raster %>%
      magick::image_read() %>%
      magick::image_flip()

    background_add_on <-list(
      ggplot2::geom_point(data = dplyr::filter(sample_coords, !barcodes %in% trajectory_bc),
                          mapping = ggplot2::aes(x = x, y = y),
                          color = "white", size = pt_size, alpha = 0.01),
      ggplot2::annotation_raster(raster = st_image,
                                 xmin = 0, ymin = 0,
                                 xmax = img_info$width,
                                 ymax = img_info$height)
    )


  } else {

    background_add_on <-  ggplot2::geom_point(data = dplyr::filter(sample_coords, !barcodes %in% trajectory_bc),
                                              mapping = ggplot2::aes(x = x, y = y),
                                              color = "lightgrey", size = pt_size)

  }


  # 5. Plotting -------------------------------------------------------------

  ggplot2::ggplot() +
    background_add_on +
    ggplot_add_on +
    ggplot2::geom_segment(data = trajectory_sgmt_df,
                          mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                          color = "black", size = sgmt_size,
                          arrow = ggplot2::arrow(length = ggplot2::unit(x = 0.125, "inches"))) +
    ggplot2::theme_void() +
    ggplot2::coord_equal()

}


#' @title Trajectory line plots
#'
#' @description Displays values along a trajectory direction.
#'
#' @param object A valid spata-object.
#' @param trajectory_name The trajectory to plot specified as a character vector
#' of length one.
#' @param of_sample The sample from which the trajectory derives specfified as a
#' character vector of length one.
#' @param features The features you want to display specified as a character
#' vector (only numeric features can be plotted with \code{plotTrajectoryFeatures()}).
#' @param genes The genes to be displayed specified as a character vector.
#' @param average_genes Logical value. If set to TRUE coords_df will be joined with
#' the mean expression values of all genes specified.
#' @param gene_sets The gene sets to be displayed specified as a character
#' vector.
#' @param method_gs The method according to which gene sets will be handled
#' specified as a character of length one. This can be either \emph{mean} or one
#' of \emph{gsva, ssgsea, zscore, or plage}. The latter four will be given to
#' \code{gsva::GSVA()}.
#' @param smooth_method The smoothing method that will be used specified as a
#' character vector of length one (e.g. \emph{"lm", "glm", "gam", "loess"}).
#' @param smooth_se Logical value. If set to TRUE the confidence interval will be
#' displayed.
#' @param smooth_span The amount of smoothing specified as a numeric vector of
#' length one given to \code{stats::loess()}.
#' @param verbose Logical value. If set to TRUE informative messages with respect
#' to the computational progress made will be printed.
#'
#' @return Returns a ggplot-object that can be additionally customized according
#' to the rules of the ggplot2-framework.
#'
#' @export

plotTrajectoryFeatures <- function(object,
                                   trajectory_name,
                                   of_sample,
                                   features = "percent.mt",
                                   smooth_method = "loess",
                                   smooth_span = 0.2,
                                   smooth_se = T,
                                   verbose = TRUE){


  # 1. Control --------------------------------------------------------------

  validation(object)

  of_sample <- check_sample(object = object, sample_input = of_sample, desired_length = 1)

  features <-
    check_features(object, features = features, valid_classes = c("numeric", "integer")) %>%
    base::unname()


  # 2. Data wrangling -------------------------------------------------------

  t_object <-
    trajectory(object = object, trajectory = trajectory_name, of_sample = of_sample)

  coords_with_feature <-
    t_object@compiled_trajectory_df %>%
    dplyr::mutate(order_binned = plyr::round_any(projection_length, accuracy = 5, f = floor)) %>%
    joinWithFeatures(object = object,
                     coords_df = .,
                     features = features,
                     smooth = F,
                     verbose = verbose)

  result_df <-
    coords_with_feature %>%
    dplyr::group_by(trajectory_part, order_binned) %>%
    dplyr::summarise(dplyr::across(.cols = dplyr::all_of(x = {{features}}), ~ mean(., na.rm = TRUE)), .groups = "drop_last") %>%
    dplyr::mutate(trajectory_part_order = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(trajectory_order = dplyr::row_number()) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(features),
                        names_to = "Features",
                        values_to = "Values")


  vline_df <-
    result_df %>%
    dplyr::group_by(trajectory_part) %>%
    dplyr::filter(trajectory_order %in% c(base::min(trajectory_order), base::max(trajectory_order)) &
                    trajectory_part_order == 1 &
                    Features == features[1])


  # 3. Plotting -------------------------------------------------------------

  ggplot2::ggplot(data = result_df, mapping = ggplot2::aes(x = trajectory_order, y = Values)) +
    ggplot2::geom_vline(data = vline_df[-1,],
                        mapping = ggplot2::aes(xintercept = trajectory_order), linetype = "dashed", color = "grey") +
    ggplot2::geom_smooth(size = 1.5, span = smooth_span, method = smooth_method, formula = y ~ x,
                         mapping = ggplot2::aes(color = Features), se = smooth_se) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "inches"))),
      axis.line.y = ggplot2::element_line()
    ) +
    ggplot2::labs(x = "Direction", y = NULL)


}



#' @rdname plotTrajectoryFeatures
#'
#' @export
plotTrajectoryGenes <- function(object,
                                trajectory_name,
                                of_sample,
                                genes,
                                average_genes = F,
                                smooth_method = "loess",
                                smooth_span = 0.2,
                                smooth_se = T,
                                verbose = T){


  # 1. Control --------------------------------------------------------------

  validation(object)

  of_sample <- check_sample(object = object, sample_input = of_sample, desired_length = 1)


  if(base::length(genes) > 5 && base::isFALSE(average_genes) && base::isTRUE(verbose)){

    base::message("In order to plot more than 5 genes we recommend 'plotTrajectoryHeatmap()'.")

  }

  if(average_genes){

    columns <- "mean_genes"

    y_title <- "Mean expression score"

    rna_assay <- exprMtr(object = object, of_sample = of_sample)
    genes <- check_genes(object, genes = genes, max_length = 10, rna_assay = rna_assay)

    if(base::length(genes) == 1){

      average_genes <- F
      plot_title <- genes
      base::warning("Can not average one gene. Treating 'average_genes' as FALSE.")

    }


  } else {

    rna_assay <- exprMtr(object = object, of_sample = of_sample)
    genes <- check_genes(object, genes = genes, max_length = 10, rna_assay = rna_assay)

    columns <- genes

    y_title <- "Expression score"

    plot_title <- NULL

  }


  # 2. Data wrangling -------------------------------------------------------

  t_object <-
    trajectory(object = object, trajectory_name = trajectory_name, of_sample = of_sample)

  coords_with_genes <-
    t_object@compiled_trajectory_df %>%
    dplyr::mutate(order_binned = plyr::round_any(projection_length, accuracy = 5, f = floor)) %>%
    joinWithGenes(object = object,
                  coords_df = .,
                  genes = genes,
                  average_genes = average_genes,
                  verbose = verbose)

  result_df <-
    coords_with_genes %>%
    dplyr::group_by(trajectory_part, order_binned) %>%
    dplyr::summarise(dplyr::across(.cols = dplyr::all_of(x = {{columns}}), ~ mean(., na.rm = TRUE)), .groups = "drop_last") %>%
    dplyr::mutate(trajectory_part_order = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(trajectory_order = dplyr::row_number())


  if(!average_genes){

    result_df <-
      tidyr::pivot_longer(data = result_df,
                          cols = dplyr::all_of(genes),
                          names_to = "genes",
                          values_to = "expr_score")

  } else {

    result_df <-
      dplyr::select(result_df, expr_score = mean_genes, genes = mean_genes, dplyr::everything())

  }

  vline_df <-
    result_df %>%
    dplyr::group_by(trajectory_part) %>%
    dplyr::filter(trajectory_order %in% c(base::min(trajectory_order), base::max(trajectory_order)) &
                    trajectory_part_order == 1 &
                    genes == genes[1])


  # 3. Plotting -------------------------------------------------------------

  ggplot2::ggplot(data = result_df, mapping = ggplot2::aes(x = trajectory_order, y = expr_score)) +
    ggplot2::geom_vline(data = vline_df[-1,],
                        mapping = ggplot2::aes(xintercept = trajectory_order), linetype = "dashed", color = "grey") +
    ggplot2::geom_smooth(size = 1.5, span = smooth_span, method = smooth_method, formula = y ~ x,
                         mapping = ggplot2::aes(color = genes), se = smooth_se) +
    ggplot2::scale_y_continuous(breaks = base::seq(0 , 1, 0.2), labels = base::seq(0 , 1, 0.2)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "inches"))),
      axis.line.y = ggplot2::element_line()
    ) +
    ggplot2::labs(x = "Direction", y = y_title, color = "Genes", title = plot_title)


}


#' @rdname plotTrajectoryFeatures
#'
#' @export
plotTrajectoryGeneSets <- function(object,
                                   trajectory_name,
                                   of_sample,
                                   gene_sets = "Neftel_NPC_Comb",
                                   method_gs = "mean",
                                   smooth_method = "loess",
                                   smooth_span = 0.2,
                                   smooth_se = T,
                                   verbose = T){


  # 1. Control --------------------------------------------------------------

  validation(object)

  of_sample <- check_sample(object = object, sample_input = of_sample, desired_length = 1)

  gene_sets <- check_gene_sets(object, gene_sets = gene_sets, max_length = 10)


  # 2. Data wrangling -------------------------------------------------------

  t_object <-
    trajectory(object = object, trajectory = trajectory_name, of_sample = of_sample)

  coords_with_feature <-
    t_object@compiled_trajectory_df %>%
    dplyr::mutate(order_binned = plyr::round_any(projection_length, accuracy = 5, f = floor)) %>%
    joinWithGeneSets(object = object,
                     coords_df = .,
                     gene_sets = gene_sets,
                     method_gs = method_gs,
                     verbose = verbose,
                     smooth = FALSE)

  result_df <-
    coords_with_feature %>%
    dplyr::group_by(trajectory_part, order_binned) %>%
    dplyr::summarise(dplyr::across(.cols = dplyr::all_of(x = {{gene_sets}}), ~ mean(., na.rm = TRUE)), .groups = "drop_last") %>%
    dplyr::mutate(trajectory_part_order = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(trajectory_order = dplyr::row_number()) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(gene_sets),
                        names_to = "gsets",
                        values_to = "expr_score")


  vline_df <-
    result_df %>%
    dplyr::group_by(trajectory_part) %>%
    dplyr::filter(trajectory_order %in% c(base::min(trajectory_order), base::max(trajectory_order)) &
                    trajectory_part_order == 1 &
                    gsets == gene_sets[1])

  # 3. Plotting -------------------------------------------------------------


   ggplot2::ggplot(data = result_df, mapping = ggplot2::aes(x = trajectory_order, y = expr_score)) +
    ggplot2::geom_vline(data = vline_df[-1,],
                        mapping = ggplot2::aes(xintercept = trajectory_order), linetype = "dashed", color = "grey") +
    ggplot2::geom_smooth(size = 1.5, span = smooth_span, method = smooth_method, formula = y ~ x,
                         mapping = ggplot2::aes(color = gsets), se = smooth_se) +
    ggplot2::scale_y_continuous(breaks = base::seq(0 , 1, 0.2), labels = base::seq(0 , 1, 0.2)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "inches"))),
      axis.line.y = ggplot2::element_line()
    ) +
    ggplot2::labs(x = "Direction", y = "Expression score", color = "Gene sets")


}











