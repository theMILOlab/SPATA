
# Dimensional reduction related -------------------------------------------
plotDimRed <- function(object,
                       method_dr,
                       of_sample = "",
                       color_to = NULL,
                       method_gs = "mean",
                       pt_size = 2,
                       pt_alpha = 1,
                       pt_clrsp = "inferno",
                       pt_clrp = "milo",
                       verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_pt(pt_size = pt_size, pt_alpha = pt_alpha, pt_clrsp = pt_clrsp)

  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample)
  color_to <- check_color_to(color_to = color_to,
                             max_length = 1,
                             all_genes = getGenes(object, in_sample = of_sample),
                             all_gene_sets = getGeneSets(object),
                             all_features = getFeatureNames(object))

  # -----


  # 2. Extract dimensional reduction ----------------------------------------

  dimRed_df <- getDimRedData(object, method_dr = method_dr, of_sample = of_sample)

  # -----

  # 3. Join data and prepare ggplot add-ons ---------------------------------

  # if of length one and feature
  if("features" %in% base::names(color_to)){

    dimRed_df <- joinWithFeatures(object = object,
                                  coords_df = dimRed_df,
                                  features = color_to$features,
                                  smooth = FALSE,
                                  verbose = verbose)

    # assemble ggplot add on
    ggplot_add_on <- list(
      ggplot2::geom_point(data = dimRed_df, size = pt_size, alpha = pt_alpha,
                          mapping = ggplot2::aes_string(x = stringr::str_c(base::tolower(method_dr), 1, sep = ""),
                                                        y = stringr::str_c(base::tolower(method_dr), 2, sep = ""),
                                                        color = color_to$features)),
      confuns::scale_color_add_on(aes = "color", clrsp = pt_clrsp, clrp = pt_clrp, variable = dplyr::pull(dimRed_df, color_to$features)),
      ggplot2::labs(color = color_to)
    )

    # if of length one and gene set
  } else if("gene_sets" %in% base::names(color_to)){

    dimRed_df <- joinWithGeneSets(object = object,
                                  coords_df = dimRed_df,
                                  gene_sets = color_to$gene_sets,
                                  method_gs = method_gs,
                                  smooth = FALSE,
                                  verbose = verbose)

    # assemble ggplot add-on
    ggplot_add_on <- list(
      ggplot2::geom_point(data = dimRed_df, size = pt_size, alpha = pt_alpha,
                          mapping = ggplot2::aes_string(x = stringr::str_c(base::tolower(method_dr), 1, sep = ""),
                                                        y = stringr::str_c(base::tolower(method_dr), 2, sep = ""),
                                                        color = color_to$gene_sets)),
      confuns::scale_color_add_on(aes = "color", clrsp = pt_clrsp),
      ggplot2::labs(color = "Expr.\nscore", title = stringr::str_c("Gene set: ", color_to, " (", method_gs, ")", sep = ""))
    )

  } else if("genes" %in% base::names(color_to)){

    rna_assay <- exprMtr(object, of_sample = of_sample)

    dimRed_df <- joinWithGenes(object = object,
                               coords_df = dimRed_df,
                               genes = color_to$genes,
                               average_genes = TRUE,
                               smooth = FALSE,
                               verbose = verbose)

    # assemble ggplot add-on
    ggplot_add_on <- list(
      ggplot2::geom_point(data = dimRed_df, size = pt_size, alpha = pt_alpha,
                          mapping = ggplot2::aes_string(x = stringr::str_c(base::tolower(method_dr), 1, sep = ""),
                                                        y = stringr::str_c(base::tolower(method_dr), 2, sep = ""),
                                                        color = "mean_genes")),
      confuns::scale_color_add_on(aes = "color", clrsp = pt_clrsp),
      ggplot2::labs(color = "Mean expr.\nscore")
    )


  } else {

    base::warning("Could not map color to specified argument 'color_to'! (Hint: Features and Gene sets need to be of length one.)")

    ggplot_add_on <- NULL

  }

  # -----

  # 4. Plotting -------------------------------------------------------------

  if(base::is.list(ggplot_add_on)){

    ggplot2::ggplot(data = dimRed_df) +
      ggplot_add_on +
      ggplot2::theme_classic() +
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

  # -----

}


#' @title Plot dimensional reduction
#'
#' @inherit check_sample params
#' @inherit check_color_to params
#' @inherit check_method params
#' @inherit check_pt params
#' @inherit verbose params
#'
#' @return Returns a ggplot-object that can be additionally customized according
#' to the rules of the ggplot2-framework.
#'
#' @export
#'

plotUMAP <- function(object,
                     of_sample = "",
                     color_to = NULL,
                     method_gs = "mean",
                     pt_size = 2,
                     pt_alpha = 1,
                     pt_clrsp = "inferno",
                     pt_clrp = "milo",
                     verbose = TRUE){

  plotDimRed(object = object,
             method_dr = "UMAP",
             of_sample = of_sample,
             color_to = color_to,
             method_gs = method_gs,
             pt_size = pt_size,
             pt_alpha = pt_alpha,
             pt_clrsp = pt_clrsp,
             pt_clrp = pt_clrp,
             verbose = verbose)

}

#' @rdname plotUMAP
#'
#' @export
plotTSNE <- function(object,
                     of_sample = "",
                     color_to = NULL,
                     method_gs = "mean",
                     pt_size = 2,
                     pt_alpha = 1,
                     pt_clrsp = "inferno",
                     pt_clrp = "milo",
                     verbose = TRUE){

  plotDimRed(object = object,
             method_dr = "TSNE",
             of_sample = of_sample,
             color_to = color_to,
             method_gs = method_gs,
             pt_size = pt_size,
             pt_alpha = pt_alpha,
             pt_clrsp = pt_clrsp,
             pt_clrp = pt_clrp,
             verbose = verbose)

}

# -----


# State plots -------------------------------------------------------------


#' @title Gene set state plot
#'
#' @description Takes four gene sets and visualizes the relative
#' expression of these four gene sets for every barcode by computing it's respective
#' x- and y- coordinates in the state plot. (See details.)
#'
#' \itemize{
#'  \item{ \code{plotFourStates()} Takes the spata-object as the starting point and creates
#'  the necessary data.frame from scratch according to additional parameters.}
#'  \item{ \code{plotFourStates2()} Takes a data.frame as input.}
#'  }
#'
#' @param data A data.frame containing at least the variables \emph{barcodes, \code{states.}}.
#' Whereby the states-variables contain the respective expression values of the specified
#' gene sets.
#' @param states The gene sets defining the four states specified as a character vector
#' of length 4.
#' @inherit check_color_to params
#' @inherit check_pt params
#' @inherit check_display params
#' @inherit verbose params
#'
#' @inherit plot_family return
#'
#' @export

plotFourStates <- function(object,
                           of_sample = "",
                           states,
                           color_to = NULL,
                           method_gs = "mean",
                           average_genes = FALSE,
                           pt_size = 1.5,
                           pt_alpha = 0.9,
                           pt_clrsp = "inferno",
                           pt_clrp = "milo",
                           display_labels = TRUE,
                           assign = FALSE,
                           assign_name,
                           verbose = TRUE){


  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)

  check_pt(pt_size, pt_alpha, pt_clrsp)
  check_assign(assign, assign_name)
  check_method(method_gs = method_gs)

  # adjusting check
  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)
  states <- check_gene_sets(object, gene_sets = states, max_length = 4)

  if(base::length(states) != 4){

    base::stop(stringr::str_c(base::length(states), "valid gene sets provided.",
                              "Need four.",sep = " "))

  }

  all_genes <- getGenes(object, in_sample = of_sample)
  all_gene_sets <- getGeneSets(object)
  all_features <- getFeatureNames(object)

  if(!base::is.null(color_to)){
    color_to <- check_color_to(color_to = color_to,
                               all_features = all_features,
                               all_gene_sets = all_gene_sets,
                               all_genes = all_genes,
                               max_length = 1)
  }

  # -----

  # 2. Data extraction ------------------------------------------------------

  data <-
    getCoordinates(object = object,
                   of_sample = of_sample) %>%
    joinWithGeneSets(object,
                     coords_df = .,
                     gene_sets = states,
                     normalize = TRUE,
                     method_gs = method_gs,
                     verbose = verbose)

  if(!base::is.null(color_to)){

    if("genes" %in% base::names(color_to)){

      data <-
        joinWithGenes(object,
                      coords_df = data,
                      genes = color_to$genes,
                      average_genes = FALSE,
                      normalize = TRUE,
                      verbose = verbose)

    } else if("gene_sets" %in% base::names(color_to)){

      data <-
        joinWithGeneSets(object,
                         coords_df = data,
                         gene_sets = color_to$gene_sets,
                         method_gs = method_gs,
                         normalize = TRUE,
                         verbose = verbose)

    } else if("features" %in% base::names(color_to)){

      data <-
        joinWithFeatures(object,
                         coords_df = data,
                         features = color_to$features,
                         verbose = verbose)

    }

  }

  # -----

  # 3. Plotting -------------------------------------------------------------

  hlpr_assign(assign = assign,
              object = list("point" = data),
              name = assign_name)

  plotFourStates2(data = data,
                 states = states,
                 color_to = base::unlist(color_to, use.names = FALSE),
                 pt_size = pt_size,
                 pt_alpha = pt_alpha,
                 pt_clrsp = pt_clrsp,
                 pt_clrp = pt_clrp,
                 display_labels = display_labels)

  # -----

}

#' @rdname plotFourStates
#' @export
#'
plotFourStates2 <- function(data,
                            states,
                            color_to = NULL,
                            pt_size = 1.5,
                            pt_alpha = 0.9,
                            pt_clrsp = "inferno",
                            pt_clrp = "milo",
                            display_labels = TRUE){


  # 1. Control --------------------------------------------------------------

  # lazy check
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

  if(!base::length(states) == 4){

    base::stop("Argument 'states' needs to be of length 4.")

  }
  if(!base::all(states %in% base::colnames(data))){

    base::stop("All elements of argument 'states' must be variables of data.frame 'data'.")

  }

  check_pt(pt_size = pt_size, pt_alpha = pt_alpha, pt_clrsp = pt_clrsp)

  # -----


  # 2. Data wrangling -------------------------------------------------------

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

  # -----



  # 3. Additional add ons ---------------------------------------------------

  states <- hlpr_gene_set_name(states)
  color_to_lab <- hlpr_gene_set_name(color_to)

  xlab <- base::bquote(paste("log2(GSV-Score "[.(states[3])]*" - GSV-Score "[.(states[4])]*")"))
  ylab <- base::bquote(paste("log2(GSV-Score "[.(states[2])]*" - GSV-Score "[.(states[1])]*")"))


  if(!base::is.null(color_to)){

    variable <- dplyr::pull(plot_df, var = {{color_to}})

  } else {

    variable <- "discrete"

  }

  # -----

  max <- base::max(plot_df$pos_x, plot_df$pos_y)

  ggplot2::ggplot(data = plot_df) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "lightgrey") +
    ggplot2::geom_hline(yintercept = 0,  linetype = "dashed", color = "lightgrey") +
    ggplot2::geom_point(mapping = ggplot2::aes_string(x = "pos_x", y = "pos_y", color = color_to),
                        size = pt_size, alpha = pt_alpha, data = plot_df) +
    ggplot2::scale_x_continuous(limits = c(-max*1.1, max*1.1), expand = c(0,0)) +
    ggplot2::scale_y_continuous(limits = c(-max*1.1, max*1.1), expand = c(0,0)) +
    confuns::scale_color_add_on(clrp = pt_clrp, clrsp = pt_clrsp, variable = variable) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = xlab, y = ylab, color = color_to_lab)

}

# -----



# Plot distribution -------------------------------------------------------


#' @title Visualize variable distribution
#'
#' @description Visualizes the distribution of values of a set of variables for the
#' whole sample across specific subgroups.
#'
#' \itemize{
#'  \item{ \code{plotDistribution()} Takes the spata-object as the starting point and creates the
#'  necessary data.frame from scratch according to additional parameters.}
#'  \item{ \code{plotDistribution2()} Takes a data.frame as the starting point.}
#'  \item{ \code{plotDistributionAcross()} Takes the spata-object as the starting point and creates the
#'  necessary data.frame from scratch according to additional parameters. It takes an additional
#'  argument which allows to display the variables value-distribution across subgroups.}
#' }
#'
#'
#' @param data The data.frame containing numeric variables.
#' @param variables The numeric variables whose distribution you want to display. Specified as
#' a character vector.
#' @param across A categorical feature across which the value-distribution is displayed. Must be
#' a \emph{character or factor}-element of \code{SPATA::getFeatureNames()}.
#' @param plot_type One of \emph{'histogram', 'density', 'violin', 'boxplot' and 'ridgeplot'}.
#' @param binwidth The binwidth to use if \code{plot_type} is specified as \emph{'histogram'}.
#' @param ... additional arguments to \code{ggplot2::facet_wrap()}
#'
#' @inherit check_sample params
#' @inherit check_method params
#' @inherit verbose params
#' @inherit normalize params
#' @inherit check_assign params
#' @inherit clrp params
#'
#' @export

plotDistribution <- function(object,
                             of_sample = "",
                             variables,
                             method_gs = "mean",
                             plot_type = "histogram",
                             clrp = "milo",
                             binwidth = 0.05,
                             ... ,
                             normalize = TRUE,
                             assign = FALSE,
                             assign_name,
                             verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  if(!plot_type %in% c("histogram", "density", "ridgeplot", "boxplot", "violin")){

    base::stop("Argument 'plot_type' needs to be one of 'histogram', 'density', 'ridgeplot', 'boxplot', 'violin'.")

  }

  if(plot_type %in% c("violin", "ridgeplot", "boxplot")){

    max_length = 10

  } else {

    max_length = 25

  }

  validation(object)
  check_assign(assign, assign_name)
  confuns::is_value(clrp, "character", "clrp")

  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)

  all_features <- getFeatureNames(object)
  all_genes <- getGenes(object)
  all_gene_sets <- getGeneSets(object)

  variables <-
    check_variables(variables = variables,
                    all_features = all_features,
                    all_gene_sets = all_gene_sets,
                    all_genes = all_genes,
                    max_length = max_length,
                    simplify = TRUE)

  # -----

  # 2. Extract and wrangle with data ----------------------------------------

  data <-
    getCoordinates(object = object,
                   of_sample = of_sample)


  if(base::any(variables %in% all_features)){

    data <-
      joinWithFeatures(object = object,
                       coords_df = data,
                       features = variables[variables %in% all_features],
                       smooth = FALSE,
                       verbose = verbose
                       )

  }

  if(base::any(variables %in% all_gene_sets)){

    data <-
      joinWithGeneSets(object = object,
                       coords_df = data,
                       gene_sets = variables[variables %in% all_gene_sets],
                       method_gs = method_gs,
                       smooth = FALSE,
                       verbose = verbose)

  }

  if(base::any(variables %in% all_genes)){

    data <-
      joinWithGenes(object = object,
                    coords_df = data,
                    genes = variables[variables %in% all_genes],
                    average_genes = FALSE,
                    verbose = verbose)

  }

  data <-
    tidyr::pivot_longer(
      data = data,
      cols = dplyr::all_of(x = variables),
      names_to = "variables",
      values_to = "values"
    )

  # -----



  # 3. Display add on -------------------------------------------------------

  data$variables <- hlpr_gene_set_name(string = data$variables)

  if(plot_type == "histogram"){

    display_add_on <-
      list(
        ggplot2::geom_histogram(mapping = ggplot2::aes(x = values, fill = variables),
                                color = "black", binwidth = binwidth,
                                data = data),
        ggplot2::theme_classic(),
        ggplot2::labs(y = NULL)
      )

  } else if(plot_type == "density"){

    display_add_on <-
      list(
        ggplot2::geom_density(mapping = ggplot2::aes(x = values, fill = variables),
                              color = "black", data = data),
        ggplot2::theme_classic(),
        ggplot2::labs(y = "Density")
      )

  } else if(plot_type == "ridgeplot"){

    display_add_on <-
      list(
        ggridges::geom_density_ridges(mapping = ggplot2::aes(x = values, y = variables, fill = variables),
                                      color = "black", alpha = 0.825, data = data),
        ggridges::theme_ridges(),
        ggplot2::scale_fill_discrete(labels = base::rev(base::unique(variables))),
        ggplot2::labs(y = NULL)
      )

  } else if(plot_type == "violin"){


    display_add_on <-
      list(
        ggplot2::geom_violin(mapping = ggplot2::aes(x = values, y = variables, fill = variables),
                             color = "black", data = data),
        ggplot2::theme_classic(),
        ggplot2::labs(y = NULL),
        ggplot2::coord_flip()
      )

  } else if(plot_type == "boxplot"){

    display_add_on <-
      list(
        ggplot2::geom_boxplot(mapping = ggplot2::aes(x = values, y = variables, fill = variables),
                              color = "black", data = data),
        ggplot2::theme_classic(),
        ggplot2::labs(y = NULL),
        ggplot2::coord_flip()
      )

  }

  if(base::length(variables) > 1 && !plot_type  %in% c("ridgeplot", "violin", "boxplot")){

    facet_add_on <-
      list(ggplot2::facet_wrap(facets = . ~ variables, ...))

  } else {

    facet_add_on <- NULL

  }

  # -----

  ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = values)) +
    display_add_on +
    facet_add_on +
    confuns::scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(color = "black"),
      axis.text.x = ggplot2::element_text(color = "black"),
      strip.text.y = ggplot2::element_text(angle = 0, face = "italic", size = 14),
      strip.placement = "outside",
      strip.background = ggplot2::element_rect(color = "white", fill = "white"),
      panel.spacing.y = ggplot2::unit(10, "pt"),
      legend.position = "none"
    ) +
    ggplot2::labs(x = NULL)

}


#' @rdname plotDistribution
#' @export
plotDistribution2 <- function(data,
                              variables = "all",
                              plot_type = "histogram",
                              clrp = "milo",
                              binwidth = 0.05,
                              ... ){

  # 1. Control --------------------------------------------------------------

  # lazy check
  confuns::is_value(clrp, "character", "clrp")

  stopifnot(base::is.data.frame(data))
  if(!base::is.null(variables)){confuns::is_vec(variables, "character", "variables")}

  if(!plot_type %in% c("histogram", "density", "ridgeplot", "boxplot", "violin")){

    base::stop("Argument 'plot_type' needs to be one of 'histogram', 'density', 'ridgeplot', 'boxplot', 'violin'.")

  }

  # check variable input
  confuns::is_vec(variables, "character", "variables")

  if(base::all(variables == "all")){

    if(base::isTRUE(verbose)){base::message("Argument 'variables' set to 'all'. Plotting all numeric variables.")}

    cnames <- base::colnames(dplyr::select_if(.tbl = data, .predicate = base::is.numeric))

    variables <- cnames[!cnames %in% c("x", "y", "umap1", "umap2")]

  } else {

    check_list <-
      purrr::map(variables, function(i){c("numeric", "integer")}) %>%
      magrittr::set_names(value = variables)

    confuns::check_data_frame(
      df = data,
      var.class = check_list,
      ref = "data"
    )

    if(base::isTRUE(verbose)){"All specified variables found."}

  }

  # -----

  # 2. Shift data -----------------------------------------------------------

  expr_data <-
    tidyr::pivot_longer(
      data = data[, variables],
      cols = dplyr::all_of(x = variables),
      names_to = "variables",
      values_to = "values"
    )

  # -----


  # 3. Display add on -------------------------------------------------------

  expr_data$variables <- hlpr_gene_set_name(string = expr_data$variables)

  if(plot_type == "histogram"){

    display_add_on <-
      list(
        ggplot2::geom_histogram(mapping = ggplot2::aes(x = values, fill = variables),
                                color = "black", binwidth = binwidth,
                                data = expr_data),
        ggplot2::labs(y = NULL)
      )

  } else if(plot_type == "density"){

    display_add_on <-
      list(
        ggplot2::geom_density(mapping = ggplot2::aes(x = values, fill = variables),
                              color = "black", data = expr_data),
        ggplot2::labs(y = "Density")
      )

  } else if(plot_type == "ridgeplot"){

    display_add_on <-
      list(
        ggridges::geom_density_ridges(mapping = ggplot2::aes(x = values, y = variables, fill = variables),
                                      color = "black", alpha = 0.825, data = expr_data),
        #ggridges::theme_ridges(),
        ggplot2::labs(y = NULL)
      )

  } else if(plot_type == "violin"){

    display_add_on <-
      list(
        ggplot2::geom_violin(mapping = ggplot2::aes(x = values, y = variables, fill = variables),
                             color = "black", data = expr_data),
        ggplot2::labs(y = NULL),
        ggplot2::coord_flip()
      )

  } else if(plot_type == "boxplot"){

    display_add_on <-
      list(
        ggplot2::geom_boxplot(mapping = ggplot2::aes(x = values, y = variables, fill = variables),
                              color = "black", data = expr_data),
        ggplot2::labs(y = NULL),
        ggplot2::coord_flip()
      )

  }

  if(base::length(variables) > 1 && !plot_type  %in% c("ridgeplot", "violin", "boxplot")){

    facet_add_on <-
      list(ggplot2::facet_wrap(facets = . ~ variables, ...))

  } else {

    facet_add_on <- NULL

  }

  theme_add_on <- base::ifelse(test = plot_type == "ridgeplot",
                               yes = list(ggridges::theme_ridges()),
                               no = list(ggplot2::theme_classic()))

  # 4. Plotting -------------------------------------------------------------

  ggplot2::ggplot(data = expr_data, mapping = ggplot2::aes(x = values)) +
    display_add_on +
    facet_add_on +
    theme_add_on +
    confuns::scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(color = "black"),
      axis.text.x = ggplot2::element_text(color = "black"),
      strip.text.y = ggplot2::element_text(angle = 0, face = "italic", size = 14),
      strip.placement = "outside",
      strip.background = ggplot2::element_rect(color = "white", fill = "white"),
      panel.spacing.y = ggplot2::unit(10, "pt"),
      legend.position = "none"
    ) +
    ggplot2::labs(x = NULL)

}

#' @rdname plotDistribution
#' @export
plotDistributionAcross <- function(object,
                                   of_sample = "",
                                   variables,
                                   across,
                                   method_gs = "mean",
                                   plot_type = "violin",
                                   binwidth = 0.05,
                                   clrp = "milo",
                                   ... ,
                                   normalize = TRUE,
                                   assign = FALSE,
                                   assign_name,
                                   verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  if(!plot_type %in% c("histogram", "density", "ridgeplot", "boxplot", "violin")){

    base::stop("Argument 'plot_type' needs to be one of 'histogram', 'density', 'ridgeplot', 'boxplot', 'violin'.")

  }

  validation(object)
  check_assign(assign, assign_name)

  confuns::is_value(clrp, "character", "clrp")

  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)
  across <- check_features(object, feature = across, valid_classes = c("character", "factor"), max_length = 1)

  all_features <- getFeatureNames(object, of_class = c("integer", "numeric"))
  all_genes <- getGenes(object, in_sample = of_sample)
  all_gene_sets <- getGeneSets(object)

  variables <-
    check_variables(variables = variables,
                    all_features = all_features,
                    all_gene_sets = all_gene_sets,
                    all_genes = all_genes,
                    simplify = TRUE)

  # -----

  # 2. Extract and wrangle with data ----------------------------------------

  data <-
    getCoordinates(object = object,
                   of_sample = of_sample)


  if(across %in% getFeatureNames(object)){

    data <-
      joinWithFeatures(object = object,
                       coords_df = data,
                       features = across,
                       smooth = FALSE,
                       verbose = verbose
      )

  }


  if(base::any(variables %in% all_features)){

    data <-
      joinWithFeatures(object = object,
                       coords_df = data,
                       features = c(variables[variables %in% all_features]),
                       smooth = FALSE,
                       verbose = verbose
      )

  }

  if(base::any(variables %in% all_gene_sets)){

    data <-
      joinWithGeneSets(object = object,
                       coords_df = data,
                       gene_sets = variables[variables %in% all_gene_sets],
                       method_gs = method_gs,
                       smooth = FALSE,
                       verbose = verbose)

  }

  if(base::any(variables %in% all_genes)){

    data <-
      joinWithGenes(object = object,
                    coords_df = data,
                    genes = variables[variables %in% all_genes],
                    average_genes = FALSE,
                    verbose = verbose)

  }

  data <-
    tidyr::pivot_longer(
      data = data,
      cols = dplyr::all_of(x = variables),
      names_to = "variables",
      values_to = "values"
    )

  # -----

  # 3. Display add on -------------------------------------------------------

  if(plot_type == "histogram"){

    display_add_on <-
      list(
        ggplot2::geom_histogram(mapping = ggplot2::aes(x = values, fill = !!rlang::sym(across)),
                                color = "black", binwidth = binwidth,
                                data = data),
        ggplot2::labs(y = NULL)
      )

  } else if(plot_type == "density"){

    display_add_on <-
      list(
        ggplot2::geom_density(mapping = ggplot2::aes(x = values, fill = !!rlang::sym(across)),
                              color = "black", data = data,alpha = 0.825),
        ggplot2::labs(y = "Density")
      )

  } else if(plot_type == "ridgeplot"){

    display_add_on <-
      list(
        ggridges::geom_density_ridges(mapping = ggplot2::aes(x = values, y = as.factor(!!rlang::sym(across)), fill = !!rlang::sym(across)),
                                      color = "black", data = data, alpha = 0.825),
        ggplot2::scale_fill_discrete(guide = ggplot2::guide_legend(reverse = TRUE)),
        ggplot2::labs(y = across, x = NULL)
      )

  } else if(plot_type == "violin"){

    display_add_on <-
      list(
        ggplot2::geom_violin(mapping = ggplot2::aes(x = !!rlang::sym(across), y = values, fill = !!rlang::sym(across)),
                             color = "black", data = data),
        ggplot2::labs(y = NULL, x = across)
      )

  } else if(plot_type == "boxplot"){

    display_add_on <-
      list(
        ggplot2::geom_boxplot(mapping = ggplot2::aes(x = !!rlang::sym(across), y = values, fill = !!rlang::sym(across)),
                              color = "black", data = data),
        ggplot2::labs(y = NULL, x = across)
      )

  }

  if(base::length(variables) > 1){

    facet_add_on <-
      list(ggplot2::facet_wrap(facets = . ~ variables, ...))

  } else {

    facet_add_on <- NULL

  }

  # -----


  # 4. Plotting -------------------------------------------------------------

  ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = values)) +
    display_add_on +
    facet_add_on +
    confuns::scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(color = "black"),
      axis.text.x = ggplot2::element_text(color = "black"),
      strip.text.y = ggplot2::element_text(angle = 0, face = "italic", size = 14),
      strip.placement = "outside",
      strip.background = ggplot2::element_rect(color = "white", fill = "white"),
      panel.spacing.y = ggplot2::unit(10, "pt")
    ) +
    ggplot2::labs(x = NULL)

}


#' @rdname plotDistribution
#' @export
plotDistributionAcross2 <- function(data,
                                    variables = "all",
                                    across,
                                    plot_type = "violin",
                                    binwidth = 0.05,
                                    clrp = "milo",
                                    ... ,
                                    normalize = TRUE,
                                    assign = FALSE,
                                    assign_name,
                                    verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  if(!plot_type %in% c("histogram", "density", "ridgeplot", "boxplot", "violin")){

    base::stop("Argument 'plot_type' needs to be one of 'histogram', 'density', 'ridgeplot', 'boxplot', 'violin'.")

  }
  if(plot_type %in% c("violin", "ridgeplot", "boxplot")){

    max_length = 10

  } else {

    max_length = 25

  }


  confuns::is_value(clrp, "character", "clrp")

  # check across input
  confuns::is_value(across, "character", "across")
  confuns::check_data_frame(
    df = data,
    var.class = list(c("character", "factor")) %>% magrittr::set_names(across),
    ref = "data"
  )

  # check variable input
  confuns::is_vec(variables, "character", "variables")

  if(base::all(variables == "all")){

    if(base::isTRUE(verbose)){base::message("Argument 'variables' set to 'all'. Plotting all numeric variables.")}

    cnames <- base::colnames(dplyr::select_if(.tbl = data, .predicate = base::is.numeric))

    variables <- cnames[!cnames %in% c("x", "y", "umap1", "umap2")]

  } else {

    check_list <-
      purrr::map(variables, function(i){c("numeric", "integer")}) %>%
      magrittr::set_names(value = variables)

    confuns::check_data_frame(
      df = data,
      var.class = check_list,
      ref = "data"
    )

    if(base::isTRUE(verbose)){"All specified variables found."}

  }

  # -----

  # 2. Data extraction ------------------------------------------------------

  data <-
    tidyr::pivot_longer(
      data = data,
      cols = dplyr::all_of(x = variables),
      names_to = "variables",
      values_to = "values"
    )

  # -----

  # 3. Display add on -------------------------------------------------------

  if(plot_type == "histogram"){

    display_add_on <-
      list(
        ggplot2::geom_histogram(mapping = ggplot2::aes(x = values, fill = !!rlang::sym(across)),
                                color = "black", binwidth = binwidth,
                                data = data),
        ggplot2::labs(y = NULL)
      )

  } else if(plot_type == "density"){

    display_add_on <-
      list(
        ggplot2::geom_density(mapping = ggplot2::aes(x = values, fill = !!rlang::sym(across)),
                              color = "black", data = data,alpha = 0.825),
        ggplot2::labs(y = "Density")
      )

  } else if(plot_type == "ridgeplot"){

    display_add_on <-
      list(
        ggridges::geom_density_ridges(mapping = ggplot2::aes(x = values, y = as.factor(!!rlang::sym(across)), fill = !!rlang::sym(across)),
                                      color = "black", data = data, alpha = 0.825),
        ggplot2::scale_fill_discrete(guide = ggplot2::guide_legend(reverse = TRUE)),
        ggplot2::labs(y = across, x = NULL)

      )

  } else if(plot_type == "violin"){

    display_add_on <-
      list(
        ggplot2::geom_violin(mapping = ggplot2::aes(x = !!rlang::sym(across), y = values, fill = !!rlang::sym(across)),
                             color = "black", data = data),
        ggplot2::labs(y = NULL, x = across)
      )

  } else if(plot_type == "boxplot"){

    display_add_on <-
      list(
        ggplot2::geom_boxplot(mapping = ggplot2::aes(x = !!rlang::sym(across), y = values, fill = !!rlang::sym(across)),
                              color = "black", data = data),
        ggplot2::labs(y = NULL, x = across)
      )

  }

  if(base::length(variables) > 1){

    facet_add_on <-
      list(ggplot2::facet_wrap(facets = . ~ variables, ...))

  } else {

    facet_add_on <- NULL

  }

  # -----

  # 4. Plotting -------------------------------------------------------------

  ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = values)) +
    display_add_on +
    facet_add_on +
    confuns::scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(color = "black"),
      axis.text.x = ggplot2::element_text(color = "black"),
      strip.text.y = ggplot2::element_text(angle = 0, face = "italic", size = 14),
      strip.placement = "outside",
      strip.background = ggplot2::element_rect(color = "white", fill = "white"),
      panel.spacing.y = ggplot2::unit(10, "pt")
    ) +
    ggplot2::labs(x = NULL)

}


#' @title Monocle3 Pseudotime
#'
#' @description A wrapper around \code{monocle3::plot_cells()}.
#'
#' @param object A valid spata-object.
#' @param use_cds_file A directory leading to a .rds file containing a valid
#' cell_data_set-object previously calculated for the specified object. Specified
#' as a character value. If set to FALSE the cell_data_set object will be created
#' from scratch.
#' @param save_cds_file A filename/directory (that does not already exists) under which the used or created cell_data_set-object
#' is going to be stored specified as a character value. Should end with .rds.
#' @param preprocess_method Given to \code{monocle3::preprocess_cds()} if \code{use_cds_file} isn't a character string.
#' @param cluster_method Given to \code{monocle3::cluster_cells()} if \code{use_cds_file} isn't a character string.
#' @param ... Additional arguments given to \code{monocle3::plot_cells()}.
#' @param verbose Logical value. If set to TRUE informative messages with respect
#' to the computational progress made will be printed.
#'
#' (Warning messages will always be printed.)
#'
#' @return Returns a list of one or two ggplot-objects that can be additionally customized according
#' to the rules of the ggplot2-framework.
#'
#' \itemize{
#'  \item{\emph{"pseudotime"}: Monocle3-Umap colored by feature 'pseudotime'. }
#'  \item{\emph{\code{color_to}}: Monocle3-Umap colored by input-character-value. (if specified)}
#' }
#'
#' @export
#'

plotPseudotime <- function(object,
                           use_cds_file = FALSE,
                           save_cds_file = FALSE,
                           preprocess_method = "PCA",
                           cluster_method = "leiden",
                           color_to = NULL,
                           ...,
                           verbose = TRUE){

  check_object(object)
  if(!base::is.null(color_to)){confuns::is_value(color_to, "character", "color_to")}

  cds <-
    hlpr_compile_cds(object = object,
                     use_cds_file = use_cds_file,
                     save_cds_file = save_cds_file,
                     preprocess_method = preprocess_method,
                     cluster_method = cluster_method,
                     verbose = verbose)

  plot_list <- list()

  plot_list[["pseudotime"]] <-
    monocle3::plot_cells(cds = cds,
                         color_cells_by = "pseudotime",
                         label_cell_groups = FALSE,
                         label_groups_by_cluster = FALSE,
                         label_branch_points = FALSE,
                         label_leaves = FALSE,
                         graph_label_size = 0,
                         ...)

  if(base::is.null(color_to)){
    plot_list[[color_to]] <-
      monocle3::plot_cells(cds = cds,
                           color_cells_by = color_to,
                           ...)
  }

  base::return(plot_list)

}


#' @title Plot segmentation
#'
#' @description Displays the segmentation of a specified sample that was drawn with
#' \code{SPATA::createSegmentation()}.
#'
#' @inherit check_sample params
#' @inherit check_pt params
#'
#' @inherit plot_family return
#'
#' @export

plotSegmentation <- function(object,
                             of_sample = "",
                             pt_size = 2){

  # control
  check_object(object)
  of_sample <- check_sample(object, of_sample, desired_length = 1)
  check_pt(pt_size = pt_size)

  # data extraction
  plot_df <-
    getCoordinates(object, of_sample = of_sample) %>%
    joinWithFeatures(object, coords_df = ., features = "segment", verbose = FALSE)

  segment_df <- dplyr::filter(plot_df, segment != "")

  if(base::nrow(segment_df) == 0){base::stop(glue::glue("Sample {of_sample} has not been segmented yet."))}

  # plotting
  ggplot2::ggplot() +
    ggplot2::geom_point(data = plot_df, mapping = ggplot2::aes(x = x, y = y), size = pt_size, color = "lightgrey") +
    ggplot2::geom_point(data = segment_df, mapping = ggplot2::aes(x = x, y = y, color = segment)) +
    ggplot2::theme_void() +
    ggplot2::labs(color = "Segments") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5)))

}



