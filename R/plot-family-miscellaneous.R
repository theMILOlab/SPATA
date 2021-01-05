

# Clustering --------------------------------------------------------------


#' @title Visualize clustering results
#'
#' @description Plots a dendrogram of the distance matrix calculated via \code{runSpatialCorrelation()}.
#'
#' @inherit check_sample params
#' @inherit check_method params
#' @param ... Additional arguments given to \code{ggdendro::ggdendrogram()}
#'
#' @return plot_family return
#' @export

plotGeneDendrogram <- function(object,
                               of_sample = "",
                               method_hclust = "complete",
                               ...){

  check_object(object)
  check_method(method_hclust = method_hclust)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  sp_cor <- getSpatialCorrelationResults(object, of_sample = of_sample)

  hcluster_out <-
    stats::hclust(d = sp_cor$dist_mtr, method = method_hclust)

  ggdendro::ggdendrogram(data = hcluster_out, labels = FALSE, ...)

}


# -----


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
                       pt_clr = "black",
                       normalize = TRUE,
                       verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_pt(pt_size = pt_size, pt_alpha = pt_alpha, pt_clrsp = pt_clrsp)

  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample)

  if(!base::is.null(color_to)){

    color_to <- check_color_to(color_to = color_to,
                               all_genes = getGenes(object, in_sample = of_sample),
                               all_gene_sets = getGeneSets(object),
                               all_features = getFeatureNames(object))

  } else {

    color_to <- list("color" = pt_clr)

  }

  # -----


  # 2. Data extraction and plot preparation ---------------------------------

  dim_red_df <-
    getDimRedDf(object, method_dr = method_dr, of_sample = of_sample)

  plot_list <-
    hlpr_scatterplot(objcet = object,
                     spata_df = dim_red_df,
                     color_to = color_to,
                     pt_size = pt_size,
                     pt_alpha = pt_alpha,
                     pt_clrp = pt_clrp,
                     pt_clrsp = pt_clrsp,
                     method_gs = method_gs,
                     normalize = normalize,
                     verbose = verbose)

  # -----

  # 3. Plotting -------------------------------------------------------------

  x <- stringr::str_c(method_dr, 1, sep = "")
  y <- stringr::str_c(method_dr, 2, sep = "")

  ggplot2::ggplot(data = plot_list$data, mapping = ggplot2::aes_string(x = x, y = y)) +
    plot_list$add_on +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  # -----

}


#' @title Plot dimensional reduction
#'
#' @description Displays the dimensional reduction and maps gene, gene-set
#' or feature information onto the color-aesthetic.
#'
#' @inherit check_sample params
#' @inherit check_color_to params
#' @inherit check_method params
#' @param n_pcs Numeric value. Determines the number of principal components to be plotted.
#' Must be an even number.
#' @inherit check_pt params
#' @inherit verbose params
#'
#' @return Returns a ggplot-object that can be additionally customized according
#' to the rules of the ggplot2-framework.
#'
#' @export
#'

plotUmap <- function(object,
                     of_sample = "",
                     color_to = NULL,
                     method_gs = "mean",
                     pt_size = 2,
                     pt_alpha = 1,
                     pt_clrsp = "inferno",
                     pt_clrp = "milo",
                     pt_clr = "black",
                     normalize = TRUE,
                     verbose = TRUE){

  plotDimRed(object = object,
             method_dr = "umap",
             of_sample = of_sample,
             color_to = color_to,
             method_gs = method_gs,
             pt_size = pt_size,
             pt_alpha = pt_alpha,
             pt_clrsp = pt_clrsp,
             pt_clrp = pt_clrp,
             pt_clr = pt_clr,
             normalize = normalize,
             verbose = verbose)

}

#' @rdname plotUmap
#' @export
plotTsne <- function(object,
                     of_sample = "",
                     color_to = NULL,
                     method_gs = "mean",
                     pt_size = 2,
                     pt_alpha = 1,
                     pt_clrsp = "inferno",
                     pt_clrp = "milo",
                     normalize = TRUE,
                     verbose = TRUE){

  plotDimRed(object = object,
             method_dr = "tsne",
             of_sample = of_sample,
             color_to = color_to,
             method_gs = method_gs,
             pt_size = pt_size,
             pt_alpha = pt_alpha,
             pt_clrsp = pt_clrsp,
             pt_clrp = pt_clrp,
             pt_clr = pt_clr,
             normalize = normalize,
             verbose = verbose)

}



#' @rdname plotUmap
#' @export
plotPca <- function(object,
                    of_sample = "",
                    color_to = NULL,
                    n_pcs = 10,
                    method_gs = "mean",
                    pt_size = 1.5,
                    pt_alpha = 0.9,
                    pt_clrp = "milo",
                    pt_clrsp = "inferno",
                    pt_clr = "black",
                    normalize = TRUE,
                    verbose = TRUE,
                    ...){

  # 1. Control --------------------------------------------------------------

  confuns::make_available(..., verbose = verbose)

  # check input
  check_object(object)

  of_sample <-
    check_sample(object = object, of_sample = of_sample, of.length = 1)

  if(!base::is.null(color_to)){

    color_to <-
      check_color_to(color_to = color_to,
                     all_features = getFeatureNames(object),
                     all_genes = getGenes(object),
                     all_gene_sets = getGeneSets(object))

  } else {

    color_to <- list("color" = pt_clr)
  }

  # get data
  pca_df <- getPcaDf(object, of_sample = of_sample)

  # check principal component input
  confuns::is_value(x = n_pcs, mode = "numeric")

  n_pcs <- 1:n_pcs
  total_pcs <- (base::ncol(pca_df)-2)
  length_n_pcs <- base::length(n_pcs)

  if(length_n_pcs > total_pcs){

    base::stop(glue::glue("Input for argument 'n_pcs' ({length_n_pcs}) exceeds total number of pcs ({total_pcs})."))

  } else if(length_n_pcs %% 2 != 0){

    base::stop(glue::glue("Input for argument 'n_pcs' must be an even number."))

  }

  uneven_pcs <- stringr::str_c("PC", n_pcs[n_pcs %% 2 != 0], sep = "")
  even_pcs <- stringr::str_c("PC", n_pcs[n_pcs %% 2 == 0], sep = "")

  # -----


  # 2. Data wrangling ------------------------------------------------------

  selected_df <-
    dplyr::select(.data= pca_df,
                  barcodes, sample, dplyr::all_of(x = c(even_pcs, uneven_pcs))
    )

  even_df <-
    tidyr::pivot_longer(
      data = selected_df,
      cols = dplyr::all_of(even_pcs),
      names_to = "even_pcs",
      values_to = "y"
    ) %>% dplyr::select(-dplyr::all_of(x = c(uneven_pcs))) %>%
    dplyr::mutate(pc_numeric_even = stringr::str_remove(even_pcs, pattern = "PC") %>%
                    base::as.numeric(),
                  pc_partner = pc_numeric_even - 1 )

  uneven_df <-
    tidyr::pivot_longer(
      data = selected_df,
      cols = dplyr::all_of(uneven_pcs),
      names_to = "uneven_pcs",
      values_to = "x"
    ) %>% dplyr::select(-dplyr::all_of(even_pcs)) %>%
    dplyr::mutate(pc_numeric_uneven = stringr::str_remove(uneven_pcs, pattern = "PC") %>%
                    base::as.numeric(),
                  pc_partner = pc_numeric_uneven)

  joined_df <-
    dplyr::left_join(x = even_df, y = uneven_df, by = c("barcodes", "sample", "pc_partner")) %>%
    dplyr::mutate(pc_pairs = stringr::str_c("PC ", pc_numeric_uneven, " & ", "PC ", pc_numeric_even, sep = ""))

  unique_pc_pairs <- dplyr::pull(joined_df, var = "pc_pairs") %>% base::unique()

  dim_red_df <-
    dplyr::mutate(.data = joined_df, pc_pairs = base::factor(pc_pairs, levels = unique_pc_pairs))

  # -----

  plot_list <-
    hlpr_scatterplot(object = object,
                     spata_df = dim_red_df,
                     color_to = color_to,
                     pt_size = pt_size,
                     pt_alpha = pt_alpha,
                     pt_clrp = pt_clrp,
                     pt_clrsp = pt_clrsp,
                     method_gs = method_gs,
                     normalize = normalize,
                     verbose = verbose)

  ggplot(data = plot_list$data, mapping = ggplot2::aes(x = x, y = y)) +
    plot_list$add_on +
    confuns::call_flexibly(fn = "facet_wrap", fn.ns = "ggplot2", verbose = verbose, v.fail = ggplot2::facet_wrap(facets = . ~ pc_pairs),
                           default = list(facets = stats::as.formula(. ~ pc_pairs))) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    )

}


#' @title Plot Pca Variation
#'
#' @description Displays a scree plot of the current pca data stored.
#'
#' @inherit check_sample params
#' @inherit getPcaMtr params
#'
#' @inherit plot_family return
#' @export

plotPcaVariation <- function(object,
                             of_sample = "",
                             n_pcs = 15){


  # 1. Control --------------------------------------------------------------

  check_object(object)

  confuns::is_value(x = n_pcs, mode = "numeric")

  of_sample <- check_sample(object, of_sample = of_sample)

  # 2. Data extraction ------------------------------------------------------

  pca_mtr <- getPcaMtr(object, n_pcs = n_pcs, of_sample = of_sample)

  plot_df <-
    base::apply(X = pca_mtr, MARGIN = 2, FUN = stats::sd) %>%
    base::as.data.frame() %>%
    magrittr::set_colnames(value = "sdev") %>%
    tibble::rownames_to_column(var = "pcs") %>%
    dplyr::mutate(
      pcs = base::factor(x = pcs, levels = base::colnames(pca_mtr)),
      group = "group"
    )

  # 3. Plotting -------------------------------------------------------------

  ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(x = pcs, y = sdev)) +
    ggplot2::geom_col(mapping = ggplot2::aes(y = sdev), fill = "steelblue") +
    ggplot2::geom_path(mapping = ggplot2::aes(group = group), size = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line() ,
      panel.grid.minor.y = ggplot2::element_line()
    ) +
    ggplot2::labs(y = "Standard Deviation", x = NULL)

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
                               all_genes = all_genes)
  }

  # -----

  # 2. Data extraction ------------------------------------------------------

  data <-
    getCoordsDf(object = object,
                   of_sample = of_sample) %>%
    joinWithGeneSets(object,
                     spata_df = .,
                     gene_sets = states,
                     normalize = TRUE,
                     method_gs = method_gs,
                     verbose = verbose)

  if(!base::is.null(color_to)){

    if("genes" %in% base::names(color_to)){

      data <-
        joinWithGenes(object,
                      spata_df = data,
                      genes = color_to$genes,
                      average_genes = FALSE,
                      normalize = TRUE,
                      verbose = verbose)

    } else if("gene_sets" %in% base::names(color_to)){

      data <-
        joinWithGeneSets(object,
                         spata_df = data,
                         gene_sets = color_to$gene_sets,
                         method_gs = method_gs,
                         normalize = TRUE,
                         verbose = verbose)

    } else if("features" %in% base::names(color_to)){

      data <-
        joinWithFeatures(object,
                         spata_df = data,
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

    confuns::is_value(color_to, "character", "color_to")

    ref.input <- base::as.character(glue::glue("'color_to'-input: '{color_to}'"))

    ref.against <- base::as.character(glue::glue("'data'-variables"))

    color_to <- confuns::check_vector(
      input = color_to,
      against = base::colnames(data),
      verbose = TRUE,
      ref.input = ref.input,
      ref.against = ref.against)

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

  shifted_df <-
    tidyr::pivot_longer(
      data = data,
      cols = dplyr::all_of(states),
      names_to = "gene_set",
      values_to = "gene_set_expr"
    )

  # figure out which of the four states is a barcode's maximum
  # by filtering for it groupwise
  max_localisation <-
    dplyr::group_by(shifted_df, barcodes) %>%
    dplyr::filter(gene_set_expr == max(gene_set_expr)) %>%
    dplyr::ungroup() %>%
    # rename the remaining gene sets to 'max_gene_set'
    dplyr::select(barcodes, max_gene_set = gene_set, max_expr = gene_set_expr) %>%
    # assign the vertical localistion of the state plot depending on where the maximum occured
    dplyr::mutate(max_loc = dplyr::if_else(max_gene_set %in% states[1:2], true = "top", false = "bottom"))

  # calculate the x-position
  with_x_positions <-
    dplyr::left_join(x = data, y = max_localisation, by = "barcodes") %>%
      dplyr::mutate(
        pos_x = dplyr::case_when(
          max_loc == "top" & !!sym(states[1]) > !!sym(states[2]) ~ (log2(abs((!!sym(states[1]) - !!sym(states[2])) + 1)) * -1),
          max_loc == "top" & !!sym(states[2]) > !!sym(states[1]) ~ log2(abs((!!sym(states[2]) - !!sym(states[1])) + 1)),
          max_loc == "bottom" & !!sym(states[3]) > !!sym(states[4]) ~ (log2(abs((!!sym(states[3]) - !!sym(states[4])) + 1)) * -1),
          max_loc == "bottom" & !!sym(states[4]) > !!sym(states[3]) ~ log2(abs((!!sym(states[4]) - !!sym(states[3])) + 1)))
      )

  # calculate the y-position
  plot_df <-
    dplyr::group_by(with_x_positions, barcodes) %>%
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

  max <- base::max(base::abs(plot_df$pos_x), base::abs(plot_df$pos_y))

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



# Distribution -------------------------------------------------------


#' Title
#'
#' @param spata_df
#'
#' @return
#' @export
#'
#' @examples
#'
plotStatisticsInteractive <- function(spata_df){

  spata_df <- dplyr::select(spata_df, -dplyr::all_of(x = c("sample", "barcodes")))

  confuns::plot_statistics(df = spata_df, 25)

}


#' @title Distribution of continuous values
#'
#' @description Visualizes the distribution of values of a set of variables.
#'
#' \itemize{
#'  \item{ \code{plotDistribution()} Takes the spata-object as the starting point and creates the
#'  necessary data.frame from scratch according to additional parameters.}
#'  \item{ \code{plotDistribution2()} Takes a data.frame as the starting point.}
#'  }
#'
#'
#' @param df A data.frame that contains the numeric variables specified in \code{variables}.
#' @inherit variables_num params
#' @inherit across params
#' @param plot_type Character value. One of \emph{'histogram', 'density', 'violin', 'boxplot' and 'ridgeplot'}.
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
  all_genes <- getGenes(object = object, in_sample = of_sample)
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
                       spata_df = data,
                       features = variables[variables %in% all_features],
                       smooth = FALSE,
                       verbose = verbose
                       )

  }

  if(base::any(variables %in% all_gene_sets)){

    data <-
      joinWithGeneSets(object = object,
                       spata_df = data,
                       gene_sets = variables[variables %in% all_gene_sets],
                       method_gs = method_gs,
                       smooth = FALSE,
                       verbose = verbose)

  }

  if(base::any(variables %in% all_genes)){

    data <-
      joinWithGenes(object = object,
                    spata_df = data,
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
plotDistribution2 <- function(df,
                              variables = "all",
                              plot_type = "histogram",
                              clrp = "milo",
                              binwidth = 0.05,
                              verbose = TRUE,
                              ... ){

  # 1. Control --------------------------------------------------------------

  # lazy check
  confuns::is_value(clrp, "character", "clrp")

  stopifnot(base::is.data.frame(df))
  if(!base::is.null(variables)){confuns::is_vec(variables, "character", "variables")}

  if(!plot_type %in% c("histogram", "density", "ridgeplot", "boxplot", "violin")){

    base::stop("Argument 'plot_type' needs to be one of 'histogram', 'density', 'ridgeplot', 'boxplot', 'violin'.")

  }

  # check variable input
  confuns::is_vec(variables, "character", "variables")

  if(base::all(variables == "all")){

    if(base::isTRUE(verbose)){base::message("Argument 'variables' set to 'all'. Extracting all valid, numeric variables.")}

    cnames <- base::colnames(dplyr::select_if(.tbl = df, .predicate = base::is.numeric))

    valid_variables <- cnames[!cnames %in% c("x", "y", "umap1", "umap2", "tsne1", "tsne2")]

  } else {

    check_list <-
      purrr::map(variables, function(i){c("numeric", "integer", "double")}) %>%
      magrittr::set_names(value = variables)

    confuns::check_data_frame(
      df = df,
      var.class = check_list,
      ref = "df"
    )

    valid_variables <- variables

    if(base::isTRUE(verbose)){"All specified variables found."}

  }

  n_valid_variables <- base::length(valid_variables)
  ref <- base::ifelse(n_valid_variables > 1,
                      yes = "different variables. (This can take a few seconds.)",
                      no = "variable.")
  if(base::isTRUE(verbose)){base::message(glue::glue("Plotting {n_valid_variables} {ref}"))}

  # -----

  # 2. Shift data -----------------------------------------------------------

  expr_data <-
    tidyr::pivot_longer(
      data = dplyr::select(.data = df, dplyr::all_of(x = valid_variables)),
      cols = dplyr::all_of(x = valid_variables),
      names_to = "valid_variables",
      values_to = "values"
    )

  # -----


  # 3. Display add on -------------------------------------------------------

  expr_data$valid_variables <- hlpr_gene_set_name(string = expr_data$valid_variables)

  if(plot_type == "histogram"){

    display_add_on <-
      list(
        ggplot2::geom_histogram(mapping = ggplot2::aes(x = values, fill = valid_variables),
                                color = "black", binwidth = binwidth,
                                data = expr_data),
        ggplot2::labs(y = NULL)
      )

  } else if(plot_type == "density"){

    display_add_on <-
      list(
        ggplot2::geom_density(mapping = ggplot2::aes(x = values, fill = valid_variables),
                              color = "black", data = expr_data),
        ggplot2::labs(y = "Density")
      )

  } else if(plot_type == "ridgeplot"){

    display_add_on <-
      list(
        ggridges::geom_density_ridges(mapping = ggplot2::aes(x = values, y = valid_variables, fill = valid_variables),
                                      color = "black", alpha = 0.825, data = expr_data),
        #ggridges::theme_ridges(),
        ggplot2::labs(y = NULL)
      )

  } else if(plot_type == "violin"){

    display_add_on <-
      list(
        ggplot2::geom_violin(mapping = ggplot2::aes(x = values, y = valid_variables, fill = valid_variables),
                             color = "black", data = expr_data),
        ggplot2::labs(y = NULL),
        ggplot2::coord_flip()
      )

  } else if(plot_type == "boxplot"){

    display_add_on <-
      list(
        ggplot2::geom_boxplot(mapping = ggplot2::aes(x = values, y = valid_variables, fill = valid_variables),
                              color = "black", data = expr_data),
        ggplot2::labs(y = NULL),
        ggplot2::coord_flip()
      )

  }

  if(base::length(valid_variables) > 1 && !plot_type  %in% c("ridgeplot", "violin", "boxplot")){

    facet_add_on <-
      list(ggplot2::facet_wrap(facets = . ~ valid_variables, ...))

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


#' @title Distribution of continuous values
#'
#' @description Visualizes the distribution of values of a set of variables for the
#' whole sample across specific subgroups.
#'
#' \itemize{
#'  \item{ \code{plotDistributionAcross()} Takes the spata-object as the starting point and creates the
#'  necessary data.frame from scratch according to additional parameters.}
#'  \item{ \code{plotDistributionAcross2()} Takes a data.frame as the starting point.}
#'  }
#'
#' @inherit plotDistribution params return
#' @inherit across params
#'
#' @return
#' @export
#'

plotDistributionAcross <- function(object,
                                   of_sample = "",
                                   variables,
                                   across,
                                   across_subset = NULL,
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
                       spata_df = data,
                       features = across,
                       smooth = FALSE,
                       verbose = verbose
      )

  }


  if(base::any(variables %in% all_features)){

    data <-
      joinWithFeatures(object = object,
                       spata_df = data,
                       features = c(variables[variables %in% all_features]),
                       smooth = FALSE,
                       verbose = verbose
      )

  }

  if(base::any(variables %in% all_gene_sets)){

    data <-
      joinWithGeneSets(object = object,
                       spata_df = data,
                       gene_sets = variables[variables %in% all_gene_sets],
                       method_gs = method_gs,
                       smooth = FALSE,
                       verbose = verbose)

  }

  if(base::any(variables %in% all_genes)){

    data <-
      joinWithGenes(object = object,
                    spata_df = data,
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

  data <- hlpr_subset_across(data, across, across_subset)


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


#' @rdname plotDistributionAcross
#' @export
plotDistributionAcross2 <- function(df,
                                    variables = "all",
                                    across,
                                    across_subset = NULL,
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
    df = df,
    var.class = list(c("character", "factor")) %>% magrittr::set_names(across),
    ref = "df"
  )

  # check variable input
  confuns::is_vec(variables, "character", "variables")

  if(base::all(variables == "all")){

    if(base::isTRUE(verbose)){base::message("Argument 'variables' set to 'all'. Extracting all valid, numeric variables.")}

    cnames <- base::colnames(dplyr::select_if(.tbl = df, .predicate = base::is.numeric))

    variables <- cnames[!cnames %in% c("x", "y", "umap1", "umap2", "tsne1", "tsne2")]

  } else {

    check_list <-
      purrr::map(variables, function(i){c("numeric", "integer")}) %>%
      magrittr::set_names(value = variables)

    confuns::check_data_frame(
      df = df,
      var.class = check_list,
      ref = "df"
    )

    if(base::isTRUE(verbose)){"All specified variables found."}

  }

  # -----

  # 2. Data extraction ------------------------------------------------------

  data <-
    tidyr::pivot_longer(
      data = df,
      cols = dplyr::all_of(x = variables),
      names_to = "variables",
      values_to = "values"
    )

  data <- hlpr_subset_across(data, across, across_subset)

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


#' @title Distribution of discrete features
#'
#' @description Visualize the distribution of discrete features.
#'
#' @inherit check_sample params
#' @inherit check_features params
#' @param feature_compare Character vector or NULL. The discrete feature you want to compare the
#' features of \code{features} to.
#' @param clrp clrp params
#' @param position Character value. Given to \code{position} of \code{ggplot2::geom_bar()}. One of
#' \emph{'stack', 'dodge'} or \emph{'fill'}.
#' @param ... Additional parameters given to \code{ggplot2::facet_wrap()}.
#' @inherit plotDistribution params return
#'
#' @export

plotDistributionDiscrete <- function(object,
                                     of_sample = "",
                                     features,
                                     feature_compare = NULL,
                                     clrp = "milo",
                                     position = "fill",
                                     ...){

  # 1. Control --------------------------------------------------------------

  check_object(object)
  confuns::check_one_of(input = position,
                        against = c("fill", "dodge", "stack"),
                        ref.input = "argument 'position'")

  of_sample <- check_sample(object, of_sample = of_sample)
  features <- check_features(object, features = features, c("character", "factor"))

  if(!base::is.null(feature_compare)){

    feature_compare <- check_features(object, features = feature_compare, c("character", "factor"), 1)

    if(feature_compare %in% features){

      base::stop("Input of argument 'feature_compare' must not be in input of argument 'features'.")

    }
  }

  # ----


  # Additional checks and data extraction -----------------------------------

  if(base::is.character(feature_compare)){

    all_features <- c(features, feature_compare)
    facet_add_on <- list(ggplot2::facet_wrap(facets = . ~ features, scales = "free_x"))
    fill <- feature_compare
    theme_add_on <- list()


  } else {

    all_features <- features

    facet_add_on <- list(ggplot2::facet_wrap(facets = . ~ features, scales = "free_x", ...))

    if(base::length(all_features) > 1){

      fill = "features"

    } else {

      fill = "values"

    }

    theme_add_on <- list(ggplot2::theme(legend.position = "none"))

    if(position == "fill" & base::length(all_features) > 1){

      position <- "stack"

      base::warning("Argument 'feature_compare' is NULL. Using 'stack' for argument 'position'.")

    }

  }


  plot_df <-
    joinWithFeatures(object = object,
                     spata_df = getSpataDf(object),
                     features = all_features,
                     verbose = FALSE) %>%
    tidyr::pivot_longer(data = .,
                        cols = dplyr::all_of(features),
                        names_to = "features",
                        values_to = "values")

  # ----

  if(position == "fill"){

    scale_y_add_on <- ggplot2::scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "25", "50", "75", "100"))

    y_title <- "Percentage"

  } else {

    scale_y_add_on <- list()

    y_title <- "Count"

  }

  ggplot2::ggplot(data = plot_df) +
    ggplot2::geom_bar(position = position, color = "black",
                      mapping = ggplot2::aes(x = values, fill = .data[[fill]])) +
    facet_add_on +
    confuns::scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
    ggplot2::theme_classic() +
    theme_add_on +
    ggplot2::theme(strip.background = ggplot2::element_blank()) +
    scale_y_add_on +
    ggplot2::labs(y = y_title, x = NULL)

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
#'  \item{\emph{pseudotime}: Monocle3-Umap colored by feature 'pseudotime'. }
#'  \item{\emph{\code{color_to}}: Monocle3-Umap colored by input of \code{color_to}. (if specified)}
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
    compileCellDataSet(object = object,
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

  if(!base::is.null(color_to)){
    plot_list[[color_to]] <-
      monocle3::plot_cells(cds = cds,
                           color_cells_by = color_to,
                           ...)
  }

  base::return(plot_list)

}

# -----


# Differential gene expression --------------------------------------------


#' @title Plot differentially expressed genes
#'
#' @description Visualizes the expression of genes across subgroups in a heatmap. It either takes the results
#' from previously conducted de-analysis and uses the expression information of the spata-object to plot a heatmap displaying the
#' differently expressed genes of every group of the feature denoted in the argument \code{across}.
#'
#' @inherit check_sample params
#' @inherit across params
#' @inherit getDeResultsDf params details
#' @param n_barcode_spots The number of barcode-spots belonging to each cluster you want to
#' include in the matrix. Should be lower than the total number of barcode-spots of every cluster
#' and can be deployed in order to keep the heatmap clear and aesthetically pleasing.
#'
#' If set to NULL (the default) it is automatically computed according to the number of genes that
#' are displayed in the heatmap.
#' @param breaks Denotes the colorspectrum breaks. If set to NULL the breaks are set automatically. If a
#' numeric vector is specified it is taken as input. If a function is specified the expression matrix is
#' passed to it as the first argument and the length of \code{hm_colors} as the second argument.
#' @param genes Character vector or NULL. If you want to display specific genes irrespective of de-anaylsis results you
#' can specifiy them in \code{genes}. If \code{genes} is specified that way arguments referring to de-anylsis results are
#' ignored and only the genes specified are taken and displayed.
#' @inherit verbose params
#' @param hm_colors A vector of colors to be used for the continuous heatmap annotation..
#' @param hm_clrp The colorpanels used for the discrete heatmap annotation. Run \code{all_colorpanels()}
#' to obtain valid input options.
#' @param ... Additional parameters given to \code{pheatmap::pheatmap()}.
#'
#' @return A heatmap of class 'pheatmap'.
#' @export

plotDeHeatmap <- function(object,
                          of_sample = "",
                          across,
                          across_subset = NULL,
                          method_de = "wilcox",
                          max_adj_pval = 0.5,
                          n_highest_lfc = 50,
                          n_lowest_pval = 50,
                          breaks = NULL,
                          genes = NULL,
                          n_barcode_spots = NULL,
                          hm_clrp = "milo",
                          hm_colors = viridis::viridis(15),
                          verbose = TRUE,
                          ...){

  confuns::make_available(...)

  # 1. Control --------------------------------------------------------------

  #lazy check
  check_object(object)

  confuns::are_values("hm_clrp", mode = "character")

  if(!base::is.null(across_subset)){ confuns::is_vec(x = across_subset, mode = "character")}

  # adjusting check
  across <- check_features(object, features = across, valid_classes = c("character", "factor"), max_length = 1)
  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  # ------

  # 2. Data extraction and pipeline -----------------------------------------

  # extract information on the genes to plot
  if(base::is.character(genes)){

    genes <- check_genes(object, genes = genes)

    if(base::isTRUE(verbose)){"Argument 'genes' has been specified. Ignoring de-related arguments."}

    de_df <- NULL

  } else {

    de_df <- getDeResultsDf(object = object,
                           across = across,
                           across_subset = across_subset,
                           of_sample = of_sample,
                           max_adj_pval = max_adj_pval,
                           n_highest_lfc = n_highest_lfc,
                           n_lowest_pval = n_lowest_pval)

    across_subset <- base::unique(de_df[[across]])

    genes <- dplyr::pull(de_df, var = "gene")

  }

  # make sure that across subset is of type factor
  if(base::is.factor(across_subset)){ across_subset <- base::levels(across_subset)}


  # data.frame that provides barcode-spots and cluster belonging
  if(base::is.null(n_barcode_spots)){

    n_barcode_spots <- base::length(genes) / base::length(across_subset)

  } else {

    confuns::is_value(x = n_barcode_spots, mode = "numeric")

  }

  barcodes_df <-
    joinWithFeatures(object, spata_df = getSpataDf(object, of_sample), features = across, verbose = FALSE) %>%
    confuns::check_across_subset(df = ., across = across, across.subset = across_subset) %>%
    dplyr::group_by(!!rlang::sym(across)) %>%
    dplyr::slice_sample(n = n_barcode_spots)

  # heatmap gaps
  if(base::is.null(de_df)){

    gaps_row <- NULL

  } else {

    gaps_row <-
      dplyr::group_by(de_df, !!rlang::sym(across)) %>%
      dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(positions = base::cumsum(count)) %>%
      dplyr::pull(positions) %>%
      base::as.numeric()

    gaps_row <- gaps_row[1:(base::length(gaps_row)-1)]

  }

  gaps_col <-
    dplyr::group_by(barcodes_df, !!rlang::sym(across)) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(positions = base::cumsum(count)) %>%
    dplyr::pull(positions) %>%
    base::as.numeric()

  gaps_col <- gaps_col[1:(base::length(gaps_col)-1)]

  # heatmap annotation
  annotation_col <-
    dplyr::select(.data = barcodes_df, !!rlang::sym(across)) %>%
    base::as.data.frame()

  base::rownames(annotation_col) <- dplyr::pull(barcodes_df, barcodes)

  # determine colors used
  if(hm_clrp == "default"){

    color_vec <- NA

  } else {

    color_vec <- confuns::color_vector(clrp = hm_clrp)

  }

  if(!base::all(base::is.na(color_vec))){

    colors <- color_vec[base::seq_along(across_subset)]

    annotation_colors <-
      purrr::set_names(x = list(colors), nm = across) %>%
      purrr::map(.f = ~ purrr::set_names(x = .x, nm = across_subset))

  } else {

    annotation_colors <- NA

  }

  # -----

  # 3. Plotting -------------------------------------------------------------

  base::message("Plotting heatmap. This can take a few seconds.")

  expr_mtr <-
    getExpressionMatrix(object, of_sample)[genes, barcodes_df$barcodes]

  if(base::is.null(breaks)){

    breaks_input <-
      hlpr_breaks(mtr = expr_mtr,
                  length_out = base::length(hm_colors))

  } else if(base::is.numeric(breaks)){

    breaks_input <- breaks

  } else if(base::is.function(breaks)){

    breaks_input <- breaks(expr_mtr, base::length(hm_colors))

  }

  confuns::call_flexibly(
    fn = "pheatmap", fn.ns = "pheatmap",
    default = list(mat = expr_mtr,
                   scale = "row",
                   breaks = breaks_input,
                   annotation_col = annotation_col,
                   cluster_cols = FALSE,
                   cluster_rows = FALSE,
                   show_colnames = FALSE,
                   color = hm_colors,
                   annotation_names_col = FALSE,
                   annotation_colors = annotation_colors,
                   gaps_row = gaps_row,
                   gaps_col = gaps_col)
  )

}

# -----


# Plot segmentation -------------------------------------------------------

#' @title Plot segmentation
#'
#' @description Displays the segmentation of a specified sample that was drawn with
#' \code{SPATA::createSegmentation()}.
#'
#' @inherit check_sample params
#' @inherit check_pt params
#' @param encircle Logical. If set to TRUE the segments are enclosed in a polygon.
#'
#' @inherit plot_family return
#'
#' @export

plotSegmentation <- function(object,
                             of_sample = "",
                             encircle = TRUE,
                             pt_size = 2,
                             pt_clrp = "milo"){

  # control
  check_object(object)
  of_sample <- check_sample(object, of_sample, desired_length = 1)
  check_pt(pt_size = pt_size)

  # data extraction
  plot_df <-
    getCoordinates(object, of_sample = of_sample) %>%
    joinWithFeatures(object, spata_df = ., features = "segment", verbose = FALSE)

  segment_df <-
    dplyr::filter(plot_df, !segment %in% c("", "none")) %>%
    tidyr::drop_na()

  if(base::nrow(segment_df) == 0){base::stop(glue::glue("Sample {of_sample} has not been segmented yet."))}

  if(base::isTRUE(encircle)){

    encircle_add_on <-
      ggforce::geom_mark_hull(data = segment_df, mapping = ggplot2::aes(x = x, y = y, color = segment, fill = segment))


  } else {

    encircle_add_on <- list()

  }

  # plotting
  ggplot2::ggplot() +
    ggplot2::geom_point(data = plot_df, mapping = ggplot2::aes(x = x, y = y), size = pt_size, color = "lightgrey") +
    ggplot2::geom_point(data = segment_df, size = pt_size, mapping = ggplot2::aes(x = x, y = y, color = segment)) +
    encircle_add_on +
    confuns::scale_color_add_on(aes = "fill", variable = "discrete", clrp = pt_clrp) +
    confuns::scale_color_add_on(aes = "color", variable = "discrete", clrp = pt_clrp) +
    ggplot2::theme_void() +
    ggplot2::labs(fill = "Segments", color = "Segments")

}

# -----

# Autoencoder -------------------------------------------------------------


#' @title Plot total variance of different neural networks
#'
#' @description Visualizes the results of \code{assessAutoencoderOptions()} by displaying the
#' total variance of each combination of an activation function (such as \emph{relu, sigmoid})
#' and the number of bottleneck neurons.
#'
#' The results depend on further adjustments like number of layers, dropout and number of epochs.
#'
#' @inherit check_object params
#'
#' @return plot_family return
#' @export
#'

plotAutoencoderAssessment <- function(object, activation_subset = NULL, clrp = "milo", verbose = TRUE){

  check_object(object)

  assessment_list <- getAutoencoderAssessment(object)

  plot_df <- assessment_list$df

  if(base::is.character(activation_subset)){

    confuns::check_vector(input = activation_subset,
                          against = activation_fns,
                          ref.input = "input for argumetn 'activation_subset'",
                          ref.against = "valid activation functions.")

    plot_df <- dplyr::filter(plot_df, activation %in% activation_subset)

  }

  if(base::isTRUE(verbose)){

    msg <- glue::glue("Additional set up of neural network: \n\nEpochs: {epochs}\nDropout: {dropout}\nLayers: {layers}",
                      epochs = assessment_list$set_up$epochs,
                      dropout = assessment_list$set_up$dropout,
                      layers = glue::glue_collapse(x = assessment_list$set_up$layers, sep = ", ", last = " and "))

    base::writeLines(text = msg)

  }


  ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(x = bottleneck, y = total_var)) +
    ggplot2::geom_point(mapping = ggplot2::aes(color = activation), size = 3) +
    ggplot2::geom_line(mapping = ggplot2::aes(group = activation, color = activation), size = 1.5) +
    ggplot2::facet_wrap(facets = . ~ activation) +
    scale_color_add_on(aes = "color", clrp = clrp, variable = plot_df$activation) +
    ggplot2::theme_bw()  +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = "Number of Bottleneck Neurons", y = "Total Variance")

}

#' @title Plot scaled vs. denoised expression
#'
#' @description Compares the distribution of the expression levels of \code{genes}
#' between the scaled matrix and the denoised matrix.
#'
#' @inherit check_sample params
#' @inherit runAutoencoderDenoising params
#' @inherit check_pt params
#' @inherit plot_family return
#'
#' @details This function requires a denoised matrix in slot @@data generated
#' by \code{runAutoEncoderDenoising()} as well as a scaled matrix.
#'
#' @export

plotAutoencoderResults <- function(object,
                                   genes,
                                   mtr_name = "denoised",
                                   scales = "free",
                                   pt_size = 1.5,
                                   pt_alpha = 0.75,
                                   pt_clrp = "milo",
                                   verbose = TRUE,
                                   ...){

  confuns::make_available(..., verbose = verbose)

  # 1. Control --------------------------------------------------------------

  check_object(object)
  check_pt(pt_size = pt_size, pt_alpha = pt_alpha, pt_clrp = pt_clrp)

  of_sample <- check_sample(object, of_sample = "")
  genes <- check_genes(object, genes = genes, of.length = 2, fdb_fn = "stop")

  denoised <- getExpressionMatrix(object, of_sample = of_sample, mtr_name = mtr_name)
  scaled <- getExpressionMatrix(object, of_sample = of_sample, mtr_name = "scaled")

  # 2. Join data ------------------------------------------------------------

  plot_df <-
    base::rbind(
      data.frame(base::t(denoised[genes, ]), type = "Denoised"),
      data.frame(base::t(scaled[genes, ]), type = "Scaled")
    ) %>%
    dplyr::mutate(type = base::factor(x = type, levels = c("Scaled", "Denoised")))


  # 3. Plot -----------------------------------------------------------------

  ggplot2::ggplot(data = plot_df, ggplot2::aes(x = .data[[genes[1]]], y = .data[[genes[2]]], color = type)) +
    ggplot2::geom_point(alpha = pt_alpha, size = pt_size) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x) +
    confuns::call_flexibly(fn = "facet_wrap", fn.ns = "ggplot2", default = list(facets = stats::as.formula(. ~ type), scales = scales)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      legend.position = "none"
    ) +
    scale_color_add_on(variable = "discrete", clrp = pt_clrp)

}

# -----



# Hotspot plotting --------------------------------------------------------

#' Title
#'
#' @param object
#' @param of_sample
#' @param of_hotspots
#' @param of_gene_sets
#' @param clrp
#' @param top_n
#' @param with_ties
#'
#' @return
#' @export
#'
#' @examples
plotHotspotEnrichment <- function(object,
                                  of_sample = "",
                                  of_hotspots = "",
                                  of_gene_sets = "",
                                  clrp = "milo",
                                  top_n = 10,
                                  with_ties = FALSE){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  hotspots <- check_hotspots(object = object, of_sample = of_sample, of_hotspots = of_hotspots)

  gene_sets <- check_gene_sets(object = object, gene_sets = of_gene_sets)

  enrichment_df <-
    getHotspotSuggestions(object, of_sample = of_sample)$gse_df %>%
    dplyr::filter(ont %in% {{gene_sets}} & hotspots %in% {{hotspots}}) %>%
    dplyr::group_by(hotspots) %>%
    dplyr::mutate(enrichment_score = confuns::normalize(enrichment_score)) %>%
    dplyr::slice_max(order_by = enrichment_score, with_ties = with_ties, n = top_n) %>%
    dplyr::mutate(
      ont = forcats::as_factor(x = ont),
      ont = tidytext::reorder_within(x = ont, by = enrichment_score, within = hotspots)
    )

  ggplot2::ggplot(enrichment_df, mapping = ggplot2::aes(x = enrichment_score, y = ont)) +
    ggplot2::geom_col(mapping = ggplot2::aes(fill = hotspots)) +
    ggplot2::facet_wrap(facets = . ~ hotspots, nrow = 3, scales = "free_y", drop = TRUE) +
    tidytext::scale_y_reordered() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "none",
      panel.border = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank()
    ) +
    scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
    ggplot2::labs(x = "Enrichment Score", y = NULL)

}


#' Title
#'
#' @param object
#' @param of_sample
#' @param of_hotspots
#'
#' @return
#' @export
#'
#' @examples
plotIntraHotspotDistance <- function(object, of_sample = "", of_hotspots = ""){

  of_hotspots <- check_hotspots(object,
                                of_sample = of_sample,
                                of_hotspots = of_hotspots)

  getHotspotSuggestions(object, of_sample = of_sample)$df %>%
    dplyr::filter(hotspots %in% {{of_hotspots}}) %>%

    ggplot2::ggplot(mapping = ggplot2::aes(x = intra_cluster_dist)) +
    ggridges::geom_density_ridges(
      mapping = ggplot2::aes(fill = hotspots, y = hotspots),
      alpha = 0.825
    ) +
    ggplot2::theme_classic() +
    scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(color = "black"),
      # axis.text.y = ggplot2::element_blank(),
      #axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(),
      strip.text.y = ggplot2::element_text(angle = 0, face = "italic", size = 14),
      strip.placement = "outside",
      strip.background = ggplot2::element_rect(color = "white", fill = "white"),
      panel.spacing.y = ggplot2::unit(10, "pt"),
      legend.position = "none",
    ) +
    ggplot2::labs(x = "Intrahotspot Distance", y = NULL, fill = "Hotspot")

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
#' @examples
plotSurfaceHotspots <- function(object,
                                of_sample = "",
                                plot_type = "expression",
                                clrsp = "inferno",
                                pt_alpha = 0.5,
                                pt_size = 2,
                                pt_clr = "lightgrey",
                                smooth = TRUE,
                                smooth_span = 0.02,
                                display_points = TRUE,
                                display_image = FALSE,
                                display_labels = FALSE,
                                ...){

  hps <-
    getHotspotSuggestions(object = object, of_sample = of_sample)

  # 3. Set up add ons -------------------------------------------------------

  # density vs. expression

  if(plot_type == "expression"){

    plot_df <-
      purrr::map_df(.x = check_hotspots(object, of_sample = of_sample),
                    .f = function(htsp){

                      genes <-
                        dplyr::filter(.data = hps$df, hotspots == {{htsp}}) %>%
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

                    })

    hotspot_add_on <-
      list(
        ggplot2::geom_point(data = plot_df, mapping = ggplot2::aes(x = x, y = y, color = mean_genes),
                            size = pt_size, alpha = pt_alpha),
        ggplot2::facet_wrap(facets = . ~ hotspot, ...),
        scale_color_add_on(clrsp = clrsp),
        ggplot2::theme_void(),
        ggplot2::theme(legend.position = "none")
      )

    panel_grid_major <- ggplot2::element_blank()
    display_points <- FALSE


  } else{

    coords_df <-
      getCoordsDf(object = object, of_sample = of_sample) %>%
      dplyr::select(x, y) %>%
      dplyr::mutate(include = TRUE)

    plot_df <-
      dplyr::select(.data = hps$df, x = center_x, y = center_y) %>%
      dplyr::mutate(include = TRUE) %>%
      base::rbind(., coords_df)


    if(plot_type == "density_filled"){

      hotspot_add_on <-
        list(
          ggplot2::geom_density2d_filled(data = plot_df,
                                         mapping = ggplot2::aes(x = x, y = y),
                                         ...),
          ggplot2::scale_fill_viridis_d(option = clrsp)
        )

      panel_grid_major <- ggplot2::element_blank()
      display_points <- FALSE

    } else if(plot_type == "density_2d"){

      hotspot_add_on <-
        list(
          ggplot2::geom_density2d(data = plot_df,
                                  mapping = ggplot2::aes(x = x, y = y),
                                  ...),
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
        data = hps$info,
        mapping = ggplot2::aes(x = center_x, y = center_y, label = hotspots)
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
