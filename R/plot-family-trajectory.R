#' @title Plot trajectory
#'
#' @description Displays a trajectory of a specified sample that was
#' drawn with \code{SPATA::createTrajectories()}.
#'
#' @inherit check_sample params
#' @inherit check_trajectory params
#' @inherit check_method params
#' @inherit check_color_to params
#' @inherit check_smooth params
#' @inherit check_pt params
#' @inherit check_display params
#' @inherit verbose params
#' @param sgmt_size The size of the segment arrrow specified as a numeric value.
#'
#' @inherit plot_family return
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
                           pt_clr = "red",
                           pt_clrsp = "inferno",
                           sgmt_size = 1,
                           display_image = TRUE,
                           display_title = FALSE,
                           verbose = T){


  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_method(method_gs = method_gs)
  check_pt(pt_size, pt_alpha, pt_clrsp, pt_clr = pt_clr)
  check_display(display_title, display_image)
  check_smooth(smooth = smooth, smooth_span = smooth_span)

  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)

  check_trajectory(object, trajectory_name, of_sample)

  if(!base::is.null(color_to)){

    color_to <- check_color_to(color_to = color_to,
                               all_gene_sets = getGeneSets(object),
                               all_genes = getGenes(object),
                               all_features = getFeatureNames(object),
                               max_length = 1)

  }

  # -----

  # 2. Extract data ---------------------------------------------------------

  t_object <- getTrajectoryObject(object = object, trajectory_name = trajectory_name, of_sample = of_sample)

  trajectory_bc <- dplyr::pull(t_object@compiled_trajectory_df, barcodes)
  trajectory_sgmt_df <- t_object@segment_trajectory_df

  bc_traj <-
    coordsTrajectory(object, of_sample = of_sample, trajectory_name = trajectory_name) %>%
    dplyr::pull(barcodes)

  background_df <- coordsSpatial(object, of_sample = of_sample)

  # 3. Determine trajectory geom_point layer --------------------------------

  # if of length one and feature
  if("features" %in% base::names(color_to)){

    coords_df <- joinWithFeatures(object = object,
                                  coords_df = background_df,
                                  features = color_to$features,
                                  smooth = smooth,
                                  smooth_span = smooth_span,
                                  verbose = verbose) %>%
      dplyr::filter(barcodes %in% bc_traj)

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
                                  coords_df = background_df,
                                  gene_sets = color_to$gene_sets,
                                  method_gs = method_gs,
                                  smooth = smooth,
                                  smooth_span = smooth_span,
                                  verbose = verbose) %>%
      dplyr::filter(barcodes %in% bc_traj)


    # labs-add-on
    labs_add_on <- hlpr_labs_add_on(input = color_to$gene_sets,
                                    input_str = "Gene set:",
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
                               coords_df = background_df,
                               genes = color_to$genes,
                               average_genes = TRUE,
                               smooth = smooth,
                               smooth_span = smooth_span,
                               verbose = verbose) %>%
      dplyr::filter(barcodes %in% bc_traj)

    color_str <- base::ifelse(test = base::length(color_to$genes) == 1,
                              yes = color_to$genes,
                              no = "Mean expr.\nscore")

    # labs-add-on
    labs_add_on <- hlpr_labs_add_on(input = color_to,
                                    input_str = "Genes:",
                                    color_str = color_str,
                                    display_title = display_title)

    # assemble ggplot add-on
    ggplot_add_on <- list(
      ggplot2::geom_point(data = coords_df, size = pt_size, alpha = pt_alpha,
                          mapping = ggplot2::aes_string(x = "x", y = "y", color = "mean_genes")),
      ggplot2::scale_color_viridis_c(option = pt_clrsp),
      labs_add_on
    )


  } else if(base::is.null(color_to)){

    coords_df <- dplyr::filter(background_df, barcodes %in% bc_traj)

    # labs-add-on
    if(base::isTRUE(display_title)){

      labs_add_on <- ggplot2::labs(title = glue::glue("Trajectory: {trajectory_name}."))

    } else {

      labs_add_on <- NULL

    }

    ggplot_add_on <- list(ggplot2::geom_point(data = coords_df, size = pt_size, alpha = 1, color = pt_clr,
                                              mapping = ggplot2::aes(x = x, y = y)),
                          labs_add_on)

  }

  # -----

  ggplot2::ggplot() +
    hlpr_image_add_on2(object, display_image, of_sample) +
    ggplot2::geom_point(data = background_df, ggplot2::aes(x = x, y = y),
                        alpha = base::ifelse(display_image, 0.1, 0.75),
                        color = base::ifelse(display_image, "white", "lightgrey")) +
    ggplot_add_on +
    ggplot2::geom_segment(data = trajectory_sgmt_df,
                          mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                          color = "black", size = sgmt_size,
                          arrow = ggplot2::arrow(length = ggplot2::unit(x = 0.125, "inches"))) +
    ggplot2::theme_void() +
    ggplot2::coord_equal()


}



#' @title Trajectory trends distribution
#'
#' @description Visualizes the distribution of the assessment-scores
#'  (residuals-area-under-the-curve) of a trajectory via histograms.
#'
#' @param atdf An assessed trajectory data.frame returned by
#'  \code{assessTrajectoryTrends()}.
#' @param limits The minimum and maximum auc-values to include. Given to
#' \code{ggplot2::scale_x_continuous()}.
#' @param plot_type One of \emph{'histogram', 'density', and 'ridgeplot'}.
#'
#' @inherit plot_family return
#' @export
#'

plotTrajectoryAssessment <- function(atdf,
                                     limits = c(NA, NA),
                                     plot_type = "density",
                                     binwidth = 1){

  # 1. Control --------------------------------------------------------------

  confuns::is_value(x = plot_type, mode = "character", ref = "plot_type")

  confuns::check_data_frame(
    df = atdf,
    var.class = list(
      pattern = c("character", "factor"),
      auc = c("numeric", "integer", "double")
    ))

  if("genes" %in% colnames(atdf)){

    var <- "genes"

  } else if("gene_sets" %in% colnames(atdf)){

    var <- "gene_sets"

  } else {

    base::stop("Data.frame 'atdf' must have a variable named 'gene_sets' or 'genes'.")

  }

  base::stopifnot(base::is.character(dplyr::pull(atdf, {{var}})))

  # -----

  # 2. Plotting -------------------------------------------------------------

  if(plot_type == "histogram"){

    display_add_on <- list(
      ggplot2::geom_histogram(mapping = ggplot2::aes(x = auc, fill = pattern),
                              binwidth = binwidth, color = "black", data = atdf),
      ggplot2::facet_wrap(facets = . ~ pattern)
    )

  } else if(plot_type == "density"){

    display_add_on <- list(
      ggplot2::geom_density(mapping = ggplot2::aes(x = auc, fill = pattern),
                            color = "black", data = atdf),
      ggplot2::facet_wrap(facets = . ~ pattern)
    )

  } else if(plot_type == "ridgeplot"){

    display_add_on <- list(
      ggridges::geom_density_ridges(mapping = ggplot2::aes(x = auc, y = pattern, fill = pattern),
                                    color = "black", data = atdf, alpha = 0.75),
      ggridges::theme_ridges()
    )

  } else {

    base::stop("Argument 'plot_type' needs to be one of 'histogram', 'density' or 'ridgeplot'")
  }

  # -----

  ggplot2::ggplot(data = atdf) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Area under the curve [residuals]",
                  y = NULL) +
    display_add_on +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      legend.position = "none") +
    ggplot2::scale_x_continuous(limits = limits)

}



#' @title Trajectory line plots
#'
#' @description Displays values along a trajectory direction with
#' a smoothed lineplot.
#'
#' @inherit check_sample params
#' @inherit check_features params
#' @inherit check_gene_sets params
#' @inherit check_method params
#' @inherit check_genes params
#' @inherit average_genes params
#' @inherit check_smooth params
#' @inherit verbose params
#'
#' @inherit plot_family return
#'
#' @export

plotTrajectoryFeatures <- function(object,
                                   trajectory_name,
                                   of_sample,
                                   features = "percent.mt",
                                   smooth_method = "loess",
                                   smooth_span = 0.2,
                                   smooth_se = TRUE,
                                   verbose = TRUE){


  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_smooth(smooth_span = smooth_span, smooth_method = smooth_method)

  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)
  features <- check_features(object, features = features, valid_classes = c("numeric", "integer"))

  check_trajectory(object, trajectory_name, of_sample)

  # -----

  # 2. Data wrangling -------------------------------------------------------

  t_object <-
    trajectory(object = object, trajectory = trajectory_name, of_sample = of_sample)

  result_df <-
    hlpr_summarize_trajectory_df(object = object,
                                 ctdf = t_object@compiled_trajectory_df,
                                 accuracy = 5,
                                 variables = features,
                                 verbose = verbose)  %>%
    dplyr::group_by(features) %>%
    dplyr::mutate(
      values = confuns::normalize(x = values)
    ) %>%
    dplyr::ungroup()

  vline_df <-
    result_df %>%
    dplyr::group_by(trajectory_part) %>%
    dplyr::filter(trajectory_order %in% c(base::min(trajectory_order), base::max(trajectory_order)) &
                    trajectory_part_order == 1 &
                    features == features[1])

  # -----

  ggplot2::ggplot(data = result_df, mapping = ggplot2::aes(x = trajectory_order,
                                                           y = values,
                                                           color = features)) +
    ggplot2::geom_vline(data = vline_df[-1,],
                        mapping = ggplot2::aes(xintercept = trajectory_order), linetype = "dashed", color = "grey") +
    ggplot2::geom_smooth(size = 1.5, span = smooth_span, method = smooth_method, formula = y ~ x,
                         se = smooth_se) +
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
#' @export
plotTrajectoryGenes <- function(object,
                                trajectory_name,
                                of_sample,
                                genes,
                                average_genes = FALSE,
                                smooth_method = "loess",
                                smooth_span = 0.2,
                                smooth_se = T,
                                verbose = T){


  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_smooth(smooth_span = smooth_span, smooth_method = smooth_method)

  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)

  check_trajectory(object, trajectory_name, of_sample)

  if(base::length(genes) > 5 && base::isFALSE(average_genes) && base::isTRUE(verbose)){

    base::message("In order to plot more than 5 genes we recommend 'plotTrajectoryHeatmap()'.")

  }

  if(base::isTRUE(average_genes)){

    y_title <- "Mean expression score"

    rna_assay <- exprMtr(object = object, of_sample = of_sample)
    genes <- check_genes(object, genes = genes, rna_assay = rna_assay)

    if(base::length(genes) == 1){

      average_genes <- FALSE
      base::warning("Can not average one gene. Treating 'average_genes' as FALSE.")
      y_title <- "Expression score"

    }

    labs_add_on <- hlpr_labs_add_on(input = genes,
                                    input_str = "Genes: ",
                                    color_str = NULL,
                                    display_title = TRUE)

  } else {

    rna_assay <- exprMtr(object = object, of_sample = of_sample)
    genes <- check_genes(object, genes = genes, max_length = 10, rna_assay = rna_assay)
    y_title <- "Expression score"
    labs_add_on <- NULL

  }

  # -----

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

  # adapt genes in case normalization failed in some cases
  if(!base::isTRUE(average_genes)){

    genes <- genes[genes %in% base::colnames(coords_with_genes)]

  } else {

    genes <- "mean_genes"

  }

  result_df <-
    coords_with_genes %>%
    dplyr::group_by(trajectory_part, order_binned) %>%
    dplyr::summarise(dplyr::across(.cols = dplyr::all_of(x = {{genes}}), ~ mean(., na.rm = TRUE)), .groups = "drop_last") %>%
    dplyr::mutate(trajectory_part_order = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(trajectory_order = dplyr::row_number())


  if(!average_genes){

    result_df <-
      tidyr::pivot_longer(data = result_df,
                          cols = dplyr::all_of(genes),
                          names_to = "genes",
                          values_to = "values") %>%
      dplyr::group_by(genes) %>%
      dplyr::mutate(values = confuns::normalize(x = values)) %>%
      dplyr::ungroup()

  } else {

    result_df <-
      dplyr::select(result_df, values = mean_genes, genes = mean_genes, dplyr::everything()) %>%
      dplyr::mutate(values = confuns::normalize(x = values))

  }

  vline_df <-
    result_df %>%
    dplyr::group_by(trajectory_part) %>%
    dplyr::filter(trajectory_order %in% c(base::min(trajectory_order), base::max(trajectory_order)) &
                    trajectory_part_order == 1 &
                    genes == genes[1])

  # -----

  ggplot2::ggplot(data = result_df, mapping = ggplot2::aes(x = trajectory_order,
                                                           y = values,
                                                           color = genes)) +
    ggplot2::geom_vline(data = vline_df[-1,],
                        mapping = ggplot2::aes(xintercept = trajectory_order), linetype = "dashed", color = "grey") +
    ggplot2::geom_smooth(size = 1.5, span = smooth_span, method = smooth_method, formula = y ~ x,
                         se = smooth_se) +
    ggplot2::scale_y_continuous(breaks = base::seq(0 , 1, 0.2), labels = base::seq(0 , 1, 0.2)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "inches"))),
      axis.line.y = ggplot2::element_line()
    ) +
    ggplot2::labs(x = "Direction", y = y_title, color = "Genes") +
    labs_add_on

}


#' @rdname plotTrajectoryFeatures
#' @export
plotTrajectoryGeneSets <- function(object,
                                   trajectory_name,
                                   of_sample,
                                   gene_sets = "Neftel_NPC_Comb",
                                   method_gs = "mean",
                                   smooth_method = "loess",
                                   smooth_span = 0.2,
                                   smooth_se = TRUE,
                                   verbose = TRUE){


  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_smooth(smooth_span = smooth_span, smooth_method = smooth_method)
  check_method(method_gs = method_gs)

  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)
  gene_sets <- check_gene_sets(object, gene_sets = gene_sets, max_length = 10)

  check_trajectory(object, trajectory_name, of_sample)

  # -----

  # 2. Data wrangling -------------------------------------------------------

  t_object <-
    trajectory(object = object, trajectory = trajectory_name, of_sample = of_sample)

  result_df <-
    hlpr_summarize_trajectory_df(object = object,
                                 ctdf = t_object@compiled_trajectory_df,
                                 accuracy = 5,
                                 variables = gene_sets,
                                 method_gs = method_gs,
                                 verbose = verbose,
                                 normalize = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(gene_sets) %>%
    dplyr::mutate(values = confuns::normalize(x = values)) %>%
    dplyr::ungroup()


  vline_df <-
    result_df %>%
    dplyr::group_by(trajectory_part) %>%
    dplyr::filter(trajectory_order %in% c(base::min(trajectory_order), base::max(trajectory_order)) &
                    trajectory_part_order == 1 &
                    gene_sets == gene_sets[1])

  # -----

  ggplot2::ggplot(data = result_df, mapping = ggplot2::aes(x = trajectory_order,
                                                           y = values,
                                                           color = gene_sets)) +
    ggplot2::geom_vline(data = vline_df[-1,],
                        mapping = ggplot2::aes(xintercept = trajectory_order), linetype = "dashed", color = "grey") +
    ggplot2::geom_smooth(size = 1.5, span = smooth_span, method = smooth_method, formula = y ~ x,
                         se = smooth_se) +
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



#' @title Expression dynamic in heatmap
#'
#' @description Displays gene- or gene-set-expression values along a trajectory
#' direction with a smoothed heatmap.
#'
#' @inherit check_sample params
#' @inherit check_trajectory params
#' @inherit verbose params
#' @inherit check_smooth params
#' @param variables The variables of interest specified as a character vector:
#'
#' \itemize{
#'  \item{ \strong{Gene sets}: Must be in \code{getGeneSets()}}
#'  \item{ \strong{Genes}: Must be in \code{getGenes()}}
#'  }
#'
#' All elements of the specified character vector must either belong to
#' gene sets or to genes.
#'
#' @param arrange_rows Alter the way the rows of the heatmap
#' are displayed in order to highlight patterns. Currently either \emph{'maxima'}
#' or \emph{'minima'}.
#'
#' Given to \code{accross} of \code{confuns::arrange_rows()}.
#'
#' @param show_row_names Logical. If set to TRUE the variable elements
#' will be displayed at the rownames of the heatmap.
#' @param split_columns Logial. If set to TRUE the heatmap is vertically
#' splitted according to the trajectory parts.
#' @param clrsp The color spectrum to be used. Needs to be one of \emph{'inferno', 'magma', 'plasma', 'cividis' or 'viridis'}.
#'
#' @return A drawn heatmap.
#' @export
#'

plotTrajectoryHeatmap <- function(object,
                                  of_sample,
                                  trajectory_name,
                                  variables = NULL,
                                  method_gs = "mean",
                                  show_row_names = TRUE,
                                  arrange_rows = "none",
                                  verbose = TRUE,
                                  split_columns = TRUE,
                                  borders = TRUE,
                                  smooth_span = 0.5,
                                  clrsp = "inferno"){

  # 1. Control --------------------------------------------------------------

  # all checks
  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  check_trajectory(object, trajectory_name = trajectory_name, of_sample = of_sample)
  check_method(method_gs = method_gs)
  check_pt(pt_clrsp = clrsp)

  variables <- check_variables(variables = variables,
                               all_gene_sets = getGeneSets(object),
                               all_genes = getGenes(object),
                               max_slots = 1)

  var_type <- base::names(variables)

  # -----

  # 2. Data wrangling -------------------------------------------------------

  t_object <- getTrajectoryObject(object = object,
                                  trajectory_name = trajectory_name,
                                  of_sample = of_sample)

  # join ctdf with genes and pivot it
  wide_tdf <-
    hlpr_summarize_trajectory_df(
      object = object,
      ctdf = t_object@compiled_trajectory_df,
      variables = variables[[1]],
      accuracy = 5,
      verbose = verbose) %>%
    dplyr::ungroup() %>%
    dplyr::group_by({{var_type}}) %>%
    dplyr::mutate(values = confuns::normalize(x = values)) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(id_cols = dplyr::all_of(var_type),
                       names_from = c("trajectory_part", "trajectory_order"),
                       names_sep = "_",
                       values_from = "values")

  # -----

  # 3. Heatmap column split -------------------------------------------------

  # if the heatmap is to be splitted into the trajectory parts
  n_parts <- base::length(base::unique(t_object@compiled_trajectory_df$trajectory_part))

  if(base::isTRUE(split_columns) && n_parts > 1){

    parts_subset <-
      base::colnames(x = wide_tdf) %>%
      stringr::str_subset(pattern = "[^genes]") %>%
      stringr::str_extract_all(pattern = "^Part\\s\\d{1,}") %>%
      base::unlist()

    column_split_vec <-
      base::sapply(X = base::unique(parts_subset),
                   FUN = function(i){

                     n_part <-
                       base::length(parts_subset[parts_subset == i])

                     base::rep(x = i, n_part * 10)

                   }) %>%
      base::unlist(use.names = FALSE)

    column_split_vec <- base::factor(x = column_split_vec, levels = base::unique(column_split_vec))
    column_title <- base::unique(column_split_vec)

    show_column_names <- F

    num_cn <- (base::ncol(wide_tdf) - 1) * 10
    column_names <- base::vector(mode = "character", length = num_cn)

  } else {

    column_split_vec <- NULL
    column_title <- NULL

    show_column_names <- T

    num_cn <- (base::ncol(wide_tdf) - 1) * 10
    column_names <- base::vector(mode = "character", length = num_cn)

  }


  # arrange rows
  if(base::all(arrange_rows == "maxima") | base::all(arrange_rows == "minima")){

    wide_tdf <-
      confuns::arrange_rows(df = wide_tdf,
                            across = arrange_rows,
                            verbose = verbose)

  }

  # -----

  # 4. Smooth rows ----------------------------------------------------------

  row_info <- dplyr::pull(.data = wide_tdf, var_type)

  mtr <- as.matrix(dplyr::select(.data = wide_tdf, -{{var_type}}))
  mtr_smoothed <- matrix(0, nrow = nrow(mtr), ncol = ncol(mtr) * 10)

  if(verbose){

    base::message(glue::glue("Smoothing values with smoothing span: {smooth_span}."))

  }

  for(i in 1:base::nrow(mtr)){

    x <- 1:ncol(mtr)

    values <- base::as.numeric(mtr[i,])

    y <- (values - base::min(values))/(base::max(values) - base::min(values))

    model <- stats::loess(formula = y ~ x, span = smooth_span)

    mtr_smoothed[i,] <- stats::predict(model, seq(1, base::max(x) , length.out = num_cn))

  }

  plot_df <-
    base::as.data.frame(mtr_smoothed) %>%
    dplyr::mutate({{var_type}} := row_info)

  # -----


  # Plot heatmap ------------------------------------------------------------

  row_labels <- row_info
  mtr_smoothed <- as.matrix(plot_df[,base::sapply(plot_df, base::is.numeric)])

  base::colnames(mtr_smoothed) <- column_names

  if(verbose){

    base::message("Plotting....")

  }

  if(var_type == "gene_sets"){

    row_labels <- hlpr_gene_set_name(string = row_labels)

  }

  hm <-
    ComplexHeatmap::Heatmap(
      matrix = mtr_smoothed,
      name = "Expr.\nscore",
      column_title = column_title,
      column_title_side = "bottom",
      column_names_rot = 0,
      row_title_side = "left",
      row_names_side = "left",
      cluster_rows = F,
      cluster_columns = F,
      show_column_names = show_column_names,
      show_row_names = show_row_names,
      row_labels = row_labels,
      col = viridisLite::viridis(n = 1000, option = clrsp),
      border = TRUE,
      column_split = column_split_vec
    )

  # -----

  base::return(ComplexHeatmap::draw(hm))

}



#' @title Display trajectory fit
#'
#' @description Displayes the trend of a trajectory in comparison to a variety
#' of models / mathematical curves and thus how well it fits to each of them.
#'
#' @param rtdf A ranked trajectory data.frame.
#' @param gene_set The gene set of interest specified as a character value.
#' @param gene The gene of interest specified as a character value.
#' @param display_residuals Logical. If set to TRUE the residuals are displayed
#' via a ret line.
#' @param ... Additional parameters given to \code{ggplot2::facet_wrap()}.
#'
#' @inherit plot_family return
#' @export
#'

plotTrajectoryFit <- function(rtdf,
                              variable,
                              display_residuals = FALSE,
                              ...){

  # 1. Control --------------------------------------------------------------

  check_rtdf(rtdf = rtdf,
             variable = variable)

  # -----

  # 2. Data wrangling -------------------------------------------------------

  column <- colnames(rtdf)[1]

  # extract and shift
  rtdf <- dplyr::filter(.data = rtdf, !!rlang::sym(column) == {{variable}})

  data <- dplyr::select(.data = rtdf$data[[1]], trajectory_order, values_Expression = values)

  models <-
    tidyr::pivot_longer(
      data = rtdf$models[[1]],
      cols = dplyr::starts_with("p_"),
      values_to = "values_Fitted curve",
      names_to = "pattern",
      names_prefix = "p_"
    )

  joined_df <-
    dplyr::left_join(x = models, y = data, key = "trajectory_order")

  # add residuals
  if(base::isTRUE(display_residuals)){

    residuals <-
      tidyr::pivot_longer(
        data = rtdf$residuals[[1]],
        cols = dplyr::starts_with("p_"),
        values_to = "values_Residuals",
        names_to = "pattern",
        names_prefix = "p_"
      )

    joined_df <-
      dplyr::left_join(x = joined_df,
                       y = residuals,
                       by = c("trajectory_order", "pattern"))

    legend_position = "right"

  } else {

    legend_position = "none"

  }

  # shift to plottable df
  plot_df <-
    tidyr::pivot_longer(
      data = joined_df,
      cols = dplyr::all_of(x = tidyselect::starts_with("values")),
      names_to = "origin",
      values_to = "all_values",
      names_prefix = "values_"
    ) %>%
    dplyr::mutate(
      pattern = hlpr_name_models(pattern)
    )

  # -----

  ggplot2::ggplot(mapping = ggplot2::aes(x = trajectory_order, y = all_values)) +
    ggplot2::geom_line(size = 1, alpha = 0.75, color = "blue4", linetype = "dotted",
                       data = dplyr::filter(plot_df, origin == "Fitted curve")
    ) +
    ggplot2::geom_path(
      mapping = ggplot2::aes(group = origin, color = origin),
      size = 1, data = dplyr::filter(plot_df, origin %in% c("Residuals", "Expression"))
    ) +
    ggplot2::facet_wrap(~ pattern, ...) +
    ggplot2::scale_color_manual(values = c("Expression" = "forestgreen",
                                           "Residuals" = "tomato")) +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"))),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(color = "black", size = 11),
                   legend.position = legend_position) +
    ggplot2::labs(x = "Trajectory direction", y = NULL, color = NULL, caption = variable)

}












