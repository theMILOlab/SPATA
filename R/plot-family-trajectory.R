#' @title Plot trajectory
#'
#' @description Displays the spatial course of spatial trajectory that was
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
#' @inherit check_uniform_genes params
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
                           of_sample = "",
                           color_to = NULL,
                           method_gs = "mean",
                           smooth = FALSE,
                           smooth_span = 0.02,
                           pt_size = 2.5,
                           pt_alpha = 0.4,
                           pt_clr = "red",
                           pt_clrp = "milo",
                           pt_clrsp = "inferno",
                           sgmt_size = 1,
                           display_image = TRUE,
                           display_title = FALSE,
                           uniform_genes = "discard",
                           verbose = TRUE){


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
                               all_features = getFeatureNames(object))

  }

  # -----

  # 2. Extract data ---------------------------------------------------------

  t_object <- getTrajectoryObject(object = object, trajectory_name = trajectory_name, of_sample = of_sample)

  trajectory_bc <- dplyr::pull(t_object@compiled_trajectory_df, barcodes)
  trajectory_sgmt_df <- t_object@segment_trajectory_df

  bc_traj <- ctdf(t_object) %>% dplyr::pull(barcodes)

  background_df <-
    getCoordinates(object, of_sample = of_sample) %>%
    dplyr::mutate(trajectory = dplyr::if_else(barcodes %in% bc_traj, "yes", "no"))


  # 3. Determine additional layers ------------------------------------------

  # if of length one and feature
  if("features" %in% base::names(color_to)){

    labs_add_on <- hlpr_labs_add_on(input = color_to, input_str = "Feature:",
                                    color_str = color_to,
                                    display_title = display_title)

    color_to_value <- base::unlist(color_to, use.names = FALSE)

    # if of length one and gene set
  } else if("gene_sets" %in% base::names(color_to)){

    # labs-add-on
    labs_add_on <- hlpr_labs_add_on(input = color_to$gene_sets,
                                    input_str = "Gene set:",
                                    color_str = hlpr_gene_set_name(color_to$gene_sets),
                                    display_title = display_title)

    color_to_value <- base::unlist(color_to, use.names = FALSE)

  } else if("genes" %in% base::names(color_to)){

    color_str <- base::ifelse(test = base::length(color_to$genes) == 1,
                              yes = color_to$genes,
                              no = "Mean expr.\nscore")

    color_to_value <- "mean_genes"

    # labs-add-on
    labs_add_on <- hlpr_labs_add_on(input = color_to,
                                    input_str = "Genes:",
                                    color_str = color_str,
                                    display_title = display_title)

  } else if(base::is.null(color_to)){

    coords_df <- dplyr::filter(background_df, barcodes %in% bc_traj)

    # labs-add-on
    if(base::isTRUE(display_title)){

      labs_add_on <- ggplot2::labs(title = glue::glue("Trajectory: {trajectory_name}."))

    } else {

      labs_add_on <- NULL

    }

    ggplot_add_on <- list(
      ggplot2::geom_point(data = background_df, size = pt_size, color = pt_clr,
                          mapping = ggplot2::aes(x = x,  y = y, alpha = trajectory)),
      ggplot2::scale_alpha_manual(values = c("yes" = 1, "no" = pt_alpha), guide = FALSE))

  }

  if(!base::is.null(color_to)){

    background_df <-
      joinWithVariables(object = object,
                        spata_df = background_df,
                        variables = color_to,
                        method_gs = method_gs,
                        average_genes = TRUE,
                        uniform_genes = uniform_genes,
                        smooth = smooth,
                        smooth_span = smooth_span,
                        verbose = verbose)

    ggplot_add_on <- list(
      ggplot2::geom_point(data = background_df, size = pt_size,
                          mapping = ggplot2::aes(x = x,  y = y, alpha = trajectory,
                                                 color = .data[[color_to_value]])),
      ggplot2::scale_alpha_manual(values = c("yes" = 1, "no" = pt_alpha), guide = FALSE),
      confuns::scale_color_add_on(aes = "color",
                                  clrsp = pt_clrsp,
                                  clrp = pt_clrp,
                                  variable = dplyr::pull(background_df, color_to_value))
    )

  }

  # -----

  ggplot2::ggplot() +
    hlpr_image_add_on2(object, display_image, of_sample) +
    ggplot_add_on +
    ggplot2::geom_segment(data = trajectory_sgmt_df,
                          mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                          color = "black", size = sgmt_size,
                          arrow = ggplot2::arrow(length = ggplot2::unit(x = 0.125, "inches"))) +
    ggplot2::theme_void() +
    ggplot2::coord_equal() +
    labs_add_on

}


#' @title Trajectory line plots
#'
#' @description Displays values along a trajectory direction with
#' a smoothed lineplot.
#'
#' @inherit check_sample params
#' @inherit check_features params
#' @param discrete_feature Character value. The discrete feature of interest.
#' @inherit check_gene_sets params
#' @inherit check_method params
#' @inherit check_genes params
#' @inherit average_genes params
#' @inherit check_smooth params
#' @inherit verbose params
#' @param display_trajectory_parts Logical. If set to TRUE the returned plot
#' visualizes the parts in which the trajectory has been partitioned while beeing
#' drawn.
#' @param clrp Character value. The color panel to be used.
#'  Run \code{all_colorpanels()} to see valid input.
#'
#' @inherit plot_family return
#'
#' @export

plotTrajectoryFeatures <- function(object,
                                   trajectory_name,
                                   of_sample = "",
                                   features = "percent.mt",
                                   smooth_method = "loess",
                                   smooth_span = 0.2,
                                   smooth_se = TRUE,
                                   clrp = "milo",
                                   verbose = TRUE){


  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_smooth(smooth_span = smooth_span, smooth_method = smooth_method, smooth_se = smooth_se)

  confuns::is_value(clrp, "character", "clrp")

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
    dplyr::group_by(variables) %>%
    dplyr::mutate(
      values = confuns::normalize(x = values)
    ) %>%
    dplyr::ungroup()

  vline_df <-
    result_df %>%
    dplyr::group_by(trajectory_part) %>%
    dplyr::filter(trajectory_order %in% c(base::min(trajectory_order), base::max(trajectory_order)) &
                  trajectory_part_order == 1 &
                  variables == features[1])

  # -----

  ggplot2::ggplot(data = result_df, mapping = ggplot2::aes(x = trajectory_order,
                                                           y = values,
                                                           color = variables)) +
    ggplot2::geom_vline(data = vline_df[-1,],
                        mapping = ggplot2::aes(xintercept = trajectory_order), linetype = "dashed", color = "grey") +
    ggplot2::geom_smooth(size = 1.5, span = smooth_span, method = smooth_method, formula = y ~ x,
                         se = smooth_se) +
    confuns::scale_color_add_on(variable = "discrete", clrp = clrp) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "inches"))),
      axis.line.y = ggplot2::element_line()
    ) +
    ggplot2::labs(x = "Trajectory Direction", y = NULL, color = "Features")

}

#' @rdname plotTrajectoryFeatures
#' @export
plotTrajectoryFeaturesDiscrete <- function(object,
                                           trajectory_name,
                                           of_sample = "",
                                           discrete_feature,
                                           clrp = "milo",
                                           display_trajectory_parts = FALSE,
                                           verbose = TRUE,
                                           ...){


# 1. Control --------------------------------------------------------------

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, 1)
  check_trajectory(object, trajectory_name = trajectory_name, of_sample = of_sample)

  feature <- check_features(object, discrete_feature, valid_classes = c("character", "factor"), 1)

  # -----


# 2. Data wrangling -------------------------------------------------------

  cns_trajectory <-
    getTrajectoryObject(object, trajectory_name = trajectory_name)

  compiled_trajectory_df <- cns_trajectory@compiled_trajectory_df

  joined_df <- joinWith(object,
                        spata_df = compiled_trajectory_df,
                        features = feature,
                        verbose = verbose)

  plot_df <-
    dplyr::mutate(.data = joined_df,
                  order_binned = plyr::round_any(x = projection_length, accuracy = 5, f = base::floor),
                  trajectory_order = stringr::str_c(trajectory_part, order_binned, sep = "_"))

  if(base::isTRUE(display_trajectory_parts)){

    facet_add_on <- list(
      ggplot2::facet_wrap(. ~ trajectory_part, scales = "free_x", ...)
    )

  } else {

    facet_add_on <- NULL

  }

  # -----

  ggplot2::ggplot(data = plot_df) +
    ggplot2::geom_bar(mapping = ggplot2::aes(x = trajectory_order, fill = .data[[discrete_feature]]), position = "fill", width = 0.9) +
    confuns::scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
    facet_add_on +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "inches"))),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = "Trajectory Direction", y = NULL)
}



#' @rdname plotTrajectoryFeatures
#' @export
plotTrajectoryGenes <- function(object,
                                trajectory_name,
                                of_sample = "",
                                genes,
                                clrp = "milo",
                                average_genes = FALSE,
                                smooth_method = "loess",
                                smooth_span = 0.2,
                                smooth_se = TRUE,
                                display_trajectory_parts = TRUE,
                                verbose = TRUE){


  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_smooth(smooth_span = smooth_span, smooth_method = smooth_method, smooth_se = smooth_se)

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
                  spata_df = .,
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


  if(!base::isTRUE(average_genes)){

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

  if(base::isTRUE(display_trajectory_parts)){

    vline_df <-
      result_df %>%
      dplyr::group_by(trajectory_part) %>%
      dplyr::filter(trajectory_order %in% c(base::min(trajectory_order), base::max(trajectory_order)) &
                      trajectory_part_order == 1 &
                      genes == genes[1])

    trajectory_part_add_on <- list(
      ggplot2::geom_vline(data = vline_df[-1,],
                          mapping = ggplot2::aes(xintercept = trajectory_order), linetype = "dashed", color = "grey")
    )

  } else {

    trajectory_part_add_on <- NULL
  }

  # -----

  ggplot2::ggplot(data = result_df, mapping = ggplot2::aes(x = trajectory_order,
                                                           y = values,
                                                           color = genes)) +
    trajectory_part_add_on +
    ggplot2::geom_smooth(size = 1.5, span = smooth_span, method = smooth_method, formula = y ~ x,
                         se = smooth_se) +
    ggplot2::scale_y_continuous(breaks = base::seq(0 , 1, 0.2), labels = base::seq(0 , 1, 0.2)) +
    confuns::scale_color_add_on(variable = "discrete", clrp = clrp) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "inches"))),
      axis.line.y = ggplot2::element_line()
    ) +
    ggplot2::labs(x = "Trajectory Direction", y = y_title, color = "Genes") +
    labs_add_on

}


#' @rdname plotTrajectoryFeatures
#' @export
plotTrajectoryGeneSets <- function(object,
                                   trajectory_name,
                                   of_sample = "",
                                   gene_sets = "Neftel_NPC_Comb",
                                   method_gs = "mean",
                                   smooth_method = "loess",
                                   smooth_span = 0.2,
                                   smooth_se = TRUE,
                                   clrp = "milo",
                                   display_trajectory_parts = TRUE,
                                   verbose = TRUE){


  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_smooth(smooth_span = smooth_span, smooth_method = smooth_method, smooth_se = smooth_se)
  check_method(method_gs = method_gs)

  confuns::is_value(clrp, "character", "clrp")

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
    dplyr::group_by(variables) %>%
    dplyr::mutate(values = confuns::normalize(x = values)) %>%
    dplyr::ungroup()

  if(base::isTRUE(display_trajectory_parts)){

    vline_df <-
      result_df %>%
      dplyr::group_by(trajectory_part) %>%
      dplyr::filter(trajectory_order %in% c(base::min(trajectory_order), base::max(trajectory_order)) &
                      trajectory_part_order == 1 &
                      variables == gene_sets[1])

    trajectory_part_add_on <- list(
      ggplot2::geom_vline(data = vline_df[-1,],
                          mapping = ggplot2::aes(xintercept = trajectory_order), linetype = "dashed", color = "grey")
    )

  } else {

    trajectory_part_add_on <- NULL
  }

  # -----

  ggplot2::ggplot(data = result_df, mapping = ggplot2::aes(x = trajectory_order,
                                                           y = values,
                                                           color = variables)) +
    trajectory_part_add_on +
    ggplot2::geom_smooth(size = 1.5, span = smooth_span, method = smooth_method, formula = y ~ x,
                         se = smooth_se) +
    confuns::scale_color_add_on(variable = "discrete", clrp = clrp) +
    ggplot2::scale_y_continuous(breaks = base::seq(0 , 1, 0.2), labels = base::seq(0 , 1, 0.2)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "inches"))),
      axis.line.y = ggplot2::element_line()
    ) +
    ggplot2::labs(x = "Trajectory Direction", y = "Expression score", color = "Gene sets")


}



#' @title Expression dynamic in heatmap
#'
#' @description Displays variable-expression values along a trajectory
#' direction with a smoothed heatmap (from left to right).
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
#' @param show_rownames Logical. If set to TRUE the variable elements
#' will be displayed at the rownames of the heatmap.
#' @param split_columns Logial. If set to TRUE the heatmap is vertically
#' splitted according to the trajectory parts.
#' @param hm_colors A vector of colors to be used.
#' @param ... Additional parameters given to \code{pheatmap::pheatmap()}
#'
#' @return A heatmap of class 'pheatmap'.
#' @export
#'

plotTrajectoryHeatmap <- function(object,
                                  trajectory_name,
                                  of_sample = "",
                                  variables,
                                  method_gs = "mean",
                                  arrange_rows = "none",
                                  show_rownames = FALSE,
                                  show_colnames = FALSE,
                                  split_columns = TRUE,
                                  smooth_span = 0.5,
                                  verbose = TRUE,
                                  hm_colors = viridis::inferno(150),
                                  ...){

  # 1. Control --------------------------------------------------------------

  # all checks
  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  check_trajectory(object, trajectory_name = trajectory_name, of_sample = of_sample)
  check_method(method_gs = method_gs)

  variables <- check_variables(variables = variables,
                               all_gene_sets = getGeneSets(object),
                               all_genes = getGenes(object),
                               max_slots = 1)

  var_type <- "variables"
  smooth <- TRUE

  # -----

  # 2. Data wrangling -------------------------------------------------------

  t_object <- getTrajectoryObject(object = object,
                                  trajectory_name = trajectory_name,
                                  of_sample = of_sample)

  # join ctdf with genes and pivot it
  stdf <-
    hlpr_summarize_trajectory_df(
      object = object,
      ctdf = t_object@compiled_trajectory_df,
      variables = variables[[1]],
      accuracy = 5,
      verbose = verbose) %>%
    dplyr::ungroup()

  wide_tdf <-
    dplyr::group_by(.data = stdf, {{var_type}}) %>%
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

    gaps_col <-
      dplyr::select(.data = stdf, trajectory_part, trajectory_part_order) %>%
      dplyr::distinct() %>%
      dplyr::group_by(trajectory_part) %>%
      dplyr::summarise(count = dplyr::n()) %>%
      dplyr::mutate(positions = base::cumsum(count) * 10) %>%
      dplyr::pull(positions) %>%
      base::as.numeric()

  } else {

    gaps_col <- NULL

  }


  # -----

  # 4. Smooth rows ----------------------------------------------------------

  mtr <- base::as.matrix(dplyr::select(.data = wide_tdf, -{{var_type}}))
  base::rownames(mtr) <- dplyr::pull(.data = wide_tdf, var_type)

  keep <- base::apply(mtr, MARGIN = 1,
                      FUN = function(x){

                        dplyr::n_distinct(x) != 1

                      })

  n_discarded <- base::sum(!keep)

  if(base::isTRUE(smooth) && n_discarded != 0){

    discarded <- base::rownames(mtr)[!keep]

    discarded_ref <- stringr::str_c(discarded, collapse = ', ')

    mtr <- mtr[keep, ]

    base::warning(glue::glue("Discarded {n_discarded} variables due to uniform expression. (Can not smooth uniform values.): '{discarded_ref}'"))

  }

  mtr_smoothed <- matrix(0, nrow = nrow(mtr), ncol = ncol(mtr) * 10)
  base::rownames(mtr_smoothed) <- base::rownames(mtr)

  if(base::isTRUE(smooth)){

    if(verbose){

      base::message(glue::glue("Smoothing values with smoothing span: {smooth_span}."))

    }

    for(i in 1:base::nrow(mtr)){

      x <- 1:base::ncol(mtr)

      values <- base::as.numeric(mtr[i,])

      y <- (values - base::min(values))/(base::max(values) - base::min(values))

      model <- stats::loess(formula = y ~ x, span = smooth_span)

      mtr_smoothed[i,] <- stats::predict(model, seq(1, base::max(x) , length.out = base::ncol(mtr)*10))

    }

  }

  # arrange rows
  if(base::all(arrange_rows == "maxima") | base::all(arrange_rows == "minima")){

    mtr_smoothed <-
      confuns::arrange_rows(df = base::as.data.frame(mtr_smoothed),
                            across = arrange_rows,
                            verbose = verbose) %>% base::as.matrix()

  }

  # -----

  # Plot heatmap ------------------------------------------------------------

  pheatmap::pheatmap(
    mat = mtr_smoothed,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    color = hm_colors,
    gaps_col = gaps_col[1:(base::length(gaps_col)-1)],
    show_colnames = show_colnames,
    show_rownames = show_rownames,
    ...
  )

  # -----


}




#' @title Display trajectory fit
#'
#' @description Displays the trend of a trajectory in comparison to a variety
#' of models / mathematical curves.
#'
#' @inherit check_sample params
#' @inherit check_trajectory params
#' @param variable The gene or gene set of interest specified as a character value.
#' @param display_residuals Logical. If set to TRUE the residuals are displayed
#' via a red line.
#' @param ... Additional parameters given to \code{ggplot2::facet_wrap()}.
#' @inherit hlpr_summarize_trajectory_df params
#'
#' @inherit plot_family return
#' @export

plotTrajectoryFit <- function(object,
                              trajectory_name,
                              of_sample = "",
                              variable,
                              method_gs = "mean",
                              accuracy = 5,
                              display_residuals = FALSE,
                              verbose = TRUE,
                              ...){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_trajectory(object, trajectory_name, of_sample)
  check_method(method_gs = method_gs)

  # adjusting check
  of_sample <- check_sample(object, of_sample, 1)

  variable <- check_variables(variables = variable,
                              all_gene_sets = getGeneSets(object),
                              all_genes = getGenes(object, in_sample = of_sample),
                              max_length = 1,
                              max_slots = 1) %>%
    base::unlist(use.names = FALSE)

  # -----

  # 2. Data wrangling -------------------------------------------------------

  stdf <- getTrajectoryDf(object = object,
                          trajectory_name = trajectory_name,
                          of_sample = of_sample,
                          variables = variable,
                          method_gs = method_gs,
                          accuracy = accuracy,
                          verbose = verbose,
                          normalize = TRUE)


  data <- dplyr::select(.data = stdf, trajectory_order, values_Expression = values)

  models <-
    tidyr::pivot_longer(
      data = hlpr_add_models(stdf),
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
        data = hlpr_add_residuals(stdf),
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
    ggplot2::geom_line(size = 1, alpha = 0.75, color = "blue4", linetype = "solid",
                       data = dplyr::filter(plot_df, origin == "Fitted curve")
    ) +
    ggplot2::geom_path(
      mapping = ggplot2::aes(group = origin, color = origin, linetype = origin),
      size = 1, data = dplyr::filter(plot_df, origin %in% c("Residuals", "Expression"))
    ) +
    ggplot2::facet_wrap(~ pattern, ...) +
    ggplot2::scale_color_manual(values = c("Expression" = "forestgreen",
                                           "Residuals" = "tomato")) +
    ggplot2::scale_linetype_discrete(c("Residuals"= "dotted", "Expression" = "solid"), guide = FALSE) +
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













