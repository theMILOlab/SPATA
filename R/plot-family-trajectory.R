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
  check_trajectory(object, trajectory_name)

  # adjusting check
  of_sample <- check_sample(object, of_sample, desired_length = 1)

  if(!base::is.null(color_to)){

    color_to <- check_color_to(color_to = color_to,
                               all_gene_sets = getGeneSets(object),
                               all_genes = getGenes(object),
                               all_features = getFeatureNames(object),
                               max_length = 1)

  }

  # -----

  # 2. Extract data ---------------------------------------------------------

  t_object <-
    trajectory(object = object, trajectory = trajectory_name, of_sample = of_sample)

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

    labs_add_on <- hlpr_labs_add_on(input = color_to, input_str = "Genes:",
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

    ggplot_add_on <- list(ggplot2::geom_point(data = coords_df, size = pt_size, alpha = 1, color = pt_clr,
                                              mapping = ggplot2::aes(x = x, y = y)))

  }

  # -----

  ggplot2::ggplot() +
    hlpr_image_add_on2(object, display_image, of_sample) +
    ggplot2::geom_point(data = background_df, ggplot2::aes(x = x, y = y), alpha = 0.1, color = "white") +
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

  # -----

  # 2. Data wrangling -------------------------------------------------------

  t_object <-
    trajectory(object = object, trajectory = trajectory_name, of_sample = of_sample)

  coords_with_feature <-
    t_object@compiled_trajectory_df %>%
    dplyr::mutate(order_binned = plyr::round_any(projection_length, accuracy = 5, f = floor)) %>%
    joinWithFeatures(object = object,
                     coords_df = .,
                     features = features,
                     smooth = FALSE,
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

  # -----

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

  # lazy check
  check_object(object)
  check_smooth(smooth_span = smooth_span, smooth_method = smooth_method)

  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)

  if(base::length(genes) > 5 && base::isFALSE(average_genes) && base::isTRUE(verbose)){

    base::message("In order to plot more than 5 genes we recommend 'plotTrajectoryHeatmap()'.")

  }

  if(average_genes){

    y_title <- "Mean expression score"

    rna_assay <- exprMtr(object = object, of_sample = of_sample)
    genes <- check_genes(object, genes = genes, max_length = 10, rna_assay = rna_assay)

    if(base::length(genes) == 1){

      average_genes <- FALSE
      base::warning("Can not average one gene. Treating 'average_genes' as FALSE.")
      y_title <- "Expression score"

    }

  } else {

    rna_assay <- exprMtr(object = object, of_sample = of_sample)
    genes <- check_genes(object, genes = genes, max_length = 10, rna_assay = rna_assay)


    y_title <- "Expression score"

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

  # -----

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
    ggplot2::labs(x = "Direction", y = y_title, color = "Genes")

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

  # -----

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

  # -----

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












