
#' Title
#'
#' @param object
#' @param of_sample
#' @param of_gene
#' @param pt_clr
#' @param pt_shape
#' @param pt_size
#' @param plotSurface
#'
#' @return
#' @export
#'
plotGeneCenter <- function(object,
                           genes,
                           plotSurfaceComparison = FALSE,
                           pt_clr = NULL,
                           pt_size = NULL,
                           of_sample = NA){

  check_object(object)
  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  pattern_df <- getPrResults(object, of_sample, method_pr = "hotspot")$df

  gene <- check_genes(object, genes = genes, valid_genes = base::unique(pattern_df$genes))

  pattern_df_flt <-
    dplyr::filter(pattern_df, genes %in% {{gene}})

  main <-
    confuns::call_flexibly(
      fn = "plotSurfaceComparison",
      fn.ns = "SPATA",
      default = list(object = object, of_sample = of_sample, variables = gene),
      v.fail = ggplot2::ggplot(),
      v.skip = list()
    )

  add_on <-
    ggplot2::geom_point(
      data = pattern_df_flt,
      mapping = ggplot2::aes(x = center_x, y = center_y, shape = genes),
      color = pt_clr,
      size = pt_size
    )

  if(base::identical(main, list())){

    base::return(add_on)

  } else {

    main + add_on + legend_right()

  }

}


#' Title
#'
#' @param object
#' @param of_sample
#' @param of_hotspots
#'
#' @return
#' @export

plotIntraPatternDistance <- function(object, of_pattern = "", clrp = NULL, of_sample = NA){

  check_object(object)
  hlpr_assign_arguments(object)

  of_pattern <-
    check_pattern(object, of_sample = of_sample, patterns = of_pattern)

  getPrSuggestion(object, of_sample = of_sample)$df %>%
    dplyr::filter(patterns %in% {{of_pattern}}) %>%

    ggplot2::ggplot(mapping = ggplot2::aes(x = intra_cluster_dist)) +
    ggridges::geom_density_ridges(
      mapping = ggplot2::aes(fill = patterns, y = patterns),
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
    ggplot2::labs(x = "Intrapattern Distance", y = NULL, fill = NULL)

}

#' Title
#'
#' @param object
#' @param of_sample
#'
#' @return
#' @export
#'
plotPrAssessment <- function(object, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  pr_res <- getPrResults(object, of_sample = of_sample)

  totwss <- pr_res$totwss

  ggplot2::ggplot(data = totwss, mapping = ggplot2::aes(x = k, y = tot_wss)) +
    ggplot2::geom_col(mapping = ggplot2::aes(y = tot_wss), fill = "steelblue") +
    ggplot2::geom_path(mapping = ggplot2::aes(group = 1), size = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(labels = totwss$k, breaks = totwss$k) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line() ,
      panel.grid.minor.y = ggplot2::element_line()
    ) +
    ggplot2::labs(y = "Total within sum of squares", x = "Number of Patterns")

}


#' Title
#'
#' @param object
#' @param of_sample
#' @param plot_type
#' @param display_points
#' @param pt_size
#' @param clrp
#'
#' @return
#' @export
#'
plotPrSummary <- function(object,
                          plot_type = "density_2d",
                          display_points = NULL,
                          pt_size = NULL,
                          clrp = NULL,
                          of_sample = NA){

  check_object(object)
  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  hsp_res <- getPrResults(object, of_sample = of_sample)

  expression <-
    plotSurfaceHotspots(object = object,
                        of_sample = of_sample,
                        plot_type = "expression",
                        pt_size = pt_size,
                        ncol = 1)

  density <-
    plotSurfaceHotspots(object = object,
                        of_sample = of_sample,
                        plot_type = plot_type,
                        display_labels = TRUE,
                        display_points = display_points,
                        pt_size = 2)

  ipd <-
    plotIntraPatternDistance(object = object,
                             of_sample = of_sample,
                             clrp = clrp)

  expression | (density / ipd)

}




