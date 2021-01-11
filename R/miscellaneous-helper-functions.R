

# Examine data.frames -----------------------------------------------------

#' @title Examine clustering results
#'
#' @description Gives an overwiew of the cluster results of e.g. `findMonocleClusters()`.
#'
#' @param cluster_df A data.frame containing the character variable \emph{barcodes}
#' as well as additional character variables representing different clustering-results.
#'
#' E.g. the return value of \code{findMonocleClusters()}
#'
#' @return A list in which every slot represents a cluster variable and it's content
#' the unique clusters (groups) it contains.
#'
#' @export
#'

examineClusterResults <- function(cluster_df){

  confuns::check_data_frame(
    df = cluster_df,
    var.class = list(
      barcodes = "character"
    ),
    ref = "cluster_df"
  )

  dplyr::select(.data = cluster_df, -barcodes) %>%
    purrr::discard(.x = ., .p = base::is.numeric) %>%
    purrr::map(.f = function(i){base::unique(i) %>% base::sort()})

}



#' @title Examine differentially expression results
#'
#' @description Visualizes the number of differentially expressed genes
#' and/or the distribution of logFC-values per cluster.
#'
#' @inherit check_de_df params
#' @inherit clrp params
#' @param return Character value. Denotes what to plot. One of \emph{'both', 'n_genes'}
#' and \emph{'d_logFC'}
#' @param ... Additional parameters given to \code{ggplot2::facet_wrap()}
#'
#' @inherit plot_family return
#'
#' @export


examineDeResults <- function(de_df, clrp = "milo", return = "both", binwidth = 10, ...){

  check_de_df(de_df)
  check_pt(pt_clrp = clrp)

  confuns::is_value(return, mode = "character", "return")
  confuns::check_one_of(input = return,
                        against = c("both", "n_genes", "d_logFC"))

  if(return == "both"){

    ggplot2::ggplot(data = de_df, mapping = ggplot2::aes(x = cluster)) +
      ggplot2::geom_bar(mapping = ggplot2::aes(fill = cluster), color = "black") +
      scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(x = "Clusters", y = NULL, title = "Number of differentially expressed genes") +
    ggplot2::ggplot(data = de_df, mapping = ggplot2::aes(x = avg_logFC)) +
      ggplot2::geom_histogram(mapping = ggplot2::aes(fill = cluster), binwidth = binwidth, color = "black") +
      scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
      ggplot2::facet_wrap(. ~ cluster, scales = "free", ...) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        legend.position = "none",
        strip.background = ggplot2::element_blank()) +
      ggplot2::labs(x = "Average logFC-values", y = NULL, title = "Distribution of logFC-values")

  } else if(return == "n_genes"){

    ggplot2::ggplot(data = de_df, mapping = ggplot2::aes(x = cluster)) +
      ggplot2::geom_bar(mapping = ggplot2::aes(fill = cluster), color = "black") +
      scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(x = "Clusters", y = NULL, title = "Number of differentially expressed genes")

  } else if(return == "d_logFC"){

    ggplot2::ggplot(data = de_df, mapping = ggplot2::aes(x = avg_logFC)) +
      ggplot2::geom_histogram(mapping = ggplot2::aes(fill = cluster), color = "black") +
      scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
      ggplot2::facet_wrap(. ~ cluster, scales = "free", ...) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        legend.position = "none",
        strip.background = ggplot2::element_blank()) +
      ggplot2::labs(x = "Average logFC-values", y = NULL, title = "Distribution of logFC-values")

  }

}


#' @title Examine trajectory-moddeling results
#'
#' @description Visualizes the distribution of the assessment-scores
#'  (residuals-area-under-the-curve) of a trajectory.
#'
#' @inherit check_atdf params
#' @param limits The minimum and maximum auc-values to include. Given to
#' \code{ggplot2::scale_x_continuous()}.
#' @param plot_type One of \emph{'histogram', 'density', and 'ridgeplot'}.
#' @param ... additional arguments given to \code{ggplot2::facet_wrap()}.
#'
#' @inherit plot_family return
#' @export
#'

examineTrajectoryAssessment <- function(atdf,
                                        limits = c(0, 10),
                                        plot_type = "histogram",
                                        binwidth = 0.5,
                                        clrp = "milo",
                                        ...){

  # 1. Control --------------------------------------------------------------

  confuns::is_value(plot_type,"character", "plot_type")
  confuns::is_value(clrp, "character", "clrp")
  check_atdf(atdf)

  var <- "variables"

  base::stopifnot(base::is.character(dplyr::pull(atdf, {{var}})))

  # -----

  # 2. Plotting -------------------------------------------------------------

  atdf <- dplyr::filter(atdf, dplyr::between(auc, left = limits[1], right = limits[2]))

  if(plot_type == "histogram"){

    display_add_on <- list(
      ggplot2::geom_histogram(mapping = ggplot2::aes(x = auc, fill = pattern),
                              binwidth = binwidth, color = "black", data = atdf),
      ggplot2::facet_wrap(facets = . ~ pattern, ...)
    )

  } else if(plot_type == "density"){

    display_add_on <- list(
      ggplot2::geom_density(mapping = ggplot2::aes(x = auc, fill = pattern),
                            color = "black", data = atdf),
      ggplot2::facet_wrap(facets = . ~ pattern, ...)
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
    confuns::scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
    display_add_on +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      legend.position = "none")

}


# -----






