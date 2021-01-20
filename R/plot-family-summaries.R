



# De-analysis -------------------------------------------------------------


#' Title
#'
#' @param object
#' @param across
#' @param across_subset
#' @param method_de
#' @param max_adj_pval
#' @param clrp
#' @param plot_type
#' @param scales
#' @param of_sample
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plotDeaSummary <- function(object,
                           across = NULL,
                           across_subset = NULL,
                           relevel = NULL,
                           method_de = NULL,
                           max_adj_pval = NULL,
                           clrp = NULL,
                           plotDeaGeneCount = list(display_title = TRUE),
                           plotDeaLogFC = list(display_title = FALSE),
                           plotDeaPvalues = list(display_title = FALSE),
                           verbose = NULL,
                           of_sample = NA){


  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  confuns::check_one_of(
    input = plot_type,
    against = validPlotTypes(fn_name = "plotDeaSummary")
  )

  default_list <-
    list("object" = object,
         "max_adj_pval" = max_adj_pval,
         "method_de" = method_de,
         "across" = across,
         "across_subset" = across_subset,
         "relevel" = relevel,
         "clrp" = clrp,
         "of_sample" = of_sample,
         "verbose" = verbose)

  # 2. Plotting -------------------------------------------------------------

  dea_gene_count <-
    confuns::call_flexibly(
      fn = "plotDeaGeneCount", fn.ns = "SPATA",
      default = default_list,
      v.fail = patchwork::plot_spacer(),
      v.skip = patchwork::plot_spacer()
    )

  dea_log_fc <-
    confuns::call_flexibly(
      fn = "plotDeaLogFC", fn.ns = "SPATA",
      default = default_list,
      v.fail = patchwork::plot_spacer(),
      v.skip = patchwork::plot_spacer()
    )

  dea_pvalues <-
    confuns::call_flexibly(
      fn = "plotDeaPvalues", fn.ns = "SPATA",
      default = default_list,
      v.fail = patchwork::plot_spacer(),
      v.skip = patchwork::plot_spacer()
    )


 (patchwork::plot_spacer() / dea_gene_count) | (dea_log_fc / dea_pvalues)

}


#' Title
#'
#' @param object
#' @param across
#' @param across_subset
#' @param relevel
#' @param method_de
#' @param binwidth
#' @param clrp
#' @param plot_type
#' @param scales
#' @param limits_x
#' @param of_sample
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plotDeaPvalues <- function(object,
                           across = NULL,
                           across_subset = NULL,
                           relevel = NULL,
                           method_de = NULL,
                           max_adj_pval = NULL,
                           binwidth = NULL,
                           clrp = NULL,
                           plot_type = "histogram",
                           scales = NULL,
                           limits_x = c(NA, NA),
                           display_title = NULL,
                           of_sample = NA,
                           ...){

  confuns::make_available(...)

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  confuns::check_one_of(
    input = plot_type,
    against = validPlotTypes(fn_name = "plotDeaPvalues")
  )

  if(plot_type == "histogram"){

    default_list <-
      c(list(mapping = ggplot2::aes(fill = .data[[across]]), color = "black"),
        "binwidth" = binwidth)

  } else {

    default_list <-
      list(mapping = ggplot2::aes(fill = .data[[across]]), color = "black")

  }

  de_df <- getDeaResultsDf(object = object,
                           across = across,
                           across_subset = across_subset,
                           relevel = relevel,
                           method_de = method_de,
                           max_adj_pval = max_adj_pval,
                           of_sample = of_sample)


  # 2. Plotting -------------------------------------------------------------

  ggplot2::ggplot(data = de_df, mapping = ggplot2::aes(x = p_val_adj)) +
    confuns::call_flexibly(
      fn = stringr::str_c("geom", plot_type, sep = "_"), fn.ns = "ggplot2",
      default = default_list,
    ) +
    scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
    confuns::call_flexibly(
      fn = "facet_wrap", fn.ns = "ggplot2",
      default = list(facets = stats::as.formula(stringr::str_c("~", across)), scales = scales, drop = TRUE)
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      strip.background = ggplot2::element_blank()
      ) +
    ggplot2::scale_x_continuous(limits = limits_x) +
    ggplot2::labs(x = "Adjusted p-values", y = NULL) +
    hlpr_display_title(display_title, title = stringr::str_c("Adj. p-value cutoff:", max_adj_pval, sep = " "))

}


#' @rdname plotDeaPvalues
#' @export
plotDeaLogFC <- function(object,
                         across = NULL,
                         across_subset = NULL,
                         relevel = NULL,
                         method_de = NULL,
                         max_adj_pval = NULL,
                         binwidth = NULL,
                         clrp = NULL,
                         plot_type = "histogram",
                         scales = NULL,
                         limits_x = c(NA, NA),
                         display_title = NULL,
                         of_sample = NA,
                         ...){


  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  confuns::check_one_of(
    input = plot_type,
    against = validPlotTypes(fn_name = "plotDeaLogFC")
  )

  if(plot_type == "histogram"){

    default_list <-
      c(list(mapping = ggplot2::aes(fill = .data[[across]]), color = "black"),
        "binwidth" = binwidth)

  } else {

    default_list <-
      list(mapping = ggplot2::aes(fill = .data[[across]]), color = "black")

  }

  de_df <- getDeaResultsDf(object = object,
                           across = across,
                           across_subset = across_subset,
                           relevel = relevel,
                           method_de = method_de,
                           max_adj_pval = max_adj_pval,
                           of_sample = of_sample)


  # 2. Plotting -------------------------------------------------------------

  ggplot2::ggplot(data = de_df, mapping = ggplot2::aes(x = avg_logFC)) +
    confuns::call_flexibly(
      fn = stringr::str_c("geom", plot_type, sep = "_"), fn.ns = "ggplot2",
      default = default_list,
    ) +
    scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
    confuns::call_flexibly(
      fn = "facet_wrap", fn.ns = "ggplot2",
      default = list(facets = stats::as.formula(stringr::str_c("~", across)), scales = scales, drop = TRUE)
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      strip.background = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(limits = limits_x) +
    ggplot2::labs(x = "Average logFC-values", y = NULL) +
    hlpr_display_title(display_title, title = stringr::str_c("Adj. p-value cutoff:", max_adj_pval, sep = " "))


}

#' Title
#'
#' @param object
#' @param across
#' @param across_subset
#' @param relevel
#' @param method_de
#' @param max_adj_pval
#' @param clrp
#' @param of_sample
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plotDeaGeneCount <- function(object,
                             across = NULL,
                             across_subset = NULL,
                             relevel = FALSE,
                             method_de = NULL,
                             max_adj_pval = NULL,
                             clrp = NULL,
                             clrp_adjust = NULL,
                             display_title = NULL,
                             sort_by_count = TRUE,
                             of_sample = NA,
                             ...){


  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  dea_df <- getDeaResultsDf(object = object,
                           across = across,
                           across_subset = across_subset,
                           relevel = relevel,
                           method_de = method_de,
                           max_adj_pval = max_adj_pval,
                           of_sample = of_sample)

  if(base::isTRUE(sort_by_count)){

    dea_df <- dplyr::mutate(dea_df, {{across}} := forcats::fct_infreq(f = !!rlang::sym(across)))

  }

  # 2. Plotting -------------------------------------------------------------

  ggplot2::ggplot(data = dea_df, mapping = ggplot2::aes(x = .data[[across]])) +
    ggplot2::geom_bar(mapping = ggplot2::aes(fill = .data[[across]]), color = "black") +
    scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp, clrp.adjust = clrp_adjust, ...) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(y = "Number of differentially expressed genes") +
    hlpr_display_title(display_title, title = stringr::str_c("Adj. p-value cutoff:", max_adj_pval, sep = " "))


}


# -----




