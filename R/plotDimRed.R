#' @title Plot dimensional reduction
#'
#' @param object A valid spata-object.
#' @param method_dr The dimensional reduction method from which to extract the
#' positions specified as a character of length one. Needs to be one of "UMAP"
#' or "TSNE".
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
#'

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
