#' @include S4-generic-functions.R
NULL

#' @title This is a text dummy
#'
#' @description A member of the \code{lazy-check_*()}-family.
#'
#' @details Members of the \code{lazy-check_*()}-family take the arguments
#' of their kind - that are used in the function they are called in - and
#' checks whether these arguments input fit the requirements. They stop and return an
#' error immediately once they stumble upon something invalid. They do not alter or adjust input
#' and return TRUE if the whole function has been executed without anything
#' invalid being found.
#'
#' @return A logical value TRUE if nothing invalid has been detected or an informative
#' error message.

lazy_check_dummy <- function(){}



#################################################################################################



#' @title Check assign input
#'
#' @param assign Logical. If set to TRUE a named list will be assigned to the global
#' environment. This list contains data and information to rebuild or additionally
#' customize the output plot of this function.
#' @param assign_name The name the assigned list is supposed to have specified as
#' a single character value.
#'
#' @inherit lazy_check_dummy description details return
#' @export

check_assign <- function(assign = FALSE,
                         assign_name = character(1)){


  if(!base::is.logical(assign)){

    base::stop("Argument 'assign' needs to be logical.")

  }

  if(base::isTRUE(assign)){

    if(!base::is.character(assign_name) | !base::length(assign_name) == 1){

      base::stop("Argument 'assign_name' needs to be a single character value.")

    }

    if(assign_name == ""){

      base::stop("Argument 'assign_name' must not be ''.")

    }

    if(base::exists(x = assign_name, where = .GlobalEnv)){

      base::stop(stringr::str_c("It already exists an object named '",
                                assign_name, "' in the global environment.",
                                sep = ""))

    }


  }

  base::return(TRUE)

}



#' @title Check compiled trajectory data.frame
#'
#' @param ctdf A compiled trajectory data.frame containing the variables
#' \emph{'barcodes', 'sample', 'x', 'y', 'projection_length', 'trajectory_part'}.
#'
#' @inherit lazy_check_dummy description details return
#' @export
#'

check_compiled_trajectory_df <- function(ctdf){

  check_spata_df(spata_df = ctdf)
  check_coordinate_variables(data = ctdf, x = "x", y = "y")

  vc <- confuns::variable_classes2(data = ctdf)

  if(!base::all(c("projection_length", "trajectory_part") %in% base::names(vc))){
    base::stop("Variables must contain 'projection_length' and 'trajectory_part'.")
  }

  if(vc["projection_length"] != "numeric"){
    base::stop("Variable 'projection_length' needs to be of class numeric.")
  }

  if(vc["trajectory_part"] != "character"){
    base::stop("Variable 'projection_length' needs to be of class character.")
  }

}



#' @title Check coordinate variables
#'
#' @param data A data.frame containing the variables of interest as well
#' as the variables needed to map onto the x- and y axis of the plot.
#' @param x The name of the numeric variable to be plotted on the x axis.
#' @param y The name of the numeric variable to be plotted on the y axis.
#'
#' @inherit lazy_check_dummy description details return
#' @export

check_coordinate_variables <- function(data, x = "x", y = "y"){

  if(!base::all(c(x, y) %in% base::colnames(data))){

    base::stop(glue::glue("Provided data.frame needs to have numeric variables '{x}' and '{y}'."))

  }

  if(base::any(!base::sapply(X = data[,c(x,y)], FUN = base::is.numeric))){

    base::stop(glue::glue("Both variables '{x}' and '{y}' of 'data' need to be numeric."))

  }

  base::return(TRUE)

}



#' @title Check coordinate data.frame
#'
#' @param spata_df A data.frame that contains at least the character variables
#' \emph{samples} and \emph{barcodes}.
#'
#' @inherit lazy_check_dummy description details return
#' @export

check_spata_df <- function(spata_df){

  if(!base::is.data.frame(spata_df)){

    base::stop("Argument 'spata_df' needs to be a data.frame.")

  } else if(!base::all(c("barcodes", "sample") %in% base::colnames(spata_df))){

    base::stop("'spata_df' needs to have 'barcodes' and 'sample' variables.")

  } else {

    classes <- base::sapply(X = spata_df[,c("barcodes", "sample")],
                            FUN = base::class)

    if(!base::all(classes == "character")){

      base::stop("Variables 'barcodes' and 'sample' need to be of class character.")

    } else {

      base::return(TRUE)

    }

  }

}

#' @rdname check_spata_df
#' @export
check_coords_df <- function(coords_df){

  warning("check_coords_df() is deprecated")

  if(!base::is.data.frame(spata_df)){

    base::stop("Argument 'spata_df' needs to be a data.frame.")

  } else if(!base::all(c("barcodes", "sample") %in% base::colnames(spata_df))){

    base::stop("'spata_df' needs to have 'barcodes' and 'sample' variables.")

  } else {

    classes <- base::sapply(X = spata_df[,c("barcodes", "sample")],
                            FUN = base::class)

    if(!base::all(classes == "character")){

      base::stop("Variables 'barcodes' and 'sample' need to be of class character.")

    } else {

      base::return(TRUE)

    }

  }

}



#' @title Check display input
#'
#' @param display_image Logical. If set to TRUE the histology image of the specified sample
#' is displayed underneath the plot.
#'
#' @param display_title Logical. If set to TRUE an informative title is displayed.
#'
#' @inherit lazy_check_dummy description details return
#' @export

check_display <- function(display_title = FALSE,
                          display_image = FALSE){

  if(!base::is.logical(display_title)){

    base::stop("Argument 'display_title' needs to be logical.")

  }

  if(!base::is.logical(display_image)){

    base::stop("Argument 'display_image' needs to be logical.")

  }

}



#' @title Check feature data.frame
#'
#' @param feature_name Character value. The name of the feature that is to be added
#'  to the object.
#' @param feature_df A data.frame that contains the feature and the key-variables.
#'
#' @inherit lazy_check_dummy description details return
#' @export

check_feature_df <- function(feature_name,
                             feature_df){

  if(!base::length(feature_name) == 1 | !base::is.character(feature_name)){

    base::stop("Argument 'feature_name' needs to be a single character value.")

  }

  if(!base::is.data.frame(feature_df)){

    base::stop("Argument 'feature_df' needs to be a data.frame.")

  } else if(!base::all(c("barcodes", "sample", feature_name) %in% base::colnames(feature_df))){

    base::stop(glue::glue("Data.frame 'feature_df' needs to have the variables 'barcodes', 'sample' and '{{feature_name}}'."))

  } else {

    classes <- base::sapply(X = feature_df[,c("barcodes", "sample")],
                            FUN = base::class)

    if(!base::all(classes == "character")){

      base::stop("Variables 'barcodes' and 'sample' need to be of class character.")

    } else {

      base::return(TRUE)

    }

  }

}



#' @title Check method input
#'
#' @param method_dr Character value. The dimensional reduction method of
#' interest specified as a single character value. (Currently
#' either \emph{'UMAP'} or \emph{'TSNE'}).
#' @param method_gs Character value. The method according to which gene sets will be handled
#' specified as a character of length one. This can be either \emph{'mean'} or one
#' of \emph{'gsva', 'ssgsea', 'zscore', or 'plage'}. The latter four will be given to
#' \code{gsva::GSVA()}.
#' @param method_padj Character value. The method according to which the adjusted p-values will
#' be calculated. Given to \code{stats::p.adjust()}.
#' @param method_de Character value. Given to argument \code{test.use} of \code{Seurat::FindAllMarkers()}.
#' @param method_ovl Character value. One of \emph{'classic', 'bayesian'}. Decides
#' according to which method the spatial overlap is calculated.
#' @inherit lazy_check_dummy description details return
#' @export

check_method <- function(method_dr = NULL,
                         method_gs = NULL,
                         method_de = NULL,
                         method_padj = NULL,
                         method_ovl = NULL){

  # dimensional reduction methods -------------------------------------------

  if(!base::is.null(method_dr)){

    if(!base::is.character(method_dr) || base::length(method_dr) != 1){

      stop("Argument 'method_dr' needs to be a single character value.")

    } else if(!method_dr %in% c("UMAP", "TSNE")) {

      stop("Argument 'method_dr' needs to be  'UMAP' or 'TSNE'.")

    }

  }

  # -----

  # gene set methods --------------------------------------------------------

  if(!base::is.null(method_gs)){

    if(!base::is.character(method_gs) || base::length(method_gs) != 1){

      stop("Argument 'method_gs' needs to be a single character value.")

    } else if(!method_gs %in% c("mean", "gsva", "ssgsea", "zscore", "plage")) {

      stop("Argument 'method_dr' needs to be  one of: 'mean', 'gsva', 'ssgsea', 'zscore', 'plage'.")

    }

  }

  # -----


  # adjusted pvalue methods -------------------------------------------------

  if(!base::is.null(method_padj)){

    if(!base::is.character(method_padj) || base::length(method_padj) != 1){

      stop("Argument 'method_padj' needs to be a single character value.")

    } else if(!method_padj %in% stats::p.adjust.methods) {

      stop("Argument 'method_padj' needs to be  one of: 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none'.")

    }

  }

  # -----


  # differential expression methods -----------------------------------------

  if(!base::is.null(method_de)){

    if(!base::is.character(method_de) || base::length(method_de) != 1){

      stop("Argument 'method_de' needs to be a single character value.")

    } else if(!method_de %in% c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2")) {

      stop("Argument 'method_de' needs to be  one of: 'wilcox', 'bimod', 'roc', 't', 'negbinom', 'poisson', 'LR', 'MAST', 'DESeq2'.")

    }

  }

  # -----

  # find overlap methods ----------------------------------------------------

  if(!base::is.null(method_ovl)){

    if(!base::is.character(method_ovl) || base::length(method_ovl) != 1){

      stop("Argument 'method_ovl' needs to be a single character value.")

    } else if(!method_ovl %in% c("classic", "bayesian")) {

      stop("Argument 'method_ovl' needs to be  one of: 'classic', 'bayesian'.")

    }

  }

}



#' @title Check monocle input
#'
#' @param preprocess_method Monocle3 - description:
#'
#' A string specifying the initial dimension method to use,
#' currently either PCA or LSI. For LSI (latent semantic
#' indexing), it converts the (sparse) expression matrix into
#' tf-idf matrix and then performs SVD to decompose the gene
#' expression / cells into certain modules / topics. Default
#' is "PCA".
#'
#' @param reduction_method Monocle3 - description:
#'
#' A character string specifying the algorithm to use for
#' dimensionality reduction. Currently "UMAP", "tSNE", "PCA"
#' and "LSI" are supported.
#'
#' @param cluster_method Monocle3 - description:
#'
#' String indicating the clustering method to use. Options are
#' "louvain" or "leiden". Default is "leiden". Resolution parameter
#' is ignored if set to "louvain".
#'
#' @param k Monocle3 - description:
#'
#' Integer number of nearest neighbors to use when creating the k
#' nearest neighbor graph for Louvain/Leiden clustering. k is related
#' to the resolution of the clustering result, a bigger k will result
#' in lower resolution and vice versa. Default is 20.
#'
#' @param num_iter Monocle3 - description:
#'
#' Integer number of iterations used for Louvain/Leiden clustering.
#' The clustering result giving the largest modularity score will be
#' used as the final clustering result. Default is 1. Note that if
#' num_iter is greater than 1, the random_seed argument will be ignored
#' for the louvain method.
#'
#' @return
#'
#' @details With respect to the arguments \code{preprocess_method},
#' \code{reduction_method} and \code{cluster_method}:
#'
#' If you provide a vector instead of a single character value (string)
#' the function will iterate over all inputs via a for-loop to compute all
#' valid combinations.
#'
#' @export
#'

check_monocle_input <- function(preprocess_method,
                                reduction_method,
                                cluster_method,
                                k = base::numeric(1),
                                num_iter = base::numeric(1)){

  confuns::is_value(k, "numeric", "k")
  confuns::is_value(num_iter, "numeric", "num_iter")

  if(!base::is.null(preprocess_method) &&
     base::any(!preprocess_method %in% c("PCA", "LSI"))){

    base::stop("Invalid input for argument 'preprocess_method'. Valid inputs are: 'PCA', 'LSI'")

  }

  if(!base::is.null(reduction_method) &&
     base::any(!reduction_method %in% c("UMAP", "tSNE", "PCA", "LSI"))){

    base::stop("Invalid input for argument 'reduction_method'. Valid inputs are: 'UMAP', 'tSNE', 'PCA', 'LSI'")

  }

  if(!base::is.null(cluster_method) &&
     base::any(!cluster_method %in% c("louvain", "leiden"))){

    base::stop("Invalid input for argument 'cluster_method'. Valid inputs are: 'louvain', 'leiden'")

  }

}


#' Check spata object input
#'
#' @param object A valid spata-object.
#'
#' @inherit lazy_check_dummy description details return
#' @export

check_object <- function(object){

  validation(object)

  base::return(TRUE)

}



#' @title Check pt input
#'
#' @param pt_size The size of the points specified as a single numeric value.
#' @param pt_alpha The transparency of the points specified as single numeric value.
#' @param pt_clrsp The color spectrum to be used if the specified variable displayed by
#' color is continuous. Run \code{all_colorspectra()} to see valid input..
#' @param pt_clrsp_dir The direction of the color spectrum specified as either \emph{1}
#' or \emph{-1}.
#' @param pt_clrp The color panel to be used if the specified variable displayed by
#' color is categorical/discrete. Run \code{all_colorpanels()} to see valid input.
#' @param pt_clr The base color of every point displayed in the plot.
#'
#' @inherit lazy_check_dummy description details return
#' @export

check_pt <- function(pt_size = NULL,
                     pt_alpha = NULL,
                     pt_clrsp = NULL,
                     pt_clrsp_dir = NULL,
                     pt_clrp = NULL,
                     pt_clr = NULL){


  if(!base::is.null(pt_size) && !base::is.numeric(pt_size)){

    base::stop("Argument 'pt_size' needs to be a single numeric value.")

  }

  if(!base::is.null(pt_alpha) && !base::is.numeric(pt_alpha)){

    base::stop("Argument 'pt_alpha' needs to be a single numeric value.")

  }

  if(!base::is.null(pt_clrsp) && !pt_clrsp %in% base::unlist(confuns::all_colorspectra(), use.names = FALSE)){

    base::stop("Invalid input for argument 'pt_clrsp'.")

  }

  if(!base::is.null(pt_clrp) && !pt_clrp %in% base::unlist(confuns::all_colorpanels(), use.names = FALSE)){

    base::stop("Invalid input for argument 'pt_clrp'.")

  }

  if(!base::is.null(pt_clrsp_dir) && !pt_clrsp_dir %in% c(1, -1)){

    base::stop("Argument 'pt_clrsp_dir' needs to be either 1 or -1")

  }

  base::return(TRUE)

}



#' @title Check ranked trajectory data.frame input
#'
#' @param atdf An assessed trajectory data.frame.
#' @param rtdf A ranked trajectory data.frame.
#' @param variable The gene or gene set of interest specified as a character
#' value.
#'
#' @inherit lazy_check_dummy description details return
#' @export
#'

check_atdf <- function(atdf){

  confuns::check_data_frame(
    df = atdf,
    var.class = list(
      variables = "character",
      pattern = "character",
      auc = c("double", "integer", "numeric")
    )
  )

}

#' @rdname check_atdf
#' @export
check_rtdf <- function(rtdf, variable = NULL){

  # check classes
  confuns::check_data_frame(df = rtdf,
                             var.class =
                               list(
                                 variables = "character",
                                 data = "list",
                                 residuals = "list",
                                 auc = "list"),
                             ref = "rtdf")

  base::return(base::invisible(TRUE))

}



#' @title Check smooth input
#'
#' @param df A data.frame that is to be smoothed spatially. That data frame must have
#' numeric \emph{x}- and \emph{y}-variables.
#' @param smooth Logical. If set to TRUE values will be smoothed according to the
#' \code{smoooth_}-parameters.
#' @param smooth_span The amount of smoothing specified as a single numeric value.
#' @param smooth_method The smoothing method that will be used specified as a
#' single character value (e.g. \emph{"lm", "glm", "gam", "loess"}).
#' @param smooth_se Logical. If set to TRUE the confidence interval will be
#' displayed.
#'
#' @inherit lazy_check_dummy description details return
#' @export
#'

check_smooth <- function(df = NULL,
                         smooth = NULL,
                         smooth_span = NULL,
                         smooth_method = NULL,
                         smooth_se = NULL){

  if(!base::is.null(smooth) &&
     !base::isTRUE(smooth) &
     !base::isFALSE(smooth)){

    base::stop("Argument 'smooth' needs to be TRUE or FALSE.")

  }

  if(!base::is.null(smooth) && base::isTRUE(smooth)){

    if(!base::is.null(df) &&
       !base::all(c("x", "y") %in% base::colnames(df))){

      base::stop("Input data.frame doesn't contain x and y variables." )

    }

  }

  if(!base::is.null(smooth_span) && !base::is.numeric(smooth_span)){

    base::stop("Argument 'smooth_span' needs to be numeric.")

  }

  if(!base::is.null(smooth_method)){

    if(!base::is.character(smooth_method) |
       !base::length(smooth_method) == 1){

      base::stop("Argument 'smooth_method' needs to be a single character value.")

    }

  }

  if(!base::is.null(smooth_se) &&
     !base::isTRUE(smooth_se) &
     !base::isFALSE(smooth_se)){

    base::stop("Argument 'smooth_se' needs to be TRUE or FALSE.")

  }

  base::return(TRUE)

}



#' @title Check stdf-input
#'
#' @param stdf A summarized trajectory data.frame
#'
#' @inherit lazy_check_dummy description details return
#' @export
#'

check_summarized_trajectory_df <- function(stdf){

  confuns::check_data_frame(
    df = stdf, var.class = list(
      trajectory_part = "character",
      trajectory_part_order = c("numeric", "integer"),
      trajectory_order = c("numeric", "integer"),
      values = c("numeric", "integer"),
      variables = "character"
    ),
    ref = "stdf")

}



#' Check trajectory name input
#'
#' @inherit check_sample params
#' @param trajectory_name The trajectory of interest specified
#' as a single character value.
#'
#'
#' @inherit lazy_check_dummy description details return
#' @export

check_trajectory <- function(object,
                             trajectory_name,
                             of_sample){

  if(!base::length(trajectory_name) == 1 || !base::is.character(trajectory_name)){

    base::stop("Argument 'trajectory_name' needs to be a single character value.")

  } else if(!trajectory_name %in% getTrajectoryNames(object, of_sample = of_sample)){

    base::stop(stringr::str_c("There is no trajectory of name '", trajectory_name,
                              "' in sample '", of_sample, "'.", sep = ""))

  }


}


#' @title Check uniform genes input
#'
#' @param uniform_genes Character value. If set to \emph{'discard'} genes that are
#' uniformly expressed across all barcode-spots of the specified coordinates
#' data.frame are discarded. If set to \emph{'keep'} they are kept.
#'
#' @return
#' @export
#'

check_uniform_genes <- function(uniform_genes){

  confuns::is_value(uniform_genes, "character", "uniform_genes")

  if(!uniform_genes %in% c("keep", "discard")){

    base::stop("Argument 'uniform genes' must be set to 'keep' or 'discard'.")

  } else {

    base::return(base::invisible(TRUE))
  }

}








