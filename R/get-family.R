#' Obtain barcodes of a sample
#'
#' @inherit check_sample params
#'
#' @return All barcodes of the specified sample(s) as a character vector.
#' @export

getBarcodes <- function(object, of_sample = "all"){

  cdf <- coordinates(object = object, of_sample = of_sample)

  return(dplyr::pull(cdf, barcodes))

}



#' @title Obtain spatial coordinates
#'
#' @inherit check_sample params
#' @param of_segment Character value. Specifies the segment of interest.
#'
#' @return A data.frame containing the variables \emph{barcods, sample, x, y}
#' (and \emph{segment} if specified).
#' @export
#'

getCoordinates <- function(object,
                           of_sample = ""){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)

  # adjusting check
  of_sample <- check_sample(object, of_sample, 1)

  # -----

  # 2. Data wrangling -------------------------------------------------------

  coords <- coordinates(object, of_sample = of_sample)

  # -----

  base::return(coords)


}

#' @rdname getCoordinates
#' @export
getCoordinatesSegment <- function(object,
                                  of_segment,
                                  of_sample = ""){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)

  # adjusting check
  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)
  bc_segm <- check_segment(object, segment_name = of_segment, of_sample = of_sample)

  # -----

  # 2. Data wrangling -------------------------------------------------------

  coords <-
    coordinates(object = object, of_sample = of_sample) %>%
    dplyr::filter(barcodes %in% bc_segm) %>%
    dplyr::mutate(segment_name = {{of_segment}})

  # -----

  base::return(coords)

}


#' Obtain count and expression matrix
#'
#' @inherit check_sample params
#'
#' @return The expression/count matrix of the specified object and sample(s).
#' @export

getExpressionMatrix <- function(object,
                                of_sample = ""){

  # lazy check
  check_object(object)

  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample)

  rna_assay <- exprMtr(object = object, of_sample = of_sample)

  base::return(rna_assay)

}

#' @rdname getExpressionMatrix
#' @export
getCountMatrix <- function(object,
                           of_sample = ""){

  # lazy check
  check_object(object)

  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample)

  rna_assay <- countMtr(object = object, of_sample = of_sample)

  base::return(rna_assay)

}



# Dimensional reduction ---------------------------------------------------

#' @title Obtain dimensional reduction data
#'
#' @inherit check_sample params
#' @inherit check_method params
#'
#' @return A data.frame that contains the unique identifiers
#' (keys): \emph{barcodes, sample} and:.
#'
#'  \itemize{
#'   \item{ \code{getTsneData()}: \emph{tsne1, tsne2}}
#'   \item{ \code{getUmapData()}: \emph{umap1, umap2}}
#'   }
#'
#' @export
#'

getDimRedData <- function(object,
                          of_sample = "",
                          method_dr = c("UMAP", "TSNE")){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_method(method_dr = method_dr)

  # adjusting check
  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  # -----

  # 2. Data extraction ------------------------------------------------------

  dr_strings <- stringr::str_c(base::tolower(x = method_dr), 1:2, sep = "")

  dim_red_df <-
    methods::slot(object = object@dim_red, name = method_dr) %>%
    dplyr::filter(sample %in% of_sample) %>%
    dplyr::select(dplyr::all_of(x = c("barcodes", "sample", dr_strings))) %>%
    tibble::remove_rownames()

  # -----

  if(base::nrow(dim_red_df) == 0){

    base::stop("There seems to be no data for method: ", method_dr)

  }

  base::return(dim_red_df)

}

#' @rdname getDimRedData
#' @export
getUmapData <- function(object,
                        of_sample = ""){


  getDimRedData(object = object,
                of_sample = of_sample,
                method_dr = "UMAP")

}

#' @rdname getDimRedData
#' @export
getTsneData <- function(object,
                        of_sample = ""){

  getDimRedData(object = object,
                of_sample = of_sample,
                method_dr = "TSNE")

}
# -----


# Genes and gene set related ----------------------------------------------

#' @title Overview about the current gene sets
#'
#' @param object A valid spata-object.
#'
#' @return A data.frame with two variables \emph{Class} and \emph{Available Gene
#' Sets} indicating the number of different gene sets the classes contain.
#'
#' @export

getGeneSetOverview <- function(object){

  # lazy check
  check_object(object)

  # main part
  gene_sets_df <- dplyr::ungroup(object@used_genesets)

  gene_sets <- object@used_genesets$ont

  if(base::nrow(gene_sets_df) == 0){

    base::message("Gene-set data.frame is empty.")
    base::return(data.frame())

  } else {

    gene_set_classes <- stringr::str_extract(string = gene_sets, pattern = "^.+?(?=_)")

    dplyr::mutate(gene_sets_df, gs_type = gene_set_classes) %>%
      dplyr::select(-gene) %>%
      dplyr::distinct() %>%
      dplyr::pull(gs_type) %>%
      base::table() %>%
      base::as.data.frame() %>%
      magrittr::set_colnames(value = c("Class", "Available Gene Sets"))

  }


}



#' @title Obtain gene set names
#'
#' @param object A valid spata-object.
#' @param of_class A character vector indicating the classes from which to obtain
#' the gene set names. (Which classes exist in the current gene set data.frame can
#' be obtained e.g. with \code{geneSetOverview()}). If set to \emph{"all"} all
#' gene sets are returned.
#' @param index A regular expression according to which the gene set names to be returned
#' are filtered again.
#' @param simplify Logical. If set to TRUE the list to be returned is simplified
#' into a character vector.
#'
#'
#' @return A list named according to the input of argument \code{of_class}. Each element of
#' the returned list is a character vector containing the names of gene sets of the specified classes.
#' The list is coalesced to an unnamed vector if \code{simplify} is set to TRUE.
#'
#' @export

getGeneSets <- function(object, of_class = "all", index = NULL, simplify = TRUE){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  stopifnot(base::is.character(index) | base::is.null(index))

  if(!base::is.character(of_class)){

    stop("Please specify 'of_class' as a character vector.")

  }

  # -----

  # 2. Main part ------------------------------------------------------------

  gene_sets_df <- object@used_genesets

  # 2.1 Extract gene sets according to 'of_class' ----------
  if(base::length(of_class) == 1 && of_class == "all"){

    res_list <- base::unique(gene_sets_df$ont)

  } else {

    # get gene sets for all elements of 'of_class' in a list
    res_list <-
      base::lapply(X = of_class, FUN = function(i){

        subset <-
          gene_sets_df$ont %>%
          stringr::str_subset(pattern = stringr::str_c("^", i, sep = "")) %>%
          base::unique()

        if(base::length(subset) == 0){

          base::warning(stringr::str_c("Could not find any gene set of class:", i, sep = " "))

          base::return(NULL)

        } else {

          base::return(subset)

        }

      })

    base::names(res_list) <- of_class

    # discard list elements if 'of_class' element wasn't found
    res_list <-
      purrr::discard(.x = res_list, .p = base::is.null)

  }

  # -----


  # 2.2 Adjust output according to 'index' ----------

  if(base::isTRUE(simplify)){

    res_list <- base::unlist(res_list) %>% base::unname()

  }


  if(!base::is.null(index) && base::is.list(res_list)){

    res_list <-
      base::lapply(X = res_list,
                   FUN = function(i){

                     i[stringr::str_detect(string = i, pattern = index)]

                   })

  } else if(!base::is.null(index) && base::is.character(res_list)){

    res_list <-
      res_list[stringr::str_detect(string = res_list, pattern = index)]

  }

  # -----
  if(base::is.null(res_list)){

    base::stop("Did not find any gene-set.")

  } else {

    base::return(res_list)

  }

}



#' @title Obtain gene names
#'
#' @param object A valid spata-object.
#' @param of_gene_sets A character vector specifying the gene sets from which to
#' return the gene names.
#' @param in_sample The sample(s) in which the genes have to be expressed in order
#' to be included.
#' @param simplify Logical. If set to TRUE the list to be returned will be simplified
#' into a character vector.
#'
#' @return A list named according to the input of \code{of_gene_sets} in which each element is
#' a character vector containing the names of genes the specific gene set is
#' composed of. Is coalesced to a vector if \code{simplify} is set to TRUE.
#'
#' @export

getGenes <- function(object,
                     of_gene_sets = "all",
                     in_sample = "",
                     simplify = TRUE){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)

  if(!is.character(of_gene_sets) | base::length(of_gene_sets) == 0){

    stop("Argument 'of_gene_sets' is empty or invalid. Has to be a character vector of length one or more.")

  }

  # adjusting check
  in_sample <- check_sample(object = object, of_sample = in_sample)

  # -----


  # 2. Main part ------------------------------------------------------------

  rna_assay <- exprMtr(object = object, of_sample = in_sample)

  # -----

  # 2.2 Return all existing genes if desired ----------

  if(base::all(of_gene_sets == "all")){

    base::return(base::rownames(rna_assay))

  }

  # -----

  # 2.3 Return a subset of genes ----------
  if(!base::all(of_gene_sets == "all")){

    gene_sets_df <- object@used_genesets
    of_gene_sets <- check_gene_sets(object, of_gene_sets)

    genes_list <-
      base::lapply(X = of_gene_sets,
                   FUN = function(i){

                     genes <-
                       dplyr::filter(gene_sets_df, ont == i) %>%
                       dplyr::pull(gene)

                     genes_in_sample <-
                       genes[genes %in% base::rownames(rna_assay)]

                       return(genes_in_sample)

                     })

    base::names(genes_list) <- of_gene_sets

    if(base::isTRUE(simplify)){

      genes_list <-
        genes_list %>%
        base::unname() %>%
        base::unlist() %>%
        base::unique()

    }

    base::return(genes_list)

  }

  # -----

}



# -----


# Feature related ---------------------------------------------------------

#' @title Obtain feature names
#'
#' @param object A valid spata-object.
#' @param of_class Character vector. Specify the classes a feature must be of for
#' it's name to be returned.
#'
#' @return A named character vector of the variables in the feature data slot.
#' @export

getFeatureNames <- function(object, of_class = NULL){

  check_object(object)
  if(!base::is.null(of_class)){confuns::is_vec(of_class, "character", "of_class")}

  feature_names <- base::colnames(object@fdata)

  classes <- base::sapply(object@fdata[,feature_names], base::class)

  base::names(feature_names) <- classes

  if(!base::is.null(of_class)){
    feature_names <- feature_names[classes %in% of_class]
  }

  base::return(feature_names[!feature_names %in% c("sample", "barcodes")])

}


#' Obtain feature data
#'
#' @inherit check_sample params
#'
#' @return The feature data data.frame of the specfied object and sample(s).
#' @export
#'

getFeatureData <- function(object, of_sample = ""){

  check_object(object)
  of_sample <- check_sample(object, of_sample)

  featureData(object = object,
              of_sample = of_sample)

}


#' @title Obtain a feature variable
#'
#' @description Extracts the specified feature variables from the
#' feature data.
#'
#' @inherit check_sample params
#' @inherit check_features params
#' @param return Character value. Determine the way the variable is returned.
#'
#' (Only relevant if input of arugment \code{features} is of length greater than 1.)
#'
#' @param unique Logical. If set to TRUE and argument \code{features} is
#' of length 1 only it's unique values are returned.
#'
#' @return A data.frame or a vector.
#' @export
#'

getFeatureVariables <- function(object,
                                features,
                                of_sample = "",
                                return = "data.frame",
                                unique = FALSE){


  # 1. Control --------------------------------------------------------------

  check_object(object)
  features <- check_features(object, features)

  confuns::is_value(return, "character", "return")
  stopifnot(return %in% c("data.frame", "vector"))

  of_sample <- check_sample(object, of_sample)

  # -----

  # 2. Extracting -----------------------------------------------------------


  if(base::length(features) == 1 && return == "vector"){

    var <-
      getFeatureData(object, of_sample) %>%
      dplyr::pull(var = {{features}})

    if(base::isTRUE(unique)){

      base::unique(var) %>%
        base::return()

    } else {

      base::return(var)

    }

  } else if(base::length(features) == 1){

    getFeatureData(object, of_sample) %>%
      dplyr::select(barcodes, sample, {{features}})

  } else {

    getFeatureData(object, of_sample) %>%
      dplyr::select(barcodes, sample, dplyr::all_of(features)) %>%
      base::return()

  }

}


# -----


# Segmentation related ----------------------------------------------------

#' @title Obtain segment names
#'
#' @inherit check_sample params
#'
#' @return A list named according to the \code{of_sample} in which each element is
#' a character vector containing the names of segments which were drawn for the
#' specific sample.
#'
#' @export

getSegmentNames <- function(object,
                            of_sample = "",
                            simplify = TRUE){

  # lazy check
  check_object(object)

  # adjusting check
  of_sample <- check_sample(object, of_sample = of_sample)

  # main part
  res_list <-
    base::lapply(X = of_sample,
                FUN = function(i){

                  segment_names <-
                    featureData(object) %>%
                    dplyr::filter(sample == i) %>%
                    dplyr::pull(segment) %>% base::unique()

                  if(base::length(segment_names) == 1 && segment_names == ""){

                     base::warning(stringr::str_c("There seems to be no segmentation for '", i, "'."))

                     base::return(NULL)

                    }

                  return(segment_names[segment_names != ""])

                })

  base::names(res_list) <- of_sample

  res_list <- purrr::discard(.x = res_list, .p = base::is.null)


  if(base::isTRUE(simplify)){

    base::unlist(res_list, use.names = FALSE) %>%
      base::return()

  } else {

    base::return(res_list)

  }


}



# -----


# Trajectory related ------------------------------------------------------

#- 'getTrajectoryComment()' is documented in 'S4_generic_functions.R' -#

#' @title Obtain trajectory names
#'
#' @inherit check_sample params
#'
#' @return A list named according to the \code{of_sample} in which each element is
#' a character vector containing the names of trajectories which were drawn for the
#' specific sample.
#'
#' @export

getTrajectoryNames <- function(object, of_sample = "all", simplify = TRUE){

  # lazy check
  check_object(object)

  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample)

  # main part
  t_names_list <-
    base::lapply(X = of_sample, FUN = function(i){

      t_names <-
        base::names(object@trajectories[[i]])

      if(base::length(t_names) == 0){

        base::message(stringr::str_c("No trajectories found in sample: ", i, sep = ""))

        base::return(NULL)

      } else {

        base::return(t_names)

      }

    })

  base::names(t_names_list) <- of_sample

  t_names_list <- purrr::discard(.x = t_names_list, .p = is.null)

  if(base::isTRUE(simplify)){

    t_names_list <- base::unlist(t_names_list) %>% base::unname()

  }

  if(!base::length(t_names_list) == 0){

    base::return(t_names_list)

  } else {

    base::return(base::invisible(NULL))

  }


}



#' @title Obtain a summarized trajectory data.frame
#'
#' @description Computes the expression trends of all specified variables
#' along the direction of the spatial trajectory.
#'
#' @inherit check_sample params
#' @inherit check_trajectory params
#' @inherit hlpr_summarize_trajectory_df params
#'
#' @return A summarized trajectory data.frame.
#'
#' @inherit hlpr_summarize_trajectory_df details
#'
#' @export

getTrajectoryDf <- function(object,
                            trajectory_name,
                            of_sample = "",
                            variables,
                            method_gs = "mean",
                            accuracy = 5,
                            normalize = TRUE,
                            verbose = TRUE){


  tobj <-
    getTrajectoryObject(object, trajectory_name, of_sample)

  stdf <-
    hlpr_summarize_trajectory_df(object,
                                 ctdf = tobj@compiled_trajectory_df,
                                 accuracy = accuracy,
                                 variables = variables,
                                 method_gs = method_gs,
                                 verbose = verbose,
                                 normalize = normalize)

  base::return(stdf)

}

#' @title Obtain trajectory object
#'
#' @inherit check_sample params
#' @inherit check_trajectory params
#'
#' @return An object of class \code{spatialTrajectory}.
#' @export

getTrajectoryObject <- function(object, trajectory_name, of_sample = ""){

  of_sample <- check_sample(object = object,
                            of_sample = of_sample,
                            desired_length = 1)

  trajectory(object = object,
             trajectory_name = trajectory_name,
             of_sample = of_sample)

}




# -----
