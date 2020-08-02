# Expression matrix related

#' Obtain expression matrix
#'
#' @inherit check_sample params
#'
#' @return The expression matrix of the specified object and sample.
#' @export

getExpressionMatrix <- function(object,
                                of_sample = "all"){

  # lazy check
  check_object(object)

  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample)

  rna_assay <- exprMtr(object = object, of_sample = of_sample)

  base::return(rna_assay)

}


#' Obtain barcodes of a sample
#'
#' @inherit check_sample params
#'
#' @return All barcodes of the specified sample(s).
#' @export

getBarcodes <- function(object, of_sample = "all"){

  cdf <- coordinates(object = object, of_sample = of_sample)

  return(dplyr::pull(cdf, barcodes))

}



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
  gene_sets_df <- object@used_genesets
  gene_sets <- object@used_genesets$ont

  gene_set_classes <- stringr::str_extract(string = gene_sets, pattern = "^.+?(?=_)")

  dplyr::mutate(gene_sets_df, gs_type = gene_set_classes) %>%
    dplyr::select(-gene) %>%
    dplyr::distinct() %>%
    dplyr::pull(gs_type) %>%
    base::table() %>%
    base::as.data.frame() %>%
    magrittr::set_colnames(value = c("Class", "Available Gene Sets"))


}



#' @title Obtain gene set names
#'
#' @param object A valid spata-object.
#' @param of_class A character vector indicating the classes from which to obtain
#' the gene set names. (Which classes exist in the current gene set data.frame can
#' be obtained e.g. with \code{geneSetOverview()}). If set to \emph{"all"} all
#' gene sets are returned.
#' @param index A regular expression according to which the gene set names to be returned
#' will be filtered again.
#' @param simplify Logical. If set to TRUE the list to be returned will be simplified
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

  base::return(res_list)

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
                     in_sample = "all",
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
        base::unlist()

    }

    base::return(genes_list)

  }

  # -----

}



# Feature related ---------------------------------------------------------

#' @title Obtain feature names
#'
#' @param object A valid spata-object.
#'
#' @return A named character vector of the variables in the feature data slot.
#' @export

getFeatureNames <- function(object){

  check_object(object)

  feature_names <- base::colnames(object@fdata)

  base::names(feature_names) <-
    base::sapply(object@fdata[,feature_names], base::class)

  base::return(feature_names[!feature_names %in% c("sample", "barcodes")])

}



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
                            of_sample = "all"){

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

  base::return(res_list)


}



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

getTrajectoryNames <- function(object,
                               of_sample = "all"){

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

  if(!base::length(t_names_list) == 0){

    base::return(t_names_list)

  } else {

    base::return(base::invisible(NULL))

  }


}



