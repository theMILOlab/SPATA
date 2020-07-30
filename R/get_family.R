# Expression matrix related

#' Obtain expression matrix
#'
#' @param object A valid spata-object.
#' @param of_sample The sample(s) from which to extract the expression data
#' specified as a character vector.
#'
#' @return The expression matrix of the specified object.
#' @export
#'

getExpressionMatrix <- function(object,
                             of_sample = "all"){

  validation(object)
  of_sample <- check_sample(object = object, sample_input = of_sample)

  return(exprMtr(object = object, of_sample = of_sample))


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
#'

getGeneSetOverview <- function(object){

  validation(x = object)

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
#' The list coalesced to an unnamed vector if \code{simplify} is set to TRUE.
#'
#' @export
#'

getGeneSets <- function(object, of_class = "all", index = NULL, simplify = TRUE){

  validation(x = object)
  stopifnot(base::is.character(index) | base::is.null(index))

  gene_sets_df <- object@used_genesets

  if(!base::is.character(of_class)){

    stop("Please specify 'of_class' as a character vector.")

  }

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


  # Adjust output -----------------------------------------------------------

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

  base::return(res_list)

}



#' @title Obtain gene names
#'
#' @param object A valid spata-object.
#' @param of_gene_sets A character vector specifying the gene sets from which to
#' return the gene names.
#' @param in_sample The sample(s) in which the genes have to be expressed in order
#' to be included.
#' @param rna_assay Old argument used to provide an Rna-assay.
#' @param simplify Logical. If set to TRUE the list to be returned will be simplified
#' into a character vector.
#'
#' @return A list named according to the input of \code{of_gene_sets} in which each element is
#' a character vector containing the names of genes the specific gene set is
#' composed of. Is simplified to a vector if the number of elements in the returned list
#' is one or \code{simplify} is set to TRUE.
#'
#' @export

getGenes <- function(object,
                     of_gene_sets = "all",
                     in_sample = "all",
                     rna_assay = NULL,
                     simplify = FALSE){

  validation(x = object)

  if(!is.null(rna_assay) && is.matrix(rna_assay)){

    warning("argument rna_assay is deprecated. do not use it while programming with getGenes()!")

  }

  # control:
  if(!is.character(of_gene_sets) | length(of_gene_sets) == 0){

    stop("Argument 'of_gene_sets' is empty or invalid. Has to be a character vector of length one or more.")

  }

  if(base::all(c(of_gene_sets, in_sample) == "all")){

    base::return(base::rownames(exprMtr(object)))

  } else {

    in_sample <- check_sample(object = object, sample_input = in_sample)
    rna_assay <- exprMtr(object = object, of_sample = in_sample)

    gene_sets_df <- object@used_genesets

    if(base::any(!c("ont", "gene") %in% base::colnames(gene_sets_df)) |
       !is.data.frame(gene_sets_df) |
       base::nrow(gene_sets_df) == 0){

      stop("Please make sure that the provided object contains a valid gene sets data.frame.")

    }


    # if only genes of specific gene sets are desired
    if(base::length(of_gene_sets) != 1){

      gene_sets_ctrl <- base::vector(mode = "logical", length = length(of_gene_sets))
      not_found <- base::vector(mode = "character")

      for(i in seq_along(of_gene_sets)){

        if(of_gene_sets[i] %in%  gene_sets_df$ont){

          gene_sets_ctrl[i] <- T

        } else{

          gene_sets_ctrl[i] <- F
          not_found[length(not_found)+1] <- of_gene_sets[i]

        }

      }


      of_gene_sets <- of_gene_sets[gene_sets_ctrl]

      if(length(of_gene_sets) >= 1 & length(not_found) >= 1){

        not_found_clpsd <- stringr::str_c(not_found, collapse = ", ")

        base::message(stringr::str_c("Could not find gene sets: ", not_found_clpsd, ".", sep = ""))

      } else if(length(of_gene_sets) == 0){

        stop("Could not find any of the provided gene sets.")

      }

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


      if(base::isTRUE(simplify) && base::length(genes_list) > 1){

        genes_list <-
          genes_list %>%
          base::unname() %>%
          base::unlist()

      }


      base::return(genes_list)

    } else if(base::all(of_gene_sets == "all")) { # if all genes are desired

      all_genes <-
        dplyr::pull(gene_sets_df, "gene") %>%
        base::unique()

      all_genes_in_sample <-
        all_genes[all_genes %in% base::rownames(rna_assay)]

      return(all_genes_in_sample)

    } else if(base::length(of_gene_sets) == 1){

      if(of_gene_sets %in% gene_sets_df$ont){

        genes <-
          gene_sets_df %>%
          dplyr::filter(ont == of_gene_sets) %>%
          dplyr::pull(gene)

        return(genes)

      } else {

        stop(stringr::str_c("Did not find gene set", of_gene_sets, sep = ": "))

      }


    }

  }

}



# Feature related ---------------------------------------------------------

#' @title Obtain feature names
#'
#' @param object A valid spata-object.
#'
#' @return A named character vector of the variables in the feature data slot.
#' @export
#'

getFeatureNames <- function(object){

  validation(x = object)

  feature_names <- base::colnames(object@fdata)

  base::names(feature_names) <-
    base::sapply(object@fdata[,feature_names], base::class)

  base::return(feature_names[!feature_names %in% c("sample", "barcodes")])

}



# Segmentation related ----------------------------------------------------

#' @title Obtain segment names
#'
#' @param object A valid spata-object.
#' @param of_sample The samples from which to obtain the segment names specified
#' as a character vector.
#'
#' @return A list named according to the \code{of_sample} in which each element is
#' a character vector containing the names of segments which were drawn for the
#' specific sample. Is simplified to a vector if the number of samples in the returned list
#' is one.
#'
#' @export
#'

getSegmentNames <- function(object,
                            of_sample){

  validation(x = object)

  of_sample <- check_sample(object, sample_input = of_sample)

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


  if(base::length(res_list) == 1){

    base::return(base::as.character(res_list[[1]]))

  } else {

    base::return(res_list)
  }


}



# Trajectory related ------------------------------------------------------

#- 'getTrajectoryComment()' is documented in 'S4_generic_functions.R' -#

#' @title Obtain trajectory names
#'
#' @param object A valid spata object.
#' @param of_sample The samples from which to obtain the trajectory names specified
#' as a character vector.
#'
#' @return A list named according to the \code{of_sample} in which each element is
#' a character vector containing the names of trajectories which were drawn for the
#' specific sample. Is simplified to a vector if the number of samples in the returned list
#' is one.
#'
#' @export
#'
getTrajectoryNames <- function(object,
                               of_sample = "all"){

  validation(x = object)

  of_sample <- check_sample(object = object, sample_input = of_sample)

  t_names_list <-
    base::lapply(X = of_sample, FUN = function(i){

      t_names <-
        base::names(object@trajectories[[i]])

      if(base::length(t_names) == 0){

        message(stringr::str_c("No trajectories found in sample: ", i, sep = ""))

        return(NULL)

      } else {

        return(t_names)

      }

    })

  base::names(t_names_list) <- of_sample

  t_names_list <- purrr::discard(.x = t_names_list, .p = is.null)

  if(length(t_names_list) == 1){

    return(base::as.character(t_names_list[[1]]))

  } else {

    return(t_names_list)

  }

}


