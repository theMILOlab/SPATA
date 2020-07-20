
# Gene set related --------------------------------------------------------

#' @title Add a new gene set
#'
#' @description Stores a new gene set in the spata object.
#'
#' @param object A valid spata-object.
#' @param gs_class The class the gene set belongs to specified as a character.
#' @param gs_name The name of the new gene set specified as a character.
#' @param gs_genes The genes specified as a character vector.
#' @param overwrite Logical. Overwrites existing gene sets with the same \code{gs_class} -
#' \code{gs_name} combination.
#'
#' @return An updated spata-object.
#'
#' @details Coerces \code{gs_class} and \code{gs_name} to the final gene set name.
#'
#' @export
#'

addGeneSet <- function(object,
                       gs_class,
                       gs_name,
                       gs_genes,
                       overwrite = FALSE){

  # control
  validation(x = object)

  gs_genes <- check_genes(object, genes = gs_genes)

  if(base::any(!base::sapply(X = list(gs_class, gs_name, gs_genes),
                             FUN = base::is.character))){

    base::stop("Arguments 'gs_*' must be of class character.")

  }

  if(base::length(gs_class) != 1 | base::length(gs_name) != 1){

    base::stop("Arguments 'gs_class' and 'gs_name' must be of length one.")

  }


  name <- stringr::str_c(gs_class, gs_name, sep = "_")

  # make sure not to overwrite if overwrite == FALSE
  if(name %in% object@used_genesets$ont && base::isFALSE(overwrite)){

    base::stop(stringr::str_c("Gene set '", name, "' already exists.",
                              " Set argument 'overwrite' to TRUE in order to overwrite existing gene set."))

  } else if(name %in% object@used_genesets$ont && base::isTRUE(overwrite)) {

    object <- discardGeneSets(object, gs_names = name)

  }

  # add gene set
  object@used_genesets <-
    dplyr::add_row(
      .data = object@used_genesets,
      ont = base::rep(name, base::length(gs_genes)),
      gene = gs_genes
    )

  base::return(object)

}


#' Discard gene sets
#'
#' @param object A valid spata-object.
#' @param gs_names The gene sets to be discarded specified as a character vector.
#'
#' @return An updated spata-object.
#' @export
#'

discardGeneSets <- function(object, gs_names){

  validation(object)

  base::stopifnot(base::is.character(gs_names))

  if(base::any(gs_names %in% object@used_genesets$ont)){

    object@used_genesets <-
      dplyr::filter(object@used_genesets,
                    !ont %in% gs_names)

  } else {

    base::message(stringr::str_c("Did not find any of the specified gene sets in object."))

  }

  return(object)

}



