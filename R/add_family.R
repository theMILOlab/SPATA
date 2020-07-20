
# Gene set related --------------------------------------------------------

#' @title Add a new gene set
#'
#' @description Stores a new gene set in the spata-object.
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
#' Gene set classes and gene set names are separated by '_' and handled like this
#' in all additional gene set related functions which is why \code{gs_class} must
#' not contain any '_'.
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

  if(stringr::str_detect(string = gs_class, pattern = "_")){

    base::stop("Invalid input for argument 'gs_class'. Must not contain '_'.")

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



# Feature related ---------------------------------------------------------


#' @title Add a new feature
#'
#' @description Adds a new variable to the objects feature data.
#'
#' @param object A valid spata-object.
#' @param feature_df A data.frame that contains the variables \emph{barcodes, samples, \code{feature_name}}.
#' @param feature_name The name of the new feature specified as a character value.
#' @param overwrite Logical. If the specified feature name already exists in the
#' current spata-object this argument must be set to TRUE in order to overwrite it.
#'
#' @details Eventually the variable will be joined via \code{dplyr::left_join()} over the
#' key-variable \emph{barcodes}. Additional steps secure the joining process.
#'
#' @return An updated spata-object.
#' @export
#'

addFeature <- function(object,
                       feature_df,
                       feature_name,
                       overwrite = FALSE){

  # control
  validation(x = object)

  stopifnot(base::is.character(feature_name))
  stopifnot(base::length(feature_name) == 1)

  feature_df <- check_feature_df(feature_df = feature_df)


  # extract data
  if(feature_name %in% getFeatureNames(object) &&
     !base::isTRUE(overwrite)){

    base::stop("Specified 'feature_name' is already present in current feature data. Set overwrite to TRUE to overwrite it.")

  } else if(feature_name %in% getFeatureNames(object) &&
            base::isTRUE(overwrite)){

    fdata <-
      object@fdata %>%
      dplyr::select(-dplyr::all_of(feature_name))

  } else {

    fdata <- featureData(object)

  }

  # make sure that feature_df contains the whole sample
  barcodes_obj <- fdata$barcodes
  barcodes_feature_df <- feature_df$barcodes

  if(!base::all(barcodes_obj %in% barcodes_feature_df)){

    not_found <- barcodes_obj[!barcodes_obj %in% barcodes_feature_df]
    n_not_found <- base::length(not_found)

    stop(stringr::str_c("Barcodes of input and of object have to match entirely. ",
                        n_not_found, " barcodes of 'feature_df' were not found in object."))

  }

  # join
  object@fdata <-
    dplyr::left_join(x = fdata, y = feature_df[,c("barcodes", feature_name)], by = "barcodes")

  base::return(object)

}





