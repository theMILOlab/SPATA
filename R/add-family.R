
# Gene set related --------------------------------------------------------

#' @title Add a new gene set
#'
#' @description Stores a new gene set in the spata-object.
#'
#' @param object A valid spata-object.
#' @param class_name The class the gene set belongs to specified as a single character value.
#' @param gs_name The name of the new gene set specified as a single character value.
#' @param overwrite Logical. Overwrites existing gene sets with the same \code{class_name} -
#' \code{gs_name} combination.
#'
#' @inherit check_genes params
#'
#' @return An updated spata-object.
#'
#' @details Combines \code{class_name} and \code{gs_name} to the final gene set name.
#' Gene set classes and gene set names are separated by '_' and handled like this
#' in all additional gene set related functions which is why \code{class_name} must
#' not contain any '_'.
#'
#' @export

addGeneSet <- function(object,
                       class_name,
                       gs_name,
                       genes,
                       overwrite = FALSE){

  # lazy control
  check_object(object)

  # adjusting control
  genes <- check_genes(object, genes = genes)

  if(base::any(!base::sapply(X = list(class_name, gs_name, genes),
                             FUN = base::is.character))){

    base::stop("Arguments 'class_name', 'gs_name' and 'genes' must be of class character.")

  }

  if(base::length(class_name) != 1 | base::length(gs_name) != 1){

    base::stop("Arguments 'class_name' and 'gs_name' must be of length one.")

  }

  if(stringr::str_detect(string = class_name, pattern = "_")){

    base::stop("Invalid input for argument 'class_name'. Must not contain '_'.")

  }

  name <- stringr::str_c(class_name, gs_name, sep = "_")

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
      ont = base::rep(name, base::length(genes)),
      gene = genes
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

discardGeneSets <- function(object, gs_names){

  # lazy control
  check_object(object)

  # adjusting control
  gs_names <- check_gene_sets(object, gene_sets = gs_names)

  # discard gene sets
  object@used_genesets <-
      dplyr::filter(object@used_genesets,
                    !ont %in% gs_names)


  return(object)

}



# Feature related ---------------------------------------------------------


#' @title Add a new feature
#'
#' @description Adds a new variable to the objects feature data.
#'
#' @param object A valid spata-object.
#' @param overwrite Logical. If the specified feature name already exists in the
#' current spata-object this argument must be set to TRUE in order to overwrite it.
#'
#' @inherit check_feature_df params
#'
#' @details Eventually the new feature will be joined via \code{dplyr::left_join()} over the
#' key-variable \emph{barcodes}. Additional steps secure the joining process.
#'
#' @return An updated spata-object.
#' @export

addFeature <- function(object,
                       feature_name,
                       feature_df,
                       overwrite = FALSE){

  # lazy control
  check_object(object)
  check_feature_df(feature_name, feature_df)

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

  # make sure that feature_df$barcdoes si complete
  barcodes_obj <- fdata$barcodes
  barcodes_feature_df <- feature_df$barcodes

  if(!base::all(barcodes_obj %in% barcodes_feature_df)){

    not_found <- barcodes_obj[!barcodes_obj %in% barcodes_feature_df]
    n_not_found <- base::length(not_found)

    stop(stringr::str_c("Barcodes of input and of object have to match entirely. ",
                        n_not_found, " barcodes of pbject were not found in 'feature_df'."))

  }

  # join
  object@fdata <-
    dplyr::left_join(x = fdata,
                     y = feature_df[,c("barcodes", feature_name)],
                     by = "barcodes")

  base::return(object)

}





