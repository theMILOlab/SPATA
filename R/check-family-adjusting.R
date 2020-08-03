
#' @title This is a text dummy
#'
#' @details Members of the \code{adjusting-check_*()}-family take their
#' arguments input and compare it to a variety of requirement settings by
#' running several logical tests. If the input turns out to be appropriate
#' for the main-function they return it the way it is supposed to be returned.
#' If not, depending on the degree of deviation from the optimum, they either adjust
#' the input in order not to interrupt the function or - if not adjustable - raise an
#' error. In both cases informative messages will be printed in order to let the user
#' know what has been adjusted or what part of the input was insufficient.
#'
#' @return The original input, an adjusted version of it or an error. In the latter two
#' cases they print an informative message about what was going on.

adjusting_check_dummy <- function(){}



#################################################################################################

#' @title Check color to
#'
#' @description A member of the \code{adjusting-check_*()}-family. Takes a character
#' vector and sorts its elements into a list depending on whether it was found in
#' the input of \code{all_features}, \code{all_genes} or \code{all_gene_sets}.
#'
#' Returns a list with three slots named \emph{features}, \emph{genes} and \emph{gene_sets}
#' containing the respective found/valid input of \code{color_to}.
#'
#' @param color_to The variable to be displayed by color specified as .
#'
#'  \itemize{
#'   \item{ \strong{Gene set} as a single character value. Must be in \code{getGeneSets()}}
#'   \item{ \strong{Genes} as a character vector. If more than one gene is specified the average
#'   expression of those genes will be calculated and displayed. Must be in \code{getGenes()}}
#'   \item{ \strong{Feature} as a single character value. Must be in \code{getFeaturenNames()}}
#'   }
#'
#' @param all_features The valid features specified as a character vector.
#' @param all_genes The valid genes specified as a character vector.
#' @param all_gene_sets The valid gene sets specified as a character vector.
#' @param max_length The maximum number of elements the resulting list can have.
#'
#' @inherit adjusting_check_dummy details return
#' @export

check_color_to <- function(color_to,
                           all_features = character(),
                           all_gene_sets = character(),
                           all_genes = character(),
                           max_length = 25){

  if(base::is.list(color_to) & !base::is.data.frame(color_to)){

    color_to <- base::unlist(color_to)

  } else if(!base::is.character(color_to)){

    stop("Argument 'color_to' needs to be of class 'character' or of class 'list'.")

  }

  if(base::length(color_to) > max_length){

    base::stop(stringr::str_c("Maximum length (", max_length,
                              ") of argument 'color_to' exceeded: ",
                              base::length(color_to) ))

  }

  if(base::any(color_to %in% all_features)){

    found_features <- all_features[all_features %in% color_to]

  } else {

    found_features <- NULL

  }

  if(base::any(color_to %in% all_gene_sets)){

    found_gene_sets <- all_gene_sets[all_gene_sets %in% color_to]

  } else {

    found_gene_sets <- NULL

  }

  if(base::any(color_to %in% all_genes)){

    found_genes <- all_genes[all_genes %in% color_to]

  } else {

    found_genes <- NULL

  }

  found_all <- list("features" = found_features,
                    "gene_sets" = found_gene_sets,
                    "genes" = found_genes)

  if(base::length(base::unlist(found_all)) != base::length(color_to)){

    not_found <- color_to[!color_to %in% base::unlist(found_all)]

    not_found_string <- stringr::str_c(not_found, collapse = "', '")

    base::warning(stringr::str_c("Did not find '", not_found_string,
                                 "' of argument 'color_to' in the provided spata object.",
                                 sep = ""))

  }

  return_list <-
    purrr::discard(.x = found_all, .p = base::is.null)

  if(base::length(return_list) == 0){

    base::stop("Could not find any element of argument 'color_to' in the provided spata-object..")

  }

  base::return(return_list)

}

#' @title Check variables
#' @inherit check_color_to description
#'
#' @param variables The variables of interest specified as a character vector:
#'
#' \itemize{
#'  \item{ \strong{Gene sets}: Must be in \code{getGeneSets()}}
#'  \item{ \strong{Genes}: Must be in \code{getGenes()}}
#'  \item{ \strong{Features}: Must be a numeric one of \code{getFeaturenNames()}}
#'  }
#'
#' @param all_features Valid features.
#' @param all_gene_sets Valid gene sets.
#' @param all_genes Valid genes.
#' @param max_length Max number of variables.
#' @param simplify If set to TRUE the \code{check_variables()}-output is a vector.
#'
#' @export

check_variables <- function(variables,
                            all_features = character(),
                            all_gene_sets = character(),
                            all_genes = character(),
                            max_length = 25,
                            simplify = FALSE){

    if(base::is.list(variables) & !base::is.data.frame(variables)){

      variables <- base::unlist(variables) %>% base::unname()

    } else if(!base::is.character(variables)){

      stop("Argument 'variables' needs to be of class 'character' or of class 'list'.")

    }


    if(base::length(variables) > max_length){

      base::stop(stringr::str_c("Maximum length (", max_length,
                                ") of argument 'variables' exceeded: ",
                                base::length(variables) ))

    }

    if(base::any(variables %in% all_features)){

      found_features <- all_features[all_features %in% variables]

    } else {

      found_features <- NULL

    }

    if(base::any(variables %in% all_gene_sets)){

      found_gene_sets <- all_gene_sets[all_gene_sets %in% variables]

    } else {

      found_gene_sets <- NULL

    }

    if(base::any(variables %in% all_genes)){

      found_genes <- all_genes[all_genes %in% variables]

    } else {

      found_genes <- NULL

    }

    found_all <- list("features" = found_features,
                      "gene_sets" = found_gene_sets,
                      "genes" = found_genes)

    if(base::length(base::unlist(found_all)) != base::length(variables)){

      not_found <- variables[!variables %in% base::unlist(found_all)]

      not_found_string <- stringr::str_c(not_found, collapse = "', '")

      base::warning(stringr::str_c("Did not find feature(s), gene set(s) and/or gene(s) '", not_found_string, "' of argument 'variables'.", sep = ""))

    }

    return_list <-
      purrr::discard(.x = found_all, .p = base::is.null)

    if(base::length(return_list) == 0){

      base::stop("Could not find any of the specified feature(s), gene set(s) and/or gene(s) of argument 'variables'.")

    } else if(max_length == 1 | base::isTRUE(simplify)) {

      return_list <- base::unlist(return_list) %>% base::unname()

    }

    base::return(return_list)

}


#' @title Check feature variables input
#'
#' @description A member of the \code{adjusting-check_*()}-family. Takes a character
#' vector of feature names, checks which of the features exist and checks additionally
#' if these features match the class requirements of \code{valid_classes}.
#'
#' Returns an adjusted features-vector or raises an error.
#'
#' @param object A valid spata-object.
#' @param features The features of interest specified as a character vector.
#' @param valid_classes The feature-classes that are allowed.
#' @param max_length The maximum number of features allowed.
#'
#' @inherit adjusting_check_dummy details return
#' @export

check_features <- function(object,
                           features,
                           valid_classes = NULL,
                           max_length = NULL){

  # 1. Control --------------------------------------------------------------

  if(base::length(features) == 0 | !base::is.character(features)){

    base::stop("Invalid input for argument 'features'. Needs to be character vector of length > 0.")

  }

  # -----

  fnames <- getFeatureNames(object = object)

  # 2. Check if/how many features actually exist  ---------------------------

  if(!base::any(features %in% fnames)){

    base::stop("Could not find any of the specified features", "'. \n  Supplied features: '",
               stringr::str_c(features, collapse = "', '"), "'.")

  } else if(base::all(features %in% fnames)){

    fnames <- getFeatureNames(object)[getFeatureNames(object) %in% features]

  } else if(base::any(features %in% fnames)){

    fnames_found <- fnames[fnames %in% features]

    not_found <- stringr::str_c(features[!features %in% fnames_found], collapse = ", ")

    base::warning(stringr::str_c("Did not find feature(s):", not_found, sep = " "))

    fnames <- fnames_found

  }

  # -----

  # 3. Check which of the provided features are of valid classes ------------

  if(!base::is.null(valid_classes)){

    fclasses <- base::names(fnames)

    valid_fnames <- fnames[fclasses %in% valid_classes]

    if(length(valid_fnames) == 0){

      valid_classes <- stringr::str_c(valid_classes, collapse = "', '")

      base::stop(glue::glue("All features are of invalid classes. Valid classes are: '",
                                valid_classes, "'. \n  Supplied features: '",
                                stringr::str_c(features, collapse = "', '"), "'."))

    } else if(base::length(fnames) != base::length(valid_fnames)){

      not_valid <- stringr::str_c(fnames[!fnames %in% valid_fnames], collapse = ", ")
      valid_classes_string <- stringr::str_c(valid_classes, collapse = "' or '")

      base::warning(stringr::str_c("Ignoring feature(s) that are not of class '" ,valid_classes_string, "': ", not_valid, sep = ""))

    }

    fnames <- valid_fnames

  }

  # -----

  # 4. Check whether fnames is of desired length ----------------------------

  if(!base::is.null(max_length) &&
     base::length(fnames) > max_length) {

    base::warning(stringr::str_c("Reducing length of feature input to required length: ", max_length))
    fnames <- fnames[1:max_length]


  }

  # -----

  base::return(base::unname(fnames))

}


#' @title Check gene input
#'
#' @description A member of the \code{adjusting-check_*()}-family. Takes a character
#' vector of gene names and checks which of the genes exist.
#'
#' Returns an adjusted genes-vector or raises an error.
#'
#' @param genes The genes of interest specified as a character vector.
#' @param rna_assay The rna-assay you want to
#' look in. If set to NULL the whole rna_assay of the provided object will be used
#' with \code{exprMtr()}.
#'
#' @inherit adjusting_check_dummy details return
#' @export
check_genes <- function(object,
                        genes,
                        rna_assay = NULL,
                        max_length = NULL){

  # 1. Control --------------------------------------------------------------

  if(base::length(genes) == 0 | !base::is.character(genes)){

    base::stop("Invalid input for argument 'genes'. Needs to be character vector of length > 0.")

  }

  if(!base::is.matrix(rna_assay) && !base::is.null(rna_assay)){

    stop("Invalid input for argument 'rna_assay'.")

  }

  if(base::is.null(rna_assay)){

    rna_assay <- exprMtr(object = object)

  }

  # -----

  # 2. Check if/how many genes actually exist -------------------------------

  if(!base::any(genes %in% base::rownames(rna_assay))){

    stop("Could not find any of the specified genes.")

  } else if(base::all(genes %in% base::rownames(rna_assay))){

    genes_found <- genes

  } else if(base::any(genes %in% base::rownames(rna_assay))){

    genes_found <- base::rownames(rna_assay)[base::rownames(rna_assay) %in% genes]

    not_found <-
      genes[!genes %in% genes_found] %>% stringr::str_c(collapse = ", ")

    base::warning(stringr::str_c("Did not find genes: ", not_found, ".", sep = ""))

  }

  # -----

  # 3. Check whether genes found is of desired length -----------------------

  if(!base::is.null(max_length) &&
     base::length(genes_found) > max_length){

    base::warning(stringr::str_c("Reducing length of gene input to required length: ", max_length))

    genes_found <- genes_found[1:max_length]

  }

  # -----

  base::return(genes_found)

}


#' @title Check gene set input
#'
#' @description A member of the \code{adjusting-check_*()}-family. Takes a character
#' vector of gene set names and checks which of these exist.
#'
#' Returns an adjusted gene-set-vector or raises an error.
#'
#' @inherit check_sample params
#' @param gene_sets The gene sets of interest specified as a character vector.
#'
#' @inherit adjusting_check_dummy return details
#' @export

check_gene_sets <- function(object,
                            gene_sets,
                            max_length = NULL){


  # 1. Control --------------------------------------------------------------

  if(base::length(gene_sets) == 0 | !base::is.character(gene_sets)){

    base::stop("Invalid input for argument 'gene_sets'. Needs to be character vector of length > 0.")

  }


  # 2. Main part ------------------------------------------------------------

  gene_set_df <- object@used_genesets

  # 1. Check if/how many gene sets actually exists ----------

  if(base::all(gene_sets == "all")){

    gene_sets_found <- gene_set_df$ont %>% base::unique()

  } else if(!any(gene_sets %in% gene_set_df$ont)){

    stop("Could not find any specified geneset.")

  } else if(base::all(gene_sets %in% gene_set_df$ont)){

    gene_sets_found <- gene_sets

  } else if(base::any(gene_sets %in% gene_set_df$ont)){

    unique_gs <-
      dplyr::pull(gene_set_df, "ont") %>%
      base::unique()

    gene_sets_found <- unique_gs[unique_gs %in% gene_sets]

    not_found <- stringr::str_c(gene_sets[!gene_sets %in% gene_sets_found], collapse = ", ")

    base::warning(stringr::str_c("Could not find gene_sets:", not_found, sep = " "))

  }

  # -----


  # 2. Check whether gene sets found is of desired length ----------

  if(!base::is.null(max_length) &&
     base::length(gene_sets_found) > max_length){

    base::warning(stringr::str_c("Reducing length of gene set input to required length: ", max_length))

    gene_sets_found <- gene_sets_found[1:max_length]

  }

  # -----

  base::return(gene_sets_found)

}



#' @title Check sample input
#'
#' @description A member of the \code{adjusting-check_*()}-family. Takes a character
#' vector of sample names and checks which of these exist.
#'
#' Returns an adjusted sample-vector or raises an error.
#'
#' @param of_sample The sample(s) of interest specified as a single character value or vector.
#' @param desired_length The length the input must have.
#'
#' @inherit adjusting_check_dummy return details
#' @export

check_sample <- function(object,
                         of_sample,
                         desired_length = NULL){

  # 1. Check which samples are in the object --------------------------------

  if(!base::is.character(of_sample) | base::length(of_sample) == 0){

    stop("Please specify the sample with its name as a character vector of length > 0.")

  } else if(base::all(of_sample == "all")){

    of_sample <- samples(object = object)

    if(!base::is.null(desired_length) && base::length(of_sample) != desired_length){

      stop(stringr::str_c("Number of samples specified needs to be: ", desired_length, ". ",
                          "Setting 'of_sample' to 'all' results in ",
                          base::length(of_sample), " samples.", sep = ""))

    }

    base::return(of_sample)

  }  else if(!base::any(of_sample %in% object@samples)){

    stop("Could not find any of the specified samples in provided object.")

  } else if(base::any(of_sample %in% samples(object))){

    samples_found <- object@samples[object@samples %in% of_sample]

    if(base::length(of_sample) > 1){

      samples_found_string <- stringr::str_c(samples_found, collapse = ", ")

    } else if(base::length(of_sample) == 1){

      samples_found_string <- samples_found

    }

    if(base::length(samples_found) != base::length(of_sample)){

      base::warning(stringr::str_c("Did only find samples: ", samples_found_string, "."))

    }

  }

  # -----

  # 2. Check if length of samples found coincides with desired length -------

  if(!base::is.null(desired_length) &&
     base::length(samples_found) != desired_length){

    base::stop(stringr::str_c("Number of samples specified needs to be:", desired_length, sep = " "))

  }

  # -----

  base::return(samples_found)

}



#' @title Check segment name
#'
#' @description A member of the \code{adjusting-check_*()}-family. Takes the
#' segment name as a single character value, check whether such a segment
#' exists in the provided spata-object. If no an error is raised. Else the
#' barcodes of spots belonging to the specified segment are returned.
#'
#' @param object A valid spata-object.
#' @param segment_name The segment of interest specified as a single character
#' value.
#'
#' @inherit adjusting_check_dummy details return
#' @export

check_segment <- function(object,
                          segment_name){

  if(!is.null(segment_name) && !is.character(segment_name) && length(segment_name) != 1){

    base::stop("Argument 'segment_name' needs to be a single character value.")

  } else if(!is.null(segment_name)){

    bc_segm <-
      featureData(object, of_sample = of_sample) %>%
      dplyr::filter(segment == segment_name) %>%
      dplyr::pull(barcodes)

    if(base::length(bc_segm) == 0){

      base::stop(stringr::str_c("There is no segment of name' ", segment_name,
                                "' in sample '", of_sample, "'.", sep = ""))

    } else {

      base::return(bc_segm)

    }

  }

}

