
#' Check input
#'
#' @param color_to character or list
#' @param all_features character. all features
#' @param all_gene_sets character. all gene sets
#' @param all_genes character. all genes
#'
#' @return a named list
#' @export
#'

check_color_to <- function(color_to,
                           all_features = character(),
                           all_gene_sets = character(),
                           all_genes = character()){



  if(base::is.list(color_to) & !base::is.data.frame(color_to)){

    color_to <- base::unlist(color_to)

  } else if(!base::is.character(color_to)){

    stop("Argument 'color_to' needs to be of class 'character' or of class 'list'.")

  }


  if(base::length(color_to) > 25){

    base::stop("Argument 'color_to' needs to be of length < 25.")

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

    base::warning(stringr::str_c("Did not find gene set(s) and/or gene(s) '", not_found_string, "' of argument 'color_to'.", sep = ""))

  }

  return_list <-
    purrr::discard(.x = found_all, .p = base::is.null)

  if(base::length(return_list) == 0){

    base::stop("Could not find any of the specified gene set(s) and/or gene(s) of argument 'color_to'.")

  } else {

    base::return(return_list)

  }

}


#' @title Check input
#'
#' @description The \code{check_()}-function family checks a functions input to
#' make sure that the function can run bugfree. If the input does not suffice a
#' helping message is printed for the user to know what needs to be adjusted in
#' order to make the function work.
#'
#' @param coords_df A data.frame that contains two character variables: \emph{
#' barcodes & sample}
#'
#' @return
#' - If the input is just fine: the input.
#' - If the input suffices for the function to work but needs slight moderation:
#' a modified version of the input along with a warning-message.
#' - If the input does not suffice for the function to work: an error
#'

check_coords_df <- function(coords_df){

  if(!base::is.data.frame(coords_df)){

    base::stop("Argument 'coords_df' needs to be a data.frame.")

  } else if(!base::all(c("barcodes", "sample") %in% base::colnames(coords_df))){

    base::stop("'coords_df' needs to have 'barcodes' and 'sample' variables.")

  } else {

    classes <- base::sapply(X = coords_df[,c("barcodes", "sample")],
                            FUN = base::class)

    if(!base::all(classes == "character")){

      base::stop("Variables 'barcodes' and 'sample' need to be of class character.")

    } else {

      return(coords_df)

    }

  }


}


#' @title Check input
#'
#' @param object A valid spata-object.
#' @param features The feature-input.
#' @param valid_classes The feature-classes that are allowed.
#' @param max_length The maximum length the input can have.
#'
#' @return
#' - If the input is just fine: the input.
#' - If the input suffices for the function to work but needs slight moderation:
#' a modified version of the input along with a warning-message.
#' - If the input does not suffice for the function to work: an error
#'
#' @description The \code{check_()}-function family checks a functions input to
#' make sure that the function can run bug-free. If the input does not suffice a
#' helping message is printed for the user to know what needs to be adjusted in
#' order to make the function work.
#'

check_features <- function(object,
                           features,
                           valid_classes = NULL,
                           max_length = NULL){

  if(base::length(features) == 0 | !base::is.character(features)){

    base::stop("Invalid input for argument 'features'. Needs to be character vector of length > 0.")

  }

  fnames <- getFeatureNames(object = object)

  # 1. Check if/how many features actually exist  ---------------------------

  if(!base::any(features %in% fnames)){

    base::stop("Could not find any of the specified features.")

  } else if(base::all(features %in% fnames)){

    fnames <- getFeatureNames(object)[getFeatureNames(object) %in% features]

  } else if(base::any(features %in% fnames)){

    fnames_found <- fnames[fnames %in% features]

    not_found <- stringr::str_c(features[!features %in% fnames_found], collapse = ", ")

    base::warning(stringr::str_c("Did not find feature(s):", not_found, sep = " "))

    fnames <- fnames_found

  }


  # 2. Check which of the provided features match the 'class' requir --------

  if(!base::is.null(valid_classes)){

    fclasses <- base::names(fnames)

    valid_fnames <- fnames[fclasses %in% valid_classes]

    if(length(valid_fnames) == 0){

      base::stop("All features are of invalid classes.")

    } else if(base::length(fnames) != base::length(valid_fnames)){

      not_valid <- stringr::str_c(fnames[!fnames %in% valid_fnames], collapse = ", ")
      valid_classes_string <- stringr::str_c(valid_classes, collapse = "' or '")

      base::warning(stringr::str_c("Ignoring feature(s) that are not of class '" ,valid_classes_string, "': ", not_valid, sep = ""))

    }

    fnames <- valid_fnames

  }


# 3. Check whether fnames is of desired length ----------------------------

  if(!base::is.null(max_length) &&
     base::length(fnames) > max_length) {

    base::warning(stringr::str_c("Reducing length of feature input to required length: ", max_length))
    fnames <- fnames[1:max_length]


  }


  base::return(base::unname(fnames))

}


#' @inheritParams check_features
#' @param genes The genes-input.
#' @param rna_assay A rna-assay (e.g. derived from \code{exprMtr()}).
#'
#' @inherit check_features description return title
#'
check_genes <- function(object,
                        genes,
                        rna_assay = NULL,
                        max_length = NULL){

  if(base::length(genes) == 0 | !base::is.character(genes)){

    base::stop("Invalid input for argument 'genes'. Needs to be character vector of length > 0.")

  }

  if(!is.matrix(rna_assay) | is.null(rna_assay)){

    stop("Invalid input for argument 'rna_assay'.")

  }

  # 1. Check if/how many genes actually exist -------------------------------

  if(!base::any(genes %in% base::rownames(rna_assay))){

    stop("Could not find any of supplied genes.")

  } else if(base::all(genes %in% base::rownames(rna_assay))){

    genes_found <- genes

  } else if(base::any(genes %in% base::rownames(rna_assay))){

    genes_found <- base::rownames(rna_assay)[base::rownames(rna_assay) %in% genes]

    not_found <-
      genes[!genes %in% genes_found] %>% stringr::str_c(collapse = ", ")

    base::warning(stringr::str_c("Did not find genes: ", not_found, ".", sep = ""))

  }


  # 2. Check whether genes found is of desired length -----------------------

  if(!base::is.null(max_length) &&
     base::length(genes_found) > max_length){

    base::warning(stringr::str_c("Reducing length of gene input to required length: ", max_length))

    genes_found <- genes_found[1:max_length]

  }


  base::return(genes_found)

}



#' @inherit check_features description return title params
#'
#' @param gene_sets The gene_sets-input.
#'

check_gene_sets <- function(object,
                            gene_sets,
                            max_length = NULL){

  if(base::length(gene_sets) == 0 | !base::is.character(gene_sets)){

    base::stop("Invalid input for argument 'gene_sets'. Needs to be character vector of length > 0.")

  }

  gene_set_df <- object@used_genesets


  # 1. Check if/how many gene sets actually exists --------------------------

  if(!any(gene_sets %in% gene_set_df$ont)){

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


  # 2. Check whether gene sets found is of desired length -------------------

  if(!base::is.null(max_length) &&
     base::length(gene_sets_found) > max_length){

    base::warning(stringr::str_c("Reducing length of gene set input to required length: ", max_length))

    gene_sets_found <- gene_sets_found[1:max_length]

  }


  base::return(gene_sets_found)


}




#' @inherit check_features description return title
#'
#' @param pt_size The pt_size-input.
#' @param pt_alpha The pt_alpha-input.
#' @param pt_clrsp The pt_clrsp-input.
#'

check_pt_input <- function(pt_size = NULL,
                           pt_alpha = NULL,
                           pt_clrsp = NULL){

  if(!base::is.null(pt_clrsp) && !pt_clrsp %in% c("inferno", "magma", "plasma", "cividis", "viridis")){

    base::stop("Argument 'pt_clrsp' needs to be one of 'inferno', 'magma', 'plasma', 'cividis' or 'viridis'.")

  }

  if(!base::is.null(pt_size) && !base::is.numeric(pt_size)){

    base::stop("Argument 'pt_size' needs to be numeric.")

  }

  if(!base::is.null(pt_alpha) && !base::is.numeric(pt_alpha)){

    base::stop("Argument 'pt_alpha' needs to be numeric.")

  }

}





#' @inherit check_features description return params title
#'
#' @param sample_input The sample input.
#' @param desired_length The length the input must have.
#'

check_sample <- function(object,
                         sample_input,
                         desired_length = NULL){

  if(!base::is.character(sample_input) | base::length(sample_input) == 0){

    stop("Please specify the sample with its name as a character vector of length > 0.")

  } else if(base::length(sample_input) == 1 && sample_input == "all"){

    sample_input <- samples(object = object)

    if(!base::is.null(desired_length) && base::length(sample_input) != desired_length){

      stop(stringr::str_c("Number of samples specified needs to be: ", desired_length, ". ",
                          "Setting 'of_sample' to 'all' results in ",
                          base::length(sample_input), " samples.", sep = ""))

    }

    return(sample_input)

  }  else if(!base::any(sample_input %in% object@samples)){

    stop("Could not find sample(s) in provided object.")

  } else if(base::any(sample_input %in% samples(object))){

    samples_found <- object@samples[object@samples %in% sample_input]

    if(length(sample_input) > 1){

      samples_found_string <- stringr::str_c(samples_found, collapse = ", ")

    } else if(base::length(sample_input) == 1){

      samples_found_string <- samples_found

    }

    if(base::length(samples_found) != base::length(sample_input)){

      base::warning(stringr::str_c("Did only find samples: ", samples_found_string, "."))

    }

  }


  # 2. Check if length of samples found coincides with desired length -------


  if(!base::is.null(desired_length) &&
     base::length(samples_found) != desired_length){

    base::stop(stringr::str_c("Number of samples specified needs to be:", desired_length, sep = " "))

  }

  base::return(samples_found)

}


#' @inherit check_features return title
#'
#' @description Checks among other things whether \code{df} is spatially 'smoothable' -
#' whether it contains x and y variables.
#'
#' @param df The data.frame that is to be smoothed.
#' @param smooth The smooth input.
#' @param verbose Logical.
#'

check_smooth <- function(df, smooth, verbose = TRUE){

  if(!base::isTRUE(smooth) & !base::isFALSE(smooth)){

    base::stop("Argument 'smooth' needs to be TRUE or FALSE.")

  }

  if(base::isTRUE(smooth)){

    if(!base::all(c("x", "y") %in% base::colnames(df))){

      base::warning("Input data.frame doesn't contain x and y variables. Skip smoothing." )

      base::return(FALSE)

    } else {

      base::return(TRUE)

    }

  } else {

    base::return(FALSE)

  }


}

