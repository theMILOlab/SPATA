# DEPRECATED FUNCTIONS ----------------------------------------------------

# alphabetically arranged

#' Deprecated functions
#'
#' @include S4-documentation.R
#'

setGeneric(name = "dimRed", def = function(object, of_sample, method){

  standardGeneric(f = "dimRed")

})
setMethod(f = "dimRed", signature = "spata", definition = function(object, of_sample, method){


  warning("'dimRed' is deprecated. Use coordsDimRed() instead.")

  if(length(method) != 1){

    stop("Argument 'method' needs to be of length one.")

  } else if(!method %in% c("UMAP", "TSNE")) {

    stop("Argument 'method' needs to be  'UMAP' or 'TSNE'.")

  }

  of_sample <- check_sample(object = object, sample_input = of_sample)

  dim_red_df <-
    methods::slot(object = object@dim_red, name = method) %>%
    dplyr::filter(sample %in% of_sample)

  dim_red_df <- dim_red_df[,c(stringr::str_c(tolower(method), 1, sep = ""),
                              stringr::str_c(tolower(method), 2, sep = ""),
                              "barcodes")]


  if(nrow(dim_red_df) == 0){

    message(stringr::str_c("No dimensional reduction of method", method, "has been performed yet. Returning empty data.frame.", sep = " "))

  }

  return(dim_red_df)

})

#' @rdname dimRed
setGeneric(name = "featureNames", valueClass = "character", def = function(object){

  standardGeneric(f = "featureNames")

})
setMethod(f = "featureNames", signature = "spata", definition = function(object){

  feature_names <- base::colnames(object@fdata)

  base::names(feature_names) <-
    base::sapply(object@fdata[,feature_names], base::class)

  base::warning("The function 'featureNames' is deprecated! Use 'getFeatureNames' instead.")

  return(feature_names[!feature_names %in% c("sample", "barcodes")])


})

#' @rdname dimRed
setGeneric(name = "getTrajectoryObject", def = function(object, trajectory_name, of_sample){

  standardGeneric(f = "getTrajectoryObject")

})
setMethod(f = "getTrajectoryObject", signature = "spata", definition = function(object, trajectory_name, of_sample){

  warning("'getTrajectoryObject' is deprecated. Use 'trajectory' instead.")

  of_sample <- check_sample(object = object, sample_input = of_sample, desired_length = 1)

  if(!is.character(trajectory_name) | length(trajectory_name) != 1){

    stop("Argument 'trajectory_name' needs to be a character vector of length 1.")

  }

  t_names <- base::names(object@trajectories[[of_sample]])

  if(trajectory_name %in% t_names){

    trajectory_object <- object@trajectories[[of_sample]][[trajectory_name]]

    return(trajectory_object)

  } else {

    stop(stringr::str_c("Could not find trajectory '", trajectory_name, "' in sample: '", of_sample, "'.", sep = ""))

  }

})
