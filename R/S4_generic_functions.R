

# Accessor functions ------------------------------------------------------

setGeneric(name = "coordinates", valueClass = "data.frame", def = function(object, of_sample = "all"){

  standardGeneric("coordinates")

})
setMethod(f = "coordinates", signature = "spata", def = function(object, of_sample = "all"){

  of_sample <- check_sample(object = object, sample_input = of_sample)

  ##----- filter for bc in sample
  coords_df <-
    object@coordinates %>%
    dplyr::filter(sample %in% of_sample)

  return(coords_df)

})

setGeneric(name = "exprMtr", def = function(object, of_sample = "all"){

  standardGeneric(f = "exprMtr")

})
setMethod(f = "exprMtr", signature = "spata", definition = function(object, of_sample = "all"){

  of_sample <- check_sample(object = object, sample_input = of_sample)

  bc_in_sample <-
    object@fdata %>%
    dplyr::filter(sample %in% of_sample) %>%
    dplyr::pull(barcodes)


  exprMtr <- object@data@norm_exp[,bc_in_sample]

  # filter for genes that were actually expressed in the sample
  rows_to_subset <- base::rowSums(exprMtr) != 0

  return(exprMtr[rows_to_subset, ])

})

setGeneric(name = "featureData", valueClass = "data.frame", def = function(object, of_sample = "all"){

  standardGeneric(f = "featureData")

})
setMethod(f = "featureData", signature = "spata", definition = function(object, of_sample = "all"){


  of_sample <- check_sample(object = object, sample_input = of_sample)

  fdata <-
    as.data.frame(object@fdata) %>%
    dplyr::filter(sample %in% of_sample)


  return(fdata)

})

setGeneric(name = "getTrajectoryComment", def = function(object, ...){

  standardGeneric(f = "getTrajectoryComment")

})
setMethod(f = "getTrajectoryComment", signature = "spata", definition = function(object, trajectory_name, of_sample){


  if(!is.character(trajectory_name) | length(trajectory_name) != 1){

    stop("Argument 'trajectory_name' needs to be a character vector of length 1.")

  }

  t_names <- base::names(object@trajectories[[of_sample]])

  if(trajectory_name %in% t_names){

    trajectory_object <- object@trajectories[[of_sample]][[trajectory_name]]

    return(trajectory_object@comment)

  } else {

    stop(stringr::str_c("Could not find trajectory '", trajectory_name, "' in sample '", of_sample, "'.", sep = ""))

  }

})
setMethod(f = "getTrajectoryComment", signature = "spatialTrajectory", definition = function(object){


  return(stringr::str_c("Comment: ", object@comment))


})


setGeneric(name = "image", def = function(object, of_sample){

  standardGeneric(f = "image")

})
setMethod(f = "image", signature = "spata", definition = function(object, of_sample){

  of_sample <- check_sample(object = object, sample_input = of_sample, desired_length = 1)

  return(object@image[[of_sample]])

})

setGeneric(name = "samples", def = function(object){

  standardGeneric(f = "samples")

})
setMethod(f = "samples", signature = "spata", definition = function(object){

  object@samples

})

setGeneric(name = "trajectory", def = function(object, trajectory_name, of_sample){

  standardGeneric(f = "trajectory")

})
setMethod(f = "trajectory", signature = "spata", definition = function(object, trajectory_name, of_sample){

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



# Setter functions --------------------------------------------------------

setGeneric(name = "coordinates<-", def = function(object, value){

  standardGeneric(f = "coordinates<-")

})
setMethod(f = "coordinates<-", signature = "spata", def = function(object, value){

  object@coordinates <- value

  return(object)

})

setGeneric(name = "featureData<-", def = function(object, value){

  standardGeneric(f = "featureData<-")

})
setMethod(f = "featureData<-", signature = "spata", definition = function(object, value){


  object@fdata <- value

  return(object)

})



# To be dealt with --------------------------------------------------------

setGeneric(name = "getRankedGenes", def = function(object, of_sample, ...){

  standardGeneric(f = "getRankedGenes")

})
setMethod(f = "getRankedGenes", signature = "spata", definition = function(object, trajectory_name, of_sample){

  of_sample <- check_sample(object, sample_input = of_sample, desired_length = 1)

  t_object <- getTrajectoryObject(object, trajectory_name = trajectory_name, of_sample = of_sample)

  return(t_object@ranked_genes_df)

})


