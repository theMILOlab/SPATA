#' @title Subset a spata-object
#'
#' @description Reduces the number of samples in the
#' specified spata-object to one.
#'
#' @param object A valid spata-object.
#' @param sample Character value. The sample of interest.
#'
#' @return An updated spata-object.
#' @export
#'

subsetSpataObject <- function(object, sample){

  check_object(object)
  sample <- check_sample(object, of_sample = sample, of.length = 1)

  if(base::length(getSampleNames(object)) == 1){

    base::stop("Spata-object contains only one sample.")

  }

  new_object <- methods::new(Class = "spata",
                             coordinates = object@coordinates[sample],
                             data = object@data[sample],
                             dim_red = object@dim_red[sample],
                             fdata = object@fdata[sample],
                             images = object@images[sample],
                             information = purrr::map(.x = object@information, .f = ~ .x[sample]),
                             dea = object@dea[sample],
                             samples = sample,
                             scvelo = object@scvelo[sample],
                             trajectories = object@trajectories[sample],
                             used_genesets = object@used_genesets,
                             version = object@version)

  base::return(new_object)

}
