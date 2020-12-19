
#' @title Update a spata object
#'
#' @description Checks if the provided spata object needs any updates regarding
#' it's structure by comparing it's @@version to the newest one.
#'
#' @inherit check_object params
#'
#' @return A valid spata object.
#' @export

updateSpataObject <- function(object){

  check_object(object)

  version <- base::tryCatch(

    expr = object@version,
    error = function(error){

      # if no version found assume that it is the earliest version
      list(major = 0, minor = 0, patch = 0, dev = 9000)

      }

  )


# Version < 1.1.0 ---------------------------------------------------------

  # adds slots '@dea = list()', '@information = list()'
  # transforms '@data from S4 to list()', '@dim_red from S4 list()'
  if(version$major <= 1 & version$minor < 1){

    # extract data from old slots

    # data.frames
    coords_df <- object@coordinates
    feature_df <- object@fdata

    umap_df <- object@dim_red@UMAP
    tsne_df <- object@dim_red@TSNE

    gs_df <- object@used_genesets

    # assays
    count_mtr <- object@data@counts
    expr_mtr <- object@data@norm_exp

    # lists
    images <- object@image
    trajectories <- object@trajectories
    scvelo <- object@scvelo

    # vectors
    samples <- object@samples

    new_object <-
      methods::new(Class = "spata",
                   coordinates = coords_df,
                   data = list("counts" = count_mtr, "scaled" = expr_mtr),
                   dim_red = list("umap" = umap_df, "tsne" = tsne_df),
                   fdata = feature_df,
                   image = images,
                   information = list("autoencoder" = list()),
                   dea = list(),
                   samples = samples,
                   scvelo = scvelo,
                   trajectories = trajectories,
                   used_genesets = gs_df,
                   version = current_spata_version)

  } else {

    base::message("According to slot 'version' the provided spata object does not need any updating.")

    new_object <- object

  }

  base::return(new_object)

}
