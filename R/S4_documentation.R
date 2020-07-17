
# gene expression data ----------------------------------------------------

#' data_counts object
#'
#' @slot counts matrix.
#' @slot norm_exp matrix.
#'
#' @return S4 object
#' @export
#'

data_counts <- setClass("data_counts",
                        slots = c(counts = "matrix",
                                  norm_exp  = "matrix"))



# dimensional reduction ---------------------------------------------------

#' dim_red object
#'
#' @slot UMAP data.frame.
#' @slot TSNE data.frame.
#'
#' @return S4 object
#' @export
#'

dim_red <- setClass("dim_red",
                    slots = c(UMAP =  "data.frame",
                              TSNE ="data.frame"))



# single cell velocity ----------------------------------------------------

#' scvelo object
#'
#' @slot scvelo data.frame.
#' @slot scvelo1 data.frame.
#'
#' @return S4 object
#' @export
#'

scvelo <- setClass("scvelo",
                   slots = c(scvelo = "data.frame",
                             scvelo1 ="data.frame"))



# spatial trajectory ------------------------------------------------------

#' spatialTrajectory object
#'
#' @slot compiled_trajectory_df data.frame.
#' @slot segment_trajectory_df data.frame.
#' @slot ranked_genes_df data.frame.
#' @slot comment character.
#' @slot name character.
#' @slot sample character.
#'
#' @return S4 object
#' @export
#'

spatialTrajectory <- setClass("spatialTrajectory",
                              slots = c(
                                compiled_trajectory_df = "data.frame",
                                segment_trajectory_df = "data.frame",
                                ranked_genes_df = "data.frame",
                                comment = "character",
                                name = "character",
                                sample = "character"
                              ))



# spata object ------------------------------------------------------------

#' spata object
#'
#' @slot coordinates data.frame.
#' @slot fdata data.frame.
#' @slot samples character.
#' @slot data data_counts.
#' @slot image list.
#' @slot description character.
#' @slot dim_red dim_red.
#' @slot scvelo scvelo.
#' @slot used_genesets data.frame.
#' @slot trajectories list.
#'
#' @return S4 object
#' @export
#'

spata <- setClass("spata",
                  slots = c(coordinates ="data.frame", #coordinates: bc, x, y, sample
                            fdata = "data.frame", #fdata : bc, ...
                            samples = "character",
                            data = "data_counts" ,
                            image = "list",
                            description = "character",
                            dim_red = "dim_red", #UMAP & TSNE: bc, umap1, umap2, sample
                            scvelo = "scvelo",
                            used_genesets = "data.frame",
                            trajectories = "list"
))
