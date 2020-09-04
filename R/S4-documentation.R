
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
                        slots = c(counts = "Matrix",
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




# spatial trajectory ------------------------------------------------------

#' spatial_trajectory object
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

spatial_trajectory <- setClass("spatial_trajectory",
                                slots = c(
                                  compiled_trajectory_df = "data.frame",
                                  segment_trajectory_df = "data.frame",
                                  assessed_trajectory_df = "data.frame",
                                  comment = "character",
                                  name = "character",
                                  sample = "character"))



# spata object ------------------------------------------------------------

#' spata object
#'
#' @slot coordinates data.frame.
#' @slot data data_counts.
#' @slot dim_red dim_red.
#' @slot fdata data.frame.
#' @slot image list.
#' @slot samples character.
#' @slot scvelo scvelo.
#' @slot trajectories list.
#' @slot used_genesets data.frame.
#' @slot version list.
#'
#' @return S4 object
#' @export
#'

spata <- setClass("spata",
                  slots = c(coordinates ="data.frame", #coordinates: bc, x, y, sample
                            data = "data_counts",
                            dim_red = "dim_red", #UMAP & TSNE: bc, umap1, umap2, sample
                            fdata = "data.frame", #fdata : bc, ...
                            image = "list",
                            samples = "character",
                            scvelo = "list",
                            used_genesets = "data.frame",
                            trajectories = "list",
                            version = "list"))
