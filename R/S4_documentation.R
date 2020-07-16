
# gene expression data ----------------------------------------------------

data_counts <- setClass("data_counts",
                        slots = c(counts = "matrix",
                                  norm_exp  = "matrix"))



# dimensional reduction ---------------------------------------------------

dim_red <- setClass("dim_red",
                    slots = c(UMAP =  "data.frame",
                              TSNE ="data.frame"))



# single cell velocity ----------------------------------------------------

scvelo <- setClass("scvelo",
                   slots = c(scvelo = "data.frame",
                             scvelo1 ="data.frame"))



# spatial trajectory ------------------------------------------------------

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
