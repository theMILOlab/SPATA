#' @title Reduce the number of samples
#'
#' @description Keeps only information of the specified sample
#'
#' @param object A valid spata-object.
#' @param sample Character value. The sample of interest.
#'
#' @return An updated spata-object.
#' @export
#'

subsetSpataObject <- function(object, sample){

  check_object(object)
  sample <- check_sample(object, sample, 1)

  if(base::length(samples(object)) == 1){

    base::stop("Spata-object contains only one sample.")

  }

  counts_new <- getCountMatrix(object, sample)
  expr_mtr_new <- getExpressionMatrix(object, sample)

  fdata_new <- getFeatureData(object, sample)

  coords_new <- coordinates(object, sample)

  trajectories_new <- object@trajectories[[sample]]

  if(base::is.null(trajectories_new)){

    trajectories_new <- list(list()) %>% magrittr::set_names(sample)

  }

  sample_image <- object@image[[sample]]
  images_new <- list(sample_image) %>% magrittr::set_names(sample)

  additional_new <- list(Seurat = list(images = list()))
  additional_new$Seurat$images[[sample]] <- object@additional$Seurat$images[[sample]]

  umap_new <- getUmapData(object, sample)
  tsne_new <- getTsneData(object, sample)

  new_object <-
    methods::new("spata",
                 coordinates = coords_new,
                 data = methods::new("data_counts",
                                     counts = counts_new,
                                     norm_exp = expr_mtr_new),
                 dim_red = methods::new("dim_red",
                                        UMAP = umap_new,
                                        TSNE = tsne_new),
                 fdata = fdata_new,
                 samples = sample,
                 used_genesets = object@used_genesets,
                 trajectories = trajectories_new,
                 version = spata_version,
                 additional = additional_new)

}
