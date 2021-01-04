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

subsetBySample <- function(object, sample){

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


#' @title Subset a spata-object by segment
#'
#' @inherit check_sample params
#' @inherit initiateSpataObject_CountMtr
#' @param sgement_name Character value. The segment according to which the spata-object is
#' to be subsetted.
#'
#' @return An updated spata-object.
#' @export
#'

subsetBySegment <- function(object,
                            segment_name,
                            of_sample = "",
                            SCTransform = FALSE,
                            NormalizeData = list(normalization.method = "LogNormalize", scale.factor = 1000),
                            FindVariableFeatures = list(selection.method = "vst", nfeatures = 2000),
                            ScaleData = TRUE,
                            RunPCA = list(npcs = 60),
                            FindNeighbors = list(dims = 1:30),
                            FindClusters = list(resolution = 0.8),
                            RunTSNE = TRUE,
                            RunUMAP = list(dims = 1:30)){

  # 1. Control --------------------------------------------------------------

  check_object(object)

  of_sample <-
    check_sample(object = object, of_sample = of_sample, of.length = 1)

  barcodes <-
    check_segment(object = object, segment_name = segment_name, of_sample = of_sample)


  # 2. Data extraction ------------------------------------------------------

  gene_set_df <-
    getGeneSetDf(object)

  coords_df <-
    getCoordsDf(object = object, of_sample = of_sample)

  segment_coords_df <-
    dplyr::filter(coords_df, barcodes %in% {{barcodes}})

  old_coords_df <-
    dplyr::filter(coords_df, !barcodes %in% {{barcodes}})

  image <-
    getImage(object = object, of_sample = of_sample)

  count_mtr <-
    getCountMatrix(object = object, of_sample = of_sample)[, barcodes]

  seurat_object <-
    process_seurat_object(
      seurat_object = Seurat::CreateSeuratObject(counts = count_mtr),
      SCTransform = SCTransform,
      NormalizeData = NormalizeData,
      FindVariableFeatures = FindVariableFeatures,
      ScaleData = ScaleData,
      RunPCA = RunPCA,
      FindNeighbors = FindNeighbors,
      FindClusters = FindClusters,
      RunTSNE = RunTSNE,
      RunUMAP = RunUMAP,
      verbose = verbose
    )

  spata_object <-
    transformSeuratToSpata(
      seurat_object = seurat_object,
      assay = "RNA",
      sample_name = of_sample,
      gene_set_path = NA,
      method = "single_cell",
      verbose = verbose
    )

  spata_object <-
    setCoordsDf(object = spata_object, coords_df = segment_coords_df, of_sample = of_sample) %>%
    setImage(object = ., image = image, of_sample = of_sample) %>%
    setGeneSetDf(object = ., gene_set_df = gene_set_df)


  spata_object@information$old_coordinates <- old_coords_df

  base::return(spata_object)

}
