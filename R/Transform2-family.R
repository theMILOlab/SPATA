#' @title From_Seurat_to_SPATA
#'
#' @description This Function will return a data frame with the sourounding spots
#'
#' @param seurat_object
#' @param sample_names if NULL sample name will be created
#' @param verbose
#'
#'
#'
#' @return A spata-object.
#'
#' @importFrom Seurat ScaleData
#'
#' @export
From_Seurat_to_SPATA <- function(seurat_object, sample_names=NULL,verbose=T){


  # Input create ------------------------------------------------------------
  print("convert a scRNA-seq data set to Spata")
  if(base::is.null(sample_names)){sample_names="sample_1"}
  coordinates=data.frame(seurat_object@reductions$umap@cell.embeddings) %>% tibble::rownames_to_column("barcodes")
  names(coordinates)[2:3] <- c("x","y")
  intensity_matrix <- as.matrix(GetAssayData(seurat_object))
  meta_data <- seurat_object@meta.data %>% tibble::rownames_to_column("barcodes")
  counts <- seurat_object@assays$RNA@counts
  geneSets=SPATA::gsdf
  string_contain=length(stringr::str_detect(rownames(intensity_matrix), "_") %in% TRUE) <= 1
  if(string_contain) {
    message(" rownames of intensity matrix contain a '_'. This is not allowed and will be changed into '-' ")
    rownames(intensity_matrix) <- str_replace(rownames(intensity_matrix), "_", "-")
  }



  #check inputs
  if (base::length(base::unique(coordinates$barcodes %in% colnames(intensity_matrix)))!=1) stop("barcodes from coordinates are not similar to colnames of the intensity_matrix")

  #Add sample name to coordinates
  coordinates <- coordinates %>% dplyr::mutate(sample=sample_names) %>% dplyr::select("barcodes","sample","x","y")

  #Change Coordinates for image in SPATA

  # Start -------------------------------------------------------------------


  #Crease empty SPATA object
  obj=methods::new(Class="spata")

  obj@coordinates <- coordinates
  obj@data@norm_exp <- intensity_matrix
  obj@data@counts=counts
  obj@samples <- unique(obj@coordinates$sample)
  obj@fdata <- obj@coordinates %>% dplyr::select(barcodes, sample) %>% dplyr::mutate(segment = "")
  obj@fdata <- obj@fdata %>% dplyr::left_join(., meta_data, by="barcodes")
  names(obj@fdata)[2] <- "sample"
  obj@used_genesets=geneSets
  obj@trajectories=list(sample=list())
  names(obj@trajectories)=obj@samples

  obj@dim_red@UMAP <- data.frame(seurat_object@reductions$umap@cell.embeddings) %>% tibble::rownames_to_column("barcodes")

  return(obj)



}
