#' @title getRNAvelocity
#'
#' @description This Function will return a SPATA object containing denoised data and re-analysis
#'
#' @param SPATAobj
#' @param Velocity_from can be "spatial" if you aim to project in the spatial space or a "UMAP")
#' @param path_to_pythonscript The py skript sc2velo.py is required and in in the data folder of the SPATA github
#' @param loom loom file
#' @param path_to_python where is you python 3.8 (with installed scvelo)
#' @param folder folder where files will be saved
#'
#'
#'
#' @return A spata-object.
#'
#' @importFrom Seurat ScaleData
#'
#' @export



getRNAvelocity <- function(object,
                           Velocity_from=c("spatial", "UMAP"),
                           path_to_pythonscript="'/Users/HenrikHeiland/Desktop/T-Cell\ Project/Revisions/sc2velo.py'",
                           loom=NULL,
                           h5ad=NULL,
                           path_to_python="/Users/HenrikHeiland/opt/anaconda3/bin/python3.7 ",
                           folder=getwd()){


# loom into seurat --------------------------------------------------------
  if(is.null(h5ad)){
  ldat <- SeuratWrappers::ReadVelocity(file = loom)

  so <- Seurat::as.Seurat(x = ldat)
  require(Seurat)
  so[["RNA"]] <- so[["spliced"]]

  # Adapting the loom and spata output
  require(tidyverse)
  feature_change <-
    so@meta.data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("rownames") %>%
    tidyr::separate(.,rownames, sep=":", into=c("sample", "bc"), remove=F) %>%
    tidyr::separate(.,bc, sep="x", into=c("barcodes", "NX") ) %>%
    dplyr::mutate(barcodes=paste0(barcodes, "-1_", object@samples))
    #dplyr::select(barcodes, nFeature_ambiguous, nCount_RNA, nFeature_RNA)
  base::rownames(feature_change) <- feature_change$barcodes
  so@meta.data <- feature_change

  #adopt assays
  base::colnames(so@assays$spliced@counts) <- feature_change$barcodes
  base::colnames(so@assays$spliced@data) <- feature_change$barcodes

  base::colnames(so@assays$unspliced@counts) <- feature_change$barcodes
  base::colnames(so@assays$unspliced@data) <- feature_change$barcodes

  base::colnames(so@assays$ambiguous@counts) <- feature_change$barcodes
  base::colnames(so@assays$ambiguous@data) <- feature_change$barcodes

  base::colnames(so@assays$RNA@counts) <- feature_change$barcodes
  base::colnames(so@assays$RNA@data) <- feature_change$barcodes

  Seurat::Idents(object=so) <- "barcodes"

  base::message("... Run seurat sc_velo wrapper .... ")

  so <- Seurat::SCTransform(so, variable.features.n=2000, verbose=F)
  so <- Seurat::RunPCA(so, verbose=F)
  so <- Seurat::RunUMAP(so, dims = 1:30, verbose=F)
  so <- Seurat::FindNeighbors(so, dims = 1:30, verbose=F)
  so <- Seurat::FindClusters(so, verbose=F)

  base::message("... Adopt SPATA and loom files .... ")
  so <- subset(so, cells=object@fdata$barcodes)

  if(base::nrow(object@fdata) != base::nrow(so@meta.data)) stop("The output of SPATA and loom is unequal different number of cells will affect the analysis")

  so@meta.data$seurat_clusters <-
    so@meta.data %>%
    as.data.frame() %>%
    dplyr::left_join(., object@fdata, by="barcodes") %>%
    dplyr::pull(seurat_clusters.y)


  base::message("... Return objects  for python.... ")

  setwd(folder)

  Seurat::DefaultAssay(so) <- "RNA"
  SeuratDisk::SaveH5Seurat(so, filename = "Export.h5Seurat")
  SeuratDisk::Convert("Export.h5Seurat", dest = "h5ad")



  #Export additional files required
  if(Velocity_from=="spatial"){
    rownames(object@coordinates) <- object@coordinates$barcodes
    UMAP <-
      object@coordinates[so@meta.data$barcodes, ] %>%
      select(x,y)
    names(UMAP) <- c("UMAP_1","UMAP_2")
    write.csv(UMAP, "UMAP.csv")
  }else{
    rownames(object@dim_red@UMAP) <- object@dim_red@UMAP$barcodes
    UMAP <-
      object@dim_red@UMAP[so@meta.data$barcodes, ] %>%
      select(umap1,umap2)
    names(UMAP) <- c("UMAP_1","UMAP_2")
    write.csv(UMAP, "UMAP.csv")
  }



  code <- paste0(path_to_python," ", path_to_pythonscript, " ",
                 "-I " , "Export.h5ad" ,
                 " -U ", "UMAP.csv",
                 " -F ", folder)

  system(code)
  }else{

    #Export additional files required
    if(Velocity_from=="spatial"){
      rownames(object@coordinates) <- object@coordinates$barcodes
      UMAP <-
        object@coordinates[so@meta.data$barcodes, ] %>%
        select(x,y)
      names(UMAP) <- c("UMAP_1","UMAP_2")
      write.csv(UMAP, "UMAP.csv")
    }else{
      rownames(object@dim_red@UMAP) <- object@dim_red@UMAP$barcodes
      UMAP <-
        object@dim_red@UMAP[so@meta.data$barcodes, ] %>%
        select(umap1,umap2)
      names(UMAP) <- c("UMAP_1","UMAP_2")
      write.csv(UMAP, "UMAP.csv")
    }



    code <- paste0(path_to_python," ", path_to_pythonscript, " ",
                   "-I " , h5ad ,
                   " -U ", "UMAP.csv",
                   " -F ", folder)

    system(code)

}

  setwd(folder)
  velo=read.csv("rank_dynamical_genes.csv")
  dynamic=read.csv("rank_velocity_genes.csv")
  obs=read.csv("observations_metadata.csv")

  #implement obs into object
  object@fdata <- object@fdata %>% left_join(.,obs, by="barcodes")
  names(object@fdata)[2] <- "sample"


  return(object)

}






