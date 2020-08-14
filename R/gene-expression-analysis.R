#' @title Find marker genes
#'
#' @description Finds the differentially expressed genes of a set of subgroups.
#'
#' @inherit check_sample params
#' @param feature Character value. A categorical feature of the feature data that
#' separates all barcodes into subgroups for which eventually the marker genes will be determined.
#'
#' (Number of subgruops must be lower than 20.)
#' @inherit check_method params
#' @param p_val_adj The maximum adjusted p-value allowed in the output.
#'
#' @return A DE-object.
#' @export
#'

findDE <- function(object,
                   of_sample,
                   feature,
                   method_de = "wilcox",
                   p_val_adj = 0.05){


  # 1. Control --------------------------------------------------------------
  check_object(object)
  check_method(method_de = method_de)
  confuns::is_value(x = feature, mode = "character", "feature")

  of_sample <- check_sample(object, of_sample)

  feature <- check_features(object = object,
                            features = feature,
                            valid_classes = c("factor", "character"))

  # -----

  # 2. Data extraction ------------------------------------------------------

  fdata <- getFeatureData(object, of_sample)

  barcodes <- dplyr::pull(.data = fdata, var = barcodes)
  groups <-
    dplyr::pull(.data = fdata, var = {{feature}}) %>%
    base::unique()

  if(!base::is.factor(groups)){

    groups <- base::as.factor(groups)

  }

  num_groups <- base::length(base::levels(groups))

  if(num_groups >= 20){

    base::stopp(glue::glue("The number of different groups is to high for DE-analysis. Is currently {num_groups}. Must be lower than 20. "))

  }

  # -----

  # 3. Perform DE according to specified method -----------------------------

  seurat <- Seurat::CreateSeuratObject(object@data@counts[, barcodes])
  seurat@assays$RNA@scale.data <- getExpressionMatrix(object, of_sample)
  seurat@meta.data$orig.ident <- groups
  seurat@active.ident <- seurat@meta.data[,"orig.ident"]
  names(seurat@active.ident) <- base::rownames(seurat@meta.data)

  de <- Seurat::FindAllMarkers(seurat, test.use = method_de)
  base::rm(seurat)

  de <- dplyr::filter(de, p_val_adj < p_val_adj)

  # -----

  base::return(de)

}




