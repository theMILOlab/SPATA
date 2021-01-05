

# Slot: coordinates -------------------------------------------------------

setCoordsDf <- function(object, coords_df, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)

  confuns::check_data_frame(
    df = coords_df,
    var.class = list("barcodes" = "character",
                     "x" = c("integer", "double", "numeric"),
                     "y" = c("integer", "double", "numeric")),
    ref = "coords_df"
  )

  coords_df <- dplyr::mutate(.data = coords_df, sample = {{of_sample}})

  object@coordinates[[of_sample]] <- coords_df

  base::return(object)

}


# Slot: fdata -------------------------------------------------------------

setFeatureDf <- function(object, feature_df, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)

  confuns::check_data_frame(
    df = feature_df,
    var.class = list("barcodes" = "character"),
    ref = "feature_df"
  )

  feature_df <-
    dplyr::mutate(.data = feature_df, sample = {{of_sample}}) %>%
    dplyr::select(barcodes, sample, dplyr::everything())

  object@fdata[[of_sample]] <- feature_df

  base::return(object)

}




# Slot: gdata -------------------------------------------------------------

#' @title Add gene meta data to the object
#'
#' @description Safely adds the output of \code{computeGeneMetaData2()}
#' to the spata-object.
#'
#' @inherit check_sample params
#' @param meta_data_list Output list of \code{computeGeneMetaData2()}. An additional
#' slot named \emph{mtr_name} needs to be added manually.
#'
#' @return An updated spata-object.
#' @export
#'

addGeneMetaData <- function(object, of_sample = "", meta_data_list){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  mtr_name <- meta_data_list$mtr_name

  object@gdata[[of_sample]][[mtr_name]] <- meta_data_list

  base::return(object)

}


# Slot: data --------------------------------------------------------------

setCountMatrix <- function(object, count_mtr, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  object@data[[of_sample]][["counts"]] <- count_mtr

  base::return(object)

}

setDenoisedMatrix <- function(object, denoised_mtr, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  object@data[[of_sample]][["denoised"]] <- denoised_mtr

  base::return(object)

}

setScaledMatrix <- function(object, scaled_mtr, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  object@data[[of_sample]][["scaled"]] <- scaled_mtr

  base::return(object)

}



#' @title Add an expression matrix
#'
#' @description Adds an expression matrix to the object's data slot and
#' makes it available for all SPATA-intern function. Use \code{setActiveExpressionMatrix()}
#' to denote it as the default to use.
#'
#' @inherit check_sample params
#' @param expr_mtr A matrix in which the rownames correspond to the gene names and the
#' column names correspond to the barcode-spots.
#' @param mtr_name A character value that denotes the name of the exprssion matrix with
#' which one can refer to it in subsequent functions.
#'
#' @return An updated spata-object.
#' @export

addExpressionMatrix <- function(object, expr_mtr, mtr_name, of_sample = ""){

  check_object(object)

  confuns::is_value(x = mtr_name, mode = "character")

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  object@data[[of_sample]][[mtr_name]] <- expr_mtr

  base::return(object)

}



# Slot: dim_red -----------------------------------------------------------

setPcaDf <- function(object, pca_df, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::check_data_frame(
    df = pca_df,
    var.class = list("barcodes" = "character",
                     "PC1" = c("integer", "double", "numeric"),
                     "PC2" = c("integer", "double", "numeric")),
    ref = "pca_df"
  )

  pca_df <-
    dplyr::mutate(.data = pca_df, sample = {{of_sample}}) %>%
    dplyr::select(barcodes, sample, dplyr::everything())

  object@dim_red[[of_sample]][["pca"]] <- pca_df

  base::return(object)

}


setTsneDf <- function(object, tsne_df, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::check_data_frame(
    df = tsne_df,
    var.class = list("barcodes" = "character",
                     "tsne1" = c("integer", "double", "numeric"),
                     "tsne2" = c("integer", "double", "numeric")),
    ref = "tsne_df"
  )

  tsne_df <-
    dplyr::mutate(.data = tsne_df, sample = {{of_sample}}) %>%
    dplyr::select(barcodes, sample, dplyr::everything())

  object@dim_red[[of_sample]][["tsne"]] <- tsne_df

  base::return(object)

}


setUmapDf <- function(object, umap_df, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::check_data_frame(
    df = umap_df,
    var.class = list("barcodes" = "character",
                     "umap1" = c("integer", "double", "numeric"),
                     "umap2" = c("integer", "double", "numeric")),
    ref = "umap_df"
  )

  umap_df <-
    dplyr::mutate(.data = umap_df, sample = {{of_sample}}) %>%
    dplyr::select(barcodes, sample, dplyr::everything())

  object@dim_red[[of_sample]][["umap"]] <- umap_df

  base::return(object)

}



# Slot: image -------------------------------------------------------------

setImage <- function(object, image, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  object@images[[of_sample]] <- image

  base::return(object)

}


# Slot: information -------------------------------------------------------

#' @title Denote the default expression matrix
#'
#' @inherit check_object params
#' @param name Character value. The name of the matrix that is to be set as
#' the active expression matrix.
#'
#' @return An updated spta-object.
#' @export

setActiveExpressionMatrix <- function(object, of_sample = "",  mtr_name){

  check_object(object)
  confuns::is_value(x = mtr_name, mode = "character")

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  # check if 'name' is slot in @data
  mtr_names <- getExpressionMatrixNames(object = object, of_sample = of_sample)

  confuns::check_one_of(input = mtr_name,
                        against = mtr_names[mtr_names != "counts"],
                        ref.input = "input for argument 'mtr_name'")

  base::message(glue::glue("Active expression matrix set to '{mtr_name}'."))

  # set name
  object@information$active_mtr[[of_sample]] <- mtr_name

  base::return(object)

}


#' Title
#'
#' @inherit check_object params
#' @param assessment_list Named list with slots \code{$df} and \code{$set_up}.
#'
#' @return A spata-object.
#' @export

setAutoencoderAssessment <- function(object, assessment_list, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::check_data_frame(
    df = assessment_list$df,
    var.class = list("activation" = c("character", "factor"),
                     "bottleneck" = c("character", "factor"),
                     "total_var" = c("numeric", "integer", "double")),
    ref = "assessment_list$df"
  )

  object@information$autoencoder[[of_sample]][["assessment"]] <- assessment_list

  base::return(object)

}





# Slot: spatial  ----------------------------------------------------------



#' Title
#'
#' @param object
#' @param of_sample
#' @param hotspot_list
#'
#' @return
#' @export
#'
#' @examples
setHotspotList <- function(object, of_sample = "", hotspot_list){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  object@spatial[[of_sample]][["hotspots"]] <- hotspot_list

  base::return(object)

}


#' Title
#'
#' @inherit check_sample params
#' @param assessment_list
#'
#' @return
#' @export

setSpCorResults <- function(object,
                                         sp_cor_list,
                                         of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  object@spatial[[of_sample]][["correlation"]] <- sp_cor_list

  base::return(object)

}


#' Title
#'
#' @param object
#' @param cluster_list
#' @param of_sample
#'
#' @return
#' @export
#'
#' @examples
addSpCorCluster <- function(object,
                                         cluster_list,
                                         of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  method <- cluster_list$method

  sp_cor <- getSpCorResults(object, of_sample = of_sample)

  if(method %in% base::names(sp_cor$clusters)){

    confuns::give_feedback(
      msg = glue::glue("Overwriting preexisting results of method '{method}'."),
      verbose = verbose
    )

  }

  sp_cor[["cluster"]][[method]] <- cluster_list

  object <- setSpCorResults(object = object,
                                         sp_cor_list = sp_cor,
                                         of_sample = of_sample)

  base::return(object)

}


# Slot: used_genesets -----------------------------------------------------

setGeneSetDf <- function(object, gene_set_df){

  check_object(object)

  confuns::check_data_frame(
    df = gene_set_df,
    var.class = list("ont" = "character", "gene" = "character"),
    ref = "gene_set_df"
  )

  object@used_genesets <- gene_set_df

  base::return(object)

}






