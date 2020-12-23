

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
#' @return
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






