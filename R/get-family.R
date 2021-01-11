


# Slot: autoencoder -------------------------------------------------------

#' @title Obtain information about the optimal neural network set up
#'
#' @description Extracts the results from \code{assessAutoencoderOptions()}.
#'
#' @inherit check_object params
#'
#' @return A data.frame containing the total variance measured by \code{irlba::prcomp_irlba()} after each
#' combination of activations/bottlenecks.
#' @export

getAutoencoderAssessment <- function(object, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  assessment <- object@information$autoencoder[[of_sample]]$assessment

  if(base::identical(assessment, list()) | base::is.null(assessment)){

    base::stop("Could not find any information. It seems as if function 'assessAutoencoderOptions()' as not been called yet.")

  }

  base::return(assessment)

}


#' @title Obtain information on neural network
#'
#' @description Returns the argument input that was chosen to construct the
#' neural network that generated the matrix denoted in \code{mtr_name}.
#'
#' @inherit getExpressionMatrix params
#'
#' @return A named list.
#' @export

getAutoencoderSetUp <- function(object, mtr_name, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  nn_set_up <-
    object@information$autoencoder[[of_sample]][["nn_set_ups"]][[mtr_name]]

  if(base::is.null(nn_set_up)){

    base::stop(glue::glue("Could not find any autoencoder information for matrix '{mtr_name}' of sample '{of_sample}'"))

  }

  base::return(nn_set_up)

}

# -----


# Slot: coordinates -------------------------------------------------------

#' @title Obtain spatial coordinates
#'
#' @inherit check_sample params
#' @param of_segment Character value. Specifies the segment of interest.
#'
#' @return A data.frame containing the variables \emph{barcods, sample, x, y}
#' (and \emph{segment} if specified).
#' @export

getCoordsDf <- function(object,
                        of_sample = ""){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)

  # adjusting check
  of_sample <- check_sample(object, of_sample, 1)

  # -----

  # 2. Data wrangling -------------------------------------------------------

  coords_df <-
    object@coordinates[[of_sample]]

  # -----

  base::return(coords_df)


}

#' @rdname getCoordsDf
#' @export
getCoordinates <- getCoordsDf

#' @rdname getCoordsDf
#' @export
getCoordinatesSegment <- function(object,
                                  of_segment,
                                  of_sample = ""){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)

  # adjusting check
  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)
  bc_segm <- check_segment(object, segment_name = of_segment, of_sample = of_sample)

  # -----

  # 2. Data wrangling -------------------------------------------------------

  coords_df <-
    getCoordsDf(object = object, of_sample = of_sample) %>%
    dplyr::filter(barcodes %in% bc_segm) %>%
    dplyr::mutate(segment = {{of_segment}})

  # -----

  base::return(coords_df)

}



# -----
# Slot: dea -------------------------------------------------------------

#' @title Obtain info on de-analysis storage
#'
#' @inherit check_object params
#'
#' @return A summarizing list.
#' @export

getDeOverview <- function(object){

  check_object(object)

  all_results <-
    purrr::map(.x = object@dea, .f = function(sample){

      purrr::map(.x = sample, .f = ~ base::names(.x))

    })

  if(base::length(getSampleNames(object)) == 1){

    final_results <-
      purrr::flatten(.x = all_results)

  } else {

    final_results <- all_results

  }

  base::return(final_results)

}

#' @title Obtain de-analysis results
#'
#' @inherit check_sample params
#' @inherit across params
#' @inherit check_method params
#' @inherit filterDeDf params details
#'
#' @return A data.frame:
#'
#' \itemize{
#'   \item{\emph{gene}} Character. The differentially expressed genes.
#'   \item{\emph{cluster}} Character. The clusters (or experimental groups) across which the analysis was performed.
#'   \item{\emph{avg_logFC}} Numeric. The average log-fold change to which the belonging gene was differentially expressed..
#'   \item{\emph{p_val}} Numeric. The p-values.
#'   \item{\emph{p_val_adj}} Numeric. The adjusted p-values.
#'  }
#'
#' If \code{getDeGenes()} is used the \emph{gene}-variable is returned as a named character vector.
#'
#' @export

getDeResultsDf <- function(object,
                          of_sample = "",
                          across,
                          across_subset = NULL,
                          method_de = "wilcox",
                          max_adj_pval = NULL,
                          n_highest_lfc = NULL,
                          n_lowest_pval = NULL){

  # 1. Control --------------------------------------------------------------

  check_object(object)
  check_method(method_de = method_de)

  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  across <- check_features(object, features = across, valid_classes = c("character", "factor"), max_length = 1)

  # 2. Extract and filter ---------------------------------------------------

  de_result_list <- object@dea[[of_sample]][[across]][[method_de]]

  if(base::is.null(de_result_list)){

    base::stop(glue::glue("No de-analysis results found across '{across}' computed via method '{method_de}'."))

  }

  de_results <- filterDeDf(de_df = de_result_list[["data"]],
                           across_subset = across_subset,
                           max_adj_pval = max_adj_pval,
                           n_highest_lfc = n_highest_lfc,
                           n_lowest_pval = n_lowest_pval,
                           return = "data.frame")

  # 3. Return ---------------------------------------------------------------

  base::return(de_results)

}


#' @rdname getDeResultsDf
#' @export
getDeGenes <- function(object,
                       of_sample = "",
                       across,
                       across_subset = NULL,
                       method_de = "wilcox",
                       max_adj_pval = NULL,
                       n_highest_lfc = 50,
                       n_lowest_pval = 50){

  # 1. Control --------------------------------------------------------------

  check_object(object)
  check_method(method_de = method_de)

  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  across <- check_features(object, features = across, valid_classes = c("character", "factor"), max_length = 1)

  # 2. Extract and filter ---------------------------------------------------

  de_result_list <- object@dea[[of_sample]][[across]][[method_de]]

  if(base::is.null(de_result_list)){

    base::stop(glue::glue("No de-analysis results found across '{across}' computed via method '{method_de}'."))

  }

  de_results <- filterDeDf(de_df = de_result_list[["data"]],
                           across_subset = across_subset,
                           max_adj_pval = max_adj_pval,
                           n_highest_lfc = n_highest_lfc,
                           n_lowest_pval = n_lowest_pval,
                           return = "vector")

  # 3. Return ---------------------------------------------------------------

  base::return(de_results)

}

# -----

# Slot: data --------------------------------------------------------------

#' @title Obtain name of currently active expression matrix
#'
#' @inherit check_object params
#'
#' @return Character value.
#' @export

getActiveMatrixName <- function(object, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  object@information$active_mtr[[of_sample]]

}


#' @title Obtain count and expression matrix
#'
#' @inherit check_sample params
#' @param mtr_name Character value. The name of the expression matrix of interest. If set to NULL
#' the currently active matrix is chosen.
#'
#' @return The active expression or count matrix of the specified object and sample(s).
#' @export

getExpressionMatrix <- function(object,
                                of_sample = "",
                                mtr_name = NULL,
                                verbose = FALSE){

  # lazy control
  check_object(object)

  # adjusting control
  of_sample <- check_sample(object = object, of_sample = of_sample)

  if(base::is.null(mtr_name)){

    active_mtr <- getActiveMatrixName(object, of_sample = of_sample)

    if(base::is.null(active_mtr) || !active_mtr %in% getExpressionMatrixNames(object, of_sample = of_sample)){

      active_mtr <- base::ifelse(test = base::is.null(active_mtr), yes = "NULL", no = active_mtr)

      base::stop(glue::glue("Did not find active matrix '{active_mtr}' in data slot of sample '{of_sample}'. Don't know which matrix to return. Please denote a valid active expression matrix with 'setActiveExpressionMatrix()'."))

    }

  } else {

    if(!mtr_name %in% getExpressionMatrixNames(object, of_sample = of_sample)){

      base::stop(glue::glue("Could not find expression matrix '{mtr_name}' of sample '{of_sample}' in provided object."))

    }

    active_mtr <- mtr_name

  }

  if(base::isTRUE(verbose)){ base::message(glue::glue("Using expression matrix '{active_mtr}'."))}

  expr_mtr <-
    object@data[[of_sample]][[active_mtr]] %>%
    base::as.matrix()

  return(expr_mtr)

}

#' @rdname getExpressionMatrix
#' @export
getCountMatrix <- function(object,
                           of_sample = ""){

  # lazy control
  check_object(object)

  # adjusting control
  of_sample <- check_sample(object = object, of_sample = of_sample)

  count_mtr <- object@data[[of_sample]][["counts"]]

  if(base::is.null(count_mtr)){

    base::stop(glue::glue("Did not find count matrix of sample '{of_sample}' in provided spata-object."))

  }

  return(count_mtr)

}


#' @title Obtain names of stored expression matrices
#'
#' @inherit check_object params
#'
#' @return Character vector.
#' @export

getExpressionMatrixNames <- function(object, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  mtr_names <-
    object@data[[of_sample]] %>% base::names() %>%
    purrr::discard(.p = ~ .x == "counts")

  if(base::is.null(mtr_names) | base::identical(mtr_names, base::character(0))){

    base::stop("Could not find any expression matrices in the provided spata-object.")

  } else {

    base::return(mtr_names)

  }

}



# -----

#' @title Obtain a spata-data.frame
#'
#' @description This function is the most basic start if you want
#' to extract data for your individual analysis.
#'
#' (In order to extract the coordinates as well use \code{getCoordinates()}.)
#'
#' @inherit check_sample params
#'
#' @return A tidy data.frame containing the character variables \emph{barcodes}
#' and \emph{sample}.
#'
#' @seealso joinWith
#'
#' @export
#'

getSpataDf <- function(object, of_sample = ""){

  check_object(object)
  of_sample <- check_sample(object, of_sample)

  getCoordsDf(object, of_sample)[,c("barcodes", "sample")]

}


# Slot: dim_red ---------------------------------------------------

#' @title Obtain dimensional reduction data
#'
#' @inherit check_sample params
#' @inherit check_method params
#'
#' @return A data.frame that contains the unique identifiers
#' (keys): \emph{barcodes, sample} and:.
#'
#'  \itemize{
#'   \item{ \code{getTsneDf()}: \emph{tsne1, tsne2}}
#'   \item{ \code{getUmapDf()}: \emph{umap1, umap2}}
#'   \item{ \code{getPcaDf()}: \emph{PC1, PC2, PC3, ...PCn}}
#'   }
#'

getDimRedDf <- function(object,
                        of_sample = "",
                        method_dr = c("pca", "tsne", "umap")){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_method(method_dr = method_dr)

  # adjusting check
  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  # -----

  # 2. Data extraction ------------------------------------------------------

  dim_red_df <-
    object@dim_red[[of_sample]][[method_dr]]

  # -----

  if(base::is.null(dim_red_df) || base::nrow(dim_red_df) == 0){

    base::stop("There seems to be no data for method: ", method_dr)

  }

  base::return(dim_red_df)

}


#' @rdname getDimRedDf
#' @export
getPcaDf <- function(object,
                     of_sample = "",
                     n_pcs = 30){

  confuns::is_value(x = n_pcs, mode = "numeric")

  pca_df <-
  getDimRedDf(object = object,
              of_sample = of_sample,
              method_dr = "pca")

  subset_pcs <- stringr::str_c("PC", 1:n_pcs, sep = "")

  subsetted_pca_df <-
    dplyr::select(pca_df, barcodes, sample, dplyr::all_of(subset_pcs))

  base::return(subsetted_pca_df)

}

#' @rdname getDimRedDf
#' @export
getPcaMtr <- function(object,
                      of_sample = "",
                      n_pcs = 30){

  confuns::is_value(x = n_pcs, mode = "numeric")

  getPcaDf(object = object, n_pcs = n_pcs) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    dplyr::select_if(.predicate = base::is.numeric) %>%
    base::as.matrix()

}


#' @rdname getDimRedDf
#' @export
getUmapDf <- function(object,
                        of_sample = ""){

  getDimRedDf(object = object,
              of_sample = of_sample,
              method_dr = "umap")

}


#' @rdname getDimRedDf
#' @export
getTsneDf <- function(object,
                        of_sample = ""){

  getDimRedDf(object = object,
              of_sample = of_sample,
              method_dr = "tsne")

}


# -----


# Slot: fdata & samples ---------------------------------------------------

#' @title Obtain feature names
#'
#' @description An easy way to obtain all features of interest along with their
#' class.
#'
#' @param object A valid spata-object.
#' @param of_class Character vector. Specify the classes a feature must be of for
#' it's name to be returned.
#'
#' @return A named character vector of the variables in the feature data slot.
#' @export

getFeatureNames <- function(object, of_class = NULL, of_sample = ""){

  check_object(object)
  confuns::is_vec(x = of_class, mode = "character", skip.allow = TRUE, skip.val = NULL)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  feature_df <- getFeatureDf(object = object, of_sample = of_sample)

  feature_names <- base::colnames(feature_df)

  classes <- base::sapply(feature_df[,feature_names], base::class)

  base::names(feature_names) <- classes

  if(!base::is.null(of_class)){
    feature_names <- feature_names[classes %in% of_class]
  }

  base::return(feature_names[!feature_names %in% c("barcodes", "sample")])

}


#' Obtain feature data
#'
#' @inherit check_sample params
#'
#' @return The feature data data.frame of the specfied object and sample(s).
#' @export

getFeatureDf <- function(object, of_sample = ""){

  check_object(object)
  of_sample <- check_sample(object, of_sample)

  fdata <- object@fdata[[of_sample]]

  if(base::is.null(fdata) | base::nrow(fdata) == 0){

    base::stop(glue::glue("Could not find feature data for sample '{of_sample}'."))

  }

  base::return(fdata)

}


#' @title Obtain a feature variable
#'
#' @description Extracts the specified feature variables from the
#' feature data.
#'
#' @inherit check_sample params
#' @inherit check_features params
#' @param return Character value. One of \emph{'vector', 'data.frame'} or
#' \emph{'list'}. In order to return a vector input of \code{features} must
#' be of length one.
#' @param unique Deprecated.
#'
#' @return A data.frame or a vector.
#' @export

getFeatureVariables <- function(object,
                                features,
                                of_sample = "",
                                return = "data.frame",
                                unique = "deprecated"){

  if(unique != "deprecated"){
    base::warning("Argument 'unique' is deprecated.")
  }

  # 1. Control --------------------------------------------------------------

  check_object(object)
  features <- check_features(object, features)

  confuns::is_value(x = return, mode = "character")
  confuns::check_one_of(input = return,
                        against = c("data.frame", "vector"),
                        ref.input = "return")

  of_sample <- check_sample(object, of_sample)

  # -----

  # 2. Extracting -----------------------------------------------------------


  if(base::length(features) == 1 && return == "vector"){

    res <-
      getFeatureDf(object, of_sample = of_sample) %>%
      dplyr::pull(var = {{features}})

  } else if(return == "data.frame"){

    res <-
      getFeatureDf(object, of_sample = of_sample) %>%
      dplyr::select(barcodes, sample, dplyr::all_of(features))

  } else if(return == "list"){

    res <-
      purrr::map(.x = features,
                 .f = function(f){

                   getFeatureDf(object, of_sample) %>%
                     dplyr::pull(var = {{f}})

                 }) %>%
      magrittr::set_names(value = features)

  }

  base::return(res)

}


#' @title Obtain unique categorical feature values
#'
#' @description Extracts the unique values of discrete features.
#'
#' @inherit check_sample params
#' @inherit check_features params
#'
#' @return A vector or a named list according to the length of \code{features}.
#' @export

getFeatureValues <- function(object, of_sample = "", features){

  # 1. Control --------------------------------------------------------------

  check_object(object)
  features <- check_features(object, features, valid_classes = c("character", "factor"))

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  # -----

  # 2. Main part ------------------------------------------------------------

  if(base::length(features) == 1){

    getFeatureDf(object, of_sample = of_sample) %>%
      dplyr::pull(var = {{features}}) %>%
      base::unique() %>%
      base::return()

  } else {

    purrr::map(.x = features,
               .f = function(f){

                 getFeatureDf(object, of_sample = of_sample) %>%
                   dplyr::pull(var = {{f}}) %>%
                   base::unique() %>%
                   base::return()

               }) %>%
      magrittr::set_names(features) %>%
      base::return()
  }


}


#' @title Obtain segment names
#'
#' @inherit check_sample params
#'
#' @return A list named according to the \code{of_sample} in which each element is
#' a character vector containing the names of segments which were drawn for the
#' specific sample.
#'
#' @export

getSegmentNames <- function(object,
                            of_sample = "",
                            simplify = TRUE){

  # lazy check
  check_object(object)

  # adjusting check
  of_sample <- check_sample(object, of_sample = of_sample)

  # main part
  res_list <-
    purrr::map(.x = of_sample,
               .f = function(i){

                 segment_names <-
                   getFeatureDf(object, of_sample = of_sample) %>%
                   dplyr::pull(segment) %>%
                   base::unique()

                 if(base::length(segment_names) == 1 && segment_names %in% c("none", "")){

                   base::warning(stringr::str_c("There seems to be no segmentation for '", i, "'."))

                   base::return(NULL)

                 } else {

                   return(segment_names[!segment_names %in% c("none", "")])

                 }

               })

  base::names(res_list) <- of_sample

  res_list <- purrr::discard(.x = res_list, .p = base::is.null)

  if(base::isTRUE(simplify)){

    base::unlist(res_list, use.names = FALSE) %>%
      base::return()

  } else {

    base::return(res_list)

  }


}

#' @title Obtain sample names
#'
#' @inherit check_object params
#'
#' @return A character vector.
#'
#' @export

getSampleNames <- function(object){

  check_object(object)

  object@samples

}

#' @rdname getSampleNames
getSamples <- function(object){

  warning("getSamples is deprecated. Use getSampleNames")

  object@samples

}

# -----


# Slot: gdata -------------------------------------------------------------


#' @title Obtain gene meta data
#'
#' @inherit check_sample params
#' @param only_df Logical. If set to TRUE only the data.frame is returned.
#' If set to FALSE (the default) the whole list is returned.
#' @inherit getExpressionMatrix params
#'
#' @return A data.frame from \code{getMetaDataDf()} or a list from \code{getGeneMetaData()}.
#' @export

getGeneMetaData <- function(object, of_sample = "", mtr_name = NULL, only_df = FALSE){

  check_object(object)
  of_sample <- check_sample(object = object, of_sample = of_sample)

  if(base::is.null(mtr_name)){

    mtr_name <- getActiveMatrixName(object, of_sample = of_sample)

  }

  gdata <- object@gdata[[of_sample]][[mtr_name]]

  check_availability(
    test = (base::is.list(gdata) & !base::identical(gdata, list())),
    ref_x = glue::glue("gene meta data for expression matrix '{mtr_name}' of sample '{of_sample}'"),
    ref_fns = "computeGeneMetaData() or addGeneMetaData()"
  )

  if(base::isTRUE(only_df)){

    base::return(gdata$df)

  } else {

    base::return(gdata)

  }

}

#' @rdname getGeneMetaData
#' @export
getGeneMetaDf <- function(object, of_sample = "", mtr_name = NULL){

  getGeneMetaData(object = object, of_sample = of_sample, mtr_name = mtr_name, only_df = TRUE)

}

# -----



# Slot: images ------------------------------------------------------------

#' Title
#'
#' @param object
#' @param of_sample
#'
#' @return
#' @export
#'
#' @examples
getImage <- function(object, of_sample = ""){

  check_object(object)

  of_sample = check_sample(object, of_sample = of_sample, of.length = 1)

  object@images[[of_sample]]

}

# -----
# Slot: information -------------------------------------------------------

#' Obtain barcodes of a sample
#'
#' @inherit check_sample params
#'
#' @return All barcodes of the specified sample(s) as a character vector.
#' @export

getBarcodes <- function(object, of_sample = ""){

  check_object(object)
  of_sample <- check_sample(object = object, of_sample = of_sample)

  object@information$barcodes[[of_sample]]

}

# -----



# Slot: spatial -----------------------------------------------------------

#' @title Obtain distance measurements of spatially correlated genes
#'
#' @inherit check_sample params
#'
#' @return A data.frame or a distance matrix.
#' @export

getGeneDistMtr <- function(object, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  sp_cor <- getSpCorResults(object, of_sample = of_sample)

  base::return(sp_cor$dist_mtr)

}

getGeneDistDf <- function(object, of_sample = ""){

  getGeneDistMtr(object = object, of_sample = of_sample) %>%
    hlpr_dist_mtr_to_df()

}


#' Title
#'
#' @param object
#' @param of_sample
#'
#' @return
#' @export

getPrResults <- function(object, of_sample = "", method_pr = "hotspot"){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  pr_list <-
    object@spatial[[of_sample]][[method_pr]]

  check_availability(
    test = base::is.list(pr_list) & confuns::is_named(pr_list),
    ref_x = "pattern recognition results",
    ref_fns = glue::glue("/ running runPatternRecognition(..., method_pr = {'method_pr'})")
  )

  base::return(pr_list)

}

#' @rdname getPrResults
#' @export
getPrSuggestion <- function(object, of_sample = "", method_pr = "hotspot"){

  pr_list <-
    getPrResults(object = object, of_sample = of_sample, method_pr = method_pr)

  base::return(pr_list$suggestion)

}

#' @rdname getPrResults
#' @export
getPatternNames <- function(object, of_sample = "", method_pr = "hotspot"){

  getPrSuggestion(object, of_sample = of_sample, method_pr = method_pr)$info %>%
    dplyr::pull(var = {{method_pr}}) %>%
    base::levels()

}



#' Title
#'
#' @inherit check_sample params
#' @inherit method_hclust params
#'
#' @return
#' @export

getSpCorCluster <- function(object,
                                         method_hclust = "complete",
                                         of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  sp_cor <-
    getSpCorResults(object, of_sample = of_sample)

  cor_clusters <-
    sp_cor$clusters

  if(base::is.null(cor_clusters) | base::identical(list(), cor_clusters)){

    base::stop("Could not find any correlation clusters. It seems as if function 'clusterSpCorResults()' has not been run yet.")

  } else if(!method_hclust %in% base::names(cor_clusters)) {

    base::stop("Could not find correlation clusters for to method '{method_hclust}'.")

  } else {

    base::return(cor_clusters[[method_hclust]])

  }

}


#' Title
#'
#' @param object
#' @param of_sample
#'
#' @return
#' @export

getSpCorClusterNames <- function(object, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  sp_cor <- getSpCorResults(object, of_sample = of_sample)

  cluster_names <- base::names(sp_cor$clusters)

  if(base::is.null(cluster_names) | base::length(cluster_names) == 0){

    base::stop(glue::glue("Could not find any spatial correlation clusters for sample '{of_sample}'. It seems as if 'clusterSpCorResults()' has not been run yet."))

  } else {

    base::return(cluster_names)

  }

}


#' Title
#'
#' @inherit check_sample params
#'
#' @return
#' @export

getSpCorResults <- function(object, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  corr_assessment <-
    object@spatial[[of_sample]]$correlation

  if(base::is.null(corr_assessment)){

    base::stop(glue::glue("Could not find any correlation assessment for sample '{of_sample}'. It seems as if function 'assessSpCor()' has not been run yet."))

  } else {

    base::return(corr_assessment)

  }

}


# -----

# Slot: trajectories ------------------------------------------------------

#- 'getTrajectoryComment()' is documented in 'S4_generic_functions.R' -#



#' @title Obtain the length of a trajectory
#'
#' @description This function returns the length (the number of bins) of a trajectory
#' depending on the chosen \code{binwidth}.
#'
#' @inherit check_trajectory params
#' @inherit check_trajectory_binwidth params
#'
#' @return Numeric value.
#' @export
#'

getTrajectoryLength <- function(object,
                                trajectory_name,
                                of_sample = "",
                                binwidth = 5){


  # 1. Control --------------------------------------------------------------

  check_object(object)
  check_trajectory(object = object, trajectory_name = trajectory_name, of_sample = of_sample)

  confuns::is_value(x = binwidth, mode = "numeric")

  # -----

  # 2. Extraction -----------------------------------------------------------

  t_object <-
    getTrajectoryObject(object = object,
                        trajectory_name = trajectory_name,
                        of_sample = of_sample)

  t_object@compiled_trajectory_df %>%
    dplyr::mutate(pl_binned = plyr::round_any(x = projection_length, accuracy = binwidth, f = base::floor)) %>%
    dplyr::group_by(pl_binned, trajectory_part) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
    base::nrow()


}


#' @title Obtain trajectory names
#'
#' @inherit check_sample params
#'
#' @return A list named according to the \code{of_sample} in which each element is
#' a character vector containing the names of trajectories which were drawn for the
#' specific sample.
#'
#' @export

getTrajectoryNames <- function(object, of_sample = "all", simplify = TRUE){

  # lazy check
  check_object(object)

  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample)

  # main part
  t_names_list <-
    purrr::map(.x = of_sample, .f = function(i){

      t_names <-
        base::names(object@trajectories[[i]])

      if(base::length(t_names) == 0){

        base::message(stringr::str_c("No trajectories found in sample: ", i, sep = ""))

        base::return(NULL)

      } else {

        base::return(t_names)

      }

    })

  base::names(t_names_list) <- of_sample

  t_names_list <- purrr::discard(.x = t_names_list, .p = is.null)

  if(base::isTRUE(simplify)){

    t_names_list <- base::unlist(t_names_list) %>% base::unname()

  }

  if(!base::length(t_names_list) == 0){

    base::return(t_names_list)

  } else {

    base::return(base::invisible(NULL))

  }


}



#' @title Obtain a summarized trajectory data.frame
#'
#' @description Computes the expression trends of all specified variables
#' along the direction of the spatial trajectory.
#'
#' @inherit check_sample params
#' @inherit check_trajectory params
#' @inherit hlpr_summarize_trajectory_df params
#' @param shift_wider Logical. If set to TRUE the trajectory data.frame is
#' shifted to it's wider format. Formats can be changed via \code{shiftTrajectoryDf()}.
#'
#' @return A summarized trajectory data.frame.
#'
#' @inherit hlpr_summarize_trajectory_df details
#'
#' @export

getTrajectoryDf <- function(object,
                            trajectory_name,
                            of_sample = "",
                            variables,
                            method_gs = "mean",
                            binwidth = 5,
                            normalize = TRUE,
                            shift_wider = FALSE,
                            verbose = TRUE){


  confuns::are_values(c("normalize", "shift_wider", "verbose"), mode = "logical")

  tobj <-
    getTrajectoryObject(object, trajectory_name, of_sample)

  stdf <-
    hlpr_summarize_trajectory_df(object,
                                 ctdf = tobj@compiled_trajectory_df,
                                 binwidth = binwidth,
                                 variables = variables,
                                 method_gs = method_gs,
                                 verbose = verbose,
                                 normalize = normalize)

  if(base::isTRUE(shift_wider)){

    stdf <- shiftTrajectoryDf(stdf = stdf, shift = "wider")

  }

  base::return(stdf)

}


#' @title Obtain trajectory object
#'
#' @inherit check_sample params
#' @inherit check_trajectory params
#'
#' @return An object of class \code{spatialTrajectory}.
#' @export

getTrajectoryObject <- function(object, trajectory_name, of_sample = ""){

  check_trajectory(object = object,
                   trajectory_name = trajectory_name,
                   of_sample = of_sample)

  of_sample <- check_sample(object = object,
                            of_sample = of_sample,
                            desired_length = 1)

  object@trajectories[[of_sample]][[trajectory_name]]

}

# -----



# Slot: used_genesets ----------------------------------------------

#' @title Overview about the current gene sets
#'
#' @param object A valid spata-object.
#'
#' @return A data.frame with two variables \emph{Class} and \emph{Available Gene
#' Sets} indicating the number of different gene sets the classes contain.
#'
#' @export

getGeneSetOverview <- function(object){

  # lazy check
  check_object(object)

  # main part
  gene_sets_df <- dplyr::ungroup(object@used_genesets)

  gene_sets <- object@used_genesets$ont

  if(base::nrow(gene_sets_df) == 0){

    base::message("Gene-set data.frame is empty.")
    base::return(data.frame())

  } else {

    gene_set_classes <- stringr::str_extract(string = gene_sets, pattern = "^.+?(?=_)")

    dplyr::mutate(gene_sets_df, gs_type = gene_set_classes) %>%
      dplyr::select(-gene) %>%
      dplyr::distinct() %>%
      dplyr::pull(gs_type) %>%
      base::table() %>%
      base::as.data.frame() %>%
      magrittr::set_colnames(value = c("Class", "Available Gene Sets"))

  }

}



#' @title Obtain gene set names
#'
#' @inherit check_object params
#' @param of_class A character vector indicating the classes from which to obtain
#' the gene set names. (Which classes exist in the current gene set data.frame can
#' be obtained e.g. with \code{geneSetOverview()}). If set to \emph{"all"} all
#' gene sets are returned.
#' @param index A regular expression according to which the gene set names to be returned
#' are filtered again.
#' @param simplify Logical. If set to TRUE the list to be returned is simplified
#' into a character vector.
#'
#'
#' @return A list named according to the input of argument \code{of_class}. Each element of
#' the returned list is a character vector containing the names of gene sets of the specified classes.
#' The list is coalesced to an unnamed vector if \code{simplify} is set to TRUE.
#'
#' @export

getGeneSets <- function(object, of_class = "all", index = NULL, simplify = TRUE){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)

  confuns::is_vec(x = of_class, mode = "character")
  confuns::is_value(x = index, mode = "character", skip.allow = TRUE, skip.val = NULL)

  # -----

  # 2. Main part ------------------------------------------------------------

  gene_sets_df <- object@used_genesets

  # 2.1 Extract gene sets according to 'of_class' ----------
  if(base::length(of_class) == 1 && of_class == "all"){

    res_list <- base::unique(gene_sets_df$ont)

  } else {

    # get gene sets for all elements of 'of_class' in a list
    res_list <-
      base::lapply(X = of_class, FUN = function(i){

        subset <-
          gene_sets_df$ont %>%
          stringr::str_subset(pattern = stringr::str_c("^", i, sep = "")) %>%
          base::unique()

        if(base::length(subset) == 0){

          base::warning(stringr::str_c("Could not find any gene set of class:", i, sep = " "))

          base::return(NULL)

        } else {

          base::return(subset)

        }

      })

    base::names(res_list) <- of_class

    # discard list elements if 'of_class' element wasn't found
    res_list <-
      purrr::discard(.x = res_list, .p = base::is.null)

  }

  # -----


  # 2.2 Adjust output according to 'index' ----------

  if(base::isTRUE(simplify)){

    res_list <- base::unlist(res_list) %>% base::unname()

  }


  if(!base::is.null(index) && base::is.list(res_list)){

    res_list <-
      base::lapply(X = res_list,
                   FUN = function(i){

                     i[stringr::str_detect(string = i, pattern = index)]

                   })

  } else if(!base::is.null(index) && base::is.character(res_list)){

    res_list <-
      res_list[stringr::str_detect(string = res_list, pattern = index)]

  }

  # -----
  if(base::is.null(res_list)){

    base::stop("Did not find any gene-set.")

  } else {

    base::return(res_list)

  }

}

#' @rdname getGeneSets
#' @export
getGeneSetsInteractive <- function(object){

  check_object(object)

  gene_sets <-
    shiny::runGadget(
      shiny::shinyApp(
        ui = {shiny::fluidPage(

          shiny::fluidRow(

            shiny::HTML("<br><br><br>"),

            shiny::fluidRow(
              shiny::column(width = 6,
                            shiny::tags$h5(shiny::strong("Chosen gene-sets:")),
                            shiny::verbatimTextOutput("display_gene_sets"),
                            shiny::actionButton("return_gene_sets", "Return gene-sets")),
              shiny::column(width = 6,
                            shiny::tags$h5(shiny::strong("Choose gene-sets:")),
                            shiny::uiOutput("select_gene_sets"))
            )

          ),



        )},
        server = function(input, output, session){


          output$select_gene_sets <- shiny::renderUI({

            shinyWidgets::pickerInput("select_gene_sets",
                                      label = NULL ,
                                      choices = getGeneSets(object),
                                      selected = NULL,
                                      options = list(`live-search` = TRUE),
                                      inline = FALSE,
                                      multiple = TRUE)

          })

          output$display_gene_sets <- shiny::renderPrint({

            input$select_gene_sets

          })

          oe <- shiny::observeEvent(input$return_gene_sets, {

            shiny::stopApp(returnValue = input$select_gene_sets)

          })

        }
      )
    )

  base::return(gene_sets)

}


#' @title Obtain gene set data.frame
#'
#' @inherit check_object params
#'
#' @return A data.frame.
#' @export

getGeneSetDf <- function(object){

  check_object(object)

  object@used_genesets

}


#' @title Obtain gene names
#'
#' @inherit check_object params
#' @param of_gene_sets A character vector specifying the gene sets from which to
#' return the gene names.
#' @param of_pattern A character vector specifiying the patterns from which to return
#' the gene names. If denoted as \emph{""} the genes of all patterns are returned.
#' Set \code{simplify} to FALSE in order to return a named list.
#' @param in_sample The sample(s) in which the genes have to be expressed in order
#' to be included.
#' @param simplify Logical. If set to TRUE the list to be returned will be simplified
#' into a character vector.
#'
#' @return A list named according to the input of \code{of_gene_sets} in which each element is
#' a character vector containing the names of genes the specific gene set is
#' composed of. Is coalesced to a vector if \code{simplify} is set to TRUE.
#'
#' @export

getGenes <- function(object,
                     of_gene_sets = NULL,
                     of_hotspots = NULL,
                     of_sample = "",
                     in_sample = "",
                     simplify = TRUE){

  if(!in_sample == ""){warning("change in_sample to of_sample")}

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)

  confuns::are_vectors(c("of_gene_sets", "of_hotspots"), mode = "character",
                       skip.allow = TRUE, skip.val = NULL)


  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample)

  # -----


  # 2. Main part ------------------------------------------------------------

  # -----

  # 2.1 Return all existing genes if desired ----------

  if(!base::is.null(of_gene_sets) && base::all(of_gene_sets == "all")){warning("change of_gene_sets to NULL")}

  if(base::all(base::is.null(of_gene_sets), base::is.null(of_hotspots))){

    expr_mtr <- getExpressionMatrix(object = object, of_sample = of_sample)

    base::return(base::rownames(expr_mtr))

  }

  # -----

  # 2.2 Return a subset of genes ----------
  if(!base::is.null(of_gene_sets)){

    gene_set_df <- getGeneSetDf(object)

    of_gene_sets <- check_gene_sets(object, gene_sets = of_gene_sets)
    expr_mtr <- getExpressionMatrix(object = object, of_sample = of_sample)

    genes_list <-
      purrr::map(.x = of_gene_sets,
                 .f = function(i){

                     genes <-
                       dplyr::filter(gene_set_df, ont == i) %>%
                       dplyr::pull(gene)

                     genes_in_sample <-
                       genes[genes %in% base::rownames(expr_mtr)]

                     return(genes_in_sample)

                   }) %>%
      purrr::set_names(nm = of_gene_sets)

  } else if(!base::is.null(of_hotspots)){

    hotspots <-
      check_pattern(object = object,
                    of_sample = of_sample,
                    patterns = of_hotspots,
                    method_pr = "hotspot")

    hotspot_df <- getPrSuggestion(object,
                                  of_sample = of_sample,
                                  method_pr = "hotspot")$df

    genes_list <-
      purrr::map(.x = hotspots,
                 .f = function(hotspot){

                   dplyr::filter(hotspot_df, hotspot %in% {{hotspots}}) %>%
                     dplyr::pull(var = "genes")


                 }) %>%
      purrr::set_names(nm = hotspots)

  }

  # simplify output if specifed
  if(base::isTRUE(simplify)){

    genes_list <-
      genes_list %>%
      base::unname() %>%
      base::unlist() %>%
      base::unique()

  }

  base::return(genes_list)

  # -----

}

#' @rdname getGenes
#' @export
getSpCorGenes <- function(object,
                          of_sample = "",
                          similar_to = NULL,
                          distinct_to = NULL,
                          top_n = 25,
                          method_hclust = NULL,
                          cluster_names = NULL,
                          simplify = TRUE){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  if(base::is.character(similar_to)){

    dist_df <- getGeneDistDf(object, of_sample = of_sample)

    confuns::is_value(x = top_n, mode = "numeric")

    confuns::check_one_of(
      input = similar_to,
      against = base::unique(dist_df$gene2),
      ref.input = "input for argument 'similar_to'"
    )

    res_genes <-
      dplyr::filter(.data = dist_df, gene2 == {{similar_to}}) %>%
      dplyr::slice_min(order_by = distance, n = top_n) %>%
      dplyr::pull(var = "gene1") %>%
      base::as.character()

  } else if(base::is.character(distinct_to)){

    dist_df <- getGeneDistDf(object, of_sample = of_sample)

    confuns::is_value(x = top_n, mode = "numeric")

    confuns::check_one_of(
      input = distinct_to,
      against = base::unique(dist_df$gene2),
      ref.input = "input for argument 'distinct_to'"
    )

    res_genes <-
      dplyr::filter(.data = dist_df, gene2 == {{distinct_to}}) %>%
      dplyr::slice_max(order_by = distance, n = top_n) %>%
      dplyr::pull(var = "gene1") %>%
      base::as.character()

  } else if(base::is.character(method_hclust)){

    sp_cor_cluster <-
      getSpCorCluster(object = object, method_hclust = method_hclust)

    if(base::is.numeric(cluster_names)){

      cluster_names <- stringr::str_c("cluster", cluster_names, sep = "_")

    } else if(!base::is.character(cluster_names)){

      base::stop("Input for argument 'cluster_names' must be a character or a numeric vector.")

    }

    gene_names_list <- sp_cor_cluster$gene_names_list

    cluster_names <-
      confuns::check_vector(
        input = cluster_names,
        against = base::names(sp_cor_cluster$gene_names_list),
        ref.input = "cluster names",
        ref.against = "valid valid cluster names."
      )

    res_genes <-
      purrr::map(.x = cluster_names, .f = ~ gene_names_list[[.x]]) %>%
      purrr::set_names(nm = cluster_names)

    if(base::isTRUE(simplify) && base::length(res_genes) >= 2){

      res_genes <-
        purrr::imap(.x = res_genes, .f = function(x, cluster_name){

          base::names(x) <-
            base::rep(x = cluster_name, base::length(x))

          base::return(x)

        }) %>%
        purrr::flatten_chr()

    }

  } else {

    base::stop("Either argument 'similar_to' or 'method_hclust' must be specified.")

  }

  base::return(res_genes)

}

#' @rdname getGenes
#' @export
getGenesInteractive <- function(object){

  check_object(object)

  genes <-
    shiny::runGadget(
      shiny::shinyApp(
        ui = {shiny::fluidPage(

          shiny::fluidRow(

            shiny::HTML("<br><br><br>"),

            shiny::fluidRow(
              shiny::column(width = 6,
                            shiny::tags$h5(shiny::strong("Chosen genes:")),
                            shiny::verbatimTextOutput("display_genes"),
                            shiny::actionButton("return_genes", "Return genes")),
              shiny::column(width = 6,
                            shiny::tags$h5(shiny::strong("Choose genes:")),
                            shiny::uiOutput("select_genes"))
            )

          )

        )},
        server = function(input, output, session){

          output$select_genes <- shiny::renderUI({

            shinyWidgets::pickerInput("select_genes",
                                      label = NULL ,
                                      choices = getGenes(object),
                                      selected = NULL,
                                      options = list(`live-search` = TRUE),
                                      inline = FALSE,
                                      multiple = TRUE)

          })

          output$display_genes <- shiny::renderPrint({

            input$select_genes

          })

          oe <- shiny::observeEvent(input$return_genes, {

            shiny::stopApp(returnValue = input$select_genes)

          })

        }
      )
    )

  base::return(genes)

}



# -----





