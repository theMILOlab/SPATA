#' Obtain barcodes of a sample
#'
#' @inherit check_sample params
#'
#' @return All barcodes of the specified sample(s) as a character vector.
#' @export

getBarcodes <- function(object, of_sample = "all"){

  coords_df <- getCoordinates(object = object, of_sample = of_sample)

  return(dplyr::pull(coords_df, barcodes))

}

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
    object@coordinates %>%
    dplyr::filter(sample %in% {{of_sample}})

  # -----

  base::return(coords_df)


}

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

  coords <-
    coordinates(object = object, of_sample = of_sample) %>%
    dplyr::filter(barcodes %in% bc_segm) %>%
    dplyr::mutate(segment = {{of_segment}})

  # -----

  base::return(coords)

}



# DE-analysis -------------------------------------------------------------

#' @title Obtain info on de-analysis storage
#'
#' @inherit check_object params
#'
#' @return A summarizing list.
#' @export

getAvailableDeResults <- function(object){

  check_object(object)

  purrr::map(.x = object@dea, .f = function(sample){

    purrr::map(.x = sample, .f = ~ base::names(.x))

  })

}

#' @title Obtain de-analysis results
#'
#' @inherit check_sample params
#' @inherit across params
#' @inherit check_method params
#' @param p_val_adj Numeric value. Denotes the max. adjusted p-value a gene can have
#' in order to be kept in the results.
#' @param n_threshold Numeric value or vector of length two. See details for more.
#'
#' @details De-analysis results are stored in a data.frame in which each gene that turned out
#' to be a marker gene for a specific group is evaluated by it's p-value, p-adjusted-value,
#' logfold-change etc.
#'
#' If the arguments \code{p_val_adj} and \code{n_threshold} are set to NULL the results are returned
#' as they are which however can result in several hundred genes per group. The results can be additionally
#' subsetted. \code{getDeResults()} does that by
#'
#' \enumerate{
#'  \item{discarding genes with \emph{avg_logFC}-values that are either infinite or negative}
#'  \item{discarding genes with a p-adjusted-value higher than the threshold denoted in \code{p_val_adj}}
#'  \item{slicing the data.frame in order that for every unique group of the feature denoted in \code{across}}:
#'  \enumerate{
#'   \item{the n genes with the highest \emph{avg_logFC}-values are kept where n = \code{n_threshold[1]}}
#'   \item{the n genes with the lowest \emph{p_val_adj}-values are kept where n = \code{n_threshold[2]}}
#'   }
#'  \item{arranging the genes according to the highest \emph{avg_logFC}-values}
#'  }
#'
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

getDeResults <- function(object,
                         of_sample = "",
                         across,
                         method_de,
                         p_val_adj = 0.5,
                         n_threshold = c(50, 50)){

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

  de_results <- de_result_list[["data"]]

  # filter according to p val adjusted limit
  if(!base::is.null(p_val_adj)){

    confuns::is_value(x = p_val_adj, mode = "numeric")

    de_results <- dplyr::filter(.data = de_results, p_val_adj <= {{p_val_adj}})

  }

  # filter according to n threshold
  if(!base::is.null(n_threshold)){

    confuns::is_vec(x = n_threshold, mode = "numeric", min.length = 1, max.length = 2)

    if(base::length(n_threshold) == 1){

      de_results <- filterDE(de_df = de_results,
                             n_highest_FC = n_threshold,
                             n_lowest_pvalue = n_threshold)

    } else if(base::length(n_threshold) == 2){

      de_results <- filterDE(de_df = de_results,
                             n_highest_FC = n_threshold[1],
                             n_lowest_pvalue = n_threshold[2])

    }

  }


  # 3. Return ---------------------------------------------------------------

  base::return(de_results)

}

#' @rdname getDeResults
#' @export
getDeGenes <- function(object,
                       of_sample = "",
                       across,
                       method_de = "wilcox",
                       p_val_adj = NULL,
                       n_threshold = NULL){

  de_results <- getDeResults(object = object,
                             across = across,
                             method_de = method_de,
                             p_val_adj = p_val_adj,
                             n_threshold = n_threshold)

  # extract gene names and name according to cluster belonging
  gene_names <-
    dplyr::pull(de_results, var = "gene") %>%
    purrr::set_names(nm = dplyr::pull(de_results, var = {{across}}))

  base::return(gene_names)

}


#' @title Obtain count and expression matrix
#'
#' @inherit check_sample params
#'
#' @return The expression/count matrix of the specified object and sample(s).
#' @export

getExpressionMatrix <- function(object,
                                of_sample = ""){

  # lazy control
  check_object(object)

  # adjusitng control
  of_sample <- check_sample(object = object, of_sample = of_sample)

  bc_in_sample <-
    object@coordinates %>%
    dplyr::filter(sample %in% {{of_sample}}) %>%
    dplyr::pull(barcodes)

  expr_mtr <- object@data$scaled[,bc_in_sample]

  return(base::as.matrix(expr_mtr))

}

#' @rdname getExpressionMatrix
#' @export
getCountMatrix <- function(object,
                           of_sample = ""){

  # lazy control
  check_object(object)

  # adjusting control
  of_sample <- check_sample(object = object, of_sample = of_sample)

  bc_in_sample <-
    object@coordinates %>%
    dplyr::filter(sample %in% of_sample) %>%
    dplyr::pull(barcodes)

  count_mtr <- object@data$counts[,bc_in_sample]

  return(count_mtr)

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

  getCoordinates(object, of_sample)[,c("barcodes", "sample")]

}


# Dimensional reduction ---------------------------------------------------

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
#'   \item{ \code{getPcaDf()}: \emph{PC_1, PC_2, PC_3, ...PC_n}}
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
    object@dim_red[[method_dr]] %>%
    dplyr::filter(sample %in% of_sample)

  # -----

  if(base::nrow(dim_red_df) == 0){

    base::stop("There seems to be no data for method: ", method_dr)

  }

  base::return(dim_red_df)

}

getDimRedData <- getDimRedDf

#' @rdname getDimRedData
#' @export
getPcaDf <- function(object,
                       of_sample = ""){

  getDimRedDf(object = object,
                of_sample = of_sample,
                method_dr = "pca")

}

#' @rdname getDimRedDf
#' @export
getUmapDf <- function(object,
                        of_sample = ""){

  getDimRedDf(object = object,
                of_sample = of_sample,
                method_dr = "uamp")

}

getUmapData <- getUmapDf

#' @rdname getDimRedDf
#' @export
getTsneDf <- function(object,
                        of_sample = ""){

  getDimRedDf(object = object,
                of_sample = of_sample,
                method_dr = "tsne")

}

getTsneData <- getTsneDf

# -----


# Genes and gene set related ----------------------------------------------

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
  stopifnot(base::is.character(index) | base::is.null(index))

  if(!base::is.character(of_class)){

    stop("Please specify 'of_class' as a character vector.")

  }

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
                     of_gene_sets = "all",
                     in_sample = "",
                     simplify = TRUE){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)

  if(!is.character(of_gene_sets) | base::length(of_gene_sets) == 0){

    stop("Argument 'of_gene_sets' is empty or invalid. Has to be a character vector of length one or more.")

  }

  # adjusting check
  in_sample <- check_sample(object = object, of_sample = in_sample)

  # -----


  # 2. Main part ------------------------------------------------------------

  rna_assay <- getExpressionMatrix(object = object, of_sample = in_sample)

  # -----

  # 2.2 Return all existing genes if desired ----------

  if(base::all(of_gene_sets == "all")){

    base::return(base::rownames(rna_assay))

  }

  # -----

  # 2.3 Return a subset of genes ----------
  if(!base::all(of_gene_sets == "all")){

    gene_sets_df <- object@used_genesets
    of_gene_sets <- check_gene_sets(object, of_gene_sets)

    genes_list <-
      base::lapply(X = of_gene_sets,
                   FUN = function(i){

                     genes <-
                       dplyr::filter(gene_sets_df, ont == i) %>%
                       dplyr::pull(gene)

                     genes_in_sample <-
                       genes[genes %in% base::rownames(rna_assay)]

                       return(genes_in_sample)

                     })

    base::names(genes_list) <- of_gene_sets

    if(base::isTRUE(simplify)){

      genes_list <-
        genes_list %>%
        base::unname() %>%
        base::unlist() %>%
        base::unique()

    }

    base::return(genes_list)

  }

  # -----

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


# Feature related ---------------------------------------------------------

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

getFeatureNames <- function(object, of_class = NULL){

  check_object(object)
  if(!base::is.null(of_class)){confuns::is_vec(of_class, "character", "of_class")}

  feature_names <- base::colnames(object@fdata)

  classes <- base::sapply(object@fdata[,feature_names], base::class)

  base::names(feature_names) <- classes

  if(!base::is.null(of_class)){
    feature_names <- feature_names[classes %in% of_class]
  }

  if(base::length(getSampleNames(object)) > 1){

    base::return(feature_names[feature_names != c("barcodes")])

  } else {

    base::return(feature_names[!feature_names %in% c("barcodes", "sample")])

  }

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

  object@fdata %>%
    dplyr::filter(sample %in% {{of_sample}})

}

getFeatureData <- getFeatureDf


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

  confuns::is_value(return, "character", "return")
  confuns::check_one_of(input = return,
                        against = c("data.frame", "vector"),
                        ref.input = "return")

  of_sample <- check_sample(object, of_sample)

  # -----

  # 2. Extracting -----------------------------------------------------------


  if(base::length(features) == 1 && return == "vector"){

    res <-
      getFeatureData(object, of_sample) %>%
        dplyr::pull(var = {{features}})

  } else if(return == "data.frame"){

    res <-
      getFeatureData(object, of_sample) %>%
        dplyr::select(barcodes, sample, dplyr::all_of(features))

  } else if(return == "list"){

    res <-
      purrr::map(.x = features,
                 .f = function(f){

                   getFeatureData(object, of_sample) %>%
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

  of_sample <- check_sample(object, of_sample)

  # -----

  # 2. Main part ------------------------------------------------------------

  if(base::length(features) == 1){

    getFeatureData(object, of_sample) %>%
      dplyr::pull(var = {{features}}) %>%
      base::unique() %>%
      base::return()

  } else {

    purrr::map(.x = features,
               .f = function(f){

                 getFeatureData(object, of_sample) %>%
                   dplyr::pull(var = {{f}}) %>%
                   base::unique() %>%
                   base::return()

               }) %>%
      magrittr::set_names(features) %>%
      base::return()
  }


}

# -----


#' @title Obtain sample image
#'
#' @inherit check_sample params
#'
#' @return Object of class\emph{image}.
#' @export

getImage <- function(object, of_sample = ""){

  check_object(object)
  of_sample <- check_sample(object, of_sample, desired_length = 1)

  object@image[[of_sample]]

}

# Segmentation related ----------------------------------------------------

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
    base::lapply(X = of_sample,
                FUN = function(i){

                  segment_names <-
                    getFeatureData(object) %>%
                    dplyr::filter(sample == i) %>%
                    dplyr::pull(segment) %>% base::unique()

                  if(base::length(segment_names) == 1 && segment_names == ""){

                     base::warning(stringr::str_c("There seems to be no segmentation for '", i, "'."))

                     base::return(NULL)

                    }

                  return(segment_names[segment_names != ""])

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


# Trajectory related ------------------------------------------------------

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
    dplyr::group_by(pl_binned) %>%
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
    base::lapply(X = of_sample, FUN = function(i){

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
