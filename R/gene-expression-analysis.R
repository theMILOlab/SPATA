#' @title Find differently expressed genes
#'
#' @description This function makes use of \code{Seurat::FindAllMarkers()} to compute
#' the differently expressed genes across the groups denoted in the argument \code{across}.
#' See details for more.
#'
#' @inherit check_sample params
#' @inherit across params
#' @inherit check_method params
#' @param ... Additional arguments given to \code{Seurat::FindAllMarkers()}
#'
#' @details If \code{across} and/or \code{method_de} are vectors instead of single
#' values \code{findeDeGenes()} iterates over all combinations in a for-loop and
#' stores the results in the respective slots. (e.g.: If \code{across} = \emph{'seurat_clusters'}
#' and \code{method_de} = \emph{c('wilcox', 'bimod')} the function computes the differently expressed genes
#' across all groups found in the feature variable \emph{seurat_clusters} according to method \emph{wilcox} and
#' stores the results in the respective slot. Then it does the same according to method \emph{bimod}.)
#'
#' The results are obtainable via \code{getDeResults()} and \code{getDeGenes()}.
#'
#' @return A spata-object containing the results in slot @@dea.
#' @export

findDeGenes <- function(object,
                        of_sample = "",
                        across,
                        method_de = "wilcox",
                        verbose = TRUE,
                        ...){

  # 1. Control --------------------------------------------------------------

  # lazy
  check_object(object)

  purrr::walk(.x = method_de, .f = ~ check_method(method_de = .x))

  valid_across <-
    check_features(object = object, valid_classes = c("character", "factor"), features = across)

  # adjusting
  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  for(across in valid_across){

    for(method in method_de){

      if(base::isTRUE(verbose)){base::message(glue::glue("Calculating differently expressed genes across '{across}' with method '{method}'."))}

      object <-
        base::tryCatch({

          # make sure across-input is valid
          groups <- getFeatureVariables(object = object,
                                        features = across,
                                        of_sample = of_sample,
                                        return = "vector")

          # make sure that across-input is passed as a factor

          if(!base::is.factor(groups)){

            groups <- base::factor(x = groups, levels = base::unique(groups))

          }

          n_groups <- dplyr::n_distinct(groups)

          if(n_groups >= 20){

            base::stop(glue::glue("The number of different groups is to high for DE-analysis. Is currently {n_groups}. Must be lower than 20. "))

          } else if(n_groups < 2){

            base::stop(glue::glue("There is only one unique group in the object's '{across}'-variable. findDeGenes() needs a minimum of two different groups."))

          } else {

            base::message(glue::glue("Number of groups/clusters: {n_groups}"))

          }

          # -----

          # 2. De analysis ----------------------------------------------------------

          # prepare seurat object
          seurat_object <- Seurat::CreateSeuratObject(counts = getCountMatrix(object, of_sample))

          seurat_object@assays$RNA@scale.data <- getExpressionMatrix(object, of_sample)

          seurat_object@meta.data$orig.ident <- groups

          seurat_object@active.ident <- seurat_object@meta.data[,"orig.ident"]

          base::names(seurat_object@active.ident) <- base::rownames(seurat_object@meta.data)

          # perform analysis and remove seurat object afterwards
          de_results <-
            Seurat::FindAllMarkers(object = seurat_object, test.use = method, ...)

          base::rm(seurat_object)

          # save results in spata object
          object@dea[[of_sample]][[across]][[method]][["data"]] <-
            tibble::remove_rownames(.data = de_results) %>%
            dplyr::rename({{across}} := "cluster")

          object@dea[[of_sample]][[across]][[method_de]][["adjustments"]] <- list(...)

          object

        },

        error = function(error){

          base::message(glue::glue("Skipping de-analysis on across-input '{across}' with method '{method}' as it resulted in the following error message: {error}"))

          base::return(object)

         }
        )

    }

  }

  # -----

  base::return(object)

}




#' @title Postprocess de-results
#'
#' @description Processes the results of \code{findDE()}. See details.
#'
#' @inherit check_de_df params
#' @param n_highest_FC Numeric value. Affects the number of genes that are kept. See details.
#' @param n_lowest_pvalue Numeric value. Affects the number of genes that are kept. See details.
#' @inherit across params
#' @param return Character value. Denotes the output type. One of \emph{'data.frame', 'vector'} or \emph{'list}
#'
#' @details \code{filterDE()} processes the input by grouping the data.frame according to the unique
#' values of the \emph{cluster}-variable such that the following steps are performed for every experimental
#' group. (With "genes" we refer to the rows (observations) of \code{data}.)
#'
#' \enumerate{
#'  \item{Discards genes with \emph{avg_logFC}-values that are either infinite or negative}
#'  \item{Slices the data.frame in order that for every unique cluster of the \emph{cluster}-variable}:
#'  \enumerate{
#'   \item{the n genes with the highest \emph{avg_logFC}-values are kept where n = \code{n_highest_FC}}
#'   \item{the n genes with the lowest \emph{p_val_adj}-values are kept where n = \code{n_lowest_pvalue}}
#'   }
#'  \item{Arranges the genes according to the highest \emph{avg_logFC}-values}
#'  }
#'
#'
#' @return Depends on input of arguemnt \code{return}:
#'
#'  \itemize{
#'    \item{ \code{return} = \emph{'data.frame'}: The filtered data.frame of \code{de_df} with all it's variables.}
#'    \item{ \code{return} = \emph{'vector'}: A named vector of all genes that remain. Named by the experimental
#'    group in which they were differentially expressed.}
#'    \item{ \code{return} = \emph{'list}: A list named according to the experimental groups. Every slot of that list is
#'    a character vector containing the differentially expressed genes of the respective experimental group.}
#'   }
#'
#' @export

filterDE <- function(de_df,
                     n_highest_FC = 100,
                     n_lowest_pvalue = 100,
                     across_subset = NULL,
                     return = "data.frame"){

  # 1. Control --------------------------------------------------------------

  confuns::is_value(n_highest_FC, "numeric", "n_highest_FC")
  confuns::is_value(n_lowest_pvalue, "numeric", "n_lowest_pvalue")
  confuns::is_value(return, "character", "return")

  confuns::check_one_of(input = return,
                        against = c("data.frame", "vector", "list"),
                        ref.input = "argument 'return'")
  check_de_df(de_df)

  across <-
    dplyr::select(de_df, -dplyr::all_of(x = de_df_columns)) %>%
    base::colnames()

  # -----

  # 2. Pipeline -------------------------------------------------------------

  res_df <-
    dplyr::ungroup(de_df) %>%
    confuns::check_across_subset(df = ., across = across, across.subset = across_subset) %>%
    dplyr::filter(!avg_logFC %in% c(Inf, -Inf) &
                  !avg_logFC < 0) %>%
    dplyr::group_by(!!rlang::sym(across)) %>%
    dplyr::slice_max(avg_logFC, n = n_highest_FC, with_ties = FALSE) %>%
    dplyr::slice_min(p_val_adj, n = n_lowest_pvalue, with_ties = FALSE) %>%
    dplyr::arrange(dplyr::desc(avg_logFC), .by_group = TRUE) %>%
    dplyr::ungroup()

  # -----

  if(return == "vector"){

    dplyr::pull(res_df, gene) %>%
      magrittr::set_names(value = dplyr::pull(res_df, cluster)) %>%
      base::return()

  } else if(return == "data.frame") {

    base::return(res_df)

  } else if(return == "list"){

    purrr::map(.x = across_subset, .f = function(i){

      dplyr::filter(.data = res_df, cluster == {{i}}) %>%
        dplyr::pull(gene)

    }) %>%
      magrittr::set_names(value = across_subset) %>%
      base::return()

  }

}
