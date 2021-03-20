#' @title Find marker genes
#'
#' @description Finds the differentially expressed genes across a set of subgroups.
#'
#' @inherit check_sample params
#' @inherit across params
#' @inherit check_method params
#' @param p_val_adj The maximum adjusted p-value allowed in the output.
#'
#' @return A data.frame containing the variables \emph{'p_val', 'avg_logFC', 'pct.1', 'pct.2', 'p_val_adj', 'cluster', 'gene'}.
#' @export
#'

findDE <- function(object,
                   of_sample = "",
                   across,
                   across_subset = NULL,
                   method_de = "wilcox",
                   p_val_adj = 0.05){

  # 1. Control --------------------------------------------------------------
  check_object(object)
  check_method(method_de = method_de)
  confuns::is_value(x = across, mode = "character", "across")

  of_sample <- check_sample(object, of_sample)

  across <- check_features(object = object,
                           features = across,
                           valid_classes = c("factor", "character"))

  # -----

  # 2. Data extraction ------------------------------------------------------

  fdata <- getFeatureData(object, of_sample)

  groups <- dplyr::pull(.data = fdata, var = {{across}})

  if(!base::is.factor(groups)){

    groups <- base::as.factor(groups)

  }

  levels_groups <- base::levels(groups)

  if(!is.null(across_subset)){

    confuns::is_vec(across_subset, "character", ref = "across_subset")

    ref.against <-
      glue::glue("variable '{across}' of feature data in the specified spata-object") %>%
      base::as.character()

    across_subset <-
      confuns::check_vector(
        input = across_subset,
        against = levels_groups,
        verbose = TRUE,
        ref.input = "input 'across_subset'",
        ref.against = ref.against
      )

    # update feature data
    fdata <- dplyr::filter(fdata, !!rlang::sym(across) %in% {{across_subset}})
    groups <- dplyr::pull(.data = fdata, var = {{ across }})

    if(!base::is.factor(groups)){

      groups <- base::as.factor(groups)

    }

  }

  num_groups <- dplyr::n_distinct(groups)

  if(num_groups >= 20){

    base::stop(glue::glue("The number of different groups is to high for DE-analysis. Is currently {num_groups}. Must be lower than 20. "))

  } else if(!num_groups > 1){

    base::stop(glue::glue("There is only one unique group in the object's '{across}'-variable. findDE() needs a minimum of two different groups."))

  } else {

    base::message(glue::glue("Number of groups/clusters: {num_groups}"))

  }

  #


  # -----

  # 3. Perform DE according to specified method -----------------------------

  barcodes <- dplyr::pull(fdata, barcodes)

  seurat <- Seurat::CreateSeuratObject(object@data@counts[, barcodes])
  seurat@assays$RNA@scale.data <- getExpressionMatrix(object, of_sample)
  seurat@meta.data$orig.ident <- groups
  seurat@active.ident <- seurat@meta.data[,"orig.ident"]
  names(seurat@active.ident) <- base::rownames(seurat@meta.data)

  de <- Seurat::FindAllMarkers(seurat, test.use = method_de)
  base::rm(seurat)

  de <-
    dplyr::filter(de, p_val_adj < {{p_val_adj}}) %>%
    dplyr::select(-pct.1, -pct.2)

  if("avg_log2FC" %in% base::colnames(de)){

    de <- dplyr::rename(de, avg_logFC = avg_log2FC)

  }

  # -----

  base::return(de)

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

  if(base::is.factor(de_df$cluster)){

    de_df$cluster <- S4Vectors::unfactor(de_df$cluster)

  }

  if(!is.null(across_subset)){

    confuns::is_vec(across_subset, "character", ref = "across_subset")

    across_subset <-
      confuns::check_vector(
        input = across_subset,
        against = base::unique(de_df$cluster),
        verbose = TRUE,
        ref.input = " input 'across_subset'",
        ref.against = "variable 'clusters' of input 'de_df'"
      )

  } else if(base::is.null(across_subset)){

    across_subset <- base::unique(de_df$cluster)

  }

  # -----


  # 2. Pipeline -------------------------------------------------------------

  res_df <-
    dplyr::ungroup(de_df) %>%
    dplyr::filter(cluster %in% {{ across_subset }} &
                    !avg_logFC %in% c(Inf, -Inf) &
                    !avg_logFC < 0) %>%
    dplyr::group_by(cluster) %>%
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
