#' @title Find marker genes
#'
#' @description Finds the differentially expressed genes of a set of subgroups.
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
                   across_subset,
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

  groups <- dplyr::pull(.data = fdata, var = {{ across }})

  if(!base::is.factor(groups)){

    groups <- base::as.factor(groups)

  }

  levels_groups <- base::levels(groups)


  if(!is.null(across_subset)){

    confuns::is_vec(across_subset, "character", ref = "across_subset")

    across_subset <-
      confuns::check_vector(
        input = across_subset,
        against = levels_groups,
        verbose = TRUE,
        ref.input = " input 'across_subset'",
        ref.against = "variable '{across}' of feature data in the specified spata-object."
      )

    # update feature data
    fdata <- dplyr::filter(fdata, !!rlang::sym(across) %in% {{ across_subset }})
    groups <- dplyr::pull(.data = fdata, var = {{ across }})

    if(!base::is.factor(groups)){

      groups <- base::as.factor(groups)

    }

  }

  num_groups <- base::length(base::levels(groups))

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

  # -----

  base::return(de)

}



#' @title Filter differential genes
#'
#' @description Post-process the results of \code{findDE()} according to your needs.
#'
#' @param data A data.frame containing the numeric variables \emph{p_val, avg_logFC, p_val_adj} and
#' the character variables \emph{cluster, gene}.
#' @param max_FC Numeric value. Denotes the number of genes to keep per \emph{cluster}-group. See details.
#' @param min_pvalue Numeric value. Affects the number of genes kept per \emph{cluster}-group. See details.
#' @inherit across params
#' @param genes_only Logical. If set to TRUE the resulting gene-variable is returned as a character vector.
#'
#' @return A filtered data.frame or a character vector of gene names.
#' @export
#'

filterDE <- function(data,
                     max_FC = 100,
                     min_pvalue = 100,
                     across_subset = NULL,
                     genes_only = FALSE){


  # 1. Control --------------------------------------------------------------

  confuns::check_data_frame(
    df = data,
    var.class = list(
      p_val = "numeric",
      avg_logFC = "numeric",
      p_val_adj = "numeric",
      cluster = c("character", "factor"),
      gene = "character"
    ),
    ref = "data"
  )

  if(base::is.factor(data$cluster)){

    data$cluster <- S4Vectors::unfactor(data$cluster)

  }

  if(!is.null(across_subset)){

    confuns::is_vec(across_subset, "character", ref = "across_subset")

    across_subset <-
      confuns::check_vector(
        input = across_subset,
        against = base::unique(data$cluster),
        verbose = TRUE,
        ref.input = " input 'across_subset'",
        ref.against = "variable 'clusters' of input 'data'"
      )

  } else if(base::is.null(across_subset)){

    across_subset <- base::unique(data$cluster)

  }

  # -----


  # 2. Pipeline -------------------------------------------------------------

  res_df <-
    dplyr::ungroup(data) %>%
    dplyr::filter(cluster %in% {{ across_subset }} &
                    !avg_logFC %in% c(Inf, -Inf) &
                    !avg_logFC < 0) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(avg_logFC, n = max_FC, with_ties = FALSE) %>%
    dplyr::slice_min(p_val_adj, n = min_pvalue, with_ties = FALSE) %>%
    dplyr::arrange(dplyr::desc(avg_logFC), .by_group = TRUE) %>%
    dplyr::ungroup()

  # -----

  if(base::isTRUE(genes_only)){

    dplyr::pull(res_df, gene) %>%
      base::return()

  } else {

    base::return(res_df)

  }


}




