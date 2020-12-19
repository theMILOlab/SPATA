#' @title Cluster sample via monocle3
#'
#' @description Assign barcode spots to clusters according to different clustering
#' algorithms.
#'
#' @inherit check_object params
#' @inherit check_monocle_input params details
#' @param prefix Character value. Clustering algorithms often return only numbers as
#' names for the clusters they generate. If you want to these numbers to have a certain
#' prefix (like \emph{'Cluster'}, the default) you can specify it with this argument.
#'
#' @details This functions is a wrapper around all monocle3-cluster algorithms which
#' take several options for dimensional reduction upon which the subsequent clustering bases.
#' It iterates over all specified methods and returns a tidy data.frame in which each row represents
#' one barcode-spot uniquely identified by the variable \emph{barcodes} and in which every other variable
#' about the cluster belonging the specified combination of methods returned. E.g.:
#'
#' A call to `findMonocleClusters()` with
#'
#' \itemize{
#'  \item{\code{preprocess_method} set to \emph{'PCA'} }
#'  \item{\code{reduction_method} set to \emph{c('UMAP', 'PCA')}}
#'  \item{\code{'leiden'}, \code{k} set to \emph{5}}
#'  }
#'
#' will return a data.frame of the following variables:
#'
#' \itemize{
#'  \item{\emph{barcodes}}
#'  \item{\emph{mncl_cluster_UMAP_leiden_k5}}
#'  \item{\emph{mncl_cluster_PCA_leiden_k5}}
#'  }
#'
#' Due to the \emph{barcodes}-variable it can be easily joined to your-spata object via `addFeature()`.
#' and thus be made available for all spata-functions.
#'
#' @return A tidy spata-data.frame
#' @export
#'

findMonocleClusters <- function(object,
                                of_sample = "",
                                preprocess_method = c("PCA", "LSI"),
                                reduction_method = c("UMAP", "tSNE", "PCA", "LSI"),
                                cluster_method = c("leiden", "louvain"),
                                k = 20,
                                num_iter = 5,
                                prefix = "Cluster ",
                                verbose = TRUE){

  check_object(object)

  check_monocle_input(preprocess_method = preprocess_method,
                      reduction_method = reduction_method,
                      cluster_method = cluster_method,
                      k = k,
                      num_iter = num_iter)

  if(base::isTRUE(verbose)){base::message("Creating 'cell_data_set'-object.")}

  expression_matrix <- base::as.matrix(getCountMatrix(object, of_sample = of_sample))

  gene_metadata <- data.frame(gene_short_name = base::rownames(expression_matrix))
  base::rownames(gene_metadata) <- base::rownames(expression_matrix)

  cell_metadata <-
    getFeatureDf(object, of_sample = of_sample) %>%
    tibble::column_to_rownames(var = "barcodes")

  cds <- monocle3::new_cell_data_set(
    expression_data = expression_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata)

  # preprocess
  for(p in base::seq_along(preprocess_method)){

    if(base::isTRUE(verbose)){

      base::message(glue::glue("Preprocessing cells with method {p}/{base::length(preprocess_method)} '{preprocess_method[p]}'"))

    }

    cds <- monocle3::preprocess_cds(cds, method = preprocess_method[p])

  }

  # align

  if(base::length(of_sample) > 1){

  if(base::isTRUE(verbose)){ base::message(glue::glue("Aligning for {base::length(of_sample)} samples belonging"))}

    cds <- monocle3::align_cds(cds = cds, alignment_group = "sample")

  }


  for(p in base::seq_along(preprocess_method)){

    base::message(glue::glue("Using preprocess method '{preprocess_method[p]}':"))

    for(r in base::seq_along(reduction_method)){

      base::message(glue::glue("Reducing dimensions with reduction method {r}/{base::length(reduction_method)}: '{reduction_method[r]}' "))

      if(reduction_method[r] == "LSI" && preprocess_method[p] != "LSI"){

        base::message(glue::glue("Ignoring invalid combination. reduction-method: '{reduction_method[r]}' &  preprocess-method: '{preprocess_method[p]}'"))

      } else if(reduction_method[r] == "PCA" && preprocess_method[p] != "PCA") {

        base::message(glue::glue("Ignoring invalid combination. reduction-method: '{reduction_method[r]}' &  preprocess-method: '{preprocess_method[p]}'"))

      } else {

        cds <- monocle3::reduce_dimension(cds = cds, reduction_method = reduction_method[r], preprocess_method = preprocess_method[p], verbose = FALSE)

      }

    }

  }

  cluster_df <- data.frame(barcodes = getBarcodes(object = object))

  for(r in base::seq_along(reduction_method)){

    if(base::isTRUE(verbose)){

      base::message(glue::glue("Using reduction method {reduction_method[r]}:"))

    }

    for(c in base::seq_along(cluster_method)){

      if(base::isTRUE(verbose)){

        base::message(glue::glue("Clustering barcode-spots with method {c}/{base::length(cluster_method)}: {cluster_method[c]}"))

      }

      cds <- monocle3::cluster_cells(cds = cds,
                                     reduction_method = reduction_method[r],
                                     k = k,
                                     num_iter = num_iter,
                                     cluster_method = cluster_method[c],
                                     verbose = FALSE)

      cluster_name <- stringr::str_c("cluster", cluster_method[c], reduction_method[r],base::paste0("k", k), sep = "_")

      cluster_df <-
        monocle3::clusters(x = cds, reduction_method = reduction_method[r]) %>%
        base::as.data.frame() %>%
        tibble::rownames_to_column(var = "barcodes") %>%
        magrittr::set_colnames(value = c("barcodes", cluster_name)) %>%
        dplyr::left_join(x = cluster_df, y = ., by = "barcodes") %>%
        tibble::as_tibble()

    }

  }

  cluster_df <- purrr::map_df(.x = dplyr::select(cluster_df, -barcodes),
                              .f = function(i){

                                i <- stringr::str_c(prefix, i, sep = "")

                                if(base::is.factor(i)){

                                  S4Vectors::unfactor(i) %>%
                                    base::as.character()

                                } else {

                                  base::as.character(i)

                                }

                              }) %>%
    dplyr::mutate(barcodes = cluster_df$barcodes)

  if(base::isTRUE(verbose)){base::message("Done.")}

  base::return(cluster_df)

}

