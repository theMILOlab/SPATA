


#' @title Merge Spata Objects
#'
#' @description Takes an arbitrary number of spata-objects and merges them into one.
#'
#' @param objects A list of valid spata-objects. All sample names must be unique.
#' @param gsdf_input Determines the final input for slot @@used_genesets:
#'
#' If set to \emph{'merge'} all gene-set data.frames of all objects are joined
#' and unique gene-sets are kept. Gene-sets with the same name but different
#' genes are merged!
#'
#' If a directory is specified the directory is given to \code{loadGSDF()}.
#'
#' If a data.frame is specified that data.frame is used.
#'
#' If set to NULL the standard \code{SPATA::gsdf} is used.
#'
#' @return A merged spata object.
#' @export
#'

mergeSpataObjects <- function(objects, gsdf_input = NULL, verbose = TRUE){

  object_list <-
    purrr::keep(.x = objects,
                .p = ~ methods::is(object = .x, class2 = "spata"))

  base::message(glue::glue("Merging {base::length(object_list)} spata-objects."))

  new_object <- methods::new(Class = "spata")

  new_object@samples <-
    purrr::map(.x = object_list, .f = ~ .x@samples) %>%
    purrr::flatten_chr()

  if(base::length(new_object@samples) != dplyr::n_distinct(new_object@samples)){

    base::stop("Sample names overlap. Provided all samples of all spata-objects need to be unique.")

  }

  new_object@coordinates <-
    purrr::map_df(.x = object_list, .f = ~ .x@coordinates)

  new_object@fdata <-
    purrr::map_df(.x = object_list, .f = ~ .x@fdata)

  new_object@image <-
    purrr::map(.x = object_list, .f = ~ .x@image) %>%
    purrr::flatten() %>%
    purrr::set_names(nm = new_object@samples)

  new_object@trajectories <-
    purrr::map(.x = object_list, .f = ~ .x@trajectories) %>%
    purrr::flatten() %>%
    purrr::set_names(nm = new_object@samples)


  if(base::is.null(gsdf_input)){

    if(base::isTRUE(verbose)){base::message("Using SPATA's default gene set data.frame.")}

    new_object@used_genesets <- gsdf

  } else if(base::all(gsdf_input == "merge")){

    if(base::isTRUE(verbose)){base::message("Merging gene-set data.frames.")}

    new_object@used_genesets <-
      purrr::map_df(.x = object_list, .f = ~ .x@used_genesets) %>%
      dplyr::distinct()

  } else if(base::is.data.frame(gsdf_input)){

    if(base::isTRUE(verbose)){base::message("Using 'gsdf_input' as gene-set data.frame.")}

    new_object@used_genesets <- gsdf_input

  } else if(base::is.character(gsdf_input) & base::length(gsdf_input) == 1){

    print("x")

    new_object@used_genesets <-
      loadGSDF(gene_set_path = gsdf_input,
               verbose = verbose)

  }

# Merging matrices --------------------------------------------------------

  # S4 class
  data <- methods::new(Class = "data_counts")

  # extract counts
  count_matrices <- purrr::map(.x = object_list, .f = ~ .x@data@counts)

  all_genes <-
    purrr::map(.x = count_matrices, .f = ~ base::rownames(x = .x)) %>%
    purrr::flatten_chr() %>%
    base::unique()

  all_barcodes <-
    purrr::map(.x = count_matrices, .f = ~ base::colnames(.x)) %>%
    purrr::flatten_chr() %>%
    base::unique()

  count_matrix <-
    base::matrix(data = 0, nrow = base::length(all_genes), ncol = base::length(all_barcodes))

  base::rownames(count_matrix) <- all_genes
  base::colnames(count_matrix) <- all_barcodes

  for(i in base::seq_along(count_matrices)){

    rows <- base::rownames(count_matrices[[i]])
    cols <- base::colnames(count_matrices[[i]])

    count_matrix[rows, cols] <- base::as.matrix(count_matrices[[i]])

  }

  data@counts <- Matrix::Matrix(data = count_matrix)

  count_matrix <- NULL


  # extract scaled data

  expr_matrices <- purrr::map(.x = object_list, .f = ~ .x@data@norm_exp)

  all_genes <-
    purrr::map(.x = expr_matrices, .f = ~ base::rownames(x = .x)) %>%
    purrr::flatten_chr() %>%
    base::unique()

  all_barcodes <-
    purrr::map(.x = expr_matrices, .f = ~ base::colnames(.x)) %>%
    purrr::flatten_chr() %>%
    base::unique()

  expr_matrix <-
    base::matrix(data = 0, nrow = base::length(all_genes), ncol = base::length(all_barcodes))

  base::rownames(expr_matrix) <- all_genes
  base::colnames(expr_matrix) <- all_barcodes

  for(i in base::seq_along(expr_matrices)){

    rows <- base::rownames(expr_matrices[[i]])
    cols <- base::colnames(expr_matrices[[i]])

    expr_matrix[rows, cols] <- base::as.matrix(expr_matrices[[i]])

  }

  data@norm_exp <- expr_matrix

  expr_matrix <- NULL

  # pass to spata

  new_object@data <- data

# Merging dimensional reduction  ------------------------------------------

  umap_df <-
    purrr::map_df(.x = object_list, .f = ~ .x@dim_red@UMAP)

  tsne_df <-
    purrr::map_df(.x = object_list, .f = ~ .x@dim_red@TSNE)

  new_object@dim_red <-
    methods::new(Class = "dim_red",
                 UMAP = umap_df,
                 TSNE = tsne_df)


# Return merged object ----------------------------------------------------

  base::return(new_object)

}
