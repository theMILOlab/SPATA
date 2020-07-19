#' @title Join barcodes with additional variables
#'
#' @description Each member of the \code{joinWith()}-family takes a data.frame as
#' input that contains at least the variables \emph{barcodes} and \emph{sample}.
#' (Easily obtained with the \code{spata::coords()}-family.) It joins them over
#'  the \emph{barcodes}-variable with a specified aspect.
#'
#' @param object A valid spata-object.
#' @param coords_df A data.frame that contains variables \emph{barcodes, sample}
#' @param features The features you want to join \code{coords_df} with specified
#' as a character vector.
#' @param genes The genes you want to join \code{coords_df} with specified
#' as a character vector.
#' @param average_genes Logical value. If set to TRUE coords_df will be joined with
#' the mean expression values of all genes specified.
#' @param gene_sets The gene sets you want to join \code{coords_df} with specified
#' as a character vector.
#' @param method_gs The method according to which gene sets will be handled
#' specified as a character of length one. This can be either \emph{mean} or one
#' of \emph{gsva, ssgsea, zscore, or plage}. The latter four will be given to
#' \code{gsva::GSVA()}.
#' @param smooth Logical value. If set to TRUE values will be
#' smoothed with respect to their local neighbors using \code{stats::loess()}.
#' @param smooth_span Numeric value, given to \code{stats::loess()} if
#'  \code{smooth} is set to TRUE.
#' @param verbose Logical value. If set to TRUE informative messages with respect
#' to the computational progress made will be printed.
#'
#' (Warning messages will always be printed.)
#'
#' @return The input data.frame of \code{coords_df} joined with all the aspects
#' specified variables specified in \code{features} by the key \emph{barcodes}.
#'
#' @export

joinWithFeatures <- function(object,
                             coords_df,
                             features,
                             smooth = FALSE,
                             smooth_span = 0.02,
                             verbose = TRUE){


  # 1. Control valid input --------------------------------------------------

  validation(x = object)

  coords_df <- check_coords_df(coords_df)

  features <- check_features(object, features = features)

  smooth <- check_smooth(df = coords_df, smooth = smooth, verbose = verbose)


  # 2. Join data ------------------------------------------------------------

  fdata <-
    featureData(object, of_sample = base::unique(coords_df$sample)) %>%
    dplyr::select(dplyr::all_of(x = c("barcodes", features)))

  joined_df <-
    dplyr::left_join(x = coords_df, y = fdata, by = "barcodes")



  # 3. Smooth if specified -------------------------------------------------

  if(base::isTRUE(smooth)){

    joined_df <-
      purrr::imap_dfr(.x = joined_df,
                      .f = hlpr_smooth,
                      coords_df = joined_df,
                      verbose = verbose,
                      smooth_span = smooth_span,
                      aspect = "feature",
                      subset = features)

  }


  # 4. Return final data frame ----------------------------------------------

  base::return(joined_df)

}


#' @rdname joinWithFeatures
#' @export
joinWithGenes <- function(object,
                          coords_df,
                          genes,
                          average_genes = FALSE,
                          smooth = FALSE,
                          smooth_span = 0.02,
                          normalize = TRUE,
                          verbose = TRUE){

  # 1. Control valid input --------------------------------------------------

  validation(x = object)

  coords_df <- check_coords_df(coords_df)

  rna_assay <- exprMtr(object, of_sample = base::unique(coords_df$sample))

  genes <- check_genes(object, genes = genes, rna_assay = rna_assay)

  smooth <- check_smooth(df = coords_df, smooth = smooth, verbose = verbose)


  # 2. Extract expression values  -------------------------------------------

  rows_to_subset <- genes
  columns_to_subset <- coords_df$barcodes

  # compute mean if necessary
  if(base::length(genes) > 1 && average_genes){

    rna_assay <- base::colMeans(rna_assay[rows_to_subset, columns_to_subset])
    col_names <- "mean_genes"

  } else if(base::length(genes) > 1){

    rna_assay <- t(rna_assay[rows_to_subset, columns_to_subset])
    col_names <- genes

  } else if(base::length(genes) == 1){

    rna_assay <- rna_assay[rows_to_subset, columns_to_subset]

    if(base::isTRUE(average_genes)){
      col_names <- "mean_genes"
    } else {
      col_names <- genes
    }


  }

  # convert results to data frame with appropriate column names
  gene_vls <-
    base::as.data.frame(rna_assay, row.names = NULL) %>%
    magrittr::set_colnames(value = col_names) %>%
    dplyr::mutate(
      barcodes = base::rownames(rna_assay)
    )


  # join both
  joined_df <-
    dplyr::left_join(x = coords_df, y = gene_vls, by = "barcodes")


  # 3. Smooth if specified --------------------------------------------------

  if(base::isTRUE(smooth)){

    joined_df <-
      purrr::imap_dfr(.x = joined_df,
                      .f = hlpr_smooth,
                      coords_df = joined_df,
                      verbose = verbose,
                      smooth_span = smooth_span,
                      aspect = "gene",
                      subset = col_names)


  }

  if(base::isTRUE(normalize)){

    joined_df <-
      purrr::imap_dfr(.x = joined_df,
                      .f = hlpr_normalize_imap,
                      aspect = "Gene",
                      verbose = verbose,
                      subset = col_names
      )

  }


  # 4. Return final data frame ----------------------------------------------

  base::return(joined_df)

}


#' @rdname joinWithFeatures
#' @export
joinWithGeneSets <- function(object,
                             coords_df,
                             gene_sets,
                             method_gs = "mean",
                             smooth = FALSE,
                             smooth_span = 0.02,
                             normalize = TRUE,
                             verbose = TRUE){

  # 1. Control valid input --------------------------------------------------

  validation(x = object)

  coords_df <- check_coords_df(coords_df)

  gene_sets <- check_gene_sets(object, gene_sets = gene_sets)

  smooth <- check_smooth(df = coords_df, smooth = smooth, verbose = verbose)

  # 2. Extract gene set data and join with coords_df ------------------------

  rna_assay <- exprMtr(object = object, of_sample = base::unique(coords_df$sample))
  gene_set_df <- object@used_genesets
  joined_df <- coords_df

  for(i in base::seq_along(gene_sets)){

    # get gene names of specified gene sets
    genes <-
      gene_set_df %>%
      dplyr::filter(ont %in% gene_sets[i]) %>%
      dplyr::filter(gene %in% base::rownames(rna_assay)) %>%
      dplyr::pull(gene)

    # apply specified method to handle gene sets
    if(method_gs == "mean"){

      geneset_vls <-
        base::colMeans(rna_assay[genes, ]) %>%
        base::as.data.frame() %>%
        magrittr::set_colnames(value = gene_sets[i]) %>%
        tibble::rownames_to_column(var = "barcodes")

      if(verbose){

        base::message(stringr::str_c(
          "Calculating expression score for gene set '",
          gene_sets[i],
          "' according to method: '",
          method_gs,
          "'.",
          sep = ""))

      }


    } else if(method_gs %in% c("gsva", "ssgsea", "zscore", "plage")) {

      if(verbose){

        base::message(stringr::str_c(
          "Calculating expression score for gene set '",
          gene_sets[i],
          "' according to method: '",
          method_gs,
          "'. This might take a few moments.",
          sep = ""))

      }

      geneset_vls <-
        GSVA::gsva(expr = rna_assay[genes,],
                   gset.idx.list = gene_set_df,
                   mx.diff = 1,
                   parallel.sz = 2,
                   method = method_gs,
                   verbose = FALSE) %>%
        t() %>%
        as.data.frame() %>%
        magrittr::set_colnames(value = gene_sets[i]) %>%
        tibble::rownames_to_column(var = "barcodes")

    } else {

      stop("Please enter valid method to handle gene sets.")

    }

    # gradually add gene_set columns to joined_df
    joined_df <-
      dplyr::left_join(x = joined_df, y = geneset_vls, by = "barcodes")

  }


  # 3. Smooth and normalize if specified ------------------------------------

  if(base::isTRUE(smooth)){

    joined_df <-
      purrr::imap_dfr(.x = joined_df,
                      .f = hlpr_smooth,
                      coords_df = joined_df,
                      verbose = verbose,
                      smooth_span = smooth_span,
                      aspect = "gene set",
                      subset = gene_sets)

  }

  if(base::isTRUE(normalize)){

    joined_df <-
      purrr::imap_dfr(.x = joined_df,
                      .f = hlpr_normalize_imap,
                      aspect = "Gene set",
                      verbose = verbose,
                      subset = gene_sets)

  }


  # 4. Return final data frame ----------------------------------------------

  base::return(joined_df)

}


