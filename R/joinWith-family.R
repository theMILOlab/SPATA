#' @title Join barcodes with additional variables
#'
#' @description Each member of the \code{joinWith()}-family takes a data.frame as
#' input that contains at least the variables \emph{barcodes} and \emph{sample}.
#' (Easily obtained with the \code{get*()}-family.) It extracts the specified
#' variables and joins them over the barcode variable with the provided data.frame.
#'
#' @param object A valid spata-object.
#' @inherit check_coords_df params
#' @inherit check_features params
#' @inherit check_gene_sets params
#' @inherit check_genes params
#' @inherit check_smooth params
#' @inherit check_method params
#' @inherit verbose params
#' @inherit check_variables params
#' @inherit normalize params
#' @param average_genes Logical. If set to TRUE the average expression of the
#' specified genes is calculated and saved under one variable named 'mean_genes'.
#'
#' @details Hint: Variables of the specified data.frame \code{coords_df} that have equal names as
#' the specified features, genes and gene-sets are overwritten!
#'
#' @return The input data.frame of \code{coords_df} joined with all the
#' specified variable-elements (by the key-variable \emph{barcodes}).
#'
#' @export


joinWith <- function(object,
                     coords_df,
                     features = NULL,
                     method_gs = "mean",
                     gene_sets = NULL,
                     genes = NULL,
                     average_genes = FALSE,
                     smooth = FALSE,
                     smooth_span = 0.02,
                     verbose = TRUE,
                     normalize = TRUE){


  confuns::check_data_frame(
    df = coords_df,
    var.class = list(
      "barcodes" = "character",
      "sample" = "character"
    ),
    ref = "coords_df")

  input_list <-
    purrr::map2(.x = list(gene_sets, genes, features),
                .y = c("gene_sets", "genes", "features"),
                .f = function(x, y){

                  if(!base::is.null(x)){

                    confuns::is_vec(x = x, mode = "character", ref = y)

                    return(x)

                  } else {

                    return(x)

                  }

                })


  variables <- check_variables(
    variables = input_list,
    all_features = getFeatureNames(object),
    all_gene_sets = getGeneSets(object),
    all_genes = getGenes(object, in_sample = base::unique(coords_df$sample))
  )

  output_df <-
  joinWithVariables(
    object = object,
    coords_df = coords_df,
    variables = variables,
    method_gs = method_gs,
    average_genes = average_genes,
    smooth = smooth,
    smooth_span = smooth_span,
    verbose = verbose,
    normalize = normalize
  )

  base::return(output_df)

}

#' @rdname joinWith
#' @export
joinWithFeatures <- function(object,
                             coords_df,
                             features,
                             smooth = FALSE,
                             smooth_span = 0.02,
                             verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  # lazy check

  check_object(object)
  check_coords_df(coords_df)
  check_smooth(df = coords_df, smooth = smooth, smooth_span = smooth_span)

  # adjusting check
  features <- check_features(object, features = features)

  smooth_string <- base::ifelse(test = base::isTRUE(verbose), yes = " and smoothing ", no = " " )

  if(base::isTRUE(verbose)){

    feature_string <- base::ifelse(base::length(features) == 1, "feature", "features")

    base::message(glue::glue("Joining{smooth_string}{base::length(features)} {feature_string}."))

    }

  # overwrite check
  discard <- features[features %in% base::colnames(coords_df)]
  n_discard <- base::length(discard)

  if(n_discard > 0){

    variable_string <- base::ifelse(n_discard == 1, "variable", "variables")

    base::message(glue::glue("Overwriting {n_discard} feature-{variable_string}."))

    coords_df <- dplyr::select(.data = coords_df, -dplyr::all_of(discard))

  }

  # -----

  # 2. Join data ------------------------------------------------------------

  fdata <-
    featureData(object, of_sample = base::unique(coords_df$sample)) %>%
    dplyr::select(dplyr::all_of(x = c("barcodes", features)))

  joined_df <-
    dplyr::left_join(x = coords_df, y = fdata, by = "barcodes")

  # -----

  # 3. Smooth if specified -------------------------------------------------

  if(base::isTRUE(smooth)){

    joined_df <-
      purrr::imap_dfr(.x = joined_df,
                      .f = hlpr_smooth,
                      coords_df = joined_df,
                      smooth_span = smooth_span,
                      aspect = "feature",
                      subset = features)

  }

  if(base::isTRUE(verbose)){base::message("Done.")}

  # -----

  base::return(joined_df)

}


#' @rdname joinWith
#' @export
joinWithGenes <- function(object,
                          coords_df,
                          genes,
                          average_genes = FALSE,
                          smooth = FALSE,
                          smooth_span = 0.02,
                          normalize = TRUE,
                          verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  # lazy check

  check_object(object)
  check_coords_df(coords_df)
  check_smooth(df = coords_df, smooth = smooth, smooth_span = smooth_span)

  # adjusting check
  rna_assay <- exprMtr(object, of_sample = base::unique(coords_df$sample))
  genes <- check_genes(object, genes = genes, rna_assay = rna_assay)

  # -----

  # 2. Extract expression values and join with coords_df ---------------------

  barcodes <- coords_df$barcodes
  num_genes <- base::length(genes)

  # compute mean if necessary
  if(num_genes > 1 && average_genes){

    if(base::isTRUE(verbose) && base::isTRUE(smooth)){

      gene_string <- base::ifelse(num_genes == 1, "gene", "genes")

      base::message(glue::glue("Averaging, joining and smoothing {num_genes} {gene_string}."))

    } else if(base::isTRUE(verbose)){

      base::message(glue::glue("Averaging and joining {num_genes} genes."))

    }

    rna_assay <- base::colMeans(rna_assay[genes, barcodes])
    col_names <- "mean_genes"
    num_genes <- "averaged"

  } else if(num_genes > 1){

    if(base::isTRUE(verbose) && base::isTRUE(smooth)){

      base::message(glue::glue("Joining and smoothing {num_genes} genes."))

    } else if(base::isTRUE(verbose)){

      base::message(glue::glue("Joining {num_genes} genes."))

    }

    rna_assay <- t(rna_assay[genes, barcodes])
    col_names <- genes

  } else if(num_genes == 1){

    rna_assay <- rna_assay[genes, barcodes]

    if(base::isTRUE(average_genes)){
      col_names <- "mean_genes"
      num_genes <- "averaged"
    } else {
      col_names <- genes
    }


  }

  # convert results to data frame with appropriate column names
  gene_vls <-
    base::as.data.frame(rna_assay, row.names = NULL) %>%
    magrittr::set_colnames(value = col_names) %>%
    dplyr::mutate(
      barcodes = coords_df$barcodes
    )


  # overwrite check
  discard <- col_names[col_names %in% base::colnames(coords_df)]
  n_discard <- base::length(discard)

  if(n_discard > 0){

    variable_string <- base::ifelse(n_discard == 1, "variable", "variables")

    base::message(glue::glue("Overwriting {n_discard} gene-{variable_string}."))
    coords_df <- dplyr::select(.data = coords_df, -dplyr::all_of(discard))

  }


  # join both
  joined_df <-
    dplyr::left_join(x = coords_df, y = gene_vls, by = "barcodes")

  # -----

  # 3. Smooth and normalize if specified ------------------------------------

  if(base::isTRUE(smooth)){

    if(base::isTRUE(verbose)){

      pb <- progress::progress_bar$new(
        format = "Progress: [:bar] :percent eta: :eta",
        total = base::ncol(joined_df), clear = FALSE, width = 100)

    } else {

      pb <- NULL

    }

    joined_df <-
      purrr::imap_dfr(.x = joined_df,
                      .f = hlpr_smooth,
                      coords_df = joined_df,
                      smooth_span = smooth_span,
                      aspect = "gene",
                      subset = col_names,
                      pb = pb)


  }

  if(base::isTRUE(normalize)){

    if(base::isTRUE(verbose)){base::message(glue::glue("Normalizing values."))}

    joined_df <-
      purrr::imap_dfr(.x = joined_df,
                      .f = hlpr_normalize_imap,
                      aspect = "Gene",
                      subset = col_names
      )

  }

  if(base::isTRUE(verbose)){

    base::message("Done.")

  }

  # -----

  base::return(joined_df)

}


#' @rdname joinWith
#' @export
joinWithGeneSets <- function(object,
                             coords_df,
                             gene_sets,
                             method_gs = "mean",
                             smooth = FALSE,
                             smooth_span = 0.02,
                             normalize = TRUE,
                             verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  # lazy check

  check_object(object)
  check_coords_df(coords_df)
  check_smooth(df = coords_df, smooth = smooth, smooth_span = smooth_span)
  check_method(method_gs = method_gs)

  # adjusting check
  gene_sets <- check_gene_sets(object, gene_sets = gene_sets)

  # overwrite check
  discard <- gene_sets[gene_sets %in% base::colnames(coords_df)]
  n_discard <- base::length(discard)

  if(n_discard > 0){

    variable_string <- base::ifelse(n_discard == 1, "variable", "variables")

    base::message(glue::glue("Overwriting {n_discard} gene-set-{variable_string}."))
    coords_df <- dplyr::select(.data = coords_df, -dplyr::all_of(discard))

  }
  # -----

  # 2. Extract gene set data and join with coords_df ------------------------

  rna_assay <- exprMtr(object = object, of_sample = base::unique(coords_df$sample))
  gene_set_df <- object@used_genesets
  joined_df <- coords_df

  if(base::isTRUE(smooth)){

    x <- dplyr::pull(coords_df, var = x)
    y <- dplyr::pull(coords_df, var = y)
    smooth_string <- glue::glue(" and smoothing ")

  } else {

    smooth_string <- " "

  }


  #feedback vectors
  filter_gs <- 0.25
  ignored_gs <- glue::glue("\nIgnored gene-sets due to insufficient gene representation (less then {filter_gs*100}%) in expression matrix:")
  skipped_gs <- base::character()

  num_gs <- base::length(gene_sets)


  if(base::isTRUE(verbose)){

    gene_set_string <- base::ifelse(num_gs == 1, "gene-set", "gene-sets")

    base::message(glue::glue("Calculating{smooth_string}expression score for {base::length(gene_sets)} {gene_set_string} according to method '{method_gs}'."))

  }

  if(base::isTRUE(verbose)){

    pb <- progress::progress_bar$new(
      format = "Progress: [:bar] :percent eta: :eta",
      total = num_gs, clear = FALSE, width = 100)

  }


  for(i in base::seq_along(gene_sets)){

    if(base::isTRUE(verbose)){pb$tick()}

    # get genes of gene set
    gs_df <- dplyr::filter(gene_set_df, ont %in% gene_sets[i])

    n_genes <- base::nrow(gs_df)

    # get genes found in expression matrix
    genes <-
      dplyr::filter(gs_df, gene %in% base::rownames(rna_assay)) %>%
      dplyr::pull(gene)

    n_found_genes <- base::length(genes)

    # calculate percentage of genes found
    p_found_genes <- base::round(n_found_genes/n_genes, digits = 2)

    # make sure that percentage is equal to or higher than the threshold
    if(p_found_genes >= filter_gs){

      # apply specified method to handle gene sets
      if(method_gs == "mean"){

        geneset_vls <-
          base::colMeans(rna_assay[genes, ]) %>%
          base::as.data.frame() %>%
          magrittr::set_colnames(value = gene_sets[i]) %>%
          tibble::rownames_to_column(var = "barcodes")

      } else if(method_gs %in% c("gsva", "ssgsea", "zscore", "plage")) {

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

      }

      # smoothing
      if(base::isTRUE(smooth)){

        variable <- dplyr::pull(.data = geneset_vls, var = gene_sets[i])
        model <- stats::loess(formula = variable ~ x*y, span = smooth_span)

        geneset_vls[, gene_sets[i]] <- stats::predict(model)

      }

      # gradually add gene_set columns to joined_df
      joined_df <-
        dplyr::left_join(x = joined_df, y = geneset_vls, by = "barcodes")

    } else {

      skipped_gs <- base::append(x = skipped_gs, values = gene_sets[i])

      ignored_gs <-
        base::append(x = ignored_gs,
                     values = glue::glue("\n- '{gene_sets[i]}'. Percentage of genes found: {p_found_genes}"))

    }

  }


  # -----


  # 3. Normalize if specified -----------------------------------------------

  if(base::isTRUE(normalize)){

    if(base::isTRUE(verbose)){base::message(glue::glue("Normalizing values."))}

    # normalize
    joined_df <-
      purrr::imap_dfr(.x = joined_df,
                    .f = hlpr_normalize_imap,
                    aspect = "Gene set",
                    subset = gene_sets)

  }

  if(base::isTRUE(verbose)){base::message("Done.")}

  if(base::length(ignored_gs) > 1){

    base::message("Warning:")
    base::append(x = ignored_gs,
                 values = "\n(Run 'adjustGeneSetDf()' in order to avoid this warning message.)") %>%
    base::writeLines()

  }

  # -----

  base::return(joined_df)

}

#' @rdname joinWith
#' @export
joinWithVariables <- function(object,
                              coords_df,
                              variables,
                              method_gs = "mean",
                              average_genes = FALSE,
                              smooth = FALSE,
                              smooth_span = 0.02,
                              normalize = TRUE,
                              verbose = TRUE){

  stopifnot(base::is.list(variables))
  stopifnot(base::any(c("features", "genes", "gene_sets") %in% base::names(variables)))

  if("features" %in% base::names(variables)){

    coords_df <-
      joinWithFeatures(object = object,
                       features = variables$features,
                       coords_df = coords_df,
                       smooth = smooth,
                       smooth_span = smooth_span,
                       verbose = verbose)

  }

  if("genes" %in% base::names(variables)){

    coords_df <-
      joinWithGenes(object = object,
                    coords_df = coords_df,
                    genes = variables$genes,
                    average_genes = average_genes,
                    smooth = smooth,
                    smooth_span = smooth_span,
                    normalize = normalize,
                    verbose = verbose)

  }

  if("gene_sets" %in% base::names(variables)){

    coords_df <-
      joinWithGeneSets(object = object,
                       coords_df = coords_df,
                       gene_sets = variables$gene_sets,
                       method_gs = method_gs,
                       smooth = smooth,
                       smooth_span = smooth_span,
                       normalize = normalize,
                       verbose = verbose)
  }

  base::return(coords_df)

}


