
# Helper functions --------------------------------------------------------

#' @title Safe extraction
#'
#' @description A wrapper around \code{base::tryCatch()} with predefined error handling
#' messages if extraction from seurat object failed.
#'
#'
#' @param return_value Whatever needs to be extracted.
#' @param error_handling Either \emph{'warning} or \emph{'stop'}.
#' @param error_value What is supposed to be returned if extractino fails.
#' @param error_ref The reference for the feedback message.
#'


getFromSeurat <- function(return_value, error_handling, error_value, error_ref){

  result <-
    base::tryCatch(

      return_value,

      error = function(error){

        if(error_handling == "warning"){

          base::warning(glue::glue("Could not find {error_ref} in specified seurat object. Did you choose the correct method?"))

        } else if(error_handling == "stop"){

          base::stop(glue::glue("Could not find {error_ref} in specified seurat object. Did you choose the correct method?"))

        }

        base::return(error_value)


        })


  base::return(result)

}


# Transform functions -----------------------------------------------------

#' @title Transform seurat-object to spata-object
#'
#' @param seurat_object A valid seurat object.
#' @param method Character value. Determines the data slots from which to compile the spata-object.
#'
#'  \describe{
#'   \item{\emph{'spatial'}}{Denotes that the data to be used derived from spatial experiments. }
#'   \item{\emph{'single_cell'}}{Denotes that the data to be used derived from single cell experiments. }
#'  }
#'
#' @param coords_from Character value. Either \emph{'umap'} or \emph{'tsne'}.
#'
#'  Only relevant if \code{method} was set to \emph{'single_cell'}. Denotes the slot from which to
#'  take the surrogate coordinates.
#'
#' @param sample_name Character value. Future input for SPATA's \code{of_sample}-argument.
#' @inherit loadGSDF params
#' @inherit verbose params
#'
#' @return A spata object.
#' @export
#'

transformSeuratToSpata <- function(seurat_object,
                                   sample_name,
                                   assay = "Spatial",
                                   slice_name = "slice1",
                                   method = "spatial",
                                   coords_from = "umap",
                                   gene_set_path = NULL,
                                   verbose = TRUE){

# 1. Control --------------------------------------------------------------

  methods::is(object = seurat_object, class2 = "Seurata")

  confuns::check_one_of(input = method, against = seurat_methods, ref.input = "input for argument 'method'")

  confuns::are_values("sample_name", "method", mode = "character")

  # spata object
  spata_object <- methods::new(Class = "spata", samples = sample_name)

  if(base::is.null(gene_set_path) | base::is.character(gene_set_path)){

    spata_object@used_genesets <-
      loadGSDF(gene_set_path = gene_set_path, verbose = verbose)

  }


# 2. Extract data ---------------------------------------------------------

  if(method == "spatial"){

    # get coordinates
    coords_df <-
      getFromSeurat(
        return_value = Seurat::GetTissueCoordinates(seurat_object),
        error_handling = "stop",
        error_ref = "coordinates",
        error_value = NULL
      )

    coords_df <-
      tibble::rownames_to_column(.data = coords_df, var = "barcodes") %>%
      dplyr::rename("x" = "imagecol", "y" = "imagerow") %>%
      dplyr::mutate(sample = {{sample_name}}) %>%
      dplyr::select(barcodes, sample, x, y)

    slice <-
      getFromSeurat(
        return_value = seurat_object@images[[slice_name]],
        error_handling = "stop",
        error_ref = glue::glue("slice-object '{slice_name}'"),
        error_value = NULL
      )

    spata_object@compatibility <- list("Seurat" = list("slice" = slice))

    # get scaled matrix
    scaled_mtr <-
      getFromSeurat(
        return_value = seurat_object@assays[[assay]]@scale.data,
        error_handling = "stop",
        error_ref = "scaled matrix",
        error_value = NULL
      )

    # get count matrix
    count_matrix <-
      getFromSeurat(
        return_value = seurat_object@assays[[assay]]@counts,
        error_handling = "warning",
        error_value = base::matrix(),
        error_ref = "count matrix"
      )

    # get image
    image <-
      getFromSeurat(
        return_value = seurat_object@images$slice1[1]@image,
        error_handling = "warning",
        error_value = NULL,
        error_ref = "image"
      )

    if(!base::is.null(image)){

      image <-
        EBImage::Image(image, colormode = "Color") %>%
        EBImage::transpose()

    }

  } else if(method == "single_cell") {

    confuns::is_value(x = coords_from, mode = "character", ref = "coords_from")
    confuns::check_one_of(input = coords_from, against = seurat_coords_from_opts, ref.input = "input for argument 'coords_from'")

    first_choice <- coords_from
    second_choice <- seurat_coords_from_opts[seurat_coords_from_opts != coords_from]

    # get coordinates/ umap cell embedding
    coords_df <-
      getFromSeurat(
        return_value = base::as.data.frame(seurat_object@reductions[[first_choice]]@cell.embeddings),
        error_handling = "warning",
        error_value = NULL,
        error_ref = glue::glue("coordinates/{first_choice} cell embedding")
      )

    # try tsne if umap did not work
    if(base::is.null(coords_df)){

      base::warning(glue::glue("Trying to extract surrogate coordinates from slot {second_choice}."))

      coords_df <-
        getFromSeurat(
          return_value = base::as.data.frame(seurat_object@reductions[[second_choice]]@cell.embeddings),
          error_handling = "stop",
          error_value = NULL,
          error_ref = glue::glue("coordinates/{second_choice} cell embedding")
        )

    }

    coords_df <-
      tibble::rownames_to_column(.data = coords_df, var = "barcodes") %>%
      magrittr::set_colnames(value = c("barcodes", "x", "y")) %>%
      dplyr::mutate(sample = {{sample_name}}) %>%
      dplyr::select(barcodes, sample, x, y)

    # get expression matrix
    scaled_mtr <-
      getFromSeurat(
        return_value = seurat_object@assays[[assay]]@scale.data,
        error_handling = "stop",
        error_value = NULL,
        error_ref = "scaled matrix"
      )

    # get count matrix
    count_matrix <-
      getFromSeurat(
        return_value = seurat_object@assays[[assay]]@counts,
        error_handling = "warning",
        error_value = base::matrix(),
        error_ref = "count matrix"
      )

    # no image
    image <- NULL

  }


# 3. Postprocess ----------------------------------------------------------

  # check if barcodes are identical
  barcodes_matrix <- base::colnames(scaled_mtr) %>% base::sort()
  barcodes_coordinates <- dplyr::pull(coords_df, var = "barcodes") %>% base::sort()

  if(!base::identical(barcodes_matrix, barcodes_coordinates)){

    base::stop("The barcodes of the coordinate system and the column names of the assay must be identical. Please check the seurat object for integrity.")

  }

  # feature data

  seurat_object@meta.data$barcodes <- NULL

  fdata <-
    tibble::rownames_to_column(.data = seurat_object@meta.data, var = "barcodes") %>%
    dplyr::mutate(segment = "none") %>%
    dplyr::select(barcodes, dplyr::everything())

  # savely discard colum 'orig.ident'
  fdata <- base::tryCatch(

    dplyr::select(fdata, -orig.ident),

    error = function(error){ fdata }

  )

  spata_object <- setFeatureDf(object = spata_object, feature_df = fdata)

# 4. Pass to Spata --------------------------------------------------------


  # dimensional reduction: pca

  pca_df <- base::tryCatch({

    pca_df <-
      base::as.data.frame(seurat_object@reductions$pca@cell.embeddings) %>%
      tibble::rownames_to_column(var = "barcodes") %>%
      dplyr::select(barcodes, dplyr::everything())

    base::colnames(pca_df) <- stringr::str_remove_all(base::colnames(pca_df), pattern = "_")

    pca_df

    },

    error = function(error){

      base::warning("Could not find or transfer PCA-data. Did you process the seurat-object correctly?")

      base::return(data.frame())

    }

  )

  spata_object <- setPcaDf(object = spata_object, pca_df = pca_df)


  # dimensional reduction: umap

  umap_df <- base::tryCatch({

    base::data.frame(
      barcodes = base::rownames(seurat_object@reductions$umap@cell.embeddings),
      umap1 = seurat_object@reductions$umap@cell.embeddings[,1],
      umap2 = seurat_object@reductions$umap@cell.embeddings[,2],
      stringsAsFactors = FALSE
    ) %>% tibble::remove_rownames()

    }, error = function(error){

      base::warning("Could not find or transfer UMAP-data. Did you process the seurat-object correctly?")

      base::return(data.frame())

    }

  )

  spata_object <- setUmapDf(object = spata_object, umap_df = umap_df)


  # dimensional reduction: tsne

  tsne_df <- base::tryCatch({

    base::data.frame(
      barcodes = base::rownames(seurat_object@reductions$tsne@cell.embeddings),
      tsne1 = seurat_object@reductions$tsne@cell.embeddings[,1],
      tsne2 = seurat_object@reductions$tsne@cell.embeddings[,2],
      stringsAsFactors = FALSE
    ) %>% tibble::remove_rownames()

    }, error = function(error){

      base::warning("Could not find or transfer TSNE-data. Did you process the seurat-object correctly?")

      base::return(data.frame())

    }

  )

  spata_object <- setTsneDf(object = spata_object, tsne_df = tsne_df)


  # data matrices

  spata_object <-
    setCountMatrix(
      object = spata_object,
      count_mtr = count_matrix[base::rowSums(base::as.matrix(count_matrix)) != 0, ]
      )

  spata_object <-
    setScaledMatrix(
      object = spata_object,
      scaled_mtr = scaled_mtr[base::rowSums(base::as.matrix(scaled_mtr)) != 0, ]
      )


  # coordinates & image

  spata_object <-
    setCoordsDf(object = spata_object, coords_df = coords_df) %>%
    setImage(object = ., image = image)


  # other lists
  spata_object@information <-
    list("active_mtr" = magrittr::set_names(x = list("scaled"), value = sample_name),
         "autoencoder" = magrittr::set_names(x = list(list()), value = sample_name),
         "barcodes" = magrittr::set_names(x = list(barcodes_matrix), value = sample_name))

  spata_object@trajectories <-
    magrittr::set_names(x = list(list()), value = sample_name)

  spata_object@version <- current_spata_version


# 5. Return spata object ---------------------------------------------------

  base::return(spata_object)

}



#' @title Transform spata-object to cell-data-set (Monocle3)
#'
#' @description Takes the count matrix of your spata-object and creates a
#' cell_data_set-object with it. See details for more information on how to use
#' the arguments.
#'
#' @inherit check_object params
#' @inherit check_monocle_input params details
#' @param estimate_size_factors_args A list of arguments given to \code{monocle3::estimate_size_factors()}.
#' @param preprocess_cds_args A list of arguments given to \code{monocle3::preprocess_cds()}.
#' @param reduce_dimension_args A list of arguments given to \code{monocle3::reduce_dimension()}.
#' @param cluster_cells_args A list of arguments given to \code{monocle3::cluster_cells()}.
#' @param learn_graph_args A list of arguments given to \code{monocle3::learn_graph()}.
#' @param order_cells_args A list of arguments given to \code{monocle3::order_cells()}.
#' @param save_cds_file Character value or NULL. A file-directory (that does not already exists) under which created cell_data_set-object
#' is saved. Should end with \emph{'.RDS'}.
#' @inherit verbose params
#'
#' @details \code{compileCellDataSet()} is a convenient wrapper around all pre processing functions
#' monocle3 provides to handle it's core object - the cell_data_set - after it's initiation. Apart from unique
#' arguments this function has two argument families, denoted with \code{_method} and \code{_args}.
#'
#' Handling \code{_method}-arguments:
#'
#' Monocle3 allows to use different methods for dimensional-reduction or clustering which depend
#' on each other. These arguments take a character vector of all valid inputs. \code{compileCellDataSet()} iterates
#' over all valid combinations and returns the cell_data_set with the computed information inside.
#'
#' Handling \code{_args}-arguments.
#'
#' These arguments take named lists of arguments that are given to the respective function. The \code{_method}-arguments
#' as well as the argument \code{cds} are automatically defined and must not be included in the given lists!!! Empty lists - the default -
#' result in running the function with it's default parameters.
#'
#' The spata-objects feature data (@@fdata) is passed to the cell_data_set for it's slot \code{cell_meta_data}
#'
#' @return A monocle3::cell_data_set object.
#' @export

transformSpataToCDS <- function(object,
                                preprocess_method = "PCA",
                                reduction_method = c("PCA", "UMAP"),
                                cluster_method = "leiden",
                                estimate_size_factors_args = list(),
                                preprocess_cds_args = list(),
                                reduce_dimension_args = list(),
                                cluster_cells_args = list(),
                                learn_graph_args = list(),
                                order_cells_args = list(),
                                save_cds_file = NULL,
                                verbose = TRUE){

  check_object(object)
  confuns::is_value(preprocess_method, "character", "preprocess_method")
  confuns::is_value(cluster_method, mode = "character", "cluster_method")

  monocle_funs <-
    rlang::fn_fmls_names(fn = compileCellDataSet) %>%
    stringr::str_subset(pattern = "args$")

  for(mf in monocle_funs){

    input <- base::parse(text = mf) %>% base::eval()

    if(!base::is.list(input) | base::is.data.frame(input)){

      base::stop(glue::glue("Invalid input for argument '{mf}'. Must be a named list of arguments."))

    }

  }

  check_monocle_input(preprocess_method = preprocess_method,
                      reduction_method = reduction_method,
                      cluster_method = cluster_method)

  if(!base::is.null(save_cds_file)){

    confuns::is_value(save_cds_file, "character", "save_cds_file")
    if(base::file.exists(save_cds_file)){

      base::stop(glue::glue("Directory '{save_cds_file}' already exists. "))

    }

  }

  # check if valid cds files

  # Step 1 - Create cds -----------------------------------------------------

  if(base::isTRUE(verbose)){base::message("No cds-file specified. Performing monocle anylsis from scratch.")}

  base::stopifnot(preprocess_method %in% c("PCA", "LSI"))
  base::stopifnot(cluster_method %in% c("leiden", "louvain"))

  if(base::isTRUE(verbose)){base::message("Step 1/7 Creating 'cell data set'-object.")}

  expression_matrix <- base::as.matrix(object@data@counts)

  gene_metadata <- data.frame(gene_short_name = base::rownames(expression_matrix))
  base::rownames(gene_metadata) <- base::rownames(expression_matrix)

  cell_metadata <- data.frame(object@fdata)
  base::rownames(cell_metadata) <- object@fdata$barcodes

  cds <- monocle3::new_cell_data_set(
    expression_data = expression_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata)

  cds <- cds[,Matrix::colSums(monocle3::exprs(cds)) != 0]

  # -----



  # Step 2-4 Estimate size factors, preprocess, reduce dimensions -----------

  if(base::isTRUE(verbose)){base::message("Step 2/7 Estimating size factors.")}

  estimate_size_factors_args <- purrr::prepend(x = estimate_size_factors_args,
                                               values = list("cds" = cds))

  cds <- rlang::invoke(.fn = base::eval(base::parse(text = "monocle3::estimate_size_factors")), estimate_size_factors_args)

  if(base::isTRUE(verbose)){base::message("Step 3/7 Preprocessing cell data set.")}

  for(p in base::seq_along(preprocess_method)){

    if(base::isTRUE(verbose)){

      base::message(glue::glue("Preprocessing cells with method {p}/{base::length(preprocess_method)} '{preprocess_method[p]}'"))

    }

    preprocess_cds_args_p <- purrr::prepend(x = preprocess_cds_args,
                                            values = list("cds" = cds, "preprocess_method" = preprocess_method[p]))

    cds <- rlang::invoke(.fn = base::eval(base::parse(text = "monocle3::preprocess_cds")), preprocess_cds_args_p)

  }

  if(base::isTRUE(verbose)){base::message("Step 4/7 Reducing dimensions.")}

  for(p in base::seq_along(preprocess_method)){

    base::message(glue::glue("Using preprocess method '{preprocess_method[p]}':"))

    for(r in base::seq_along(reduction_method)){

      base::message(glue::glue("Reducing dimensions with reduction method {r}/{base::length(reduction_method)}: '{reduction_method[r]}' "))

      if(reduction_method[r] == "LSI" && preprocess_method[p] != "LSI"){

        base::message(glue::glue("Ignoring invalid combination. reduction-method: '{reduction_method[r]}' &  preprocess-method: '{preprocess}'"))

      } else if(reduction_method[r] == "PCA" && preprocess_method[p] != "PCA") {

        base::message(glue::glue("Ignoring invalid combination. reduction-method: '{reduction_method[r]}' &  preprocess-method: '{preprocess}'"))

      } else {

        reduce_dimension_args_r <- purrr::prepend(x = reduce_dimension_args,
                                                  values = list("cds" = cds,
                                                                reduction_method = reduction_method[r],
                                                                preprocess_method = preprocess_method[p]))

        cds <- base::tryCatch(

          rlang::invoke(.fn = base::eval(base::parse(text = "monocle3::reduce_dimension")), reduce_dimension_args_r),

          error = function(error){

            base::message(glue::glue("Attempting to call 'reduce_dimensions()' resulted in an error: {error$message}.
                                       Skipping."))

            base::return(cds)

          })

      }

    }

  }

  # -----

  # Step 5 Cluster cells ----------------------------------------------------

  if(base::isTRUE(verbose)){base::message("Step 5/7 Clustering cells.")}

  for(r in base::seq_along(reduction_method)){

    if(base::isTRUE(verbose)){

      base::message(glue::glue("Using reduction method {reduction_method[r]}:"))

    }

    for(c in base::seq_along(cluster_method)){

      if(base::isTRUE(verbose)){

        base::message(glue::glue("Clustering barcode-spots with method {c}/{base::length(cluster_method)}: {cluster_method[c]}"))

      }

      cluster_cells_args_c <- purrr::prepend(x = cluster_cells_args,
                                             values = list("cds" = cds,
                                                           "reduction_method" = reduction_method[r],
                                                           "cluster_method" = cluster_method[c]))

      cds <- base::tryCatch(

        rlang::invoke(.fn = base::eval(base::parse(text = "monocle3::cluster_cells")), cluster_cells_args_c),

        error = function(error){

          base::message(glue::glue("Attempting to call 'cluster_cells()' resulted in an error: {error$message}.
                                     Skipping."))

          base::return(cds)

        })

    }

  }

  if(base::isTRUE(verbose)){base::message("Step 6/7 Learning trajectory.")}

  learn_graph_args <- purrr::prepend(x = learn_graph_args, values = list(cds = cds))

  cds <- base::tryCatch(

    rlang::invoke(.fn = base::eval(base::parse(text = "monocle3::learn_graph")), learn_graph_args),

    error = function(error){

      base::message(glue::glue("Attempting to call 'learn_graph()' resulted in an error: {error$message}.
                               Skipping step 6/7."))

      base::return(cds)

    })

  if(base::isTRUE(verbose)){base::message("Step 7/7 Ordering cells.")}

  order_cells_args <- purrr::prepend(x = order_cells_args, values = list(cds = cds))

  cds <- base::tryCatch(

    rlang::invoke(.fn = base::eval(base::parse(text = "monocle3::order_cells")), order_cells_args),

    error = function(error){

      base::message(glue::glue("Attempting to call 'order_cells()' resulted in an error: {error$message}.
                               Skipping step 7/7."))

      base::return(cds)

    })

  # -----


  # Save cds-file and return ------------------------------------------------

  # save cds file if save_cds_file is specified as a character
  if(base::is.character(save_cds_file)){

    if(base::isTRUE(verbose)){

      base::message(stringr::str_c("Saving cell data set object 'cds' under directory: '", save_cds_file, "'"))

    }

    base::tryCatch(

      base::saveRDS(cds, file = save_cds_file),

      error = function(error){

        base::warning(glue::glue("Attempting to save the cell_data_set resulted in an error: {error}.
                                 Skip saving."))

      })

  }

  base::return(cds)

}


#' @title Transform spata-object to a seurat-object
#'
#' @description Takes the count matrix of your spata-object and creates a
#' Seurat-object with it. The spata-object's feature-data is passed as input
#' for the \code{meta.data}-argument of \code{Seurat::CreateSeuratObject()}.
#' If specified as TRUE or named list of arguments the respective functions are called in
#' order to pre process the object.
#'
#' The specified spata-object must contain only one sample! (use \code{subsetSpataObject()} to reduce
#' the number of samples). If you want to analyze several samples with Seurat please compile the objects one by one and
#' consider using \code{Seurat::merge()}.
#'
#' @inherit check_object params
#' @param assay Character value. The name under which the count- and expression matrix is to be saved.
#' @param ... Additional parameters given to \code{Seurat::CreateSeuratObject()}.
#' @param SCTransform A named list of arguments given to \code{Seurat::SCTransform()}, TRUE or FALSE.
#' @param NormalizeData A named list of arguments given to \code{Seurat::NormalizeData()}, TRUE or FALSE.
#' @param FindVariableFeatures A named list of arguments given to \code{Seurat::FindVariableFeatures()}, TRUE or FALSE.
#' @param ScaleData A named list of arguments given to \code{Seurat::ScaleData()}, TRUE or FALSE.
#'
#' Hint: If set to TRUE or the argument-list provided does not specify the argument \code{features} input
#' for argument \code{features} is set to \code{base::rownames(seurat_object)}.
#'
#' @param RunPCA A named list of arguments given to \code{Seurat::RunPCA()}, TRUE or FALSE.
#' @param FindNeighbors A named list of arguments given to \code{Seurat::FindNeighbors()}, TRUE or FALSE.
#' @param FindClusters A named list of arguments given to \code{Seurat::FindClusters()}, TRUE or FALSE.
#' @param RunTSNE A named list of arguments given to \code{Seurat::RunTSNE()}, TRUE or FALSE.
#' @param RunUMAP A named list of arguments given to \code{Seurat::RunUMAP()}, TRUE or FALSE.
#'
#' @details `compileSeuratObject()` is a convenient wrapper around all functions that preprocess a seurat-object
#' after it's initiation. The object is initiated by passing the spata-objects count-matrix and feature data to it whereupon
#' the functions are called in the order they are presented in this documentation. For all
#' pre processing functions apply the following instructions:
#'
#'  \itemize{
#'   \item{If set to FALSE the processing function is skipped.}
#'   \item{If set to TRUE the respective function is called with it's default argument settings. Note: \code{RunUMAP()} needs
#'   additional input!}
#'   \item{If a named list is provided the respective function called whereby the named list will provide the argument-input (the seurat-object to be constructed must not be named and will be
#'   passed to every function as the first argument!!!.)}
#'   }
#'
#' Note that certain listed functions require previous functions! E.g. if \code{RunPCA} is set to FALSE \code{RunTSNE()}
#' will result in an error. (\code{base::tryCatch()} will prevent the function from crashing.)
#'
#' @return
#' @export
#'

transformSpataToSeurat <- function(object,
                                   assay = "Spatial",
                                   ...,
                                   SCTransform = FALSE,
                                   NormalizeData = list(normalization.method = "LogNormalize", scale.factor = 1000),
                                   FindVariableFeatures = list(selection.method = "vst", nfeatures = 2000),
                                   ScaleData = TRUE,
                                   RunPCA = list(npcs = 60),
                                   FindNeighbors = list(dims = 1:30),
                                   FindClusters = list(resolution = 0.8),
                                   RunTSNE = TRUE,
                                   RunUMAP = list(dims = 1:30),
                                   verbose = TRUE){


  # 1. Control --------------------------------------------------------------

  check_object(object)
  sample <- getSampleNames(object)

  if(dplyr::n_distinct(sample) > 1){

    base::stop("The specified spata-object contains more than one sample. Please subset the object with 'subsetSpataObject()'.")

  }

  # -----

  # 2. Passing data ---------------------------------------------------------

  counts <- getCountMatrix(object)
  cnames_counts <- base::colnames(counts)

  pattern <- stringr::str_c("_", sample, "$", sep = "")
  cnames_new <- stringr::str_remove_all(string = cnames_counts, pattern = pattern)

  base::colnames(counts) <- cnames_new

  meta_data <-
    getFeatureDf(object) %>%
    dplyr::mutate(barcodes = stringr::str_remove_all(string = barcodes, pattern = pattern)) %>%
    tibble::column_to_rownames(var = "barcodes")

  seurat_object <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta_data, assay = assay, ...)

  seurat_object <- base::tryCatch({

    base::stopifnot(methods::is(object@compatibility$Seurat$slice, "VisiumV1"))

    seurat_object@images$slice1 <-
      object@compatibility$Seurat$slice

    seurat_object

    }, error = function(error){

      base::warning("The provided spata-object does not contain a valid VisiumV1-object. To use spatial features of the Seurat package you need to add that manually.")

      base::return(seurat_object)

    }
  )



  # 3. Processing seurat object ---------------------------------------------

  seurat_object <-
    process_seurat_object(
      seurat_object = seurat_object,
      assay = assay,
      SCTransform = SCTransform,
      NormalizeData = NormalizeData,
      FindVariableFeatures = FindVariableFeatures,
      ScaleData = ScaleData,
      RunPCA = RunPCA,
      FindNeighbors = FindNeighbors,
      FindClusters = FindClusters,
      RunTSNE = RunTSNE,
      RunUMAP = RunUMAP,
      verbose = verbose)


  # Passing features and images ---------------------------------------------


  if(base::isTRUE(verbose)){base::message("Done.")}

  return(seurat_object)

}








