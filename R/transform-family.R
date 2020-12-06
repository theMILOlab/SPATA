
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
                                   method,
                                   sample_name,
                                   coords_from = "umap",
                                   gene_set_path = NULL,
                                   verbose = TRUE){

# 1. Control --------------------------------------------------------------

  methods::is(object = seurat_object, class2 = "Seurata")

  confuns::is_value(x = method, mode = "character", ref = "method")
  confuns::check_one_of(input = method, against = seurat_methods, ref.input = "input for argument 'method'")

  confuns::is_value(x = sample_name, mode = "character", ref = "sample_name")


  # spata object
  spata_obj <- methods::new(Class = "spata")

  if(base::is.null(gene_set_path) | base::is.character(gene_set_path)){

    spata_obj@used_genesets <-
      loadGSDF(gene_set_path = gene_set_path, verbose = verbose)

  }


# 2. Extract data ---------------------------------------------------------

  if(method == "spatial"){

    # get coordinates
    coordinates <-
      getFromSeurat(
        return_value = Seurat::GetTissueCoordinates(seurat_object),
        error_handling = "stop",
        error_ref = "coordinates",
        error_value = NULL
      )

    coordinates <-
      tibble::rownames_to_column(.data = coordinates, var = "barcodes") %>%
      dplyr::rename("x" = "imagecol", "y" = "imagerow") %>%
      dplyr::mutate(
        sample = {{sample_name}},
        barcodes = stringr::str_c(barcodes, sample, sep = "_")) %>%
      dplyr::select(barcodes, sample, x, y)

    # get scaled matrix
    expr_matrix <-
      getFromSeurat(
        return_value = seurat_object@assays$Spatial@scale.data,
        error_handling = "stop",
        error_ref = "scaled matrix",
        error_value = NULL
      )

    # get count matrix
    count_matrix <-
      getFromSeurat(
        return_value = seurat_object@assays$Spatial@counts,
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
    coordinates <-
      getFromSeurat(
        return_value = base::as.data.frame(seurat_object@reductions[[first_choice]]@cell.embeddings),
        error_handling = "warning",
        error_value = NULL,
        error_ref = glue::glue("coordinates/{first_choice} cell embedding")
      )

    # try tsne if umap did not work
    if(base::is.null(coordinates)){

      base::warning(glue::glue("Trying to extract surrogate coordinates from slot {second_choice}."))

      coordinates <-
        getFromSeurat(
          return_value = base::as.data.frame(seurat_object@reductions[[second_choice]]@cell.embeddings),
          error_handling = "stop",
          error_value = NULL,
          error_ref = glue::glue("coordinates/{second_choice} cell embedding")
        )

    }

    coordinates <-
      tibble::rownames_to_column(.data = coordinates, var = "barcodes") %>%
      magrittr::set_colnames(value = c("barcodes", "x", "y")) %>%
      dplyr::mutate(sample = {{sample_name}}) %>%
      dplyr::select(barcodes, sample, x, y)

    # get expression matrix
    expr_matrix <-
      getFromSeurat(
        return_value = seurat_object@assays$RNA@scale.data,
        error_handling = "stop",
        error_value = NULL,
        error_ref = "scaled matrix"
      )

    # get count matrix
    count_matrix <-
      getFromSeurat(
        return_value = seurat_object@assays$RNA@counts,
        error_handling = "warning",
        error_value = base::matrix(),
        error_ref = "count matrix"
      )

    # no image
    image <- NULL

  }


# 3. Postprocess ----------------------------------------------------------

  # barcode check

  barcode_pattern <- stringr::str_c("_", sample_name, "$", sep = "")

  if(!base::all(stringr::str_detect(coordinates$barcodes, pattern = barcode_pattern))){

    coordinates <-
      dplyr::mutate(coordinates, barcodes = stringr::str_c(barcodes, {{sample_name}}, sep = "_"))

  }

  if(!base::all(stringr::str_detect(base::colnames(expr_matrix), pattern = barcode_pattern))){

    base::colnames(expr_matrix) <-
      stringr::str_c(base::colnames(expr_matrix), sample_name, sep = "_")

  }

  if(!base::all(stringr::str_detect(base::colnames(count_matrix), pattern = barcode_pattern))){

    base::colnames(count_matrix) <-
      stringr::str_c(base::colnames(count_matrix), sample_name, sep = "_")

  }


  # check if barcodes are identical
  barcodes_matrix <- base::colnames(expr_matrix) %>% base::sort()
  barcodes_coordinates <- dplyr::pull(coordinates, var = "barcodes") %>% base::sort()

  if(!base::identical(barcodes_matrix, barcodes_coordinates)){

    base::stop("The barcodes of the coordinate system and the column names of the assay must be identical. Please check the seurat object for integrity.")

  }

  # meta data
  fdata <-
    seurat_object@meta.data %>%
    tibble::rownames_to_column(var = "barcodes") %>%
    dplyr::mutate(
      sample = {{sample_name}},
      barcodes = stringr::str_c(barcodes, sample, sep = "_"),
      segment = ""
      ) %>%
    dplyr::select(barcodes, sample, dplyr::everything())

  # savely discard colum 'orig.ident'
  fdata <- base::tryCatch(

    dplyr::select(fdata, -orig.ident),

    error = function(error){ fdata }

  )

# 4. Pass to Spata --------------------------------------------------------

  # dimensional reduction
  dim_red_new <- methods::new("dim_red")

  dim_red_new@UMAP <- base::tryCatch(
    base::data.frame(
      barcodes = coordinates$barcodes,
      sample = coordinates$sample,
      umap1 = seurat_object@reductions$umap@cell.embeddings[,1],
      umap2 = seurat_object@reductions$umap@cell.embeddings[,2],
      stringsAsFactors = FALSE
    ) %>% tibble::remove_rownames() ,

    error = function(error){

      base::warning("Could not find or transfer UMAP-data. Did you process the seurat-object correctly?")

      base::return(data.frame())

    }

  )

  dim_red_new@TSNE <- base::tryCatch(

    base::data.frame(
      barcodes = coordinates$barcodes,
      sample = coordinates$sample,
      tsne1 = seurat_object@reductions$tsne@cell.embeddings[,1],
      tsne2 = seurat_object@reductions$tsne@cell.embeddings[,2],
      stringsAsFactors = F
    ) %>% tibble::remove_rownames()
    ,

    error = function(error){

      base::warning("Could not find or transfer TSNE-data. Did you process the seurat-object correctly?")

      base::return(data.frame())

    }

  )

  # data
  data <- methods::new(Class = "data_counts")

  data@counts <- count_matrix[base::rowSums(base::as.matrix(count_matrix)) != 0, ]
  data@norm_exp <- expr_matrix[base::rowSums(base::as.matrix(expr_matrix)) != 0, ]


  # passing
  spata_obj@coordinates <- coordinates

  spata_obj@samples <- sample_name

  spata_obj@fdata <- fdata

  spata_obj@data <- data

  spata_obj@dim_red <- dim_red_new

  spata_obj@trajectories <-
    magrittr::set_names(x = list(list()), value = sample_name)

  spata_obj@image <-
    magrittr::set_names(x = list(image), value = sample_name)


# 5. Return spata object ---------------------------------------------------

  base::return(spata_obj)

}
