#' @include S4-documentation.R
#'
NULL




#' @title The current spata-version
#' @export
spata_version <- base::list(major = 0,
                            minor = 0,
                            patch = 1,
                            dev = 0)


#' @title Initiate a spata-object
#'
#' @description Creates, saves and returns an object of class spata
#' from scratch. Several samples can be stored in one object. See details for more.
#'
#' @param input_paths Character vector. Specifies the 10X visium-folders from
#' which to load the information. This folder must contain the following sub directories:
#'
#' \itemize{
#'  \item{\emph{'/outs/filtered_feature_bc_matrix.h5'}}
#'  \item{\emph{'/outs/spatial/*.jpg}}
#'  }
#'
#' @param output_path Character vector or NULL. Specifies the folder in which to store
#' the object if the directory is valid.
#' @param gene_set_path Character value (or NULL). Specifies the path to a
#' .RDS-file containing a data.frame that is to be used as input for slot @@used_genesets.
#'
#'  Must have the character-variables
#'
#'    \itemize{
#'     \item{\emph{'ont'}: The gene set names.}
#'     \item{\emph{'gene'}: The belonging gene names.}
#'     }
#'
#' If set to NULL the default gene-set data.frame will used. Run \code{?gsdf} to get more information.
#'
#' @param sample_names Character vector. The sample name with which to refer to the
#' respective sample. Should start with a letter.
#'
#' @param file_name Character value. The name-suffix for the file name under which the
#' spata-object is stored. Is prefixed with \emph{'spata-obj-'}.
#' @inherit compileSeuratObject params
#' @inherit verbose params
#'
#' @details The loading and preprocessing of the spata-object  currently relies on the Seurat-package. For more advanced users the arguments
#' above starting with a capital letter allow to manipulate the way the spata-object is processed. For all of these arguments apply
#' the following instructions:
#'
#' \itemize{
#'   \item{If set to FALSE the processing function is skipped.}
#'   \item{If set to TRUE the respective function is called with it's default argument settings. Note: \code{RunUMAP()} needs
#'   additional input!}
#'   \item{If a named list is provided the respective function is called whereby the named list will provide the argument-input (the seurat-object to be constructed must not be named and will be
#'   passed to every function as the first argument!!!.)}
#'   }
#'
#' Note that certain listed functions require previous functions! E.g. if \code{RunPCA} is set to FALSE \code{RunTSNE()}
#' will result in an error. (\code{base::tryCatch()} will prevent the function from crashing but the respective slot
#' is going to be empty.) Skipping functions might result in an incomplete spata-object. Use \code{validateSpataObject()} after
#' initiating it in order to see which slots are valid and which are not.
#'
#' Handling more than one sample:
#'
#' Several samples can be stored in one object. If so, the count-matrices will be combined to one matrix which is given to the seurat-object that is temporarily
#' initiated in order to perform the pre processing steps. Sample related unambiguity with respect to the barcode's belonging is maintained
#' by suffixing the barcode-sequences with the respective sample name specified in \code{sample_names}. The meta.data data.frame of the
#' seurat-object is joined with a variable called \emph{sample} denoting the sample-belonging of every barcode which can be used as input
#' for pre processing functions.
#'
#' @return A spata-object.
#'
#' @importFrom Seurat ScaleData
#'
#' @export

initiateSpataObject_10X <- function(input_paths,
                                     sample_names,
                                     gene_set_path = NULL,
                                     output_path = NULL,
                                     file_name = NULL,
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

  confuns::is_vec(input_paths, "character", ref = "input_paths")
  confuns::is_vec(sample_names, "character", ref = "sample_names")

  if(base::any(sample_names %in% c("", "all"))){

    base::stop(glue::glue("' ' and 'all' are invalid sample names."))

  }

  if(base::length(sample_names) != base::length(input_paths)){

    base::stop("Lengths of arguments 'input_paths' and 'sample_names' must be equal.")

  } else if(base::length(sample_names) != base::length(base::unique(sample_names))){

    base::stop("Specified sample names must be unique.")

  }

  # saving if desired
  if(!base::is.null(output_path)){

    confuns::is_value(output_path, "character", ref = "output_path")
    confuns::check_directories(directories = input_paths, ref = "input_paths", type = "folders")

  }

  if(!base::is.null(output_path)){

    confuns::is_value(file_name, "character", ref = "file_name")
    confuns::check_directories(directories = output_path, ref = "output_path", type = "folders")

  }


  object_file <- base::paste0(output_path, "/spata-obj-", file_name, ".RDS")

  if(base::file.exists(object_file)){

    base::stop(glue::glue("It already exists a .RDS-file ending with '{file_name}' in the directory '{output_path}'."))

  }


  # gene set data.frame
  if(base::is.null(gene_set_path)){

    if(base::isTRUE(verbose)){base::message("No gene-set data.frame path specified.")}

  } else {

    confuns::check_directories(directories = gene_set_path, ref = "gene_set_path", type = "files")

  }

  # seurat processing
  for(fn in c("SCTransform", "NormalizeData", "FindVariableFeatures", "ScaleData",
              "RunPCA", "FindNeighbors", "FindClusters", "RunTSNE",  "RunUMAP")){

    input <- base::parse(text = fn) %>% base::eval()

    if(base::is.data.frame(input) | (!base::isTRUE(input) && !base::is.list(input) &&!base::isFALSE(input))){

      base::stop(glue::glue("Invalid input for argument '{fn}'. Must either be TRUE, FALSE or a named list."))

    }

  }


  # -----


  # 2. Read in data ---------------------------------------------------------

  if(base::isTRUE(verbose)){base::message("Step 1/6 : Reading in .h5 file.")}

  list_seurat_objects <-
    purrr::map(.x = input_paths,
               .f = function(path){

                   data_dir <- base::paste0(path, "/outs")
                   file_dir <- base::paste0(path, "/outs/filtered_feature_bc_matrix.h5")

                   if(base::file.exists(paths = file_dir)){

                     base::message(glue::glue("Loading from directory: '{data_dir}'"))

                     seurat <- Seurat::Load10X_Spatial(data.dir = data_dir, filename = "filtered_feature_bc_matrix.h5")

                     base::return(seurat)

                   } else {

                     base::message(glue::glue("Directory '{data_dir}' does not exist. Skip loading."))

                     base::return(NULL)

                   }

                 })

  # -----


  # 3. Merge count matrices -------------------------------------------------

  # list of expression matrices

  if(base::isTRUE(verbose)){base::message("Step 2/6: Extracting count matrices.")}

  gene_counts <- purrr::map(.x = list_seurat_objects,
                            .f = function(so){

                              base::as.matrix(so@assays$Spatial@counts)

                            })

  # vector of all genes
  all_genes <-
    base::lapply(X = gene_counts, FUN = function(mtr){ base::rownames(mtr) }) %>%
    base::unlist() %>%
    base::unique()

  # vector of all barcodes
  all_barcodes <- base::lapply(X = base::seq_along(list_seurat_objects),
                               FUN = function(i){

                                 barcodes <- base::colnames(list_seurat_objects[[i]])

                                 bc_with_sample_name <- stringr::str_c(barcodes, sample_names[i], sep = "_")

                                 base::return(bc_with_sample_name)

                               })

  all_barcodes <- base::unlist(all_barcodes)

  # matrix of counts
  all_counts <- base::matrix(0, nrow = base::length(all_genes), ncol = base::length(all_barcodes))

  base::rownames(all_counts) <- all_genes
  base::colnames(all_counts) <- all_barcodes

  for(i in base::seq_along(gene_counts)){

    rows <- base::rownames(gene_counts[[i]])
    cols <- stringr::str_c(base::colnames(gene_counts[[i]]), sample_names[i], sep = "_")

    all_counts[rows, cols] <- gene_counts[[i]]

  }

  # list of images
  if(base::isTRUE(verbose)){base::message("Step 3/6: Extracting images and coordinates-information.")}

  list_images <-
    purrr::map(.x = list_seurat_objects,
               .f = function(so){

                 image <-
                   EBImage::Image(so@images$slice1[1]@image, colormode = "Color") %>%
                   EBImage::transpose()

                 base::return(image)

               })

  base::names(list_images) <- sample_names

  # coordinate data.frame
  fdata_list <-
    purrr::map2(.x = list_seurat_objects,
                .y = sample_names,
                .f = function(so, sample) { # so = seurat object

                  barcodes <- base::rownames(so@meta.data)

                  coordinates <- Seurat::GetTissueCoordinates(object = so)[barcodes, ]
                  coordinates2 <- so@images$slice1@coordinates[barcodes, ]

                  fdata <- data.frame(HE_Slide = coordinates2$tissue+1,
                                      x = coordinates$imagecol,
                                      y = coordinates$imagerow,
                                      x_Slide = (coordinates2$col + 1),
                                      y_Slide = (coordinates2$row + 1),
                                      barcodes = stringr::str_c(barcodes, sample, sep = "_"),
                                      sample = sample,
                                      stringsAsFactors = FALSE)

                  return(fdata)

                })


  # 4. Seurat analysis ------------------------------------------------------

  if(base::isTRUE(verbose)){base::message("Step 4/6: Performing Seurat-analysis steps.")}

  flt_counts <- all_counts[base::rowSums(all_counts) != 0, ]
  fdata <- purrr::map_df(.x = fdata_list, .f = function(f){ base::return(f) })

  seurat_object <- Seurat::CreateSeuratObject(counts = flt_counts)
  seurat_object@meta.data$sample <- fdata$sample

  seurat_object[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_object, pattern = "^MT.")
  seurat_object[["percent.RB"]] <- Seurat::PercentageFeatureSet(seurat_object, pattern = "^RPS")

  # remove mitochondrial genes and stress genes
  exclude <- c(rownames(seurat_object@assays$RNA)[base::grepl("^RPL", rownames(seurat_object@assays$RNA))],
               rownames(seurat_object@assays$RNA)[base::grepl("^RPS", rownames(seurat_object@assays$RNA))],
               rownames(seurat_object@assays$RNA)[base::grepl("^MT-", rownames(seurat_object@assays$RNA))],
               c('JUN','FOS','ZFP36','ATF3','HSPA1A","HSPA1B','DUSP1','EGR1','MALAT1'))

  feat_keep <- base::rownames(seurat_object@assays$RNA[!(base::rownames(seurat_object@assays$RNA) %in% exclude), ])

  seurat_object <- base::subset(x = seurat_object, features = feat_keep)

  # perform pre processing steps
  functions_to_call <- c("SCTransform", "NormalizeData", "FindVariableFeatures", "ScaleData",
                         "RunPCA", "FindNeighbors", "FindClusters", "RunTSNE",  "RunUMAP")

  for(fn in functions_to_call){

    input <-
      base::parse(text = fn) %>%
      base::eval()

    if(base::isTRUE(input)){

      if(base::isTRUE(verbose)){base::message(glue::glue("Running 'Seurat::{fn}()' with it's default parameters."))}

      args <- base::list("object" = seurat_object)

      if(fn == "ScaleData"){

        args <- base::append(x = args,
                             values = list("features" = base::rownames(seurat_object)))

      }

      # ensure that function is called from Seurat-namespace
      fn <- stringr::str_c("Seurat::", fn, sep = "")

      seurat_object <- base::tryCatch(

        rlang::invoke(.fn = base::eval(base::parse(text = fn)), args),

        error = function(error){

          base::message(glue::glue("Running'Seurat::{fn}()' resulted in the following error: {error$message}. Abort and continue with next function."))

          base::return(seurat_object)

        })

    } else if(base::is.list(input) &
              !base::is.data.frame(input)){

      input_content <- base::names(input)
      keep <- !input_content %in% c("", "object")

      input <- input[keep]
      named_arguments <- input_content[keep]

      if(base::isTRUE(verbose)){

        ref_named_arguments <- stringr::str_c(named_arguments, collapse = "', '")

        base::message(glue::glue("Running 'Seurat::{fn}()' with specified parameters: '{ref_named_arguments}'"))

      }

      args <- purrr::prepend(x = input, values = seurat_object)

      if(fn == "ScaleData" && !"features" %in% base::names(args)){

        args <- base::append(x = args,
                             values = list("features" = base::rownames(seurat_object)))

      }
      # ensure that function is called from Seurat-namespace
      fn <- stringr::str_c("Seurat::", fn, sep = "")

      seurat_object <- base::tryCatch(

        rlang::invoke(.fn = base::eval(base::parse(text = fn)), args),

        error = function(error){

          base::message(glue::glue("Running'Seurat::{fn}' resulted in the following error: {error$message}. Abort and continue with next function."))

          base::return(seurat_object)

        }

      )

    } else {

      if(base::isTRUE(verbose)){base::message(glue::glue("Skip running {fn} as it's argument input is neither TRUE nor a list."))}

    }

  }
  # -----


  # 5. Create SPATA-object --------------------------------------------------

  if(base::isTRUE(verbose)){base::message("Step 5/6: Initiating spata-object.")}

  mdata <-
    dplyr::left_join(x = tibble::rownames_to_column(seurat_object@meta.data, "barcodes"),
                     y = dplyr::select(fdata, barcodes),
                     by = "barcodes")

  #Transform to S4 Class
  coords_n <- dplyr::select(.data = fdata, barcodes, sample, x, y)
  fdata_o <- dplyr::select(.data = fdata, sample)
  fdata_n <-
    dplyr::select(.data = mdata, barcodes, sample, dplyr::everything(), -orig.ident) %>%
    dplyr::mutate(segment = "")

  count_matrix <- seurat_object@assays$RNA@counts
  count_mtr <- count_matrix[base::rowSums(base::as.matrix(count_matrix)) != 0, ]

  norm_exp <- seurat_object@assays$RNA@scale.data
  norm_exp <- norm_exp[base::rowSums(norm_exp) != 0, ]

  data_counts_n <- new(Class = "data_counts",
                       counts = count_matrix,
                       norm_exp = base::as.matrix(norm_exp))


  dim_red_n <- methods::new("dim_red")

  dim_red_n@UMAP <- base::tryCatch(
    base::data.frame(
      barcodes = fdata_n$barcodes,
      sample = fdata_n$sample,
      umap1 = seurat_object@reductions$umap@cell.embeddings[,1],
      umap2 = seurat_object@reductions$umap@cell.embeddings[,2],
      stringsAsFactors = F
    ) %>% tibble::remove_rownames() ,

    error = function(error){

      base::warning("Could not find or transfer UMAP-data. Did you set up the respective Seurat-functions correctly?")

      base::return(data.frame())

    }

  )

  dim_red_n@TSNE <- base::tryCatch(

    base::data.frame(
      barcodes = fdata_n$barcodes,
      sample = fdata_n$sample,
      tsne1 = seurat_object@reductions$tsne@cell.embeddings[,1],
      tsne2 = seurat_object@reductions$tsne@cell.embeddings[,2],
      stringsAsFactors = F
    ) %>% tibble::remove_rownames()
  ,

    error = function(error){

      base::warning("Could not find or transfer TSNE-data. Did you set up the respective Seurat-functions correctly?")

      base::return(data.frame())

    }

  )


  # get gene set data.frame
  if(!base::is.null(gene_set_path)){

    if(base::isTRUE(verbose)){base::message(glue::glue("Reading in specified gene-set data.frame from directory '{gene_set_path}'."))}

    gene_set_df <- base::readRDS(file = gene_set_path)

    if(!base::is.data.frame(gene_set_df)){

      gene_set_df <- data.frame()

      base::warning(glue::glue("Invalid input from directory '{gene_set_path}'. Returning empty geneset data.frame."))

    }

  } else {

    if(base::isTRUE(verbose)){base::message("Using SPATA's default gene set data.frame.")}

    gene_set_df <- gsdf

  }

  # trajectories
  trajectory_list <-vector(mode = "list", length = base::length(sample_names))
  base::names(trajectory_list) <- sample_names

  #additional
  seurat_image_objects <- purrr::map(.x = list_seurat_objects,
                                .f = function(so){

                                  so@images$slice1

                                }) %>% magrittr::set_names(sample_names)

  additional_list <- list(Seurat = list(images = seurat_image_objects))

  # final object
  new_object <- new(Class = "spata",
                    coordinates= coords_n,
                    dim_red = dim_red_n,
                    data = data_counts_n,
                    fdata = fdata_n,
                    image = list_images,
                    samples = sample_names,
                    trajectories = trajectory_list,
                    used_genesets = gene_set_df,
                    version = spata_version,
                    additional = additional_list)

  # -----


  # 6. Save and return object -----------------------------------------------

  if(!base::is.null(output_path)){

    if(base::isTRUE(verbose)){base::message("Step 6/6: Saving spata-object.")}

    base::saveRDS(new_object, file = object_file)

    if(base::isTRUE(verbose)){

      base::message(glue::glue("The spata-object has been saved under '{object_file}'."))
      base::message("Done.")

    }

  } else {

    if(base::isTRUE(verbose)){
      base::message("Skipping step 6/6 as 'output_path' was set to NULL.")
      base::message("Done.")
      }

  }

  base::return(new_object)

}



#' @title Load and save a spata-object
#'
#' @description Wrapper around \code{base::readRDS()} and \code{base::saveRDS()}.
#'
#' @param input_path Character value. The directory leading to the spata-object.
#' @inherit check_object params
#' @inherit initiateSpataObject_10X params
#' @param overwrite Logical. Needs to be set to TRUE if the resulting directory from
#' \code{output_path} and \code{object_name} already exists.
#'
#' @export

loadSpataObject <- function(input_path){

  confuns::is_value(input_path, "character", "input_path")
  confuns::check_directories(directories = input_path, ref = "input_path", type = "files")

  spata_obj <- base::readRDS(file = input_path)

  if(!methods::is(spata_obj, "spata")){

    base::warning("Return object is not of class 'spata'!")

  }

  base::return(spata_obj)

}


#' @rdname loadSpataObject
#' @export
saveSpataObject <- function(object, output_path, object_name, overwrite = FALSE){

  # 1. Control --------------------------------------------------------------

  check_object(object)
  confuns::is_value(output_path, "character", "output_path")
  confuns::is_value(object_name, "character", "object_name")
  confuns::check_directories(output_path, ref = "output_path", type = "folders")

  # -----

  filename <- stringr::str_c(output_path, "/spata-obj-", object_name, ".RDS", sep = "")

  if(base::file.exists(filename) && !base::isTRUE(overwrite)){

    base::stop(glue::glue("The file '{filename}' already exists. Set argument 'overwrite' to TRUE in order to overwrite."))

  } else if(base::file.exists(filename) && base::isTRUE(overwrite)){

    base::message(glue::glue("Argument 'overwrite' set to TRUE - overwriting {filename} with input for argument 'object'."))

    base::file.remove(filename)

    base::saveRDS(object = object, file = filename)

  } else if(!base::file.exists(filename)){

    base::message(glue::glue("Saving object under '{filename}'."))
    base::saveRDS(object = object, file = filename)

  }

  if(base::file.exists(filename)){

    base::message("Saving successful.")
    base::return(base::invisible(TRUE))

  } else {

    base::warning("Saving failed.")
    base::return(base::invisible(FALSE))

  }

}



#' @title Save a gene set data.frame
#'
#' @description Extracts the gene-set data.frame and saves it as a .RDS-file.
#'
#' @inherit check_object params
#' @param output_path Character value. A directory leading to the folder in which
#' to store the data.frame.
#' @param filename Character value. The filename. ( \emph{'.RDS'} is attached automatically.)
#'
#' @return An invisible TRUE if saved successfully or an informative error message.
#' @export
#'

saveGeneSetDf <- function(object,
                          output_path,
                          filename){

  check_object(object)

  confuns::is_value(output_path, "character", "output_path")
  confuns::is_value(filename, "character", "filename")

  confuns::check_directories(output_path, "output_path", type = "folders")

  final_path <- stringr::str_c(output_path, "/", filename, ".RDS", sep = "")

  if(base::file.exists(final_path)){

    base::stop(glue::glue("The file '{final_path}' already exists."))

  } else if(base::nrow(object@used_genesets) == 0){

    base::stop("The objects's gene-set data.frame is empty.")

  } else {

    base::saveRDS(object = object@used_genesets, file = final_path)

    if(base::file.exists(final_path)){

      object_name <- stringr::str_c("~/", filename, ".RDS", sep = "")
      base::message(glue::glue("Gene set data.frame has been saved as '{object_name}'."))
      base::return(base::invisible(TRUE))

    } else {

      base::stop("Could not save the gene-set data.frame. Unknown error.")

    }

  }

}


















