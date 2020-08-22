#' @title Initiate a spata-object
#'
#' @description Creates, saves and returns an object of class spata
#' from scratch.
#'
#' @param input_paths Character vector. Specifies the 10X visium-folders from
#' which to load the information. These folders must contain the following sub directories:
#'
#' \itemize{
#'  \item{\emph{'/outs/filtered_feature_bc_matrix.h5'}}
#'  \item{\emph{'/outs/spatial/*.jpg}}
#'  }
#'
#' @param output_path Character value. Specifies the folder in which to store
#' the object.
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
#' @param sample_names Character vector. The sample names with which to refer to the
#' respective data. Must be of the same length as \code{input_paths}.
#' @param object_name Character value. The name-suffix for the file name under which the
#' spata-object is stored. Is prefixed with \emph{'spata-obj-'}.
#' @inherit verbose params
#'
#' @return A spata-object.
#' @export

initiateSpataObject_10X <- function(input_paths,
                                    output_path,
                                    gene_set_path = NULL,
                                    sample_names,
                                    object_name,
                                    verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  confuns::is_vec(input_paths, "character", ref = "input_path")

  confuns::is_vec(sample_names, "character", ref = "sample_names")
  if(base::any(sample_names %in% c("", "all"))){

    base::stop(glue::glue("' ' and 'all' are invalid sample names."))

  }

  confuns::is_value(output_path, "character", ref = "output_path")
  confuns::is_value(object_name, "character", ref = "object_name")

  confuns::check_directories(directories = input_paths, ref = "input_paths", type = "files")
  confuns::check_directories(directories = output_path, ref = "output_path", type = "files")

  object_file <- base::paste0(output_path, "/spata-obj-", object_name, ".RDS")

  if(base::file.exists(object_file)){

    base::stop(glue::glue("It already exists a .RDS-file ending with '{object_name}' in the directory '{output_path}'."))

  }

  if(base::length(sample_names) != base::length(input_paths)){

    base::stop("Lengths of arguments 'input_paths' and 'sample_names' must be equal.")

  } else if(base::length(sample_names) != base::length(base::unique(sample_names))){

    base::stop("Specified sample names must be unique.")

  }

  if(base::is.null(gene_set_path)){

    if(base::isTRUE(verbose)){base::message("No gene-set data.frame path specified.")}

  } else {

    confuns::check_directories(directories = gene_set_path, ref = "gene_set_path", type = "files")

  }

  # -----


  # 2. Read in data ---------------------------------------------------------

  if(base::isTRUE(verbose)){base::message("Step 1/6 : Reading in .h5 files.")}

  list_seurat_objects <-
    base::lapply(X = input_paths,
                 FUN = function(path){

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

                                 bc_with_sample_name <-
                                   stringr::str_c(barcodes, sample_names[i], sep = "_")

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
                   EBImage::Image(so@images$slice1[1]@image) %>%
                   EBImage::transpose()

                 base::return(image[,,1])

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

  test <- Seurat::CreateSeuratObject(counts = flt_counts)
  test[["percent.mt"]] <- Seurat::PercentageFeatureSet(test, pattern = "^MT.")
  test[["percent.RB"]] <- Seurat::PercentageFeatureSet(test, pattern = "^RPS")

  # remove mitochondrial genes and stress genes
  exclude <- c(rownames(test@assays$RNA)[base::grepl("^RPL", rownames(test@assays$RNA))],
               rownames(test@assays$RNA)[base::grepl("^RPS", rownames(test@assays$RNA))],
               rownames(test@assays$RNA)[base::grepl("^MT-", rownames(test@assays$RNA))],
               c('JUN','FOS','ZFP36','ATF3','HSPA1A","HSPA1B','DUSP1','EGR1','MALAT1'))

  feat_keep <- base::rownames(test@assays$RNA[!(base::rownames(test@assays$RNA) %in% exclude), ])

  test <- base::subset(x = test, features = feat_keep)
  test <- Seurat::NormalizeData(test, normalization.method = "LogNormalize", scale.factor = 1000)
  test <- Seurat::FindVariableFeatures(test, selection.method = "vst", nfeatures = 2000)

  test <- Seurat::ScaleData(test, features = base::rownames(test)) #"percent.mt", "percent.RB" vars.to.regress = c("Samples")

  if(base::isTRUE(verbose)){base::message("Running PCA...")}

  test <- Seurat::RunPCA(test, npcs = 60, features = Seurat::VariableFeatures(object = test))

  # clustering
  test <- Seurat::FindNeighbors(test, dims = 1:30)
  test <- Seurat::FindClusters(test, resolution = 0.8)
  test <- Seurat::RunUMAP(test, dims = 1:30, n.components = 2, n.neighbors = 50)

  # -----


  # 5. Create SPATA-object --------------------------------------------------

  if(base::isTRUE(verbose)){base::message("Step 5/6: Creating spata-object.")}

  mdata <-
    dplyr::left_join(x = tibble::rownames_to_column(test@meta.data, "barcodes"),
                     y = dplyr::select(fdata, barcodes, sample),
                     by = "barcodes")

  #Transform to S4 Class
  coords_n <- dplyr::select(.data = fdata, barcodes, sample, x, y)
  fdata_o <- dplyr::select(.data = fdata, sample)
  fdata_n <-
    dplyr::select(.data = mdata, barcodes, sample, dplyr::everything(), -orig.ident) %>%
    dplyr::mutate(segment = "")

  count_matrix <- test@assays$RNA@counts
  count_mtr <- count_matrix[base::rowSums(base::as.matrix(count_matrix)) != 0, ]

  norm_exp <- test@assays$RNA@scale.data
  norm_exp <- norm_exp[base::rowSums(norm_exp) != 0, ]

  data_counts_n <- new(Class = "data_counts",
                       counts = count_matrix,
                       norm_exp = base::as.matrix(norm_exp))

  dim_red_n <- new("dim_red",
                   UMAP = data.frame(
                     barcodes = fdata_n$barcodes,
                     sample = fdata_n$sample,
                     umap1 = test@reductions$umap@cell.embeddings[,1],
                     umap2 = test@reductions$umap@cell.embeddings[,2],
                     stringsAsFactors = F
                   ),
                   TSNE = data.frame(
                     barcodes = character(),
                     sample = character(),
                     tsne1 = numeric(),
                     tsne2 = numeric(),
                     stringsAsFactors = F
                   ))


  # get gene set data.frame
  if(!base::is.null(gene_set_path)){

    if(base::isTRUE(verbose)){base::message("Reading in gene set data.frame.")}

    gene_set_df <- base::readRDS(file = gene_set_path)

    if(!base::is.data.frame(gene_set_df)){

      gene_set_df <- data.frame()

      base::warning(glue::glue("Invalid input from directory '{gene_set_path}'. Returning empty geneset data.frame."))

    }

  } else {

    gene_set_df <- data.frame(ont = base::character(0), gene = base::character(0))

  }

  # trajectories
  trajectory_list <-vector(mode = "list", length = base::length(sample_names))
  base::names(trajectory_list) <- sample_names


  # final object
  new_object <- new(Class = "spata",
                    coordinates= coords_n,
                    fdata = fdata_n,
                    samples = sample_names,
                    data = data_counts_n,
                    image = list_images,
                    dim_red = dim_red_n,
                    trajectories = trajectory_list,
                    used_genesets = gene_set_df)

  # -----


  # 6. Save and return object -----------------------------------------------

  if(base::isTRUE(verbose)){base::message("Step 6/6: Saving spata-object.")}

  base::saveRDS(new_object, file = object_file)

  if(base::isTRUE(verbose)){

    base::message(glue::glue("SPATA-object has been saved under '{object_file}'."))

  }

  base::return(new_object)

}



#' @title Load and save a spata-object
#'
#' @description Wrapper around \code{base::readRDS()} and \code{base::saveRDS()}.
#'
#' @param input_path Character value. The directory leading to the spata-object.
#' @inherit check_object params
#' @inherit createSpataObject_10X params
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


















