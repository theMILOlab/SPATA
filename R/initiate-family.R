


# Test functions --------------------------------------------------------

#' @title Initiate a spata-object from MALDI Experiments
#'
#' @inherit initiateSpataObject_ExprMtr params
#' @param intensity_mtr A numeric matrix to be used as the expression matrix. Rownames must
#' correspond to the genes and column names must correspond to the barcodes.
#'
#' @return A spata-object.
#'

initiateSpataObject_MALDI <- function(coords_df,
                                      intensity_mtr,
                                      sample_name,
                                      gene_set_path = NULL,
                                      output_path = NULL,
                                      file_name = NULL,
                                      pca_comp = 30,
                                      nn = 50,
                                      tsne_perplexity = 30,
                                      verbose = TRUE){

  initiateSpataObject_ExprMtr(coords_df = coords_df,
                           expr_mtr = intensity_mtr,
                           ref_expr_mtr = "intensity_mtr",
                           image = image,
                           sample_name = sample_name,
                           gene_set_path = gene_set_path,
                           output_path = output_path,
                           file_name = file_name,
                           pca_comp = pca_comp,
                           nn = nn,
                           tsne_perplexity = tsne_perplexity,
                           verbose = verbose)

}


# Exported functions ------------------------------------------------------

#' @title Initiate a spata-object from a raw count matrix
#'
#' @description Default function for any spatial related experiment whoose spata-objects are initiated with
#' a raw count matrix. See details for more information.
#'
#' @param count_mtr A numeric matrix to be used as the count matrix. Rownames must
#' correspond to the genes and column names must correspond to the barcodes.
#' @inherit initiateSpataObject_ExprMtr params return
#' @inherit compileSeuratObject params
#'
#' @details The loading and preprocessing of the spata-object  currently relies on the Seurat-package. Before any pre processing function is applied
#' mitochondrial and stress genes are discarded. For more advanced users the arguments above starting with a capital letter allow to
#' manipulate the way the spata-object is processed. For all of these arguments apply the following instructions:
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
#' is going to be empty.) Skipping functions might result in an incomplete spata-object. Use \code{validateSpataObject()}
#' to check your object for validity.
#'
#' @export

initiateSpataObject_CountMtr <- function(coords_df,
                                         count_mtr,
                                         feature_df = NULL,
                                         sample_name,
                                         image = NULL,
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

    confuns::is_value(x = sample_name, mode = "character", ref = "sample_name")

    confuns::check_data_frame(
      df = coords_df,
      var.class = list("barcodes" = "character", "x" = c("double", "integer", "numeric"), "y" = c("double", "integer", "numeric")),
      ref = "coords_df"
    )

    if(!methods::is(count_mtr, "Matrix")){

      base::stop("'count_mtr'-input needs to be of type 'Matrix'.")

    } else if(base::is.null(base::colnames(count_mtr))){

      base::stop("'count_mtr'-input needs to have column names")

    } else if(base::is.null(base::rownames(count_mtr))){

      base::stop("'count_mtr'-input needs to have row names.")

    }

    barcodes_count_mtr <- base::colnames(count_mtr) %>% base::sort()
    barcodes_coords_df <- dplyr::pull(coords_df, var = "barcodes") %>% base::sort()

    # check identical barcodes
    if(!base::identical(barcodes_count_mtr, barcodes_coords_df)){

      base::stop("Barcodes of 'coords_df'-input and column names of 'count_mtr'-input need to be identical.")

    }


    # -----

    # 2. Passing data ---------------------------------------------------------

    counts <- count_mtr

    seurat_object <-
      Seurat::CreateSeuratObject(counts = counts, meta.data = feature_df)

    seurat_object <-
      process_seurat_object(
        seurat_object = seurat_object,
        SCTransform = SCTransform,
        NormalizeData = NormalizeData,
        FindVariableFeatures = FindVariableFeatures,
        ScaleData = ScaleData,
        RunPCA = RunPCA,
        FindNeighbors = FindNeighbors,
        FindClusters = FindClusters,
        RunTSNE = RunTSNE,
        RunUMAP = RunUMAP,
        verbose = verbose
      )


    # Passing features and images ---------------------------------------------

    spata_object <-
      transformSeuratToSpata(
        seurat_object = seurat_object,
        assay = "RNA",
        sample_name = sample_name,
        gene_set_path = gene_set_path,
        method = "single_cell",
        verbose = verbose
      )

    spata_object@coordinates <- coords_df
    spata_object@image[[sample_name]] <- image

    spata_object <- check_all_barcodes(spata_object)

    # Save and return object -----------------------------------------------

    if(!base::is.null(output_path)){

      if(base::isTRUE(verbose)){base::message("Saving spata-object.")}

      base::saveRDS(spata_object, file = object_file)

      if(base::isTRUE(verbose)){

        base::message(glue::glue("The spata-object has been saved under '{object_file}'."))
        base::message("Done.")

      }

    }

    if(base::isTRUE(verbose)){base::message("Done.")}

    return(spata_object)

  }


#' @title Initiate spata object from scaled expression matrix
#'
#' @description Default function for any spatial related experiment whoose output is
#' an already processed expression/intensity matrix. See details for more information.
#'
#' @param coords_df Data.frame containing information about the positions of all
#' barcode-spots in form of a numeric \emph{x}- and \emph{y}-variable. The key-variable
#' \emph{barcodes} needs to be of type character and must be identical to the column names
#' of the input matrix (\code{expr_mtr}).
#'
#' @param expr_mtr A numeric matrix. The expression matrix to be used.
#' @param image An Image of the sample that can be displayed as the surface plot's background.
#' @inherit sample_name params
#' @inherit gene_set_path params
#' @inherit check_saving params
#' @param pca_comp Numeric value. Given to argument \code{n} of function \code{irlba::prcomp_irlba()}: Determines
#' the number of components to return.
#' @param nn Numeric value. Given to argument \code{k} of function \code{RANN::nn2()}: Determines to maximum number
#' of nearest neighbours to compute.
#' @param tsne_perplexity Numeric value. Given to argument \code{perplexity} of function \code{Rtsne::Rtsne()}.
#' @inherit verbose params
#'
#' @details After initiating the spata-object PCA is performed via \code{irlba::prcomp_irlba()} and clustering
#' is done via \code{RANN::nn2()}. (Use \code{addFeatures()} to add any clustering results of your own analysis.)
#' Additional dimensional reduction is performed via \code{Rtsne::Rtsne()} and \code{umap::umap()}.
#'
#' @return A spata-object.
#'
#' @export

initiateSpataObject_ExprMtr <- function(coords_df,
                                        expr_mtr,
                                        count_mtr = NULL,
                                        image = NULL,
                                        sample_name,
                                        gene_set_path = NULL,
                                        output_path = NULL,
                                        file_name = NULL,
                                        pca_comp = 30,
                                        nn = 50,
                                        tsne_perplexity = 30,
                                        verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  if(base::isTRUE(verbose)){base::message("Step 1/4: Checking input for validity.")}

  # check if expr matrix is a matrix
  if(!base::is.matrix(expr_mtr) | base::length(base::dim(expr_mtr)) != 2){

    base::stop(glue::glue("Input for argument 'expr_mtr' must be a matrix."))

  }

  # check if coordinate data.frame is valid
  confuns::check_data_frame(
    df = coords_df,
    var.class = list(barcodes = c("character"),
                     x = c("numeric", "double", "integer"),
                     y = c("numeric", "double", "integer")),
    ref = "coords_df"
  )

  # check value inputs
  confuns::is_value(x = sample_name, mode = "character", ref = "sample_name")

  if(!base::is.null(output_path)){

    confuns::are_values("ouptut_path", "file_name", mode = "character")

  }

  confuns::are_values("pca_comp", "nn", "tsne_perplexity", mode = "numeric")

  # check if gene names are valid
  if(base::any(stringr::str_detect(base::rownames(expr_mtr), pattern = "_"))){

    base::message("Rownames of expression matrix contain '_' which is not allowed and will be changed into '-'.")

    base::rownames(expr_mtr) <-
      stringr::str_replace_all(string = base::rownames(expr_mtr), pattern = "_", replacement = "-")

  }

  # check if barcodes are valid and idential

  expr_mtr <- hlpr_add_barcode_suffix(input = expr_mtr, sample_name = sample_name)
  coords_df <- hlpr_add_barcode_suffix(input = coords_df, sample_name = sample_name)

  barcodes_coords <- dplyr::pull(.data = coords_df, var = "barcodes") %>% base::sort()
  barcodes_expr_mtr <- base::colnames(expr_mtr) %>% base::sort()

  if(!base::identical(barcodes_coords, barcodes_expr_mtr)){

    base::stop("Barcodes of expression matrix and barcodes of the coordinate data.frame must match.")

  }

  # check is sample column exist
  if(!"sample" %in% base::colnames(coords_df)){

    base::message("Adding sample variable to coordinate data.frame.")

    coords_df <-
      dplyr::mutate(.data = coords_df,
                    sample = {{sample_name}}
      )

  }

  coords_df <- dplyr::select(.data = coords_df, barcodes, sample, x, y)

  # check if saving is possible

  object_file <- check_saving(output_path = output_path, file_name = file_name)


  # -----


  # 2. Setting up spata object ----------------------------------------------

  if(base::isTRUE(verbose)){ base::message("Step 2/4:Setting up spata-object.")}

  spata_obj <- methods::new(Class = "spata")

  spata_obj@samples <- sample_name

  # transfer assays
  data <- list()

  data$counts <- count_mtr
  data$scaled <- expr_mtr

  spata_obj@data <- data

  # transfer data.frames
  spata_obj@coordinates <- coords_df
  spata_obj@fdata <- dplyr::select(.data = coords_df, barcodes, sample) %>% dplyr::mutate(segment = "")
  spata_obj@used_genesets <- loadGSDF(gene_set_path = gene_set_path, verbose = verbose)

  # transfer lists
  spata_obj@trajectories <- purrr::set_names(x = list(list()), nm = sample_name)
  spata_obj@image <- purrr::set_names(x = list(image), nm = sample_name)
  spata_obj@information <- list("initiation" = list("timepoint" = base::Sys.time()))

  # -----


  # 3. Running analysis -----------------------------------------------------

  if(base::isTRUE(verbose)){ base::message("Step 3/4: Running analysis steps.")}

  if(base::isTRUE(verbose)){ base::message("---------- PCA Analysis 1/4 ------------- ")}

  pca_res <- irlba::prcomp_irlba((expr_mtr), n = pca_comp)

  pca <- pca_res[["rotation"]]

  base::rownames(pca) <- spata_obj@fdata$barcodes

  if(base::isTRUE(verbose)){ base::message("---------- Select Eigenvalues Analysis 2/4 ------------- ") }

  # missing?

  if(base::isTRUE(verbose)){ base::message("---------- SNN-Cluster Analysis 3/4 ------------- ") }

  nearest <- RANN::nn2(pca, k = nn, treetype = "bd", searchtype = "priority")

  edges <-
    reshape::melt(base::t(nearest$nn.idx[, 1:nn])) %>%
    dplyr::select(A = X2, B = value) %>%
    dplyr::mutate(C = 1)

  edges <-
    base::transform(edges, A = pmin(A, B), B = pmax(A, B)) %>%
    base::unique() %>%
    dplyr::rename(V1 = A, V2 = B, weight = C)

  edges$V1 <- base::rownames(pca)[edges$V1]
  edges$V2 <- base::rownames(pca)[edges$V2]

  g_df <- igraph::graph.data.frame(edges, directed = FALSE)

  graph.out <- igraph::cluster_louvain(g_df)

  clust.assign <- base::factor(x = graph.out$membership,
                               levels = base::sort(base::unique(graph.out$membership)))

  cluster_out <-
    base::data.frame(barcodes = base::rownames(pca), nn2_cluster = clust.assign) %>%
    magrittr::set_rownames(value = base::rownames(pca)) %>%
    dplyr::mutate(nn2_cluster = base::factor(nn2_cluster))

  spata_obj <- addFeatures(object = spata_obj,
                     feature_df = cluster_out,
                     feature_names = c("nn2_cluster"),
                     key_variable = "barcodes",
                     of_sample = sample_name)

  if(base::isTRUE(verbose)){ base::message("---------- Dimensional Reduction 4/4 ------------- ")}

  pca_df <-
    tibble::rownames_to_column(.data = base::as.data.frame(pca), var = "barcodes") %>%
    dplyr::mutate(sample = {{sample_name}}) %>%
    dplyr::select(barcodes, sample, dplyr::everything())

  dim_red_new <- list("pca" = pca_df)

  if(base::isTRUE(verbose)){ base::message(glue::glue("Running TSNE with argument 'perplexity' = {tsne_perplexity}."))}

  ts <- Rtsne::Rtsne(pca, perplexity = tsne_perplexity)

  dim_red_new$tsne <-
    base::data.frame(barcodes = base::rownames(pca),
                     sample = sample_name,
                     tsne1 = ts$Y[,1],
                     tsne2 = ts$Y[,2])


  if(base::isTRUE(verbose)){ base::message("Running UMAP.")}

  umap <- umap::umap(pca)

  dim_red_new$umap <-
    base::data.frame(
      barcodes = base::rownames(pca),
      sample = sample_name,
      umap1 = umap$layout[,1],
      umap2 = umap$layout[,2]
    )

  # transfer dimensional reduction
  spata_obj@dim_red <- dim_red_new

  spata_obj <- check_all_barcodes(spata_obj)

  # -----

  # 4. Saving and returning object ------------------------------------------

  hlpr_save_spata_object(object = spata_obj,
                         object_file = object_file,
                         ref_step = "4/4",
                         verbose = verbose)

  return(spata_obj)

}

#' @title Initiate a spata-object from 10X Visium
#'
#' @description Creates, saves and returns an object of class spata
#' from scratch. Several samples can be stored in one object, though we recommend to stick
#' to one. (See details for more.)
#'
#' @param input_paths Character vector. Specifies the 10X visium-folders from
#' which to load the information. This folder must contain the following sub directories:
#'
#' \itemize{
#'  \item{\emph{'/outs/filtered_feature_bc_matrix.h5'}}
#'  \item{\emph{'/outs/spatial/*.jpg}}
#'  }
#'
#' @inherit gene_set_path params
#'
#' @param sample_names Character vector. The sample name with which to refer to the
#' respective sample. Should start with a letter.
#'

#' @inherit compileSeuratObject params
#' @inherit check_saving params
#' @inherit verbose params
#'
#' @details The loading and preprocessing of the spata-object  currently relies on the Seurat-package. Before any pre processing function is applied
#' mitochondrial and stress genes are discarded. For more advanced users the arguments above starting with a capital letter allow to
#' manipulate the way the spata-object is processed. For all of these arguments apply the following instructions:
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
#' For now we recommend to stick to one sample per spata-object.
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


  confuns::check_directories(directories = input_paths, type = "folders")

  # check saving if desired
  object_file <- check_saving(output_path = output_path, file_name = file_name)


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

  if(base::isTRUE(verbose)){base::message("Step 1/4 : Reading in .h5 file.")}

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



  # 3. Seurat analysis ------------------------------------------------------

  if(base::isTRUE(verbose)){base::message("Step 2/4: Performing Seurat-analysis steps.")}

  processed_seurat_objects <-
    purrr::map(
      .x = list_seurat_objects,
      .f = process_seurat_object,
      calculate_rb_and_mt = TRUE,
      remove_stress_and_mt = TRUE,
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

  list_seurat_objects <- NULL

  # -----


  # 5. Create SPATA-object --------------------------------------------------

  if(base::isTRUE(verbose)){base::message("Step 3/4: Initiating spata-object.")}

  spata_objects <-
    purrr::map2(
      .x = processed_seurat_objects,
      .y = sample_names,
      .f = transformSeuratToSpata,
      gene_set_path = NA,
      method = "spatial",
      verbose = verbose
    )

  #processed_seurat_objects <- NULL

  if(base::length(spata_objects) > 1){

    spata_object <-
      mergeSpataObjects(spata_objects,
                        gsdf_input = gene_set_path,
                        verbose = verbose)

  } else {

    spata_object <- spata_objects[[1]]

    spata_object@used_genesets <-
      loadGSDF(gene_set_path = gene_set_path, verbose = verbose)

  }

  spata_objects <- NULL


  # -----


  # 6. Save and return object -----------------------------------------------

  hlpr_save_spata_object(object = spata_object,
                         object_file = object_file,
                         ref_step = "4/4",
                         verbose = verbose)

  base::return(spata_object)

}
