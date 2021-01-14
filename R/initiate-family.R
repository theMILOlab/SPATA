


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
                                         directory_spata = NULL,
                                         directory_seurat = NULL,
                                         combine_with_wd = FALSE,
                                         gene_set_path = NULL,
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
    confuns::give_feedback(
      msg = "Starting initiation",
      verbose = verbose
    )

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

    spata_object <-
      setCoordsDf(object = spata_object, coords_df = coords_df) %>%
      setImage(object = ., image = image)

    spata_object <- setInitiationInfo(spata_object)

    # Save and return object -----------------------------------------------

    # save seurat
    if(base::is.character(directory_seurat)){

      spata_object <-
        base::tryCatch({

          spata_object_return <-
            saveCorrespondingSeuratObject(
              seurat_object = seurat_object,
              object = spata_object,
              directory_seurat = directory_seurat,
              combine_with_wd = combine_with_wd
            )

          spata_object_return

        }, error = function(error){

          base::warning(glue::glue("Attempt to save seurat-object under '{directory_seurat}' failed with the following error message: {error}"))

          spata_object

        })

    }

    # save spata
    if(base::is.character(directory_spata)){

      spata_object <-
        base::tryCatch({

          spata_object_ret <-
            saveSpataObject(
              object = spata_object,
              directory_spata = directory_spata,
              combine_with_wd = combine_with_wd
            )

          spata_object_ret

        }, error = function(error){

          base::warning(glue::glue("Attempt to save spata-object under '{directory_spata}' failed with the following error message: {error}"))

          spata_object

        })

    }

    confuns::give_feedback(
      msg = "Initiation finished.",
      verbose = verbose
    )

    # -----

    base::return(spata_object)

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
                                        mtr_name = "scaled",
                                        sample_name = "sample1",
                                        image = NULL,
                                        directory_spata = NULL,
                                        combine_with_wd = FALSE,
                                        gene_set_path = NULL,
                                        n_pcs = 30,
                                        k = 50,
                                        tsne_perplexity = 30,
                                        ...,
                                        verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  confuns::give_feedback(
    msg = "Starting initiation.",
    verbose = verbose
  )

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
  confuns::is_value(x = sample_name, mode = "character")

  confuns::are_values("output_path", "file_name", mode = "character", skip.allow = TRUE, skip.val = NULL)

  confuns::are_values("n_pcs", "k", "tsne_perplexity", mode = "numeric")

  # check if gene names are valid
  if(base::any(stringr::str_detect(base::rownames(expr_mtr), pattern = "_"))){

    base::message("Rownames of expression matrix contain '_' which is not allowed and will be changed into '-'.")

    base::rownames(expr_mtr) <-
      stringr::str_replace_all(string = base::rownames(expr_mtr), pattern = "_", replacement = "-")

  }

  # check if barcodes are valid and idential
  barcodes_coords <- dplyr::pull(.data = coords_df, var = "barcodes") %>% base::sort()
  barcodes_expr_mtr <- base::colnames(expr_mtr) %>% base::sort()

  if(!base::identical(barcodes_coords, barcodes_expr_mtr)){

    base::stop("Barcodes of expression matrix and barcodes of the coordinate data.frame must match.")

  }

  coords_df <-
    dplyr::mutate(.data = coords_df, sample = {{sample_name}}) %>%
    dplyr::select(barcodes, sample, x, y)


  # check if saving is possible (returns NULL if output path is NULL)
  object_file <- check_saving(output_path = output_path, file_name = file_name)


  # -----


  # 2. Setting up spata object ----------------------------------------------

  if(base::isTRUE(verbose)){ base::message("Step 2/4:Setting up spata-object.")}

  spata_object <- methods::new(Class = "spata", samples = sample_name)

  # data matrices

  if(base::is.matrix(count_mtr)){

    spata_object <-
      setCountMatrix(
        object = spata_object,
        count_mtr = count_mtr[base::rowSums(base::as.matrix(count_mtr)) != 0, ]
      )

  }

  spata_object <-
    addExpressionMatrix(
      object = spata_object,
      expr_mtr = expr_mtr[base::rowSums(base::as.matrix(expr_mtr)) != 0, ],
      mtr_name = mtr_name
    )

  spata_object <- setActiveExpressionMatrix(object = spata_object, mtr_name = mtr_name)

  # transfer data.frames and image
  feature_df <-
    dplyr::mutate(.data = coords_df, segment = "none") %>%
    dplyr::select(barcodes, sample, segment)

  gene_set_df <-
    loadGSDF(gene_set_path = gene_set_path, verbose = verbose)

  spata_object <-
    setCoordsDf(object = spata_object, coords_df = coords_df) %>%
    setFeatureDf(object = ., feature_df = feature_df) %>%
    setGeneSetDf(object = ., gene_set_df = gene_set_df) %>%
    setImage(object = ., image = image)


  # transfer lists
  spata_object@trajectories <-
    purrr::set_names(x = list(list()), nm = sample_name)

  spata_object@information <-
    base::append(
      x = spata_object@information,
      values = list(
        autoencoder = magrittr::set_names(x = list(list()), nm = sample_name),
        barcodes = magrittr::set_names(x = list(base::colnames(expr_mtr)), nm = sample_name)
      ))



  # -----


  # 3. Running analysis -----------------------------------------------------

  if(base::isTRUE(verbose)){ base::message("Step 3/4: Running analysis steps.")}

  if(base::isTRUE(verbose)){ base::message("---------- PCA Analysis 1/4 ------------- ")}

  spata_object <- runPca(object = spata_object, n_pcs = n_pcs)

  if(base::isTRUE(verbose)){ base::message("---------- SNN-Cluster Analysis 2/3 ------------- ") }

  cluster_df <-
    findNearestNeighbourClusters(object = spata_object,
                                 k = k,
                                 searchtype = "priority",
                                 treetype = "bd",
                                 eps = 0,
                                 radius = 0)

  feature_names <-
    dplyr::select_if(cluster_df, .predicate = base::is.factor) %>%
    base::colnames()

  spata_object <- addFeatures(object = spata_object,
                              feature_df = cluster_df,
                              feature_names = feature_names,
                              key_variable = "barcodes")


  if(base::isTRUE(verbose)){ base::message("---------- Dimensional Reduction 3/3 ------------- ")}

  if(base::isTRUE(verbose)){ base::message("Running TSNE.") }

  spata_object <- runTsne(object = spata_object, tsne_perplexity = tsne_perplexity)


  if(base::isTRUE(verbose)){ base::message("Running UMAP.") }

  spata_object <- runUmap(object = spata_object, ...)


  # -----

  # 4. Saving and returning object ------------------------------------------

  if(base::is.character(directory_spata)){

    spata_object <-
      base::tryCatch({

        spata_object_ret <-
          saveSpataObject(
            object = spata_object,
            directory_spata = directory_spata,
            combine_with_wd = combine_with_wd
          )

        spata_object_ret

      }, error = function(error){

        base::warning(glue::glue("Attempt to save spata-object under '{directory_spata}' failed with the following error message: {error}"))

        spata_object

      })

  }

  confuns::give_feedback(
    msg = "Initiation finished.",
    verbose = verbose
  )

  # -----

  base::return(setInitiationInfo(spata_object))

}

#' @title Initiate a spata-object from 10X Visium
#'
#' @description Creates, saves and returns an object of class spata
#' from scratch. Several samples can be stored in one object, though we recommend to stick
#' to one. (See details for more.)
#'
#' @param directory_10X Character value. Specifies the 10X visium-folder from
#' which to load the information. This folder must contain the following sub directories:
#'
#' \itemize{
#'  \item{\emph{'/outs/filtered_feature_bc_matrix.h5'}}
#'  \item{\emph{'/outs/spatial/*.jpg}}
#'  }
#'
#' @inherit gene_set_path params
#'
#' @param sample_name Character value. The sample name with which to refer to the
#' respective sample. Should start with a letter.
#'

#' @inherit process_seurat_object params
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

initiateSpataObject_10X <- function(directory_10X,
                                    sample_name,
                                    directory_spata = NULL,
                                    directory_seurat = NULL,
                                    combine_with_wd = "/",
                                    gene_set_path = NULL,
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

  confuns::give_feedback(
    msg = "Starting initiation.",
    verbose = verbose
  )

  # check input for sample and directory
  confuns::check_directories(directories = directory_10X, type = "folders")

  confuns::are_values(c("directory_10X", "sample_name"), mode = "character")

  if(sample_name %in% c("", "all")){

    base::stop(glue::glue("' ' and 'all' are invalid sample names."))

  }


  # check input for gene set data.frame
  confuns::is_value(x = gene_set_path, mode = "character", skip.allow = TRUE, skip.val = NULL)

  if(base::is.null(gene_set_path)){

    if(base::isTRUE(verbose)){base::message("No gene-set data.frame path specified.")}

  } else {

    confuns::check_directories(directories = gene_set_path, ref = "gene_set_path", type = "files")

  }

  # check input for seurat processing functions
  for(fn in seurat_process_fns){

    input <- base::parse(text = fn) %>% base::eval()

    if(base::is.data.frame(input) | (!base::isTRUE(input) && !base::is.list(input) &&!base::isFALSE(input))){

      base::stop(glue::glue("Invalid input for argument '{fn}'. Must either be TRUE, FALSE or a named list."))

    }

  }


  # -----


  # 2. Read in data ---------------------------------------------------------

  if(base::isTRUE(verbose)){base::message("Step 1/4 : Reading in .h5 file.")}

  data_dir <- base::paste0(directory_10X, "/outs")
  file_dir <- base::paste0(directory_10X, "/outs/filtered_feature_bc_matrix.h5")

  if(base::file.exists(paths = file_dir)){

   base::message(glue::glue("Loading from directory: '{data_dir}'"))

   seurat_object <- Seurat::Load10X_Spatial(data.dir = data_dir, filename = "filtered_feature_bc_matrix.h5")

  } else {

   base::stop(glue::glue("Directory '{file_dir}' does not exist."))

  }

  # -----


  # 3. Seurat analysis ------------------------------------------------------

  if(base::isTRUE(verbose)){base::message("Step 2/4: Performing Seurat-analysis steps.")}

  processed_seurat_object <-
    process_seurat_object(
      seurat_object = seurat_object,
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

  # -----


  # 5. Create SPATA-object --------------------------------------------------

  if(base::isTRUE(verbose)){base::message("Step 3/4: Initiating spata-object.")}

  spata_object <-
    transformSeuratToSpata(
      seurat_object = processed_seurat_object,
      sample_name = sample_name,
      gene_set_path = gene_set_path,
      method = "spatial",
      verbose = verbose
    )

  spata_object <- setInitiationInfo(spata_object)

  # -----


  # 6. Save objects and return spata object ---------------------------------

  # save seurat
  if(base::is.character(directory_seurat)){

    spata_object <-
      base::tryCatch({

        spata_object_return <-
          saveCorrespondingSeuratObject(
            seurat_object = seurat_object,
            object = spata_object,
            directory_seurat = directory_seurat,
            combine_with_wd = combine_with_wd,
            verbose = verbose
          )

        spata_object_return

      }, error = function(error){

        base::warning(glue::glue("Attempt to save seurat-object under '{directory_seurat}' failed with the following error message: {error}"))

        spata_object

      })

  } else {

    confuns::give_feedback(
      msg = "Skip saving seurat-object.",
      verbose = verbose
    )

  }

  if(base::is.character(directory_spata)){

    spata_object <-
      base::tryCatch({

        spata_object_ret <-
          saveSpataObject(
            object = spata_object,
            directory_spata = directory_spata,
            combine_with_wd = combine_with_wd,
            verbose = verbose
          )

        spata_object_ret

      }, error = function(error){

        base::warning(glue::glue("Attempt to save spata-object under '{directory_spata}' failed with the following error message: {error}"))

        spata_object

      })

  } else {

    confuns::give_feedback(
      msg = "Skip saving spata-object.",
      verbose = verbose
    )

  }

  confuns::give_feedback(
    msg = "Initiation finished.",
    verbose = verbose
  )

  # -----

  base::return(spata_object)

}
