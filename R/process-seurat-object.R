

#' @title Wrapper around Seurat processing functions
#'
#' @param seurat_object A valid seurat-object.
#' @inherit compileSeuratObject params
#'
#' @inherit verbose params
#'
#' @return A processed seurat-object.
#'

process_seurat_object <- function(seurat_object,
                                  calculate_rb_and_mt = TRUE,
                                  remove_stress_and_mt = TRUE,
                                  SCTransform = FALSE,
                                  NormalizeData = TRUE,
                                  FindVariableFeatures = TRUE,
                                  ScaleData = TRUE,
                                  RunPCA = TRUE,
                                  FindNeighbors = TRUE,
                                  FindClusters = TRUE,
                                  RunTSNE = TRUE,
                                  RunUMAP = TRUE,
                                  verbose = TRUE){

# 1. Control --------------------------------------------------------------

  base::stopifnot(methods::is(object = seurat_object, class2 = "Seurat"))

  for(fn in seurat_process_fns){

    input <- base::parse(text = fn) %>% base::eval()

    if(base::is.data.frame(input) | (!base::isTRUE(input) && !base::is.list(input) &&!base::isFALSE(input))){

      base::stop(glue::glue("Invalid input for argument '{fn}'. Must either be TRUE, FALSE or a named list."))

    }

  }

  # calculate ribosomal and mitochondrial percentage
  if(base::isTRUE(calculate_rb_and_mt)){

    if(base::isTRUE(verbose)){

      base::message("Calculating percentage of ribosomal and mitochondrial genes.")

    }

    seurat_object[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_object, pattern = "^MT.")
    seurat_object[["percent.RB"]] <- Seurat::PercentageFeatureSet(seurat_object, pattern = "^RPS")

  }

  # remove stress and mitochondrial genes
  if(base::isTRUE(remove_stress_and_mt)){

    if(base::isTRUE(verbose)){

      base::message("Removing stress genes and mitochondrial genes.")

    }

    exclude <- c(rownames(seurat_object@assays$RNA)[base::grepl("^RPL", rownames(seurat_object@assays$RNA))],
                 rownames(seurat_object@assays$RNA)[base::grepl("^RPS", rownames(seurat_object@assays$RNA))],
                 rownames(seurat_object@assays$RNA)[base::grepl("^MT-", rownames(seurat_object@assays$RNA))],
                 c('JUN','FOS','ZFP36','ATF3','HSPA1A","HSPA1B','DUSP1','EGR1','MALAT1'))

    feat_keep <- base::rownames(seurat_object@assays$RNA[!(base::rownames(seurat_object@assays$RNA) %in% exclude), ])

    seurat_object <- base::subset(x = seurat_object, features = feat_keep)

  }



# 2. Process seurat object ------------------------------------------------

  functions_to_call <- seurat_process_fns

  for(fn in functions_to_call){

    input <-
      base::parse(text = fn) %>%
      base::eval()

    if(base::isTRUE(input)){

      if(base::isTRUE(verbose)){base::message(glue::glue("Running 'Seurat::{fn}()' with default parameters."))}

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

      if(base::isTRUE(verbose)){base::message(glue::glue("Running 'Seurat::{fn}()' with specified parameters."))}

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

          base::message(glue::glue("Running'Seurat::{fn}()' resulted in the following error: {error$message}. Abort and continue with next function."))

          base::return(seurat_object)

        }

      )

    } else {

      if(base::isTRUE(verbose)){base::message(glue::glue("Skip running '{fn}()' as it's argument input is neither TRUE nor a list."))}

    }

  }

  base::return(seurat_object)

}



seuratSingleCellToSpata <- function(seurat_object,
                                    sample_names = NULL,
                                    verbose = T){

  # Input create ------------------------------------------------------------
  base::message("convert a scRNA-seq data set to spata")

  if(!base::is.null(sample_names)){confuns::is_vec(x = sample_names, mode = "character", ref = "sample_names")}
  if(base::is.null(sample_names)){sample_names = "S1"}

  coordinates <-
    base::as.data.frame(seurat_object@reductions$umap@cell.embeddings) %>%
    tibble::rownames_to_column("barcodes") %>%
    magrittr::set_colnames(value = c("barcodes", "x", "y"))

  intensity_matrix <- base::as.matrix(Seurat::GetAssayData(seurat_object))

  meta_data <- seurat_object@meta.data %>% tibble::rownames_to_column("barcodes")

  counts <- seurat_object@assays$RNA@counts

  geneSets <- SPATA::gsdf

  genes <- base::rownames(intensity_matrix)
  string_contain <- base::length(stringr::str_detect(genes, "_") %in% TRUE) <= 1

  if(base::isTRUE(string_contain)) {

    base::message("Rownames of intensity matrix contain a '_'. This is not allowed and will be changed into '-' ")

    base::rownames(intensity_matrix) <-
      stringr::str_replace(string = base::rownames(intensity_matrix), pattern = "_", replacement = "-")

  }

  barcodes_matrix <- base::colnames(intensity_matrix) %>% base::sort()
  barcodes_coordinates <- dplyr::pull(coordinates, var = "barcodes") %>% base::sort()

  # check inputs
  if(!base::identical(barcodes_matrix, barcodes_coordinates)){

    base::stop("The barcodes of the coordinate system and the column names of the assay must be identical.")

  }

  # Add sample name to coordinates
  coordinates <-
    dplyr::mutate(.data = coordinates, sample = {{sample_names}}) %>%
    dplyr::select("barcodes","sample","x","y")

  # Change Coordinates for image in SPATA

  # Start -------------------------------------------------------------------


  #Create empty SPATA object
  obj <- methods::new(Class="spata")

  obj@coordinates <- coordinates
  obj@data@norm_exp <- intensity_matrix
  obj@data@counts <- counts
  obj@samples <- base::unique(obj@coordinates$sample)
  obj@fdata <- obj@coordinates %>% dplyr::select(barcodes, sample) %>% dplyr::mutate(segment = "")
  obj@fdata <- obj@fdata %>% dplyr::left_join(., meta_data, by="barcodes")

  names(obj@fdata)[2] <- "sample"
  obj@used_genesets=geneSets
  obj@trajectories=list(sample=list())
  names(obj@trajectories)=obj@samples

  obj@dim_red@UMAP <- data.frame(seurat_object@reductions$umap@cell.embeddings) %>% tibble::rownames_to_column("barcodes")

  return(obj)

}


seurat10XVisiumToSpata <- function(seurat_object,
                                   sample_names = NULL,
                                   verbose = TRUE){

}





