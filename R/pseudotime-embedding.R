#' @title Initiate monocle3-pseudotime analysis
#'
#' @description Loads or compiles a valid cell_data_set-object
#' that fits to the provided spata-object.
#'
#' @inherit check_object params
#' @param use_cds_file A file-directory leading to a .rds file containing a valid
#' cell_data_set-object previously calculated for the specified object. Specified
#' as a character value. If set to FALSE the cell_data_set object will be created
#' from scratch.
#' @param save_cds_file A file-directory (that does not already exists) under which the used or created cell_data_set-object
#' is going to be stored specified as a character value. Should end with .rds.
#' @param preprocess_method Given to \code{monocle3::preprocess_cds()} if \code{use_cds_file} isn't a character string.
#' @param cluster_method Given to \code{monocle3::cluster_cells()} if \code{use_cds_file} isn't a character string. Must be one of
#' \emph{'leiden', 'louvain'}.
#' @param inherit verbose
#'
#' (Warning messages will always be printed.)
#'
#' @return A monocle3::cell_data_set object.
#' @export

compileCellDataSet <- function(object,
                               use_cds_file = FALSE,
                               save_cds_file = FALSE,
                               preprocess_method = "PCA",
                               cluster_method = "leiden",
                               verbose = TRUE){


  check_object(object)
  confuns::is_value(preprocess_method, "character", "preprocess_method")
  confuns::is_value(cluster_method, mode = "character", "cluster_method")
  if(!base::isFALSE(use_cds_file)){confuns::check_directories(use_cds_file, ref = "use_cds_file", type = "files")}
  if(!base::isFALSE(save_cds_file)){

    confuns::is_value(save_cds_file, "character", "save_cds_file")
    if(base::file.exists(save_cds_file)){

      base::stop(glue::glue("Directory '{save_cds_file}' already exists. "))

    }

  }


  # check if valid cds files
  if(base::is.character(use_cds_file)){

    cds <- base::readRDS(use_cds_file)

    if(!methods::is(object = cds, class2 = "cell_data_set")){

      base::stop(stringr::str_c("File '", use_cds_file, "' is not a valid object of class 'cell_data_set'."))

    }

    barcodes_cds <- base::names(monocle3::pseudotime(cds)) %>% base::sort()
    barcodes_spata <- featureData(object)$barcodes %>% base::sort()

    if(!base::identical(barcodes_cds, barcodes_spata)){

      base::stop(stringr::str_c("The barcodes of '", use_cds_file, "' and the provided spata-object are not identical."))

    }


  } else { # start from scratch

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

    if(base::isTRUE(verbose)){base::message("Step 2/7 Estimating size factors.")}
    cds <- monocle3::estimate_size_factors(cds)

    if(base::isTRUE(verbose)){base::message("Step 3/7 Preprocessing cell data set.")}
    cds <- monocle3::preprocess_cds(cds, method = preprocess_method, num_dim = 30)

    if(base::isTRUE(verbose)){base::message("Step 4/7 Reducing dimensions.")}
    cds <- monocle3::reduce_dimension(cds, cores = 2)

    if(base::isTRUE(verbose)){base::message("Step 5/7 Clustering cells.")}
    cds <- monocle3::cluster_cells(cds, cluster_method = cluster_method)

    if(base::isTRUE(verbose)){base::message("Step 6/7 Learning trajectory.")}
    cds <- monocle3::learn_graph(cds)

  }

  if(base::isTRUE(verbose)){base::message("Step 7/7 Ordering cells.")}
  cds <- monocle3::order_cells(cds)


  # save cds file if save_cds_file is specified as a character
  if(base::is.character(save_cds_file)){

    if(base::isTRUE(verbose)){

      base::message(stringr::str_c("Saving cell data set object 'cds' under directory: '", save_cds_file, "'"))

    }

    base::saveRDS(cds, file = save_cds_file)

  }

  base::return(cds)

}
