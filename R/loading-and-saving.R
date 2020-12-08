#' @include S4-documentation.R
#'
NULL


#' @title Original load gene set data.frame
#'
#' @description Not exported due to naming issues. Kept as it is used in several
#' loading functions.

loadGSDF <- function(gene_set_path = NULL, verbose = TRUE){

  if(!base::is.null(gene_set_path)){

    confuns::is_value(x = gene_set_path, mode = "character", ref = "gene_set_path")
    confuns::check_directories(directories = gene_set_path, ref = "gene_set_path", type = "files")

    if(base::isTRUE(verbose)){base::message(glue::glue("Reading in specified gene-set data.frame from directory '{gene_set_path}'."))}

    gene_set_df <- base::readRDS(file = gene_set_path)

    if(!base::is.data.frame(gene_set_df)){

      gene_set_df <- gsdf

      base::warning(glue::glue("Input from directory '{gene_set_path}' is not a data.frame. Using SPATA's default gene set data.frame."))

    }

  } else {

    if(base::isTRUE(verbose)){base::message("Using SPATA's default gene set data.frame.")}

    gene_set_df <- gsdf

  }

  base::return(gene_set_df)

}

#' @title Load gene set data.frame
#'
#' @param gene_set_path If set to NULL the default \code{SPATA::gsdf} is used.
#' If a directory is specified the object is loaded via \code{base::readRDS()}, checked
#' and used if valid. If it is invalid the default \code{SPATA::gsdf} is used .
#' @inherit verbose params
#'
#' @return A data.frame.
#'
#' @export

loadGeneSetDf <- loadGSDF



#' @title Load and save a spata-object
#'
#' @description Wrapper around \code{base::readRDS()} and \code{base::saveRDS()}.
#'
#' @param input_path Character value. The directory leading to the spata-object.
#' @inherit check_object params
#' @inherit initiateSpataObject_10X params
#' @param overwrite Logical. Needs to be set to TRUE if the resulting directory from
#' \code{output_path} and \code{file_name} already exists.
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
saveSpataObject <- function(object, output_path, file_name, overwrite = FALSE){

  # 1. Control --------------------------------------------------------------

  check_object(object)
  confuns::is_value(output_path, "character", "output_path")
  confuns::is_value(file_name, "character", "file_name")
  confuns::check_directories(output_path, ref = "output_path", type = "folders")

  # -----

  filename <- stringr::str_c(output_path, "/spata-obj-", file_name, ".RDS", sep = "")

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
#' @param file_name Character value. The filename. ( is suffixed with \emph{'.RDS'})
#'
#' @return An invisible TRUE if saved successfully or an informative error message.
#' @export
#'

saveGeneSetDf <- function(object,
                          output_path,
                          file_name){

  check_object(object)

  confuns::is_value(output_path, "character", "output_path")
  confuns::is_value(file_name, "character", "file_name")

  confuns::check_directories(output_path, "output_path", type = "folders")

  final_path <- stringr::str_c(output_path, "/", file_name, ".RDS", sep = "")

  if(base::file.exists(final_path)){

    base::stop(glue::glue("The file '{final_path}' already exists."))

  } else if(base::nrow(object@used_genesets) == 0){

    base::stop("The objects's gene-set data.frame is empty.")

  } else {

    base::saveRDS(object = object@used_genesets, file = final_path)

    if(base::file.exists(final_path)){

      file_name <- stringr::str_c("~/", file_name, ".RDS", sep = "")
      base::message(glue::glue("Gene set data.frame has been saved as '{file_name}'."))
      base::return(base::invisible(TRUE))

    } else {

      base::stop("Could not save the gene-set data.frame. Unknown error.")

    }

  }

}


















