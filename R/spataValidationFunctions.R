#' @title Validate a spata object
#'
#' @description Takes a spata object and checks whether all slots contain suitable
#' data. If not it attempts to provide a helpful report.
#'
#' @param object A spata-object.
#'
#' @return A character string printed by \code{base::writeLines()}
#' @export
#'

validateSpataObject <- function(object){

  validation(x = object)

# 1. Examine the slot names -----------------------------------------------

  input_slots <- methods::slotNames(object) %>% sort()
  spata_slots <- c("coordinates", "data", "description", "dim_red", "fdata",
                   "image", "samples", "scvelo", "trajectories", "used_genesets")

  # check for missing input_slots
  if(!base::all(spata_slots %in% input_slots)){

    not_found <-
      stringr::str_c(spata_slots[!spata_slots %in% input_slots], collapse = "', '")

    base::message(stringr::str_c("Could not find slots: '", not_found,
                                 "'. Can not validate slots that do not exist." )
    )


  }

  # check for unknown input slots
  if(base::any(!input_slots %in% spata_slots)){

    unknown <-
      stringr::str_c(input_slots[!input_slots %in% spata_slots], collapse = "', '")

      base::message(stringr::str_c("Ignorign unknown slots: '", unknown, "'."))

    # keep only valid input_slots
    input_slots <- spata_slots[spata_slots %in% input_slots]

  }

  feedback <- base::vector(mode = "list")

  for(slot in input_slots){

    fun <- stringr::str_c("check_slot_", slot, sep = "")

    feedback[[slot]] <- base::do.call(fun, list(object))

  }

  # unlist feedback
  feedback_vec <- base::unlist(x = feedback) %>% unname()
  prefix <- stringr::str_c("Slot '", base::names(feedback), "': ", sep = "")

  # combine with descriptive prefix
  final_feedback <- stringr::str_c(prefix, feedback_vec, sep = "")

  # return results
  base::writeLines(final_feedback)

}



#' @title Check spata slots
#'
#' @description Functions that provide a report regarding the validity of
#' the respective slot.
#'
#' @param object A spata-object.
#'
#' @return A character string. (Call \code{base::writeLines() wtih that string as input.})
#' @export
#'

check_slot_coordinates <- function(object){

  coords <- object@coordinates

  messages <- base::character()

  # column names and samples
  c_colnames <- base::colnames(coords)

  if(!base::identical(c_colnames, c("barcodes", "sample", "x", "y"))){

    c_colnames <- stringr::str_c(base::colnames(coords), collapse = "', '")
    feedback <- stringr::str_c("Invalid column names.",
                               "\n Currently:  '",  c_colnames,
                               "\n Have to be: 'barcodes', 'sample', 'x', 'y'",
                               sep = "")

    messages <-
      base::append(x = messages,
                   values = feedback)

  } else {

    messages <- hlpr_compare_samples(object, df = coords, messages = messages)

  }

  # variable classes
  c_classes <- base::sapply(coords, base::class) %>% base::unname()

  if(!base::identical(c_classes, c("character", "character", "numeric", "numeric"))){

    c_classes <- stringr::str_c(c_classes, collapse = "', '")
    feedback <- stringr::str_c("Invalid column classes.",
                               "\n             'barcodes',  'sample',    'x',       'y'",
                               "\n Currently:  '", c_classes,
                               "\n Have to be: 'character', 'character', 'numeric', 'numeric'.",
                               sep = "")

    messages <-
      base::append(x = messages,
                   values = feedback)
  }

  # return

  if(base::identical(messages, base::character())){

    base::return("Valid!")

  } else {

    base::return(messages)

  }

}

#' @rdname check_slot_coordinates
#'
#' @export
check_slot_data <- function(object){

  data <- object@data

  messages <- base::character()

  if(!base::is.matrix(data@counts)){

    messages <-
      base::append(messages,
                   values = "Slot 'counts' needs to be a numeric matrix.")

  }

  if(!base::is.matrix(data@norm_exp)){

    messages <-
      base::append(messages,
                   values = "Slot 'norm_exp' needs to be a numeric matrix.")

  }

  if(base::identical(messages, base::character())){

    base::return("Valid!")

  } else {

    base::return(messages)

  }


}


#' @rdname check_slot_coordinates
#'
#' @export
check_slot_description <- function(object){

  base::return("Valid!")

}


#' @rdname check_slot_coordinates
#'
#' @export
check_slot_scvelo <- function(object){

  "Valid!"

}


#' @rdname check_slot_coordinates
#'
#' @export
check_slot_trajectories <- function(object){

  "Valid!"

}


#' @rdname check_slot_coordinates
#'
#' @export
check_slot_dim_red <- function(object){

  messages <- base::character()

  input_slots <- methods::slotNames(x = object@dim_red) %>% base::sort()
  dim_red_slots <- c("UMAP", "TSNE") %>% base::sort()

  if(!base::identical(input_slots, dim_red_slots)){

    messages <- base::append(x = messsages, values = "Invalid slot names. Have to be 'UMAP' and 'TSNE'.")

    return(messages)

  } else {

    # UMAP --------------------------------------------------------------------

    umap_df <- object@dim_red@UMAP

    # column names and samples
    u_colnames <- base::colnames(umap_df)

    if(!base::identical(u_colnames, c("barcodes", "sample", "umap1", "umap2"))){

      u_colnames <- stringr::str_c(u_colnames, collapse = "', '")
      feedback <- stringr::str_c("Invalid column names in slot 'UMAP'.",
                                 "\n Currently:  '",  u_colnames,
                                 "\n Have to be: 'barcodes', 'sample', 'umap1', 'umap2'",
                                 sep = "")

      messages <-
        base::append(x = messages,
                     values = feedback)

    } else if(base::nrow(umap_df) != 0){

      messages <- hlpr_compare_samples(object, df = umap_df, messages = messages)

    } else if(base::nrow(umap_df) == 0){

      messages <- base::append(x = messages, values = "UMAP data.frame is empty.")

    }

    # variable classes
    u_classes <- base::sapply(umap_df, base::class) %>% base::unname()

    if(!base::identical(u_classes, c("character", "character", "numeric", "numeric"))){

      u_classes <- stringr::str_c(u_classes, collapse = "', '")
      feedback <- stringr::str_c("Invalid column classes in slot 'UMAP'.",
                                 "\n             'barcodes',  'sample',    'umap1',   'umap2'",
                                 "\n Currently:  '", u_classes,
                                 "\n Have to be: 'character', 'character', 'numeric', 'numeric'.",
                                 sep = "")

      messages <-
        base::append(x = messages,
                     values = feedback)

    }


    # TSNE --------------------------------------------------------------------

    tsne_df <- object@dim_red@TSNE

    # column names and samples
    t_colnames <- base::colnames(tsne_df)

    if(!base::identical(t_colnames, c("barcodes", "sample", "tsne1", "tsne2"))){

      t_colnames <- stringr::str_c(t_colnames, collapse = "', '")
      feedback <- stringr::str_c("Invalid column names in slot 'TSNE'.",
                                 "\n Currently:  '",  t_colnames,
                                 "\n Have to be: 'barcodes', 'sample', 'tsne1', 'tsne2'",
                                 sep = "")

      messages <-
        base::append(x = messages,
                     values = feedback)

    } else if(base::nrow(tsne_df) != 0) {

      messages <- hlpr_compare_samples(object, df = tsne_df, messages = messages)

    } else if(base::nrow(tsne_df) == 0){

      messages <- base::append(x = messages, values = "TSNE data.frame is empty.")

    }

    # variable classes
    t_classes <- base::sapply(tsne_df, base::class) %>% base::unname()

    if(!base::identical(t_classes, c("character", "character", "numeric", "numeric"))){

      t_classes <- stringr::str_c(t_classes, collapse = "', '")
      feedback <- stringr::str_c("Invalid column classes in slot 'TSNE'.",
                                 "\n             'barcodes',  'sample',    'tsne1',   'tsne2'",
                                 "\n Currently:  '", t_classes,
                                 "\n Have to be: 'character', 'character', 'numeric', 'numeric'.",
                                 sep = "")

      messages <-
        base::append(x = messages,
                     values = feedback)

    }


    # Return ------------------------------------------------------------------

    if(base::identical(messages, base::character())){

      base::return("Valid!")

    } else {

      base::return(messages)

    }

  }

}


#' @rdname check_slot_coordinates
#'
#' @export
check_slot_fdata <- function(object){

  fdata <- object@fdata
  messages <- base::character()

  if(base::nrow(fdata) == 0){

    messages <- base::append(x = messages, values = "'fdata' data.frame is empty.")

  } else {

    # Column names  -----------------------------------------------------------

    f_colnames <- base::colnames(fdata)
    missing <- character(0)

    if(!c("sample") %in% f_colnames){

      missing <- base::append(x = missing, values = "sample")

    }

    if(!"barcodes" %in% f_colnames){

      missing <- base::append(x = missing, values = "barcodes")

    }

    if(!"segment" %in% f_colnames){

      missing <- base::append(x = missing, values = "segment")

    }

    if(base::length(missing) != 0){

      missing <- stringr::str_c(missing, collapse = "', '")

      messages <- base::append(x = messages,
                               values = stringr::str_c(
                                 "Missing columns in 'fdata': '",
                                 missing, "'", sep = ""
                               ))

    } else {

      # variable classes
      f_classes <- base::sapply(fdata[,c("sample", "barcodes", "segment")], base::class) %>% base::unname()

      if(!base::identical(f_classes, c("character", "character", "character"))){

        f_classes <- stringr::str_c(f_classes, collapse = "', '")
        feedback <- stringr::str_c("Invalid column classes in 'fdata'.",
                                   "\n Columns:    'barcodes',  'sample',    'segment'",
                                   "\n Classes:    '", f_classes,
                                   "\n Have to be: 'character', 'character', 'character'",
                                   sep = "")

        messages <-
          base::append(x = messages,
                       values = feedback)
      }

      # compare samples
      if(base::all(f_classes[1:2] == "character")){

        messages <- hlpr_compare_samples(object, df = fdata, messages = messages)

      }

    }

    # Return ------------------------------------------------------------------

    if(base::identical(messages, base::character())){

      base::return("Valid!")

    } else {

      base::return(messages)

    }

  }

}


#' @rdname check_slot_coordinates
#'
#' @export
check_slot_image <- function(object){

  image_list <- object@image
  messages <- base::character()

  i_samples <- base::names(image_list) %>% base::sort()
  o_samples <- samples(object) %>% base::sort()

  # if sample names match check for classes
  if(!base::identical(i_samples, o_samples )){

    i_samples <- stringr::str_c(i_samples, collapse = ", ")
    o_samples <- stringr::str_c(o_samples, collapse = ", ")

    messages <- base::append(x = messages,
                             values = stringr::str_c(
                               "Invalid name(s) in 'image-list'. Must match samples in object.",
                               "\n Image names: ", i_samples,
                               "\n In object  : ", o_samples,
                               sep = ""
                             ))

  } else {

    i_samples <- base::names(image_list)

    i_classes <- base::sapply(X = image_list, FUN = base::class) %>% base::unname()

    if(!base::all(i_classes == "Image")){

      invalid_images <- stringr::str_c(i_samples[i_classes != "Image"], collapse = ", ")

      messages <- base::append(x = messages,
                               values = stringr::str_c("Invalid class in slot(s): '",
                                                       invalid_images, "' of image-list.",
                                                       " Must be of class 'Image'.",
                                                       sep = ""))



    }

  }


  # return
  if(base::identical(messages, base::character())){

    base::return("Valid!")

  } else {

    base::return(messages)

  }






}


#' @rdname check_slot_coordinates
#'
#' @export
check_slot_used_genesets <- function(object){

  gs_df <- object@used_genesets
  messages <- base::character()

  if(base::nrow(gs_df) == 0){

    messages <- base::append(x = messages,
                             values = "'used_geneset' data.frame is empty")

  } else {

    gs_colnames <- base::colnames(gs_df)

    # check for column names and then for column classes
    if(!base::all(gs_colnames %in% c("ont", "gene"))){

      gs_colnames <-
        gs_colnames %>%
        base::sort() %>%
        stringr::str_c(collapse = "', '")

      messages <-
        base::append(x = messages,
                     values = stringr::str_c("Invalid column names in slot 'used_genesets'",
                                             "\n Currently: '", gs_colnames, "'",
                                             "\n Must be  : 'gene', 'ont'",
                                             sep = ""))


    } else {

      gs_classes <- base::sapply(X = gs_df, FUN = base::class) %>% base::unname()

      if(!base::all(gs_classes == "character")){

        gs_classes <-
          gs_classes %>%
          base::sort() %>%
          stringr::str_c(collapse = "', '")

        messages <-
          base::append(x = messages,
                       values = stringr::str_c("Invalid column classes:",
                                               "\n Currently: '", gs_classes, "'",
                                               "\n Must be  : 'character', 'character'",
                                               sep = ""))

      }

    }

  }

  # check for rownames
  if(!base::identical(base::rownames(gs_df), base::as.character(1:base::nrow(gs_df)))){

    messages <- base::append(x = messages,
                             values = "The use of row.names is discouraged.")

  }

  # return
  if(base::identical(messages, base::character())){

    base::return("Valid!")

  } else {

    base::return(messages)

  }

}


#' @rdname check_slot_coordinates
#'
#' @export
check_slot_samples <- function(object){

  samples <- object@samples

  if(base::length(samples) != 0){

    "Valid!"

  } else {

    "Slot 'samples' must not be of length zero."

  }

}



