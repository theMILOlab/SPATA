#' @title getSurroundedSpots
#'
#' @description This Function will return a data frame with the sourounding spots
#'
#' @param SPATAobj
#' @param of_sample the sample to be analyzed
#' @param assay_type Visium with 6 sourounding spots and ST with 8 spots
#'
#'
#'
#' @return A spata-object.
#'
#' @importFrom Seurat ScaleData
#'
#' @export
getSurroundedSpots <- function(object, of_sample=NULL, assay_type="Visium"){

  if(is.null(of_sample)){of_sample=object@samples}
  if(length(of_sample)>1) stop("Unclear ... of_sample, please define")


  if(assay_type=="Visium"){
    coords <- SPATA::getCoordinates(object, of_sample)
    bc_origin <- coords$barcodes
    bc_destination <- coords$barcodes

    # get grouped data.frame with all barcode combinations

    dfgr <-
      tidyr::expand_grid(bc_origin, bc_destination) %>%
      dplyr::left_join(x = ., y = dplyr::select(coords, bc_origin = barcodes, xo = x, yo = y), by = "bc_origin") %>%
      dplyr::left_join(x = ., y = dplyr::select(coords, bc_destination = barcodes, xd = x, yd = y), by = "bc_destination") %>%
      dplyr::mutate(distance = base::round(base::sqrt((xd - xo)^2 + (yd - yo)^2), digits = 0)) %>%
      dplyr::group_by(bc_origin)

    # filter barcodes that are surrounded by a total of six spots

    sufficient_bc <-
      dplyr::slice_min(.data = dfgr, order_by = distance, n = 7) %>%
      dplyr::group_by(bc_origin, distance) %>%
      dplyr::summarise(count = dplyr::n()) %>%
      dplyr::filter(distance != 0 & count == 6) %>%
      dplyr::pull(bc_origin)

    # filter for barcodes

    final_df <-
      dplyr::filter(dfgr, bc_origin %in% sufficient_bc) %>%
      dplyr::slice_min(order_by = distance, n = 7) %>%
      dplyr::mutate(sample = {{of_sample}}) %>%
      dplyr::select(sample, dplyr::everything())

    base::return(final_df)
  }else{
    if(assay_type=="ST"){

      coords <- SPATA::getCoordinates(object, of_sample)
      bc_origin <- coords$barcodes
      bc_destination <- coords$barcodes

      # get grouped data.frame with all barcode combinations

      dfgr <-
        tidyr::expand_grid(bc_origin, bc_destination) %>%
        dplyr::left_join(x = ., y = dplyr::select(coords, bc_origin = barcodes, xo = x, yo = y), by = "bc_origin") %>%
        dplyr::left_join(x = ., y = dplyr::select(coords, bc_destination = barcodes, xd = x, yd = y), by = "bc_destination") %>%
        dplyr::mutate(distance = base::round(base::sqrt((xd - xo)^2 + (yd - yo)^2), digits = 0)) %>%
        dplyr::group_by(bc_origin)

      # filter barcodes that are surrounded by a total of 8 spots

      sufficient_bc <-
        dplyr::slice_min(.data = dfgr, order_by = distance, n = 9) %>%
        dplyr::group_by(bc_origin) %>%
        dplyr::summarise(count = dplyr::n()) %>%
        dplyr::filter(count == 9) %>%
        dplyr::pull(bc_origin)

      # filter for barcodes

      final_df <-
        dplyr::filter(dfgr, bc_origin %in% sufficient_bc) %>%
        dplyr::slice_min(order_by = distance, n = 9) %>%
        dplyr::mutate(sample = {{of_sample}}) %>%
        dplyr::select(sample, dplyr::everything())

      base::return(final_df)

    }else{message("assay_type unknown")}
  }

}

