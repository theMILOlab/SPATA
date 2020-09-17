#' @title Color palette names
#' @description Returns all currently valid color panels or -spectra.
#' @return A named list.
#'
#' @details Discrete color panels derive from the ggsci-package. Continuous
#' colorspectra derive from the colorspace-package.
#'
#' @export

all_colorpanels <- confuns::all_colorpanels

#' @rdname all_colorpanels
#' @export

all_colorspectra <- confuns::all_colorspectra


#'@inherit confuns::scale_color_add_on
#'@export

scale_color_add_on <- confuns::scale_color_add_on
