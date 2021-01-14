
#' @title Ggplot add on wrapper
#' @export
legendBottom <- purrr::partial(.f = ggplot2::theme, legend.position = "bottom")

#' @rdname legend_bottom
#' @export
legendNone <- purrr::partial(.f = ggplot2::theme, legend.position = "none")

#' @rdname legend_bottom
#' @export
legendRight <- purrr::partial(.f = ggplot2::theme, legend.position = "right")

#' @rdname legend_bottom
#' @export
legendTop <- purrr::partial(.f = ggplot2::theme, legend.position = "top")
