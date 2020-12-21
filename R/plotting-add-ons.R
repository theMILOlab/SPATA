
#' @title Ggplot add on wrapper
#' @export
legend_bottom <- purrr::partial(.f = ggplot2::theme, legend.position = "bottom")

#' @rdname legend_bottom
#' @export
legend_top <- purrr::partial(.f = ggplot2::theme, legend.position = "top")


#' @rdname legend_bottom
#' @export
legend_none <- purrr::partial(.f = ggplot2::theme, legend.position = "none")
