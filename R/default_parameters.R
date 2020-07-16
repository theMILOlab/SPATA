#' Title
#'
#' @param object A valid spata-object.
#' @param coords_df A data.frame that contains variables \emph{barcodes, sample}.
#' @param features The features you want to join \code{coords_df} with specified
#' as a character vector.
#' @param genes The genes you want to join \code{coords_df} with specified
#' as a character vector.
#' @param average_genes Logical value. If set to TRUE coords_df will be joined with
#' the mean expression values of all genes specified.
#' @param gene_sets The gene sets you want to join \code{coords_df} with specified
#' as a character vector.
#' @param method_gs The method according to which gene sets will be handled
#' specified as a character of length one. This can be either \emph{mean} or one
#' of \emph{gsva, ssgsea, zscore, or plage}. The latter four will be given to
#' \code{gsva::GSVA()}.
#' @param smooth Logical value. If set to TRUE values will be
#' smoothed with respect to their local neighbors using \code{stats::loess()}.
#' @param smooth_span Numeric value, given to \code{stats::loess()} if
#'  \code{smooth} is set to TRUE.
#' @param verbose Logical value. If set to TRUE informative messages with respect
#' to the computational progress made will be printed.

#'
#' @return
#' @export
#'
#' @examples
#'
default_parameters <- function(){



}
