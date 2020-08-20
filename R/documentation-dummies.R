
#' @title Verbose
#' @param verbose Logical. If set to TRUE informative messages regarding
#' the computational progress will be printed.
#'
#' (Warning messages will always be printed.)

verbose <- function(verbose){

}



#' @title plot_family
#' @return Returns a ggplot-object that can be additionally customized according
#' to the rules of the ggplot2-framework.
#'

plot_family <- function(){}



#' @title normalize
#' @param normalize Logical. If set to TRUE values will be scaled to 0-1.
#'
#' Hint: Variables that are uniformly expressed can not be scaled and are discarded.
#'
#'

normalize <- function(normalize){}



#' @title image_dummy
#' @param image An image of class \emph{Image} to be displayed in the background.
#' Easily accessible via \code{SPATA::image()}.

image_dummy <- function(image){}



#' @title average_genes
#' @param average_genes Logical value. If set to TRUE the mean expression
#' of all specified genes will be calculated.

average_genes <- function(average_genes){}


#' @title variable
#'
#' @param variabel The variable of interest.
#'
#'  \itemize{
#'   \item{ \strong{Gene set} as a single character value. Must be in \code{getGeneSets()}}
#'   \item{ \strong{Genes} as a character vector. If more than one gene is specified the average
#'   expression of those genes will be calculated and displayed. Must be in \code{getGenes()}}
#'   \item{ \strong{Feature} as a single character value. Must be in \code{getFeaturenNames()}}
#'   }
#'
#'

variable <- function(variable){}

