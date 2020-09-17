
#' @title verbose
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
#' @param variable The variable of interest.
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


#' @title variables_num
#'
#' @param variables Character vector. The numeric variables of interest.

variables_num <- function(variables){}


#' @title clrp
#' @param clpr Character value. The color panel to be used. Run \code{all_colorpanels()} to see
#' all valid input options.

clrp <- function(clrp){}


#' @title getC_joinW_combo
#' @seealso Combine \code{getCoordinates()} and \code{joinWithGeneSets()} to obtain
#' a valid input data.frame for \code{data}.
#'
getC_joinW_combo <- function(){}


#' @title across
#' @param across Character value. The name of the discrete feature-variable of interest.
#'
#' Hint: Character- or factor- variables are discrete variables. These functionally assign the barcode spots to distinct
#' groups or clusters (e.g. \emph{segment} or \emph{seurat_clusters}) whoose properties you might want
#' to compare against each other. Use \code{getFeatureNames()} to get an overview of the
#' features variables your spata-object contains.
#'
#' @param across_subset Character vector or NULL. Specify the particular groups or clusters of interest the feature-variable
#' specified in argument \code{across} contains. If set to NULL all of them are chosen.
#'
#' Hint: Use \code{getFeatureValues()} to obtain all available groups of a certain
#' feature-variable.
#'

across <- function(across, across_subset){}



#' @title de_df
#' @param de_df A data.frame containing information about differentially expressed genes.
#' This includes the numeric variables \emph{p_val, avg_logFC, p_val_adj} and the character
#' variables \emph{cluster, gene}.

pheatmap <- function(de_df){}
