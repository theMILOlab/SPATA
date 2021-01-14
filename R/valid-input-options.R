

#' @title Obtain valid argument inputs
#'
#' @description These function simply return valid input options
#' for recurring arguments.
#'
#' @return Character vectors or named lists of such.
#' @export
#'

validActivationFunctions <- function(){

  base::return(activation_fns)

}

#' @rdname validActivationFunctions
#' @export
validColorPanels <- function(){

  confuns::all_colorpanels()

}

#' @rdname validActivationFunctions
#' @export
validColorSpectra <- function(){

  confuns::all_colorspectra()

}

#' @rdname validActivationFunctions
#' @export
validDeAnalysisMethods <- function(){

  base::return(de_methods)

}

#' @rdname validActivationFunctions
#' @export
validDefaultInstructionSlots <- function(){

  base::return(methods::slotNames(methods::new("default_instructions")))

}

#' @rdname validActivationFunctions
#' @export
validDimRedMethods <- function(){

  base::return(gene_set_emthods)

}

#' @rdname validActivationFunctions
#' @export
validDirectoryInstructionSlots <- function(){

  base::return(directory_options)

}

#' @rdname validActivationFunctions
#' @export
validHierarchicalClusterMethods <- function(){

  base::return(hclust_methods)

}

#' @rdname validActivationFunctions
#' @export
validPatternRecognitionMethods <- function(){

  base::return(pr_methods)

}


#' @rdname validActivationFunctions
#' @export
validPadjMethods <- function(){

  base::return(stats::p.adjust.methods)

}




