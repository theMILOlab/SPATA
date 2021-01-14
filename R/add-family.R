

# Autoencoder related -----------------------------------------------------

#' Title
#'
#' @inherit check_object params
#' @param set_up_list A named list with slots \code{$activation, $bottleneck, $dropout, $epochs, $layers}.
#'
#' @return A spata-object.
#' @export
#'

addAutoencoderSetUp <- function(object, mtr_name, set_up_list, of_sample = NA){

  check_object(object)
  confuns::is_list(set_up_list)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  object@information$autoencoder[[of_sample]][["nn_set_ups"]][[mtr_name]] <- set_up_list

  base::return(object)

}

# -----
# Feature related ---------------------------------------------------------


#' @title Add a new feature
#'
#' @description Adds a new variable to the objects feature data.
#'
#' @inherit check_object
#' @param overwrite Logical. If the specified feature name already exists in the
#' current spata-object this argument must be set to TRUE in order to overwrite it.
#' @param key_variable Character value. Either \emph{'barcodes'} or \emph{'coordinates'}.
#' If set to \emph{'coordinates'} the \code{feature_df}-input must contain numeric x- and
#' y- variables.
#'
#' Key variables are variables in a data.frame that uniquely identify each observation -
#' in this case each barcode-spot. In SPATA the barcode-variable is a key-variable on its own,
#' x- and y-coordinates work as key-variables if they are used combined.
#'
#' @inherit check_feature_df params
#'
#' @details Eventually the new feature will be joined via \code{dplyr::left_join()} over the
#' key-variables \emph{barcodes} or \emph{x} and \emph{y}. Additional steps secure the joining process.
#'
#' @return An updated spata-object.
#' @export

addFeatures <- function(object,
                        feature_names,
                        feature_df,
                        key_variable = "barcodes",
                        overwrite = FALSE,
                        of_sample = ""){

  # 1. Control --------------------------------------------------------------
  check_object(object)

  confuns::is_vec(x = feature_names, mode = "character")
  confuns::is_value(x = key_variable, mode = "character")

  confuns::check_one_of(input = key_variable, against = c("barcodes", "coordinates"), ref.input = "argument 'key_variable'")


  feature_names <- confuns::check_vector(
    input = feature_names,
    against = base::colnames(feature_df),
    verbose = TRUE,
    ref.input = "specified feature names",
    ref.against = "variables of provided feature data.frame")

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  if(key_variable  == "barcodes"){

    confuns::check_data_frame(df = feature_df,
                              var.class = list("barcodes" = "character"),
                              ref = "feature_df")

  } else if(key_variable == "coordinates"){

    confuns::check_data_frame(df = feature_df,
                              var.class = list(
                                "x" = c("numeric", "integer", "double"),
                                "y" = c("numeric", "integer", "double")
                              ),
                              ref = "feature_df")

  }

  # 2. Extract and compare --------------------------------------------------

  existing_fnames <- getFeatureNames(object = object, of_sample = of_sample)

  # throw error if there intersecting feature names and overwrite is FALSE
  if(base::any(feature_names %in% existing_fnames) && !base::isTRUE(overwrite)){

    found <- feature_names[feature_names %in% existing_fnames]

    if(base::length(found) > 1){

      ref <- c("are", "them")

    } else {

      ref <- c("is", "it")

    }

    found_ref <- stringr::str_c(found, collapse = "', '")

    base::stop(glue::glue("Specified feature names '{found_ref}' {ref[1]} already present in current feature data. Set overwrite to TRUE in order to overwrite {ref[2]}."))

    # discard existing, intersecting feature names if overwrite is TRUE
  } else if(base::any(feature_names %in% existing_fnames) && base::isTRUE(overwrite)){

    overwrite_features <- existing_fnames[existing_fnames %in% feature_names]

    fdata <-
      getFeatureDf(object, of_sample = of_sample) %>%
      dplyr::select(-dplyr::all_of(overwrite_features))

    #
  } else {

    fdata <- getFeatureDf(object, of_sample = of_sample)

  }

  # join over coordinates
  if(key_variable == "coordinates"){

    coords_df <-
      getCoordsDf(object, of_sample = of_sample) %>%
      purrr::map_at(.at = c("x", "y"), .f = function(i){ base::round(i, digits = 0)}) %>%
      purrr::map_df(.f = function(i){ base::return(i) })

    fdata <- dplyr::left_join(x = fdata, y = coords_df, key = "barcodes")

    feature_df <-
      purrr::map_at(.x = feature_df, .at = c("x", "y"), .f = function(i){ base::round(i, digits = 0)}) %>%
      purrr::map_df(.f = function(i){ base::return(i) }) %>%
      dplyr::left_join(y = coords_df, key = c("x", "y"))

    # feedback about how many barcode-spots can be joined
    barcodes_feature_df <- feature_df$barcodes
    barcodes_obj <- fdata$barcodes

    n_bc_feat <- base::length(barcodes_feature_df)
    n_bc_obj <- base::length(barcodes_obj)

    if(!base::all(barcodes_obj %in% barcodes_feature_df)){

      not_found <- barcodes_obj[!barcodes_obj %in% barcodes_feature_df]
      n_not_found <- base::length(not_found)

      if(n_not_found == n_bc_obj){base::stop("Did not find any barcode-spots of the specified object in input for 'feature_df'.")}

      base::warning(glue::glue("Only {n_bc_feat} barcode-spots of {n_bc_obj} were found in 'feature_df'. Not found barcode-spots obtain NAs for all features to be joined."))

    }

    new_feature_df <-
      dplyr::left_join(x = fdata,
                       y = feature_df[,c("x", "y", feature_names)],
                       by = c("x", "y")) %>%
      dplyr::select(-x, -y)

    object <- setFeatureDf(object = object, feature_df = new_feature_df, of_sample = of_sample)

    # join over coordinates
  } else if(key_variable == "barcodes") {

    # feedback about how many barcode-spots can be joined
    barcodes_feature_df <- feature_df$barcodes
    barcodes_obj <- fdata$barcodes

    n_bc_feat <- base::length(barcodes_feature_df)
    n_bc_obj <- base::length(barcodes_obj)

    if(!base::all(barcodes_obj %in% barcodes_feature_df)){

      not_found <- barcodes_obj[!barcodes_obj %in% barcodes_feature_df]
      n_not_found <- base::length(not_found)

      if(n_not_found == n_bc_obj){base::stop("Did not find any barcode-spots of the specified object in input for 'feature_df'.")}

      base::warning(glue::glue("Only {n_bc_feat} barcode-spots of {n_bc_obj} were found in 'feature_df'. Not found barcode-spots obtain NAs for all features to be joined."))

    }

    new_feature_df <-
      dplyr::left_join(x = fdata,
                       y = feature_df[,c("barcodes", feature_names)],
                       by = "barcodes")

    object <- setFeatureDf(object = object, feature_df = new_feature_df, of_sample = of_sample)

  }

  base::return(object)

}

# -----

# Gene set related --------------------------------------------------------

#' @title Add a new gene set
#'
#' @description Stores a new gene set in the spata-object.
#'
#' @inherit check_object
#' @param class_name Character value. The class the gene set belongs to..
#' @param gs_name Character value. The name of the new gene set.
#' @param overwrite Logical. Overwrites existing gene sets with the same \code{class_name} -
#' \code{gs_name} combination.
#'
#' @inherit check_genes params
#'
#' @return An updated spata-object.
#'
#' @details Combines \code{class_name} and \code{gs_name} to the final gene set name.
#' Gene set classes and gene set names are separated by '_' and handled like this
#' in all additional gene set related functions which is why \code{class_name} must
#' not contain any '_'.
#'
#' @export

addGeneSet <- function(object,
                       class_name,
                       gs_name,
                       genes,
                       overwrite = FALSE){

  # lazy control
  check_object(object)

  # adjusting control
  genes <- check_genes(object, genes = genes)

  if(base::any(!base::sapply(X = list(class_name, gs_name, genes),
                             FUN = base::is.character))){

    base::stop("Arguments 'class_name', 'gs_name' and 'genes' must be of class character.")

  }

  if(base::length(class_name) != 1 | base::length(gs_name) != 1){

    base::stop("Arguments 'class_name' and 'gs_name' must be of length one.")

  }

  if(stringr::str_detect(string = class_name, pattern = "_")){

    base::stop("Invalid input for argument 'class_name'. Must not contain '_'.")

  }

  name <- stringr::str_c(class_name, gs_name, sep = "_")

  # make sure not to overwrite if overwrite == FALSE
  if(name %in% object@used_genesets$ont && base::isFALSE(overwrite)){

    base::stop(stringr::str_c("Gene set '", name, "' already exists.",
                              " Set argument 'overwrite' to TRUE in order to overwrite existing gene set."))

  } else if(name %in% object@used_genesets$ont && base::isTRUE(overwrite)) {

    object <- discardGeneSets(object, gs_names = name)

  }

  # add gene set
  object@used_genesets <-
    dplyr::add_row(
      .data = object@used_genesets,
      ont = base::rep(name, base::length(genes)),
      gene = genes
    )

  base::return(object)

}


#' @rdname addGeneSet
#' @export
addGeneSetsInteractive <- function(object){

  check_object(object)

  new_object <-
    shiny::runApp(
      shiny::shinyApp(
        ui = function(){

          shiny::fluidPage(
            moduleAddGeneSetsUI(id = "add_gs"),
            shiny::HTML("<br><br>"),
            shiny::actionButton("close_app", label = "Close application")
          )

        },
        server = function(input, output, session){

          module_return <-
            moduleAddGeneSetsServer(id = "add_gs",
                                    object = object)


          oe <- shiny::observeEvent(input$close_app, {

            shiny::stopApp(returnValue = module_return())

          })

        }
      )
    )

  base::return(new_object)

}


#' Discard gene sets
#'
#' @inherit check_object
#' @param gs_names Character vector. The gene sets to be discarded.
#'
#' @return An updated spata-object.
#' @export

discardGeneSets <- function(object, gs_names){

  # lazy control
  check_object(object)

  # adjusting control
  gs_names <- check_gene_sets(object, gene_sets = gs_names)

  # discard gene sets
  object@used_genesets <-
      dplyr::filter(object@used_genesets,
                    !ont %in% gs_names)

  return(object)

}


# -----




# Data related ------------------------------------------------------------

#' @title Discard an expression matrix
#'
#' @description Discards the expression matrix of choice.
#'
#' @inherit getExpressionMatrix params
#'
#' @return An updated spata-object.
#' @export
#'

discardExpressionMatrix <- function(object, mtr_name, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  all_mtr_names <- getExpressionMatrixNames(object, of_sample = of_sample)

  confuns::check_one_of(
    input = mtr_name,
    against = all_mtr_names,
    ref.input = "argument 'mtr_name'"
  )

  object <- addExpressionMatrix(object = object,
                                expr_mtr = NULL,
                                mtr_name = mtr_name,
                                of_sample = of_sample)

  base::message(glue::glue("Expression matrix '{mtr_name}' discarded."))

  # feedback if discarded matrix was denoted as active matrix
  if(mtr_name == getActiveMatrixName(object, of_sample = of_sample)){

    base::warning(glue::glue("Expression matrix '{mtr_name}' was denoted as the active matrix. Make sure to denote a new one with 'setActiveExpressionMatrix()'"))

  }

  # feedback if no expression matrix left
  remaining_mtr_names <- all_mtr_names[all_mtr_names != mtr_name]

  if(base::is.null(remaining_mtr_names) | base::identical(remaining_mtr_names, base::character(0))){

    base::warning("There are no expression matrices left in the provided spata-object. Make sure to add one with 'addExpressionMatrix()'.")

  }

  base::return(object)

}


# Trajectory related ------------------------------------------------------


addTrajectoryObject <- function(object,
                                trajectory_name,
                                trajectory_object,
                                of_sample = NA){

  # 1. Control --------------------------------------------------------------

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::is_value(x = trajectory_name, mode = "character")

  base::stopifnot(methods::is(trajectory_object, class2 = "spatial_trajectory"))

  if(trajectory_name %in% getTrajectoryNames(object, of_sample = of_sample)){

    base::stop(glue::glue("Trajectory name '{trajectory_name}' is already taken."))

  } else if(trajectory_name == ""){

    base::stop("'' is not a valid trajectory name.")

  }

  # 2. Set trajectory object ------------------------------------------------

  object@trajectories[[of_sample]][[trajectory_name]] <-
    trajectory_object

  base::return(object)

}


# ------
