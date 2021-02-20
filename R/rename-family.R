#' @title Rename features
#'
#' @description Allows to rename features stored inside the @@fdata slot.
#'
#' @inherit check_sample params
#' @param ... The features to be renamed specified according to the following
#' syntax: \emph{'new_feature_name'} \code{=} \emph{'old_feature_name'}.
#'
#' @return An upated spata-object.
#' @export
#'
#' @examples #Not run:
#'
#'  object <- renameFeatures(object, "seurat_clusters_new" = "seurat_clusters")
#'

renameFeatures <- function(object, ..., of_sample = NA){

  check_object(object)
  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  rename_input <- confuns::keep_named(c(...))

  valid_feature_names <-
    check_features(object = object, features = rename_input, of_sample = of_sample)

  valid_rename_input <- rename_input[rename_input %in% valid_feature_names]

  #assign("valid_rename_input", value = valid_rename_input, envir = .GlobalEnv)

  # rename feature df
  feature_df <-
    getFeatureDf(object, of_sample = of_sample) %>%
    dplyr::rename(!!! valid_rename_input)

  # rename dea list
  dea_list <- object@dea[[of_sample]]

  dea_names <- base::names(dea_list)

  if(!base::is.null(dea_names)){

    dea_names <- valid_rename_input[valid_rename_input %in% dea_names]

    if(base::length(dea_names) >= 1){

      for(dea_name in dea_names){

        base::names(dea_list)[base::names(dea_list) == dea_name] <-
          base::names(dea_names)[dea_names == dea_name]

      }

      object@dea[[of_sample]] <- dea_list

    }

  }


  object <- setFeatureDf(object, feature_df = feature_df, of_sample = of_sample)

  base::return(object)

}



#' @title Rename cluster/group names
#'
#' @description Allows to rename groups within a discrete grouping variable (such as
#' cluster variables) of the feature data in slot @@fdata.
#'
#' @inherit check_sample params
#' @param discrete_feature Character value. The grouping variable of interest.
#' @param ... The groups to be renamed specified according to the following
#' syntax: \emph{'new_group_name'} \code{=} \emph{'old_group_name'}.
#'
#' @return An updated spata-object.
#' @export
#'
#' @examples #Not run:
#'
#'  object <-
#'     renameGroups(object = spata_object,
#'                  discrete_feature = "seurat_clusters",
#'                  "first_new_group" = "1",
#'                  "sec_new_group" = "2")
#'
#'

renameGroups <- function(object, discrete_feature, ..., of_sample = NA){

  check_object(object)
  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  discrete_feature <-
    check_features(object, features = discrete_feature, valid_classes = c("factor"))

  rename_input <- confuns::keep_named(c(...))

  if(base::length(rename_input) == 0){

    msg <- renaming_hint

    confuns::give_feedback(
      msg = msg,
      fdb.fn = "stop"
    )

  }

  feature_df <- getFeatureDf(object, of_sample = of_sample)

  valid_rename_input <-
    confuns::check_vector(
      input = base::unname(rename_input),
      against = base::levels(feature_df[[discrete_feature]]),
      fdb.fn = "warning",
      ref.input = "groups to rename",
      ref.against = glue::glue("all groups of feature '{discrete_feature}'. ({renaming_hint})")
    )

  rename_input <- rename_input[rename_input %in% valid_rename_input]

  # rename feature
  renamed_feature_df <-
    dplyr::mutate(
      .data = feature_df,
      {{discrete_feature}} := forcats::fct_recode(.f = !!rlang::sym(discrete_feature), !!!rename_input)
    )

  # rename dea list
  dea_list <- object@dea[[of_sample]][[discrete_feature]]

  if(!base::is.null(dea_list)){

    object@dea[[of_sample]][[discrete_feature]] <-
      purrr::map(
        .x = dea_list,
        .f = function(method){

          new_df <-
            dplyr::mutate(
              .data = method$data,
              {{discrete_feature}} := forcats::fct_recode(.f = !!rlang::sym(discrete_feature), !!!rename_input)
            )

          list(data = new_df, adjustments = method$adjustments)

        }
      ) %>%
      purrr::set_names(nm = base::names(dea_list))


  }

  object <- setFeatureDf(object, feature_df = renamed_feature_df, of_sample = of_sample)

  base::return(object)

}




