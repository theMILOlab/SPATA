#' @title Spatial distance computation
#'
#' @description - -
#'
#' @inherit check_sample params
#' @inherit variable params
#' @inherit check_method params
#' @param plot_type Character value. Either \emph{'heatmap'} or \emph{'lineplot'}.
#' @param bins Numeric value. If \code{variable} is a numeric variable
#' it specifies the number of bins in which the the values are divided.
#' Given to \code{dplyr::ntile()}.
#' @inherit verbose params
#'
#' @return A named list with four slots:
#'
#' \itemize{
#'  \item{\emph{plot}: The lineplot or the heatmap.}
#'  \item{\emph{data}: The computed distance data.}
#'  \item{\emph{variable}: The variable specified.}
#'  \item{\emph{sample}: The sample specified.}}
#' @export
#'

calculateDistance <- function(object,
                              of_sample = "",
                              variable,
                              plot_type,
                              method_gs = "mean",
                              bins = 10,
                              verbose = TRUE){

  confuns::is_value(plot_type, "character", "plot_type")

  if(plot_type == "heatmap"){

    calculateDistanceHeatmap(
      object = object,
      of_sample = of_sample,
      variable = variable,
      method_gs = method_gs,
      bins = bins,
      verbose = verbose
    )

  } else if(plot_type == "lineplot"){

    calculateDistanceLineplot(
      object = object,
      of_sample = of_sample,
      variable = variable,
      method_gs = method_gs,
      bins = bins,
      verbose = verbose
    )

  } else {

    base::stop("Invalid input for argument 'plot_type'")

  }


}

#' @rdname calculateDistance
#' @export
calculateDistanceLineplot <- function(object,
                                      of_sample = "",
                                      variable,
                                      method_gs = "mean",
                                      bins = 10,
                                      verbose = TRUE){


  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_method(method_gs = method_gs)

  confuns::is_value(variable, "character", "variable")
  confuns::is_value(bins, "numeric", "bins")

  # adjusting check
  of_sample <- check_sample(object, of_sample, 1)

  variable <-
    check_variables(
      variables = variable,
      all_features = getFeatureNames(object),
      all_genes = getGenes(object, in_sample = of_sample),
      all_gene_sets = getGeneSets(object),
      max_slots = 1,
      max_length = 1
    )

  # -----

  # 2. Additional merging and control steps ---------------------------------

  coords_df <- getCoordinates(object, of_sample)

  feature_df <-
    joinWithVariables(object = object,
                      coords_df = coords_df,
                      variables = variable,
                      method_gs = method_gs,
                      average_genes = TRUE,
                      verbose = verbose)

  variable <- base::unlist(variable, use.names = FALSE)

  # check if variable needs to be converted into a categorical one

  variable_vls <- dplyr::pull(feature_df, var = {{variable}})

  if(base::is.numeric(variable_vls)){

    if(base::isTRUE(verbose)){
      base::message(glue::glue("Variable '{variable}' is numeric. Binning with bins = {bins}."))
    }

    feature_df <-
      dplyr::arrange(.data = feature_df, dplyr::desc(!!rlang::sym(variable))) %>%
      dplyr::mutate(!!variable := dplyr::ntile(x = !!rlang::sym(variable), n = bins)) %>%
      dplyr::select(barcodes, {{variable}})

  }

  # -----


  # 2. Data wrangling -------------------------------------------------------

  if(base::isTRUE(verbose)){base::message("Computing distance.")}

  bc_origin <- coords_df$barcodes
  bc_destination <- coords_df$barcodes

  df_distance <-
    tidyr::expand_grid(bc_origin, bc_destination) %>%
    dplyr::left_join(x = ., y = dplyr::select(coords_df, bc_origin = barcodes, xo = x, yo = y), key = "bc_origin") %>%
    dplyr::left_join(x = ., y = dplyr::select(coords_df, bc_destination = barcodes, xd = x, yd = y), key = "bc_destination") %>%
    dplyr::left_join(x = ., y = dplyr::select(feature_df, bc_origin = barcodes, cluster_origin = !!rlang::sym(variable)), key = "cluster_origin") %>%
    dplyr::left_join(x = ., y = dplyr::select(feature_df, bc_destination = barcodes, cluster_destination = !!rlang::sym(variable)), key = "cluster_destination") %>%
    dplyr::mutate(distance = base::round(sqrt((xd - xo)^2 + (yd - yo)^2), digits = 0)) %>%
    dplyr::filter(cluster_origin == cluster_destination & distance != 0) %>%
    dplyr::select(cluster_origin, distance) %>%
    dplyr::group_by(cluster_origin) %>%
    dplyr::summarise_all(median) %>%
    dplyr::arrange(dplyr::desc(distance)) %>%
    dplyr::mutate(group = "group",
                  x = dplyr::row_number()) %>%
    base::as.data.frame()

  # -----
  p <-
    ggplot2::ggplot(data = df_distance, mapping = ggplot2::aes(x = as.factor(x), y = distance)) +
    ggplot2::geom_point() +
    ggplot2::geom_path(mapping = ggplot2::aes(group = group)) +
    ggplot2::scale_x_discrete(labels = base::unique(df_distance$cluster_origin), breaks = df_distance$x) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = variable, y = "Relative Median Distance")

  if(base::isTRUE(verbose)){base::message("Done.")}

  base::return(
    list(
      "plot" = p,
      "data" = df_distance,
      "variable" = variable,
      "sample" = of_sample
    )
  )

}

#' @rdname calculateDistance
#' @export
calculateDistanceHeatmap <- function(object,
                                     of_sample = "",
                                     variable,
                                     method_gs = "mean",
                                     bins = 10,
                                     verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_method(method_gs = method_gs)

  confuns::is_value(bins, "numeric", "bins")
  confuns::is_value(variable, "character", "variable")

  # adjusting check
  of_sample <- check_sample(object, of_sample, 1)

  variable <-
    check_variables(
      variables = variable,
      all_features = getFeatureNames(object),
      all_genes = getGenes(object, in_sample = of_sample),
      all_gene_sets = getGeneSets(object),
      max_slots = 1,
      max_length = 1
    )

  # -----

  # 2. Additional merging and control steps ---------------------------------

  coords_df <- getCoordinates(object, of_sample)

  feature_df <-
    joinWithVariables(object = object,
                      coords_df = coords_df,
                      variables = variable,
                      method_gs = method_gs,
                      average_genes = TRUE,
                      verbose = verbose)

  variable <- base::unlist(variable, use.names = FALSE)

  # check if variable needs to be converted into a categorical one

  variable_vls <- dplyr::pull(feature_df, var = {{variable}})

  if(base::is.numeric(variable_vls)){

    if(base::isTRUE(verbose)){
      base::message(glue::glue("Variable '{variable}' is numeric. Binning with bins = {bins}."))
    }

    feature_df <-
      dplyr::arrange(.data = feature_df, dplyr::desc(!!rlang::sym(variable))) %>%
      dplyr::mutate(!!variable := dplyr::ntile(x = !!rlang::sym(variable), n = bins)) %>%
      dplyr::select(barcodes, {{variable}})

  }

  # -----


  # 2. Data wrangling -------------------------------------------------------

  if(base::isTRUE(verbose)){base::message("Computing distance.")}

  bc_origin <- coords_df$barcodes
  bc_destination <- coords_df$barcodes

  df_distance <-
    tidyr::expand_grid(bc_origin, bc_destination) %>%
    dplyr::left_join(x = ., y = dplyr::select(coords_df, bc_origin = barcodes, xo = x, yo = y), key = "bc_origin") %>%
    dplyr::left_join(x = ., y = dplyr::select(coords_df, bc_destination = barcodes, xd = x, yd = y), key = "bc_destination") %>%
    dplyr::left_join(x = ., y = dplyr::select(feature_df, bc_origin = barcodes, cluster_origin = !!rlang::sym(variable)), key = "cluster_origin") %>%
    dplyr::left_join(x = ., y = dplyr::select(feature_df, bc_destination = barcodes, cluster_destination = !!rlang::sym(variable)), key = "cluster_destination") %>%
    dplyr::mutate(distance = base::round(sqrt((xd - xo)^2 + (yd - yo)^2), digits = 0)) %>%
    dplyr::filter(distance != 0) %>%
    dplyr::select(cluster_origin, cluster_destination, distance) %>%
    dplyr::group_by(cluster_origin, cluster_destination) %>%
    dplyr::summarise_all(median) %>%
    reshape2::acast(data = ., cluster_origin ~ cluster_destination, value.var = "distance")

  # -----

  p <- pheatmap::pheatmap(as.matrix(df_distance),
                          color = base::rev(viridis::inferno(50)))

  if(base::isTRUE(verbose)){base::message("Done.")}

  base::return(
    list(
      "plot" = p,
      "data" = df_distance,
      "variable" = variable,
      "sample" = of_sample
    )
  )

}




#' @title Convert a numeric variable to a discrete one
#'
#' @description A wrapper around \code{dplyr::ntile()} to bin a numeric feature
#' into a discrete one.
#'
#' @param data A data.frame containing at least the character variables \emph{barcodes}
#' and \code{sample} and the numeric variable specified in \code{num_variable}.
#' @param num_variable Character value. The name of the numeric variable that you want
#' to convert.
#' @param discr_variable Character value. The name the new discrete variable wil have.
#' @param n_bins Numeric value. The number of bins you want to distribute the
#' values of \code{num_variable} to. Given to argument \code{n} of \code{dplyr::ntile()}.
#'
#' @return The data.frame specified in \code{data} with the additional discrete variable.
#' @export
#'

convertToDiscrete <- function(data,
                              num_variable,
                              discr_variable,
                              n_bins){

  confuns::is_value(num_variable, "character", "num_variable")
  confuns::is_value(discr_variable, "character", "discr_variable")

  check_list <-
    list(c("character"),
         c("character"),
         c("numeric", "integer", "double")) %>%
    magrittr::set_names(value = c("barcodes", "sample", num_variable))

  confuns::check_data_frame(
    df = data,
    var.class = check_list,
    ref = "data")

  res_data <-
    dplyr::mutate(.data = data,
                  !!discr_variable := dplyr::ntile(x = !!rlang::sym(num_variable), n = n_bins))

  res_data[[discr_variable]] <- base::as.character(res_data[[discr_variable]])

  base::return(res_data)

}

#' @title Clustering with igraph
#'
#' @param cor_mtr A correlation matrix.
#' @param num_nn Numeric value. The maximum number of nearest neighbours to compute.
#' The default value is set to the smaller of the number of columns in data. Given
#' to \code{RANN::nn2()} as input for argument \code{k}.
#' @param color_nn Numeric value.
#' @param k Numeric value.
#'
#' @return Cluster results.
#' @export
#'

findIgraphCluster <- function(cor_mtr,
                              num_nn = 20){

  # 1. Control  -------------------------------------------------------------

  if(!base::is.matrix(cor_mtr) ||
     !base::is.numeric(cor_mtr) ||
     !base::identical(base::colnames(cor_mtr), base::rownames(cor_mtr))){

    base::stop("Input for argument 'cor_mtr' needs to be a correlation matrix.
               In order to proceed. Make sure that row- and columnames are identical and that
               all values numeric.")

  }

  # -----

  # 2. Data wrangling

  nearest <- RANN::nn2(data = cor_mtr,
                       k = num_nn,
                       treetype = "bd",
                       searchtype = "priority")

  edges_n <-
    reshape::melt(base::t(nearest$nn.idx[, 1:num_nn])) %>%
    magrittr::set_colnames(value = c("B", "A", "C")) %>%
    dplyr::select(A, B = C) %>%
    dplyr::transmute(
      V1 = base::pmin(A,B),
      V2 = base::pmax(A,B),
      weight = 1
    ) %>%
    base::unique()

  graph_out <-
    dplyr::transmute(.data = edges_n,
                     V1 = base::rownames(cor_mtr)[edges_n$V1],
                     V2 = base::rownames(cor_mtr)[edges_n$V2],
                     weight = weight
    ) %>%
    igraph::graph.data.frame(d = ., directed = FALSE) %>%
    igraph::cluster_louvain()

  cluster_assign <- base::factor(graph_out$membership, levels = base::sort(base::unique(graph_out$membership)))

  cluster_out <-
    base::data.frame(
      id = base::rownames(cor_mtr),
      cluster = cluster_assign
    ) %>%
    magrittr::set_rownames(value = base::rownames(cor_mtr))

  # -----

  base::return(cluster_out)

}


#' @rdname computeIgraphCluster
#' @export
plotIgraphHeatmap <- function(cor_mtr,
                              num_nn = 20,
                              color_nn,
                              k = 6){

  # 1. Control --------------------------------------------------------------

  confuns::is_value(x = color_nn, mode = "numeric", ref = "color_nn")

  # -----


  # 2. Data wrangling -------------------------------------------------------

  cluster <- findIgraphCluster(cor_mtr = cor_mtr, num_nn = num_nn)

  prel_breaks <- stats::quantile(x = cor_mtr, probs = base::seq(0, 1, length.out = color_nn))
  mat_breaks <- prel_breaks[!base::duplicated(prel_breaks)]

  my_breaks <- c(seq(0, 0.1, length.out = 10),
                 seq(0.11, 1, length.out = 20))

  cl <- pheatmap::pheatmap(mat = cor_mtr,
                           color = viridisLite::viridis(base::length(mat_breaks) - 1),
                           annotation_col = cluster,
                           show_colnames = F,
                           cutree_rows = base::max(base::as.numeric(cluster$cluster)),
                           cutree_cols = base::max(base::as.numeric(cluster$cluster)),
                           fontsize = 2,
                           border_color = NA,
                           breaks = my_breaks)

  cluster <- data.frame(Cluster = base::sort(stats::cutree(cl$tree_row, k = k)))
  cluster$Cluster <- paste0("Cluster_",cluster$cluster)
  cluster$ID <- base::rownames(cluster)
  cluster <- cluster[,2:1]

  # -----

  # Return heatmap ----------------------------------------------------------

  if(base::isTRUE(verbose)){base::message("Plotting heatmap. This might take a few seconds.")}

  pheatmap::pheatmap(mtr = cor_mtr,
                     color = viridisLite::viridis(length(mat_breaks) - 1),
                     annotation_col = cluster,
                     show_colnames = F,
                     cutree_rows = base::length(base::unique(cluster$Cluster)),
                     cutree_cols = base::length(base::unique(cluster$Cluster)),
                     fontsize = 4,
                     border_color = NA,
                     breaks = my_breaks)

}


#' @title Cluster sample via monocle3
#'
#' @description Assign barcode spots to clusters according to different clustering
#' algorithms.
#'
#' @inherit check_object params
#' @inherit check_monocle_input params details
#' @param prefix Character value. Clustering algorithms often return only numbers as
#' names for the clusters they generate. If you want to these numbers to have a certain
#' prefix (like \emph{'Cluster'}, the default) you can specify it with this argument.
#'
#' @details This functions is a wrapper around monocle3-cluster algorithms which
#' take several options for dimensional reduction upon which the subsequent clustering bases.
#' It iterates over all specified methods and returns a tidy data.frame in which each row represents
#' one barcode-spot uniquely identified by the variable \emph{barcodes} and in which every other variable
#' about the cluster belonging the specified combination of methods returned. E.g.:
#'
#' A call to `findMonocleClusters()` with
#'
#' \itemize{
#'  \item{\code{preprocess_method} set to \emph{'PCA'} }
#'  \item{\code{reduction_method} set to \emph{c('UMAP', 'PCA')}}
#'  \item{\code{'leiden'}, \code{k} set to \emph{5}}
#'  }
#'
#' will return a data.frame of the following variables:
#'
#' \itemize{
#'  \item{\emph{barcodes}}
#'  \item{\emph{mncl_cluster_UMAP_leiden_k5}}
#'  \item{\emph{mncl_cluster_PCA_leiden_k5}}
#'  }
#'
#' Due to the \emph{barcodes}-variable it can be easily joined to your-spata object via `addFeature()`.
#' and thus be made available for all spata-functions.
#'
#' @note With respect to the arguments \code{preprocess_method},
#' \code{reduction_method} and \code{cluster_method}:
#'
#' If a vector of character values is prvodided instead of a single character value (string)
#' the function will iterate over all inputs via a for-loop to compute all
#' valid combinations.
#'
#' @return A tidy spata-data.frame
#' @export
#'

findMonocleClusters <- function(object,
                                preprocess_method = c("PCA", "LSI"),
                                reduction_method = c("UMAP", "tSNE", "PCA", "LSI"),
                                cluster_method = c("leiden", "louvain"),
                                k = 20,
                                num_iter = 5,
                                prefix = "Cluster",
                                verbose = TRUE){

  check_object(object)

  check_monocle_input(preprocess_method = preprocess_method,
                      reduction_method = reduction_method,
                      cluster_method = cluster_method,
                      k = k,
                      num_iter = num_iter)

  if(base::isTRUE(verbose)){base::message("Creating 'cell_data_set'-object.")}

  expression_matrix <- base::as.matrix(object@data@counts)

  gene_metadata <- data.frame(gene_short_name = base::rownames(expression_matrix))
  base::rownames(gene_metadata) <- base::rownames(expression_matrix)

  cell_metadata <- data.frame(object@fdata)
  base::rownames(cell_metadata) <- object@fdata$barcodes

  cds <- monocle3::new_cell_data_set(
    expression_data = expression_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata)

  for(p in base::seq_along(preprocess_method)){

    if(base::isTRUE(verbose)){

      base::message(glue::glue("Preprocessing cells with method {p}/{base::length(preprocess_method)} '{preprocess_method[p]}'"))

    }

    cds <- monocle3::preprocess_cds(cds, method = preprocess_method[p])

  }


  for(p in base::seq_along(preprocess_method)){

    base::message(glue::glue("Using preprocess method '{preprocess_method[p]}':"))

    for(r in base::seq_along(reduction_method)){

      base::message(glue::glue("Reducing dimensions with reduction method {r}/{base::length(reduction_method)}: '{reduction_method[r]}' "))

      if(reduction_method[r] == "LSI" && preprocess_method[p] != "LSI"){

        base::message(glue::glue("Ignoring invalid combination. reduction-method: '{reduction_method[r]}' &  preprocess-method: '{preprocess}'"))

      } else if(reduction_method[r] == "PCA" && preprocess_method[p] != "PCA") {

        base::message(glue::glue("Ignoring invalid combination. reduction-method: '{reduction_method[r]}' &  preprocess-method: '{preprocess}'"))

      } else {

        cds <- monocle3::reduce_dimension(cds = cds, reduction_method = reduction_method[r], preprocess_method = preprocess_method[p], verbose = FALSE)

      }

    }

  }

  cluster_df <- data.frame(barcodes = getBarcodes(object = object))

  for(r in base::seq_along(reduction_method)){

    if(base::isTRUE(verbose)){

      base::message(glue::glue("Using reduction method {reduction_method[r]}:"))

    }

    for(c in base::seq_along(cluster_method)){

      if(base::isTRUE(verbose)){

        base::message(glue::glue("Clustering barcode-spots with method {c}/{base::length(cluster_method)}: {cluster_method[c]}"))

      }

      cds <- monocle3::cluster_cells(cds = cds,
                                     reduction_method = reduction_method[r],
                                     k = k,
                                     num_iter = num_iter,
                                     cluster_method = cluster_method[c],
                                     verbose = FALSE)

      cluster_name <- stringr::str_c("cluster", cluster_method[c], reduction_method[r],base::paste0("k", k), sep = "_")

      cluster_df <-
        monocle3::clusters(x = cds, reduction_method = reduction_method[r]) %>%
        base::as.data.frame() %>%
        tibble::rownames_to_column(var = "barcodes") %>%
        magrittr::set_colnames(value = c("barcodes", cluster_name)) %>%
        dplyr::left_join(x = cluster_df, y = ., by = "barcodes") %>%
        tibble::as_tibble()

    }

  }

  cluster_df <- purrr::map_df(.x = dplyr::select(cluster_df, -barcodes),
                              .f = function(i){

                                i <- stringr::str_c(prefix, i, sep = "")

                                if(base::is.factor(i)){

                                  S4Vectors::unfactor(i) %>%
                                    base::as.character()

                                } else {

                                  base::as.character(i)

                                }

                              }) %>%
    dplyr::mutate(barcodes = cluster_df$barcodes)

  if(base::isTRUE(verbose)){base::message("Done.")}

  base::return(cluster_df)

}

