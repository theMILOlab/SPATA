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



#' @title Clustering with igraph
#'
#' @param cor.mtr A correlation matrix.
#' @param num.nn Numeric value.
#' @param color.nn Numeric value.
#'
#' @return Cluster results.
#' @export
#'

igraph_cluster <- function(cor.mtr, num.nn = 20){

  # 1. Control  -------------------------------------------------------------

  if(!base::is.matrix(cor.mtr) ||
     !base::is.numeric(cor.mtr) ||
     !base::identical(base::colnames(cor.mtr), base::rownames(cor.mtr))){

    base::stop("Input for argument 'cor.mtr' needs to be a correlation matrix.
               In order to proceed. Make sure that row- and columnames are identical and that
               all values numeric.")

  }

  # -----

  # 2. Data wrangling

  nearest <- RANN::nn2(data = cor.mtr,
                       k = num.nn,
                       treetype = "bd",
                       searchtype = "priority")

  edges_n <-
    reshape::melt(base::t(nearest$nn.idx[, 1:num.nn])) %>%
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
                     V1 = base::rownames(cor.mtr)[edges_n$V1],
                     V2 = base::rownames(cor.mtr)[edges_n$V2],
                     weight = weight
    ) %>%
    igraph::graph.data.frame(d = ., directed = FALSE) %>%
    igraph::cluster_louvain()

  cluster_assign <- base::factor(graph_out$membership, levels = base::sort(base::unique(graph.out$membership)))

  cluster_out <-
    base::data.frame(
      id = base::rownames(cor.mtr),
      cluster = cluster_assign
    ) %>%
    magrittr::set_rownames(value = base::rownames(cor.mtr))

  # -----

  base::return(cluster_out)

}

#' @rdname igraph_cluster
#' @export
plot_igraph_heatmap <- function(cor.mtr,
                                num.nn,
                                color.nn){


  # 1. Control --------------------------------------------------------------

  confuns::is_value(x = color.nn, mode = "numeric", ref = "color.nn")

  # -----


  # 2. Data wrangling -------------------------------------------------------

  cluster <- igrap_cluster(cor.mtr = cor.mtr, num.nn = num.nn)

  prel_breaks <- stats::quantile(x = cor.mtr, probs = base::seq(0, 1, length.out = color.nn))
  mat_breaks <- prel_breaks[!base::duplicated(breaks)]

  my_breaks <- c(seq(0, 0.1, length.out = 10),
                 seq(0.11, 1, length.out = 20))


  cl <- pheatmap::pheatmap(mat = cor.mtr,
                           color = viridisLite::viridis(n, base::length(mat_breaks) - 1),
                           annotation_col = cluster,
                           show_colnames = F,
                           cutree_rows = base::max(base::as.numeric(cluster$Cluster)),
                           cutree_cols = base::max(base::as.numeric(cluster$Cluster)),
                           fontsize = 2,
                           border_color = NA,
                           breaks = my_breaks)

  cluster <- data.frame(Cluster = base::sort(stats::cutree(cl$tree_row, k = 6)))
  cluster$Cluster <- paste0("Cluster_",cluster$Cluster)
  cluster$ID - base::rownames(cluster)
  cluster <- cluster[,2:1]

  # -----

  # Return heatmap ----------------------------------------------------------

  pheatmap::pheatmap(mtr = cor.mtr,
                     color = viridisLite::viridis(length(mat_breaks) - 1),
                     annotation_col = cluster,
                     show_colnames = F,
                     cutree_rows = base::length(base::unique(cluster$Cluster)),
                     cutree_cols = base::length(base::unique(cluster$Cluster)),
                     fontsize = 4,
                     border_color = NA,
                     breaks = my_breaks)


}
