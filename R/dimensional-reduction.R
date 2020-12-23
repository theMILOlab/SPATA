
#' @title Run Principal Component Analysis
#'
#' @description Takes the expression matrix of choice and passes it to
#' \code{irlba::prcomp_irlba()}.
#'
#' @inherit check_sample params
#' @inherit getExpressionMatrix params
#' @pca_comp Numeric value. Denotes the number of principal components to be computed.
#'
#' @inherit runDimRed_dummy return
#'
#' @export

runPca <- function(object, of_sample = "", pca_comp = 30, mtr_name = NULL){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  pca_res <- runPca2(object = object,
                     pca_comp = pca_comp,
                     mtr_name = mtr_name)

  expr_mtr <- getExpressionMatrix(object, of_sample = of_sample, mtr_name = mtr_name)

  pca_df <-
    base::as.data.frame(x = pca_res[["x"]]) %>%
    dplyr::mutate(barcodes = base::colnames(expr_mtr), sample = {{of_sample}}) %>%
    dplyr::select(barcodes, sample, dplyr::everything())

  object <- setPcaDf(object = object, pca_df = pca_df)

  base::return(object)

}

#' @rdname runPca
#' @export
runPca2 <- function(object, of_sample = "", pca_comp = 30, mtr_name = NULL){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  expr_mtr <- getExpressionMatrix(object, of_sample = of_sample, mtr_name = mtr_name)

  pca_res <- irlba::prcomp_irlba(x = base::t(expr_mtr), n = pca_comp)

  base::return(pca_res)

}


#' @title Run t-Stochastic Neighbour Embedding
#'
#' @param object
#' @param of_sample
#' @param tsne_perplexity
#'
#' @return
#' @export

runTsne <- function(object, of_sample = "", tsne_perplexity = 30){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  tsne_res <- runTsne2(object = object,
                       of_sample = of_sample,
                       tsne_perplexity = tsne_perplexity)

  pca_mtr <- getPcaMtr(object = object, of_sample = of_sample)

  tsne_df <-
    base::data.frame(barcodes = base::rownames(pca_mtr),
                     sample = of_sample,
                     tsne1 = tsne_res$Y[,1],
                     tsne2 = tsne_res$Y[,2])

  object <- setTsneDf(object = object, tsne_df = tsne_df, of_sample = of_sample)

}

#' @rdname runTsne
#' @export
runTsne2 <- function(object, of_sample = "", tsne_perplexity = 30){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  pca_mtr <- getPcaMtr(object = object, of_sample = of_sample)

  tsne_res <- Rtsne::Rtsne(pca_mtr, perplexity = tsne_perplexity)

  base::return(tsne_res)

}


#' Title
#'
#' @param object
#' @param of_sample
#' @param ...
#'
#' @return
#' @export

runUmap <- function(object, of_sample = "", ...){

  check_object(object)

  of_sample <-
    check_sample(object = object, of_sample = of_sample, of.length = 1)

  umap_res <-
    runUmap2(object = object, of_sample = of_sample, ...)

  pca_mtr <-
    getPcaMtr(object = object, of_sample = of_sample)

  umap_df <-
    base::data.frame(
      barcodes = base::rownames(pca_mtr),
      sample = of_sample,
      umap1 = umap_res$layout[,1],
      umap2 = umap_res$layout[,2]
    )

  object <- setUmapDf(object = object, umap_df = umap_df, of_sample = of_sample)

  base::return(object)

}

#' @rdname runUmap
#' @export
runUmap2 <- function(object, of_sample = "", ...){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  pca_mtr <- getPcaMtr(object = object, of_sample = of_sample)

  umap_res <- umap::umap(d = pca_mtr, ...)

  base::return(umap_res)

}


