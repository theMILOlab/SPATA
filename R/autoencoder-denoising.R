#' @title Deep Autoencoder for Denoising
#'
#' @description This function constructs and uses a neural network to denoies
#' expression levels spatially.
#'
#' @inherit check_object params
#' @param layers Numeric vector of length 3. Denotes the number of neurons in the three hidden layers.
#'  (default = c(128, 64, 32))
#' @param bottleneck Numeric value. Denotes the number of bottleneck neurons. (default = 32)
#' @param dropout Numeric value. Denotes the dropout. (default = 0.1)
#' @param activation Character value. Denotes the activation function. (default = \emph{'relu'})
#' @param epochs Numeric value. Denotes the epochs of the neural network. (default = 20)
#' @param val_plot Logical. If set to TRUE a scatter plot of the result is displayed.
#' @param val_genes Character vector of length two. Denoting the genes to be used for the validation plot.
#'
#' @return A spata-object containing the denoised expression matrix.
#'
#' @importFrom Seurat ScaleData
#'
#' @export

autoEncoderDenoising <- function(object,
                                 layers = c(128, 64, 32),
                                 bottleneck = 16,
                                 dropout = 0.1,
                                 activation = "relu",
                                 epochs = 20,
                                 val_plot = FALSE,
                                 val_genes = c("CD163", "CD68"),
                                 verbose = TRUE){


# 1. Control --------------------------------------------------------------

  check_object(object)

  confuns::are_values(c("bottleneck", "dropout", "epochs"), mode = "numeric")
  confuns::are_values(c("activation"), mode = "character")
  confuns::are_values(c("val_plot", "verbose"), mode = "logical")

  confuns::is_vec(x = layers, mode = "numeric", of.length = 3)

  if(base::isTRUE(val_plot)){

    confuns::is_vec(x = val_genes, mode = "character", of.length = 2)

    val_genes <- check_genes(object, genes = val_genes, max_length = 2)

  }

  # Start Create Network ----------------------------------------------------

  ##Autoencoder

  x_train <- getExpressionMatrix(object)

  input_layer <-
    keras::layer_input(shape = c(ncol(x_train)))


  encoder <-
    input_layer %>%
    keras::layer_dense(units = layers[1], activation = activation) %>%
    keras::layer_batch_normalization() %>%
    keras::layer_dropout(rate = dropout) %>%
    keras::layer_dense(units = layers[2], activation = activation) %>%
    keras::layer_dropout(rate = dropout) %>%
    keras::layer_dense(units = layers[3], activation = activation) %>%
    keras::layer_dense(units = bottleneck)

  decoder <-
    encoder %>%
    keras::layer_dense(units = layers[3], activation = activation) %>%
    keras::layer_dropout(rate = dropout) %>%
    keras::layer_dense(units = layers[2], activation = activation) %>%
    keras::layer_dropout(rate = dropout) %>%
    keras::layer_dense(units = layers[1], activation = activation) %>%
    keras::layer_dense(units = c(ncol(x_train)))

  autoencoder_model <- keras::keras_model(inputs = input_layer, outputs = decoder)

    #summary(autoencoder_model)

  autoencoder_model %>% keras::compile(
    loss = 'mean_squared_error',
    optimizer = 'adam',
    metrics = c('accuracy')
  )

  history <-
    autoencoder_model %>%
    keras::fit(x_train, x_train, epochs = epochs, shuffle = TRUE,
               validation_data = list(x_train, x_train), verbose = verbose)

  reconstructed_points <-
    autoencoder_model %>%
    keras::predict_on_batch(x = x_train)

  base::rownames(reconstructed_points) <- base::rownames(x_train)
  base::colnames(reconstructed_points) <- base::colnames(x_train)

    #head(reconstructed_points)
  if(base::isTRUE(val_plot)){

    plot_df <-
      base::rbind(
        data.frame(base::t(reconstructed_points[c("CD163", "CD68"), ]), type = "Reconstructed / Denoised"),
        data.frame(base::t(x_train[val_genes, ]), type = "Original"))

    val_plot <-
      ggplot2::ggplot(data = plot_df, ggplot2::aes(x = CD163, y = CD68, color = type)) +
      ggplot2::geom_point(alpha = 0.75) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x) +
      ggplot2::facet_wrap(. ~ type, scales = "free") +
      ggplot2::theme_classic() +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        legend.position = "none"
      ) +
      scale_color_add_on(variable = "discrete", clrp = "milo")

    plot(val_plot)

  }

  object@data@norm_exp <- reconstructed_points

  return(object)

}
