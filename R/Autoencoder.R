#' @title Deep Autoencoder for Denoising
#'
#' @description This Function will return a SPATA object containing denoised data and re-analysis
#'
#' @param SPATAobj
#' @param Layer Set Number of Neurons in your 3 hidden layers (defauled = c("128", "64", "32"))
#' @param Bottleneck Set numbers of Bottleneck Neurons (defauled = "32")
#' @param Dropout Define Dropout (defauled = "0.1")
#' @param activation Set activation function (defauled = "relu")
#' @param epochs epochs of your NN (defauled = 20)
#' @param val_plot
#'
#'
#'
#' @return A spata-object.
#'
#' @importFrom Seurat ScaleData
#'
#' @export
#'
#'


AutoEncoderDenoising <-
  function(object,
           Layer=NULL,
           Bottleneck=NULL,
           Dropout=NULL,
           activation=NULL,
           epochs=NULL,
           val_plot=F){


# Start Create Network ----------------------------------------------------

#Hyper Parameter

if(is.null(Layer)){Layer=c(128,64,32)}
if(is.null(Bottleneck)){Bottleneck=16}
if(is.null(Dropout)){Dropout=0.1}
if(is.null(activation)){activation="relu"}
if(is.null(epochs)){epochs=20}
require(tensorflow)
require(keras)
require(reticulate)
use_implementation("tensorflow")
#keras:::implementation()
tensorflow::tf_config()
#tensorflow::install_tensorflow()
#tensorflow::tfe_enable_eager_execution(device_policy = "silent")


##Autoencoder

x_train <- getExpressionMatrix(object)
x_train[1:10,1:10]
input_layer <-
  keras::layer_input(shape = c(ncol(x_train)))


encoder <-
  input_layer %>%
  keras::layer_dense(units = Layer[1], activation = activation) %>%
  keras::layer_batch_normalization() %>%
  keras::layer_dropout(rate = Dropout) %>%
  keras::layer_dense(units = Layer[2], activation = activation) %>%
  keras::layer_dropout(rate = Dropout) %>%
  keras::layer_dense(units = Layer[3], activation = activation) %>%
  keras::layer_dense(units = Bottleneck)

decoder <-
  encoder %>%
  keras::layer_dense(units = Layer[3], activation = activation) %>%
  keras::layer_dropout(rate = Dropout) %>%
  keras::layer_dense(units = Layer[2], activation = activation) %>%
  keras::layer_dropout(rate = Dropout) %>%
  keras::layer_dense(units = Layer[1], activation = activation) %>%
  keras::layer_dense(units = c(ncol(x_train)))

autoencoder_model <- keras::keras_model(inputs = input_layer, outputs = decoder)
#summary(autoencoder_model)
autoencoder_model %>% keras::compile(
  loss='mean_squared_error',
  optimizer='adam',
  metrics = c('accuracy')
)
history <-
  autoencoder_model %>%
  keras::fit(x_train,x_train,epochs=epochs,shuffle=TRUE, validation_data= list(x_train, x_train), verbose=1
  )

reconstructed_points <-
  autoencoder_model %>%
  keras::predict_on_batch(x = x_train)
rownames(reconstructed_points) <- rownames(x_train) ; colnames(reconstructed_points) <- colnames(x_train)

#head(reconstructed_points)
if(val_plot==T){
  require(ggplot2)
  plotdf <- rbind(data.frame(t(reconstructed_points[c("CD163", "CD68"), ]), type="reco"),
                  data.frame(t(x_train[c("CD163", "CD68"), ]), type="orig"))
  p1=ggplot2::ggplot(plotdf %>%  dplyr::filter(type=="reco"), aes(x=CD163, y=CD68, col=type))+ggplot2::geom_point()+ geom_smooth(method = "lm")+ggplot2::theme_classic()
  p2=ggplot2::ggplot(plotdf %>%  dplyr::filter(type=="orig"), aes(x=CD163, y=CD68, col=type))+ggplot2::geom_point()+ geom_smooth(method = "lm")+ggplot2::theme_classic()
  p1+p2
}

object@data@norm_exp <- reconstructed_points

return(object)

}




