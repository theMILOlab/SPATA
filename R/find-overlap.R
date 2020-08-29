#' @title Spatial overlap
#'
#' @description --
#'
#' @inherit check_sample params
#' @inherit check_method params
#' @inherit check_pt params
#' @inherit verbose params
#' @param input Character vector of length 2. Specifies the genes or gene-sets
#' whose overlap you want to calculate.
#'
#' @return A named list containing the results (plots, estimated correlation value, data.frames).
#'
#' @export
#'

findOverlap <- function(object,
                        of_sample = "",
                        input,
                        method_gs = "zscore",
                        method_ovl = "bayesian",
                        pt_clrp = "milo",
                        pt_clrsp = "inferno",
                        pt_size = 2,
                        pt_alpha = 0.9,
                        verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_method(method_ovl = method_ovl, method_gs = method_gs)
  check_pt(pt_size = pt_size, pt_alpha = pt_alpha,
           pt_clrp = pt_clrp, pt_clrsp = pt_clrsp)

  # adjusting check
  of_sample <- check_sample(object, of_sample, 1)

  variables <-
    check_variables(
      variables = input,
      all_genes = getGenes(object, in_sample = of_sample),
      all_gene_sets = getGeneSets(object),
      max_length = 2,
      max_slots = 1)

  # -----


  # 2. Data extraction ------------------------------------------------------

  df <-
    joinWithVariables(
      object = object,
      coords_df = getCoordinates(object, of_sample),
      method_gs = method_gs,
      variables = variables,
      verbose = verbose,
      normalize = TRUE) %>%
    base::as.data.frame()

  # -----

  if(method_ovl == "bayesian"){

    # C - script
    model_string <- "
    model {
    for(i in 1:n) {
    x[i,1:2] ~ dmnorm(mu[], prec[ , ])
    }

    # Constructing the covariance matrix and the corresponding precision matrix.
    prec[1:2,1:2] <- inverse(cov[,])
    cov[1,1] <- sigma[1] * sigma[1]
    cov[1,2] <- sigma[1] * sigma[2] * rho
    cov[2,1] <- sigma[1] * sigma[2] * rho
    cov[2,2] <- sigma[2] * sigma[2]

    # Uninformative priors on all parameters which could, of course, be made more informative.
    sigma[1] ~ dunif(0, 1000)
    sigma[2] ~ dunif(0, 1000)
    rho ~ dunif(-1, 1)
    mu[1] ~ dnorm(0, 0.001)
    mu[2] ~ dnorm(0, 0.001)

    # Generate random draws from the estimated bivariate normal distribution
    x_rand ~ dmnorm(mu[], prec[ , ])
    }
    "

    if(base::isTRUE(verbose)){

      base::message("Starting bayesian correlation. This might take a while.")

    }

    # select value-variables
    df_vls <-
      dplyr::select(df, -barcodes, -sample, -x, -y) %>%
      base::as.data.frame()

    data_list <- list(x = df_vls, n = base::nrow(df_vls))
    inits_list <- list(mu = c(base::mean(df_vls[, 1]), base::mean(df_vls[, 2])),
                      rho = stats::cor(df_vls[, 1], df_vls[, 2]),
                      sigma = c(stats::sd(df_vls[, 1]), stats::sd(df_vls[, 1])))

    if(base::isTRUE(verbose)){

      base::message("Fitting the model.")

    }

    jags_model <-
      rjags::jags.model(base::textConnection(model_string),
                        data = data_list,
                        inits = inits_list,
                        n.adapt = 500, n.chains = 3, quiet = base::isFALSE(verbose))

    if(base::isTRUE(verbose)){

      base::message("Fitting 'Markov Chain Monte Carlo'.")

    }

    mcmc_samples <-
      rjags::coda.samples(model = jags_model,
                          variable.names = c("mu", "rho", "sigma", "x_rand"),
                          n.iter = 1000)

    if(base::isTRUE(verbose)){

      base::message("Creating output for method 'bayesian'.")

    }

    samples_mat <- as.matrix(mcmc_samples)

    estimated_correlation <-
      base::mean(samples_mat[, "rho"]) %>%
      base::round(digits = 2)

    df_plot_estimated <-
      base::as.data.frame(samples_mat[, c("x_rand[1]", "x_rand[2]")]) %>%
      magrittr::set_names(value = input) %>%
      dplyr::mutate(type = "Estimated")

    df_plot_realdata <-
      dplyr::mutate(.data = df_vls, type = "Real")

    df_bayesian <- base::rbind(df_plot_estimated, df_plot_realdata)

    # model plot
    x <- df_plot_estimated[,1]
    y <- df_plot_estimated[,2]

    model <- stats::loess(formula = y ~ x, span = 0.75)

    loess <- data.frame(x = stats::predict(model, base::seq(min(x), max(x),length.out = 3000)),
                        y = base::seq(min(x), max(x), length.out = 3000),
                        col = "Loess Fit")

    plot_bayesian <-
      ggplot2::ggplot() +
      ggplot2::geom_vline(xintercept = 0, linetype = "dotted") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
      ggplot2::geom_point(data = df_bayesian, alpha = 0.25,
                          ggplot2::aes(x = .data[[input[1]]], y = .data[[input[2]]], color = type)) +
      ggplot2::geom_line(data = loess, alpha = 0.75, size = 0.5,
                         ggplot2::aes(x = x, y = y, color = col)) +
      confuns::scale_color_add_on(variable = "discrete", clrp = pt_clrp) +
      ggplot2::theme_classic() +
      ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2)) +
      ggplot2::labs(title = base::paste0("Bayesian Spatial Co-Expression mean RHO: ",
                                         base::round(estimated_correlation, digits = 2)),
                    x = base::names(df_bayesian)[1],
                    y = base::names(df_bayesian)[2],
                    color = "Bayesian Prediction")  +
      ggplot2::guides(size = FALSE,
                      alpha = FALSE,
                      color = ggplot2::guide_legend(override.aes = list(alpha = 0.5))
                      )

    df_spatial <-
      dplyr::mutate(.data = df,
                    sum = base::rowSums(as.matrix(df_vls[,input])))

    plot_spatial <-
      ggplot2::ggplot(data = df_spatial,
                      mapping = ggplot2::aes(x = x, y = y, color = sum)) +
      ggplot2::geom_point(size = pt_size, alpha = pt_alpha) +
      confuns::scale_color_add_on(clrsp = pt_clrsp) +
      ggplot2::labs(title = base::paste0("Bayesian Spatial Co-Expression mean RHO: ",
                                         estimated_correlation),
                    x = NULL, y = NULL, color = "Bayesian Prediction") +
      cowplot::theme_map()


    ##Output data file
    output_list <- list(
      data = df,
      est_corr = estimated_correlation,
      plot_bayesian = plot_bayesian,
      df_bayesian = df_bayesian,
      fit_bayesian = loess,
      df_spatial = df_spatial,
      plot_spatial = plot_spatial
    )

    if(base::isTRUE(verbose)){base::message("Done.")}

    base::return(output_list)

  } else if(method_ovl == "classic"){

    if(base::isTRUE(verbose)){

      base::message("Creating output for method 'classic'.")

    }

    estimated_correlation <-
      stats::cor(df[, input[1]], df[, input[2]]) %>%
      base::as.numeric() %>%
      base::round(digits = 2)

    df_plot <- dplyr::mutate(.data = df, type = "Real")

    # fit a model
    x <- df_plot[,input[1]]
    y <- df_plot[,input[2]]

    model <- stats::glm(formula = y ~ x)
    glm <- data.frame(x = as.numeric(predict(model, newdata = data.frame(y = x))),
                      y = x, type = "GLM Fit")

    plot_pearson <-
    ggplot2::ggplot(data = df_plot) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dotted") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
      ggplot2::geom_point(alpha = 0.6,
                          mapping = ggplot2::aes(x = !!rlang::sym(input[1]),
                                                 y = !!rlang::sym(input[2]),
                                                 color = type)) +
      ggplot2::geom_line(data = glm, ggplot2::aes(x = x,
                                                    y = y,
                                                    color = type)) +
      confuns::scale_color_add_on(variable = "discrete", clrp = pt_clrp) +
      ggplot2::theme_classic() +
      ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.2)) +
      ggplot2::labs(title = stringr::str_c("Pearson Spatial Co-Expression mean RHO: ", estimated_correlation),
                    color = "Pearson Prediction")

    df_spatial <-
      dplyr::mutate(.data = df,
                    sum = base::rowSums(as.matrix(df[,input])))

    plot_spatial <-
      ggplot2::ggplot(data = df_spatial,
                      mapping = ggplot2::aes(x = x, y = y, color = sum)) +
      ggplot2::geom_point(size = pt_size, alpha = pt_alpha) +
      confuns::scale_color_add_on(clrsp = pt_clrsp) +
      ggplot2::labs(title = base::paste0("Pearson Spatial Co-Expression mean RHO: ",
                                         base::round(estimated_correlation, digits = 2)),
                    x = NULL, y = NULL, color = "Pearson Prediction") +
      cowplot::theme_map()

    # output list
    output_list <- list(
      data = df,
      est_corr = estimated_correlation,
      plot_pearson = plot_pearson,
      fit_glm = glm,
      plot_spatial = plot_spatial,
      df_spatial = df_spatial
    )

    if(base::isTRUE(verbose)){base::message("Done.")}

  }

  base::return(output_list)

}
