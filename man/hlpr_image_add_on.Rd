% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper-functions.R
\name{hlpr_image_add_on}
\alias{hlpr_image_add_on}
\alias{hlpr_image_add_on2}
\title{Provides the image as ggplot background}
\usage{
hlpr_image_add_on(image)

hlpr_image_add_on2(object, display_image, of_sample)
}
\arguments{
\item{image}{Image input.}

\item{object}{A valid spata-object.}

\item{display_image}{Logical value.}

\item{of_sample}{The sample(s) of interest specified as a single character value or vector.
If set to \emph{""} (the default) the first sample is chosen.}
}
\value{
Either null or a ggplot2::geom_annotation_raster
}
\description{
Provides the image as ggplot background
}
