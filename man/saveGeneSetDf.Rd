% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loading-and-saving.R
\name{saveGeneSetDf}
\alias{saveGeneSetDf}
\title{Save a gene set data.frame}
\usage{
saveGeneSetDf(object, output_path, file_name)
}
\arguments{
\item{object}{A valid spata-object.}

\item{output_path}{Character value. A directory leading to the folder in which
to store the data.frame.}

\item{file_name}{Character value. The filename. ( is suffixed with \emph{'.RDS'})}
}
\value{
An invisible TRUE if saved successfully or an informative error message.
}
\description{
Extracts the gene-set data.frame and saves it as a .RDS-file.
}
