% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pseudotime-embedding.R
\name{compileCellDataSet}
\alias{compileCellDataSet}
\title{Initiate monocle3-pseudotime analysis}
\usage{
compileCellDataSet(
  object,
  use_cds_file = FALSE,
  save_cds_file = FALSE,
  preprocess_method = "PCA",
  cluster_method = "leiden",
  verbose = TRUE
)
}
\arguments{
\item{object}{A valid spata-object.}

\item{use_cds_file}{A file-directory leading to a .rds file containing a valid
cell_data_set-object previously calculated for the specified object. Specified
as a character value. If set to FALSE the cell_data_set object will be created
from scratch.}

\item{save_cds_file}{A file-directory (that does not already exists) under which the used or created cell_data_set-object
is going to be stored specified as a character value. Should end with .rds.}

\item{preprocess_method}{Given to \code{monocle3::preprocess_cds()} if \code{use_cds_file} isn't a character string.}

\item{cluster_method}{Given to \code{monocle3::cluster_cells()} if \code{use_cds_file} isn't a character string. Must be one of
\emph{'leiden', 'louvain'}.}

\item{inherit}{verbose

(Warning messages will always be printed.)}
}
\value{
A monocle3::cell_data_set object.
}
\description{
Loads or compiles a valid cell_data_set-object
that fits to the provided spata-object.
}