% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trajectory-analysis.R
\name{assessTrajectoryTrends}
\alias{assessTrajectoryTrends}
\alias{assessTrajectoryTrends2}
\title{Trajectory trend analysis}
\usage{
assessTrajectoryTrends(
  object,
  trajectory_name,
  of_sample = "",
  variables,
  accuracy = 5,
  verbose = TRUE
)

assessTrajectoryTrends2(stdf, verbose = TRUE)
}
\arguments{
\item{object}{A valid spata-object.}

\item{trajectory_name}{The trajectory of interest specified
as a single character value.}

\item{of_sample}{The sample(s) of interest specified as a single character value or vector.
If set to \emph{""} (the default) the first sample is chosen.}

\item{variables}{Character vector. The variables of interest:

\itemize{
 \item{ \strong{Gene sets}: Must be in \code{getGeneSets()}}
 \item{ \strong{Genes}: Must be in \code{getGenes()}}
 \item{ \strong{Features}: Must be numeric ones of \code{getFeatureNames()}}
 }}

\item{verbose}{Logical. If set to TRUE informative messages regarding
the computational progress will be printed.

(Warning messages will always be printed.)}

\item{stdf}{A summarized trajectory data.frame. (e.g. obtained by
\code{getTrajectoryDf()}).}
}
\value{
A data.frame arranged by the residuals area-under-the-curve-values describing
how well a model fitted the expression trend of a gene or gene set.
}
\description{
Analyzes the trend of gene and gene-set-expressions along
trajectories by fitting a variety of mathematical models to them and
by assessing the quality of the fit.

\itemize{
 \item{\code{assessTrajectoryTrends()}: Takes a valid spata-object and constructs
 the subsequent data.frame from scratch.}
 \item{\code{assessTrajectoryTrends2()}: Takes a summarized trajectory data.frame
 returned by \code{getTrajectoryDf()}.}
 }
}
