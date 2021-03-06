% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-family.R
\name{getTrajectoryObject}
\alias{getTrajectoryObject}
\title{Obtain trajectory object}
\usage{
getTrajectoryObject(object, trajectory_name, of_sample = "")
}
\arguments{
\item{object}{A valid spata-object.}

\item{trajectory_name}{The trajectory of interest specified
as a single character value.}

\item{of_sample}{The sample(s) of interest specified as a single character value or vector.
If set to \emph{""} (the default) the first sample is chosen.}
}
\value{
An object of class \code{spatialTrajectory}.
}
\description{
Obtain trajectory object
}
