% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check-family-lazy.R
\name{check_pt}
\alias{check_pt}
\title{Check pt input}
\usage{
check_pt(
  pt_size = NULL,
  pt_alpha = NULL,
  pt_clrsp = NULL,
  pt_clrsp_dir = NULL,
  pt_clrp = NULL,
  pt_clr = NULL
)
}
\arguments{
\item{pt_size}{The size of the points specified as a single numeric value.}

\item{pt_alpha}{The transparency of the points specified as single numeric value.}

\item{pt_clrsp}{The color spectrum to be used if the specified variable displayed by
color is continuous. Run \code{all_colorspectra()} to see valid input..}

\item{pt_clrsp_dir}{The direction of the color spectrum specified as either \emph{1}
or \emph{-1}.}

\item{pt_clrp}{The color panel to be used if the specified variable displayed by
color is categorical/discrete. Run \code{all_colorpanels()} to see valid input.}

\item{pt_clr}{The base color of every point displayed in the plot.}
}
\value{
A logical value TRUE if nothing invalid has been detected or an informative
error message.
}
\description{
A member of the \code{lazy-check_*()}-family.
}
\details{
Members of the \code{lazy-check_*()}-family take the arguments
of their kind - that are used in the function they are called in - and
checks whether these arguments input fit the requirements. They stop and return an
error immediately once they stumble upon something invalid. They do not alter or adjust input
and return TRUE if the whole function has been executed without anything
invalid being found.
}
