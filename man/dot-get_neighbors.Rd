% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/array_error_metrics.R
\name{.get_neighbors}
\alias{.get_neighbors}
\title{Get keys of neighboring spots}
\usage{
.get_neighbors(i, coords)
}
\arguments{
\item{i}{An \code{integer(1)} giving a row index in \code{coords}.}

\item{coords}{A \code{tibble()} containing \code{array_row}, \code{array_col}, \code{key}, and
\code{capture_area} columns.}
}
\value{
A \code{character()} of neighboring spot keys.
}
\description{
For a given row of a \code{tibble()} containing array coordinates, find the
associated spot's neighbors (belonging to the same capture area) and return
their keys.
}
\author{
Nicholas J. Eagles
}
\keyword{internal}
