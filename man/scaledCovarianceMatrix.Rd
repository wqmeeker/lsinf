% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scaledCovarianceMatrix.R
\name{scaledCovarianceMatrix}
\alias{scaledCovarianceMatrix}
\title{Returns elements of the scaled variance-covariance martrix for a location-scale distribution.}
\usage{
scaledCovarianceMatrix(distribution, std.log.censor.time, unlist = TRUE)
}
\arguments{
\item{distribution}{The probability distribution.}

\item{std.log.censor.time}{Vector of standardized log censoring times.}

\item{unlist}{If TRUE and std.log.censor.time is a scalar, return a vector instead of a list.}
}
\value{
Silently returns  a results matrix containing the variance factors.
}
\description{
Returns elements of the scaled variance-covariance martrix for a location-scale distribution.
}
