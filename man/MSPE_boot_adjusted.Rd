% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MSPE.R
\name{MSPE_boot_adjusted}
\alias{MSPE_boot_adjusted}
\title{Bootstrap MSPE estimator of Adjusted Observed Best Predictor (Adjusted OBP).}
\usage{
MSPE_boot_adjusted(theta.OBP.adjusted, D, x, weight, L = 200)
}
\arguments{
\item{theta.OBP.adjusted}{A vector of adjusted OBPs.}

\item{D}{Scalar, Vector (or diagonal matrix) of random error variances.}

\item{x}{Design matrix without intercept term.}

\item{weight}{Weight of each small area for calculating the benchmarked OBP (generally sampling weights of each small area).}

\item{L}{is the number of bootstrap samples. the default is 200.}
}
\value{
This function will return a vector of the Bootstrap MSPE estimator of Adjusted OBP.
}
\description{
This function computes the bootstrap MSPE estimator of benchmarked adjusted OBP.
}
\references{
Bandyopadhyay R, Jiang J (2017) "Benchmarking the Observed Best Predictor"

Jiang J, Nguyen T, and Rao J. S. (2011), "Best Predictive Small Area Estimation", Journal of the American Statistical Association.
}
