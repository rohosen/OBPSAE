% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MSPE.R
\name{MSPE_naive}
\alias{MSPE_naive}
\title{Naive MSPE Estimator of Observed Best Predictor (OBP).}
\usage{
MSPE_naive(y, D, theta.OBP, A.BPE)
}
\arguments{
\item{y}{A vector of the response variable.}

\item{D}{Scalar, Vector (or diagonal matrix) of random error variances.}

\item{theta.OBP}{A vector of OBP values}

\item{A.BPE}{BPE estimate of variance of random effects or the true value, if known.}
}
\value{
This function will return the Naive MSPE estimator of the OBP.
}
\description{
This function computes the naive MSPE estimator of OBP.
}
\references{
Jiang J, Nguyen T, and Rao J. S. (2011), "Best Predictive Small Area Estimation", Journal of the American Statistical Association.
}
