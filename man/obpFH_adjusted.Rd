% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/obpFHbenchmark.R
\name{obpFH_adjusted}
\alias{obpFH_adjusted}
\title{Adjusted observed best predictor for Fay-Herriot Model.}
\usage{
obpFH_adjusted(formula, data, errorvar, weight, randvar = NULL,
  maxiter = 100, precision = 1e-04)
}
\arguments{
\item{formula}{an object of class formula (or one that can be coerced to that class): a symbolic description of the model to be fitted. The variables included in formula must have a length equal to the number of small areas. More about the model specification are given under Details.}

\item{data}{data frame containing the variable names in formula and errorvar.}

\item{errorvar}{vector containing the variances of the random errors for all small areas.}

\item{weight}{vector containing the sampling weights of small areas. If sum of the weights is not 1, the weights are normalized.}

\item{randvar}{variance of the random effect. If not supplied, BPE is estimated.}

\item{maxiter}{maximum number of iterations used in estimating randvar.}

\item{precision}{covergence tolerance limit for estimating randvar.}
}
\value{
The function will return a list containing the Adjusted OBP as follows:
\item{obpAdjusted}{a vector of adjusted OBP values.}
}
\description{
This function computes the Adjusted Observed Best Predictor (OBP) for Fay-Herriot model. The variance of the random error can be specified by the user. Otherwise the function will calculate its Best Predictive Estimator (BPE). In the process of of computing Adjusted OBP it also calculates the BPE of the regression coefficients of the fixed effect.
}
\details{
The variance of the random effect can be specified by the user. Otherwise the function will calculate its Best Predictive Estimator (BPE). In the process of of computing Adjusted OBP it also calculates the BPE of the regression coefficients of the fixed effect.\cr
\cr \code{formula} is specified in the form \code{response ~ predictors} where the predictors are separated by \code{+}. \code{formula} has an implied intercept term. To remove the intercept term, use either \code{y ~ x - 1} or \code{y ~ 0 + x}.\cr
}
\references{
Bandyopadhyay R, Jiang J (2017) "Benchmarking the Observed Best Predictor"

Jiang J, Nguyen T, and Rao J. S. (2011), "Best Predictive Small Area Estimation", Journal of the American Statistical Association.
}
