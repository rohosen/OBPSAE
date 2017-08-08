#' Hospital data
#'
#' Morris and Christiansen (1995) presented a dataset involving 23 hospitals (out of a total of 219 hospitals) that had at least 50 kidney transplants during a 27-month period.
#'
#' @format A data frame with 23 rows and 4 variables:
#' \describe{
#'   \item{Sample Size}{Number of kidney transplants in each hospital during a 27 months period (at least 50)}
#'   \item{Graft Failure Rate}{Graft failure rates for kidney transplant operations i.e. number of graft failures/sample size}
#'   \item{Severity Index}{Severity index at each hospital which is measured by the average fraction of females, blacks, children, and extremely ill kidney recipients at that hospital}
#'   \item{SE}{The standard error for the graft failure rate. The variance for the graft failure rate \eqn{D_{i}}, is approximated by \eqn{(0.2)(0.8)/n_{i}}, where 0.2 is the observed failure rate for all hospitals. Thus, \eqn{D_{i}} is assumed known.}
#' }
#'
"Hospital"

#' Income data
#'
#' This dataset is regarding the median income of families with three, four and five people for the fifty states of the U.S. and the District of Columbia.
#'
#' @format A data frame with 561 rows and 9 variables:
#' \describe{
#'   \item{State code}{State code of each of the 50 states and District of Columbia}
#'   \item{Year}{The year the data correspond to, generally ranging over 1969, 1979 to 1988 except 1980.}
#'   \item{BEA}{Per-capita income estimates produced by the Bureau of Economic Analysis (BEA).}
#'   \item{CPS/Census estimate-3persons}{Estimated median income of 3-persons families from Current Population Survey (CPS) OR Census of 1969 and 1979.}
#'   \item{CPS SE estimate 3-person}{Standard Errors of the Estimated median income of 3-persons family from Current Population Survey (CPS) OR Census, for Census the SE is set to zero.}
#'   \item{CPS/Census estimate 4persons}{Estimated median income of 3-persons families from Current Population Survey (CPS) OR Census of 1969 and 1979.}
#'   \item{CPS SE estimate 4-person}{Standard Errors of the Estimated median income of 3-persons family from Current Population Survey (CPS) OR Census, for Census the SE is set to zero.}
#'   \item{CPS/Census estimate-5persons}{Estimated median income of 3-persons families from Current Population Survey (CPS) OR Census of 1969 and 1979.}
#'   \item{CPS SE estimate 5-person}{Standard Errors of the Estimated median income of 3-persons family from Current Population Survey (CPS) OR Census, for Census the SE is set to zero.}
#' }
"Income"
