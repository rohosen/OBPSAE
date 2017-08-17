
############### 7. Jackknife function #################
####### A function of data (y,x) ###################
mcjack_bench <- function(y, x, j, D, weight, K){

  y.notj = y[-j]
  x.notj = x[-j]
  D.notj = D[-j]

  BPE.notj = bpeFH(y.notj~x.notj, errorvar = D.notj)
  beta.notj = BPE.notj$beta
  A.notj = BPE.notj$A


  ###### We assume that the above fitted spline model is the true underlying model for the data #########
  ###### We need K simulations from the spline model ######
  ###### Each simulation is a sample of size m ####

  ###### Fitting a smoothing spline model to our data #####


  fhat.notj <- function(t, p, q, knots, beta, gamma){
    Tt = matrix(1,length(t),1)
    if(p > 0){
      for(i in 1:p){
        Tt = cbind(Tt,t^i)
      }
    }

    for(i in 1:q){
      Tt = cbind(Tt, sapply(t-knots[i], function(u)max(u,0))^p)
    }

    para <- c(beta, gamma)
    t.pred <- Tt %*% para
    return(t.pred)

  }



  fit.spline.notj = fitSpline.cv(x=x.notj,y=y.notj)
  p = fit.spline.notj$p
  q = fit.spline.notj$q
  knots = fit.spline.notj$knots
  beta = fit.spline.notj$beta
  gamma = fit.spline.notj$gamma

  fittedSpline.notj <- fhat.notj(x, p, q, knots, beta, gamma)

  sim.data <- lapply( 1:K,  FUN = spline.sim, f.hat = fittedSpline.notj, A.hat = A.notj, Di = D)
  vec1 = NULL
  vec2 = NULL
  for(i in 1:K){
    vec1 = c(vec1,sim.data[[i]]$Y)
    vec2 = c(vec2,sim.data[[i]]$theta)
  }
  sim.matrix.Y = matrix(vec1,ncol=K)
  sim.matrix.theta = matrix(vec2,ncol=K)


  ##### Calculate the OBP for each of the above K simulated sample from the fitted spline model ####

  ## Adjusted MSPE Estimator
  obpAdj.wrapper <- function(y_wrapper){
    for (ii in ls(parent.frame())){
      assign(ii,get(ii,parent.frame()))
    }
    return(obpFH_adjusted(y_wrapper~x, errorvar = D, weight = weight)$obpAdjusted)
  }
  OBP.adj.sim <- apply(sim.matrix.Y, 2, FUN = obpAdj.wrapper)
  b.phi.adj <- MSPE.calc(OBP.adj.sim, sim.matrix.theta)

  ## Augmented MSPE Estimator
  obpAug.wrapper <- function(y_wrapper){
    for (ii in ls(parent.frame())){
      assign(ii,get(ii,parent.frame()))
    }
    return(obpFH_augmented(y_wrapper~x, errorvar = D, weight = weight)$obpAugmented)
  }
  OBP.aug.sim <- apply(sim.matrix.Y, 2, FUN = obpAug.wrapper)
  b.phi.aug <- MSPE.calc(OBP.aug.sim, sim.matrix.theta)


  return(list(b.phi.adj = b.phi.adj, b.phi.aug = b.phi.aug))

}





#' McSpline MSPE estimator of benchmarked OBP
#'
#' This function computes the McSpline MSPE estimator of the benchmarked Observed Best Predictor.
#' @param formula an object of class formula (or one that can be coerced to that class): a symbolic description of the model to be fitted. The variables included in formula must have a length equal to the number of small areas. More about the model specification are given under Details.
#' @param data optional data frame containing the variable names in \code{formula}.
#' @param errorvar vector containing the variances of the random errors for all small areas.
#' @param weight vector containing the sampling weights of small areas. If sum of the weights is not 1, the weights are normalized.
#' @param A.BPE optional BPE estimate of random effects variance.
#' @param K number of Monte Carlo simulations. Default is 1000.
#' @details
#' \code{formula} is specified in the form \code{response ~ predictor} where the predictor is univariate. \code{formula} has an implied intercept term. To remove the intercept term, use either \code{y ~ x - 1} or \code{y ~ 0 + x}.\cr
#' \cr If \code{A.BPE} is missing, the function computes its BPE from data. User can also provide the true A instead if that is known.
#' @return The function will return a list with the following objects.
#' \item{McSpline_Adj_bench}{ McSpline MSPE estimator of the adjusted OBP.}
#' \item{McSpline_Aug_bench}{ McSpline MSPE estimator of the augmented OBP.}
#' @references Bandyopadhyay R, Jiang J (2017) "Benchmarking the Observed Best Predictor"
#' @references Jiang J, Nguyen T, and Rao J. S. (2011), "Best Predictive Small Area Estimation", Journal of the American Statistical Association.
#' @export







MSPE_McSpline_Benchmark <- function(formula, data, errorvar, weight, A.BPE, K=1000){
  if (!missing(data)) {
    formuladata <- model.frame(formula, na.action = na.omit, data)
    X <- model.matrix(formula, data)
    #vardir <- data[, namevar]
  }
  else {
    formuladata <- model.frame(formula, na.action = na.omit)
    X <- model.matrix(formula)
  }
  if (missing(A.BPE)){
    A.BPE = obpFH(formula,data,errorvar)$A.BPE
  }
  y <- formuladata[,1]
  if (attr(attributes(formuladata)$terms, "response") == 1) textformula <- paste(formula[2], formula[1], formula[3])else textformula <- paste(formula[1], formula[2])
  if (length(na.action(formuladata)) > 0) stop("Argument formula=", textformula, " contains NA values.")
  if (any(is.na(errorvar))) stop("Argument errorvar=", namevar, " contains NA values.")
  if (dim(X)[2] > 2) stop("McSpline takes only one x.")

  x = as.vector(X[,2])

  fit.spline = fitSpline.cv(x=x,y=y)$f.hat

  if(is.matrix(errorvar)) {
    D = diag(errorvar)
  }else if(length(errorvar)==1) {
    D = c(rep(errorvar,length(y)))
  }else {
    D = errorvar
  }

  sim.data <- lapply( 1:K,  FUN = spline.sim, f.hat = fit.spline, A.hat = A.BPE, Di = D)

  vec1 = NULL
  vec2 = NULL
  for(i in 1:K){
    vec1 = c(vec1,sim.data[[i]]$Y)
    vec2 = c(vec2,sim.data[[i]]$theta)
  }
  sim.matrix.Y = matrix(vec1,ncol=K)
  sim.matrix.theta = matrix(vec2,ncol=K)

  ##### Calculate the OBP for each of the above K simulated sample from the fitted spline model ####


  ## Adjusted MSPE Estimator
  obpAdj.wrapper <- function(y_wrapper){
    for (ii in ls(parent.frame())){
      assign(ii,get(ii,parent.frame()))
    }
    return(obpFH_adjusted(y_wrapper~x, errorvar = D, weight = weight)$obpAdjusted)
  }
  OBP.adj.sim <- apply(sim.matrix.Y, 2, FUN = obpAdj.wrapper)
  b.phi.adj <- MSPE.calc(OBP.adj.sim, sim.matrix.theta)

  ## Augmented MSPE Estimator
  obpAug.wrapper <- function(y_wrapper){
    for (ii in ls(parent.frame())){
      assign(ii,get(ii,parent.frame()))
    }
    return(obpFH_augmented(y_wrapper~x, errorvar = D, weight = weight)$obpAugmented)
  }
  OBP.aug.sim <- apply(sim.matrix.Y, 2, FUN = obpAug.wrapper)
  b.phi.aug <- MSPE.calc(OBP.aug.sim, sim.matrix.theta)



  return(list(McSpline_Adj_bench = b.phi.adj, McSpline_Aug_bench = b.phi.aug))
}





#' McJack MSPE estimator of benchmarked OBP
#'
#' This function computes the McJack MSPE benchmark estimator.
#' @param formula an object of class formula (or one that can be coerced to that class): a symbolic description of the model to be fitted. The variables included in formula must have a length equal to the number of small areas. More about the model specification are given under Details.
#' @param data optional data frame containing the variable names in \code{formula}.
#' @param errorvar vector containing the variances of the random errors for all small areas.
#' @param weight vector containing the sampling weights of small areas. If sum of the weights is not 1, the weights are normalized.
#' @param A.BPE optional BPE estimate of random effects variance.
#' @param K number of Monte Carlo simulations. Default is 1000.
#' @param returnMcSpline logical. Returns McSpline MSPE benchmark estimators with McJack (default).
#' @details
#' \code{formula} is specified in the form \code{response ~ predictor} where the predictor is univariate. \code{formula} has an implied intercept term. To remove the intercept term, use either \code{y ~ x - 1} or \code{y ~ 0 + x}.\cr
#' \cr If \code{A.BPE} is missing, the function computes the BPE from data. User can also provide the true A instead if that is known.
#' @return The function will return a list with the following objects.
#' \item{McJack_Adj_bench}{McJack MSPE estimator of the adjusted OBP.}
#' \item{McJack_Aug_bench}{McJack MSPE estimator of the augmented OBP.}
#' \item{McSpline_Adj_bench}{McSpline MSPE estimator of the adjusted OBP. This is returned by default. To turn this feature off, set \code{returnMcSpline = FALSE}}
#' \item{McSpline_Aug_bench}{McSpline MSPE estimator of the augmented OBP. This is returned by default. To turn this feature off, set \code{returnMcSpline = FALSE}}
#' @references Bandyopadhyay R, Jiang J (2017) "Benchmarking the Observed Best Predictor"
#' @references Jiang J, Nguyen T, and Rao J. S. (2011), "Best Predictive Small Area Estimation", Journal of the American Statistical Association.
#' @export



MSPE_McJack_Benchmark <- function(formula, data, errorvar, weight, A.BPE, K=1000, returnMcSpline=TRUE){
  if (!missing(data)) {
    formuladata <- model.frame(formula, na.action = na.omit, data)
    X <- model.matrix(formula, data)
    #vardir <- data[, namevar]
  }
  else {
    formuladata <- model.frame(formula, na.action = na.omit)
    X <- model.matrix(formula)
  }
  if (missing(A.BPE)){
    A.BPE = obpFH(formula,data,errorvar)$A.BPE
  }
  y <- formuladata[,1]
  if (attr(attributes(formuladata)$terms, "response") == 1) textformula <- paste(formula[2], formula[1], formula[3])else textformula <- paste(formula[1], formula[2])
  if (length(na.action(formuladata)) > 0) stop("Argument formula=", textformula, " contains NA values.")
  if (any(is.na(errorvar))) stop("Argument errorvar=", namevar, " contains NA values.")
  if (dim(X)[2] > 2) stop("McSpline takes only one x.")

  x = as.vector(X[,2])
  m = length(y)

  McSpline_bench = MSPE_McSpline_Benchmark(y~x, errorvar=errorvar, weight = weight, A.BPE=A.BPE, K=K)
  b.phi.adj = McSpline_bench$McSpline_Adj_bench
  b.phi.aug = McSpline_bench$McSpline_Aug_bench

  b.phi.adj.notj = matrix(NA,m,m)
  b.phi.aug.notj = matrix(NA,m,m)

  if(is.matrix(errorvar)) {
    D = diag(errorvar)
  }else if(length(errorvar)==1) {
    D = c(rep(errorvar,length(y)))
  }else {
    D = errorvar
  }

  for(j in 1:m){
    BenchMcJack = mcjack_bench(y, x, j, D, weight, K)
    b.phi.adj.notj[,j] = BenchMcJack$b.phi.adj
    b.phi.aug.notj[,j] = BenchMcJack$b.phi.aug
  }


  #### Calculate the Jacknife MSPE estimator

  #### Adjusted OBP
  b.phi.adj.mat = matrix(b.phi.adj,m,1)%*%matrix(1,nrow=1,ncol=m)
  Adj.factor.mat1 = b.phi.adj.notj - b.phi.adj.mat
  Adj.factor1 = rowSums(Adj.factor.mat1)
  b.phi.adj.jackknife = b.phi.adj - ((m-1)/m)*Adj.factor1


  #### Augmented OBP
  b.phi.aug.mat = matrix(b.phi.aug,m,1)%*%matrix(1,nrow=1,ncol=m)
  Adj.factor.mat2 = b.phi.aug.notj - b.phi.aug.mat
  Adj.factor2 = rowSums(Adj.factor.mat2)
  b.phi.aug.jackknife = b.phi.aug - ((m-1)/m)*Adj.factor2

  if(returnMcSpline)return(list(McJack_Adj_bench = b.phi.adj.jackknife, McJack_Aug_bench = b.phi.aug.jackknife, McSpline_Adj_bench = b.phi.adj, McSpline_Aug_bench = b.phi.aug))
  else return(list(McJack_Adj_bench = b.phi.adj.jackknife, McJack_Aug_bench = b.phi.aug.jackknife))


}
