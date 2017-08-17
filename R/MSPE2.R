
#' @import mvtnorm
spline.sim <- function(i, f.hat, A.hat, Di){
  m = length(f.hat)
  ### Error generated from N(0, D)
  e = c(mvtnorm::rmvnorm(1, rep(0,m), diag(Di)))
  ### Random effect generated from N(0, A)
  v = rnorm(m, 0, sqrt(A.hat))
  ### True small area mean
  theta = f.hat + v ### smspline.y is the fitted y from the smoothing spline
  ### generate data
  Y = theta + e
  return(list(Y=Y, theta= theta))
}

##### 4. MSPE Function for calculating MSPE from the original and estimated matrices of theta
MSPE.calc<-function(theta.OBP.mat,theta.mat){
  m = nrow(theta.OBP.mat)
  MSPE<-rep(0,m)
  sq.error.mat = (theta.OBP.mat - theta.mat)^2
  MSPE = rowMeans(sq.error.mat)
  return(MSPE)
}

############### 7. Jackknife function #################
####### A function of data (y,x) ###################
mcjack <- function(y, x, j, D, K){

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

  obpFH.wrapper <- function(y_wrapper){
    for (ii in ls(parent.frame())){
      assign(ii,get(ii,parent.frame()))
    }
    return(obpFH(y_wrapper~x, errorvar = D)$theta.OBP)
  }
  theta.OBP.sim <- apply(sim.matrix.Y, 2, FUN = obpFH.wrapper)

  ##### the Naive MSPE estimator
  bphi <- MSPE.calc(theta.OBP.sim, sim.matrix.theta)

  return(bphi)

}



#' McSpline MSPE estimator
#'
#' This function computes the McSPline MSPE estimator of the Observed Best Predictor (OBP).
#' @param formula an object of class formula (or one that can be coerced to that class): a symbolic description of the model to be fitted. The variables included in formula must have a length equal to the number of small areas. More about the model specification are given under Details.
#' @param data optional data frame containing the variable names in \code{formula}.
#' @param errorvar vector containing the variances of the random error for each small area.
#' @param A.BPE optional BPE estimate of random effects variance.
#' @param K number of Monte Carlo simulations. Default is 1000.
#' @details
#' \code{formula} is specified in the form \code{response ~ predictor} where the predictor is univariate. \code{formula} has an implied intercept term. To remove the intercept term, use either \code{y ~ x - 1} or \code{y ~ 0 + x}.\cr
#' \cr If \code{A.BPE} is missing, the function computes the BPE from data. User can also provide the true A instead if that is known.
#' @return The function will return a list with the following object.
#' \item{McSpline}{ McSpline estimator of the MSPE of OBP.}
#' @references Bandyopadhyay R, Jiang J (2017) "Benchmarking the Observed Best Predictor"
#' @references Jiang J, Nguyen T, and Rao J. S. (2011), "Best Predictive Small Area Estimation", Journal of the American Statistical Association.
#' @export


MSPE_McSpline <- function(formula, data, errorvar, A.BPE, K=1000){
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

  obpFH.wrapper <- function(y_wrapper){
    for (ii in ls(parent.frame())){
      assign(ii,get(ii,parent.frame()))
    }
    return(obpFH(y_wrapper~x, errorvar = D)$theta.OBP)
  }
  theta.OBP.sim <- apply(sim.matrix.Y, 2, FUN = obpFH.wrapper)

  ##### the Naive MSPE estimator
  b.phi <- MSPE.calc(theta.OBP.sim, sim.matrix.theta)

  return(list(McSpline=b.phi))
}





#' McJack MSPE estimator
#'
#' This function computes the McJack MSPE estimator of the Observed Best Predictor (OBP).
#' @param formula an object of class formula (or one that can be coerced to that class): a symbolic description of the model to be fitted. The variables included in formula must have a length equal to the number of small areas. More about the model specification are given under Details.
#' @param data optional data frame containing the variable names in \code{formula}.
#' @param errorvar vector containing the variances of the random error for each small area.
#' @param A.BPE optional BPE estimate of random effects variance.
#' @param K number of Monte Carlo simulations. Default is 1000.
#' @param returnMcSpline logical. Returns McSpline estimator with McJack (default).
#' @details
#' \code{formula} is specified in the form \code{response ~ predictor} where the predictor is univariate. \code{formula} has an implied intercept term. To remove the intercept term, use either \code{y ~ x - 1} or \code{y ~ 0 + x}.\cr
#' \cr If \code{A.BPE} is missing, the function computes the BPE from data. User can also provide the true A instead if that is known.
#' @return The function will return a list with the following objects.
#' \item{McJack}{ McJack estimator of the MSPE of OBP.}
#' \item{McSpline}{ McSpline estimator of the MSPE of OBP. This is returned by default. To turn this feature off, set \code{returnMcSpline = FALSE}.}
#' @references Bandyopadhyay R, Jiang J (2017) "Benchmarking the Observed Best Predictor"
#' @references Jiang J, Nguyen T, and Rao J. S. (2011), "Best Predictive Small Area Estimation", Journal of the American Statistical Association.
#' @export






MSPE_McJack <- function(formula, data, errorvar, A.BPE, K=1000, returnMcSpline=TRUE){
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


  b.phi = MSPE_McSpline(y~x, errorvar=errorvar, A.BPE=A.BPE, K=K)$McSpline

  b.phi.notj = matrix(NA,m,m)

  if(is.matrix(errorvar)) {
    D = diag(errorvar)
  }else if(length(errorvar)==1) {
    D = c(rep(errorvar,length(y)))
  }else {
    D = errorvar
  }

  for(j in 1:m){
    b.phi.notj[,j] = mcjack(y, x, j, D, K)
  }


  #### Calculate the Jacknife MSPE estimator
  b.phi.mat = matrix(b.phi,m,1)%*%matrix(1,nrow=1,ncol=m)
  Adj.factor.mat = b.phi.notj - b.phi.mat
  Adj.factor = rowSums(Adj.factor.mat)
  b.phi.jackknife = b.phi - ((m-1)/m)*Adj.factor

  if (returnMcSpline)return(list(McJack = b.phi.jackknife, McSpline = b.phi))
  else return(list(McJack = b.phi.jackknife))


}
