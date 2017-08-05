#' Observed best predictor for Fay-Herriot model.
#' 
#' This function computes the observed best predictor (OBP) for Fay-Herriot model. The variance of the random error can be specified by the user. Otherwise the function will calculate its Best Predictive Estimator (BPE). In the process of of computing OBP it also calculates the BPE of the regression coefficients of the fixed effect 
#' @param formula an object of class formula (or one that can be coerced to that class): a symbolic description of the model to be fitted. The variables included in formula must have a length equal to the number of domains D. Details of model specification are given under Details.
#' @param data data frame containing the variable names in formula and errorvar.
#' @param errorvar vector containing the D sampling variances of direct estimators for each domain. The values must be sorted as the variables in formula.
#' @param randvar varinace of the random effect. If not supplied, BPE is estimated.
#' @param maxiter maximum number of iterations used in estimating randvar.
#' @param precision covergence tolerance limit for estimating randvar.
#' @return The function will return a list consisting of the OBP of the small area mean, BPE of the regression coefficient of the fixed effect and BPE of variance of the random effect (if not specified by the user).
#' @references Jiang J, Nguyen T, and Rao J. S. (2011), "Best Predictive Small Area Estimation", Journal of the American Statistical Association.
#' @import sae
#' @export

obpFH <- function(formula, data, errorvar, randvar=NULL, maxiter=100, precision=1e-04){
  #namevar <- deparse(substitute(errorvar))
  if (!missing(data)) {
    formuladata <- model.frame(formula, na.action = na.omit, data)
    X <- model.matrix(formula, data)
    #vardir <- data[, namevar]
  }
  else {
    formuladata <- model.frame(formula, na.action = na.omit)
    X <- model.matrix(formula)
  }
  y <- formuladata[,1]
  if (attr(attributes(formuladata)$terms, "response") == 1) textformula <- paste(formula[2], formula[1], formula[3])
  else textformula <- paste(formula[1], formula[2])
  if (length(na.action(formuladata)) > 0) stop("Argument formula=", textformula, " contains NA values.")
  if (any(is.na(errorvar))) stop("Argument errorvar=", namevar, " contains NA values.")

  D = diag(errorvar) #diagonalize error variance
  I = diag(1,length(y)) #identity matrix for Fay-Herriot model
  
  ##### Q(A) function #####
  Q.A <- function(A){
    V <- D + A * I #%*% t(Z)
    B <- A * ginv(V)
    T <- diag(1,nrow(B)) - B
    P.TX <- T %*% X %*% ginv(t(X) %*% T^2 %*% X) %*% t(X) %*% T
    Q <- t(y) %*% T %*% (diag(1,nrow(P.TX)) - P.TX) %*% T %*% y + 2 * A * (sum(diag(T)))
    return(Q)
  }
  
  if (!is.null(randvar)){
    A <- randvar
  }
  else{
    A <- optimize(f = Q.A, interval = c(0,10000), tol = precision)$minimum
  }
  
  # Variance of y
  V <- D + A * I #%*% t(Z)
  
  # Calculation for beta.BPE
  B <- A * ginv(V)
  T <- diag(1,nrow(B)) - B
  K <- t(X) %*% T^2 %*% X
  beta.BPE <- ginv(K) %*% t(X) %*% T^2 %*% y
  
  theta <- y - T %*% (y - X %*% beta.BPE)
  return(list(theta.OBP = theta, A.BPE = A, beta.BPE = beta.BPE))
}


#' Best predictive estimator for Fay-Herriot model.
#' 
#' This function computes the best predictive estimator (BPE) for Fay-Herriot model. The variance of the random error can be specified by the user. Otherwise the function will calculate its Best Predictive Estimator (BPE). In the process of of computing OBP it also calculates the BPE of the regression coefficients of the fixed effect 
#' @param formula an object of class formula (or one that can be coerced to that class): a symbolic description of the model to be fitted. The variables included in formula must have a length equal to the number of domains D. Details of model specification are given under Details.
#' @param data data frame containing the variable names in formula and errorvar.
#' @param errorvar vector containing the D sampling variances of direct estimators for each domain. The values must be sorted as the variables in formula.
#' @param randvar varinace of the random effect. If not supplied, BPE is estimated.
#' @param maxiter maximum number of iterations used in estimating randvar.
#' @param precision covergence tolerance limit for estimating randvar.
#' @return The function will return a list consisting of the OBP of the small area mean, BPE of the regression coefficient of the fixed effect and BPE of variance of the random effect (if not specified by the user).
#' @references Jiang J, Nguyen T, and Rao J. S. (2011), "Best Predictive Small Area Estimation", Journal of the American Statistical Association.
#' @import sae
#' @export


bpeFH <- function(formula, data, errorvar, randvar=NULL, maxiter=100, precision=1e-04){
  #namevar <- deparse(substitute(errorvar))
  if (!missing(data)) {
    formuladata <- model.frame(formula, na.action = na.omit, data)
    X <- model.matrix(formula, data)
    #vardir <- data[, namevar]
  }
  else {
    formuladata <- model.frame(formula, na.action = na.omit)
    X <- model.matrix(formula)
  }
  y <- formuladata[,1]
  if (attr(attributes(formuladata)$terms, "response") == 1) textformula <- paste(formula[2], formula[1], formula[3])
  else textformula <- paste(formula[1], formula[2])
  if (length(na.action(formuladata)) > 0) stop("Argument formula=", textformula, " contains NA values.")
  if (any(is.na(errorvar))) stop("Argument errorvar=", namevar, " contains NA values.")
  
  D = diag(errorvar) #diagonalize error variance
  I = diag(1,length(y)) #identity matrix for Fay-Herriot model
  
  ##### Q(A) function #####
  Q.A <- function(A){
    V <- D + A * I #%*% t(Z)
    B <- A * ginv(V)
    T <- diag(1,nrow(B)) - B
    P.TX <- T %*% X %*% ginv(t(X) %*% T^2 %*% X) %*% t(X) %*% T
    Q <- t(y) %*% T %*% (diag(1,nrow(P.TX)) - P.TX) %*% T %*% y + 2 * A * (sum(diag(T)))
    return(Q)
  }
  
  if (!is.null(randvar)){
    A <- randvar
  }
  else{
    A <- optimize(f = Q.A, interval = c(0,10000), tol = precision)$minimum
  }
  
  # Variance of y
  V <- D + A * I #%*% t(Z)
  
  # Calculation for beta.BPE
  B <- A * ginv(V)
  T <- diag(1,nrow(B)) - B
  K <- t(X) %*% T^2 %*% X
  beta.BPE <- ginv(K) %*% t(X) %*% T^2 %*% y
  return(list(A.BPE = A, beta.BPE = beta.BPE))
}




