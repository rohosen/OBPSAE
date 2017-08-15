#' Benchmarked observed best predictor for Fay-Herriot model. (agbranch)
#'
#' This function computes the benchmarked observed best predictor (OBP) for Fay-Herriot model. The variance of the random error can be specified by the user. Otherwise the function will calculate its Best Predictive Estimator (BPE). In the process of of computing OBP it also calculates the BPE of the regression coefficients of the fixed effect
#' @param formula an object of class formula (or one that can be coerced to that class): a symbolic description of the model to be fitted. The variables included in formula must have a length equal to the number of small areas. Details of model specification are given under Details.
#' @param data optional data frame containing the variable names in \code{formula}.
#' @param errorvar vector containing the variances of the random error for each small area.
#' @param weight vector containing the sampling weights of small areas. Default is uniform. If sum of the weights is not 1, the weights are normalized.
#' @param method string specifying the benchmarking method. Options are "adjusted" and "augmented". Computes both if not specified. See Details for more usage information.
#' @param randvar varinace of the random effect. If not supplied, BPE is estimated.
#' @param maxiter maximum number of iterations used in estimating randvar.
#' @param precision covergence tolerance limit for estimating randvar.
#' @details
#' \code{formula} is specified in the form \code{response ~ predictors} where the predictors are separated by \code{+}. \code{formula} has an implied intercept term. To remove the intercept term, use either \code{y ~ x - 1} or \code{y ~ 0 + x}.\cr
#' \cr If \code{method} is set to "adjusted", only \code{obpAdjusted} is returned.\cr
#' \cr If \code{method} is set to "augmented", \code{obpAugmented, A.BPE.aug} and \code{beta.BPE.aug} are returned.\cr
#'
#' @return The function will return a list with all the following objects by default.
#' \item{obpAdjusted}{a vector of adjusted OBP values.}
#' \item{obpAugmented}{a vector of augmented OBP values.}
#' \item{A.BPE.aug}{ BPE of variance of random effects.}
#' \item{beta.BPE.aug}{ BPE of fixed effects regression coefficients.}
#' @references Bandyopadhyay R, Jiang J (2017) "Benchmarking of Observed Best Predictor"
#' @references Jiang J, Nguyen T, and Rao J. S. (2011), "Best Predictive Small Area Estimation", Journal of the American Statistical Association.
#' @import sae
#' @export

obpFHbenchmark <- function(formula, data, errorvar, weight, method=c("adjusted","augmented"), randvar=NULL, maxiter=100, precision=1e-04){
  result <- list(obp = NA, fit = list(betaBPE=NA, randvarBPE=NA, convergence = TRUE, iterations = 0))
  #namevar <- deparse(substitute(errorvar))
  if (length(method) >= 3) stop("Supports only adjusted and augmented method of benchmarking.")
  result = list()
  for(mm in method){
    if(!(tolower(mm) %in% c("adjusted","augmented"))) stop("Error! Use 'adjusted' or 'augmented'")
    if(tolower(mm) == "adjusted") result = c(result, obpFH_adjusted(formula, data, errorvar, weight, randvar, maxiter, precision))
    if(tolower(mm) == "augmented") result = c(result, obpFH_augmented(formula, data, errorvar, weight, randvar, maxiter, precision))
  }
  return(result)
}


#' Adjusted observed best predictor for Fay-Herriot Model
#'
#' This function computes the Adjusted observed best predictor (OBP) for Fay-Herriot model. The variance of the random error can be specified by the user. Otherwise the function will calculate its Best Predictive Estimator (BPE). In the process of of computing OBP it also calculates the BPE of the regression coefficients of the fixed effect
#' @param formula an object of class formula (or one that can be coerced to that class): a symbolic description of the model to be fitted. The variables included in formula must have a length equal to the number of domains D. Details of model specification are given under Details.
#' @param data data frame containing the variable names in formula and errorvar.
#' @param errorvar vector containing the D sampling variances of direct estimators for each domain. The values must be sorted as the variables in formula.
#' @param weight vector containing the sampling weights of small areas. Default is uniform. If sum of the weights is not 1, the weights are normalized.
#' @param randvar varinace of the random effect. If not supplied, BPE is estimated.
#' @param maxiter maximum number of iterations used in estimating randvar.
#' @param precision covergence tolerance limit for estimating randvar.
#' @return The function will return a list containing the Adjusted OBP.
#' @references Bandyopadhyay R, Jiang J (2017) "Benchmarking of Observed Best Predictor"
#' @references Jiang J, Nguyen T, and Rao J. S. (2011), "Best Predictive Small Area Estimation", Journal of the American Statistical Association.
#' @import sae
#' @export

obpFH_adjusted <- function(formula, data, errorvar, weight, randvar=NULL, maxiter=100, precision=1e-04){

####################################
######### Adjusted method ##########
####################################
  if (!missing(data)) {
    formuladata <- model.frame(formula, na.action = na.omit, data)
    #X <- model.matrix(formula, data)
    #vardir <- data[, namevar]
  }
  else {
    formuladata <- model.frame(formula, na.action = na.omit)
    #X <- model.matrix(formula)
  }
  y <- formuladata[,1]
  if (attr(attributes(formuladata)$terms, "response") == 1) textformula <- paste(formula[2], formula[1], formula[3])
  else textformula <- paste(formula[1], formula[2])
  if (length(na.action(formuladata)) > 0) stop(paste0("Argument formula=", textformula, " contains NA values."))
  if (any(is.na(errorvar))) stop(paste0("Argument errorvar=", namevar, " contains NA values."))
  if (sum(weight) == 0 | any(weight < 0)) stop("Weights should be positive.")

### Step 1: Calculate original OBP by FH model
  obpOriginal = obpFH(formula,data,errorvar,randvar,maxiter,precision)

### Step 2: Adjust for benchmarking
  w = weight/sum(weight)
  adjFactor = w / sum(w^2)
  e = rep(1,length(y)) ## necessary unit vector
  adjTerm = sum(w * y) - sum(w * obpOriginal$theta.OBP)
  obpAdjusted = obpOriginal$theta.OBP + adjFactor * adjTerm
  return(list(obpAdjusted = obpAdjusted))
}


#' Augmented observed best predictor for Fay-Herriot model.
#'
#' This function computes the Augmented observed best predictor (OBP) for Fay-Herriot model. The variance of the random error can be specified by the user. Otherwise the function will calculate its Best Predictive Estimator (BPE). In the process of of computing OBP it also calculates the BPE of the regression coefficients of the fixed effect
#' @param formula an object of class formula (or one that can be coerced to that class): a symbolic description of the model to be fitted. The variables included in formula must have a length equal to the number of domains D. Details of model specification are given under Details.
#' @param data data frame containing the variable names in formula and errorvar.
#' @param errorvar vector containing the D sampling variances of direct estimators for each domain. The values must be sorted as the variables in formula.
#' @param weight vector containing the sampling weights of small areas. Default is uniform. If sum of the weights is not 1, the weights are normalized.
#' @param randvar varinace of the random effect. If not supplied, BPE is estimated.
#' @param maxiter maximum number of iterations used in estimating randvar.
#' @param precision covergence tolerance limit for estimating randvar.
#' @return The function will return a list containing the Augmented OBP.
#' @references Bandyopadhyay R, Jiang J (2017) "Benchmarking of Observed Best Predictor"
#' @references Jiang J, Nguyen T, and Rao J. S. (2011), "Best Predictive Small Area Estimation", Journal of the American Statistical Association.
#' @import sae
#' @export

obpFH_augmented <- function(formula, data, errorvar, weight, randvar=NULL, maxiter=100, precision=1e-04){

  ####################################
  ######### Augmented method ##########
  ####################################
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
  if (sum(weight) == 0 | any(weight < 0)) stop("Weights should be positive.")

  D = diag(errorvar) #diagonalize error variance
  I = diag(1,length(y)) #identity matrix for Fay-Herriot model
  w = weight/sum(weight)


  ##### Q(A) function #####
  Q.A <- function(A){
    for (ii in ls(parent.frame())){
      assign(ii,get(ii,parent.frame()))
    }
    D.vec <- errorvar
    gamma = A / (A + D.vec)
    gamma.inv = 1/(1 - gamma)
    gamma.inv.w = w * gamma.inv

    ### Create augmented X
    X.a = cbind(X, gamma.inv.w)

    V <- D + A * I #%*% t(Z)
    B <- A * ginv(V)
    T <- diag(1,nrow(B)) - B
    P.TX <- T %*% X.a %*% ginv(t(X.a) %*% T^2 %*% X.a) %*% t(X.a) %*% T
    Q <- t(y) %*% T %*% (diag(1,nrow(P.TX)) - P.TX) %*% T %*% y + 2 * A * (sum(diag(T)))
    return(Q)
  }

  if (!is.null(randvar)){
    A <- randvar
  }
  else{
    A <- optimize(f = Q.A, interval = c(0,1000), tol = precision)$minimum
  }

  D.vec <- errorvar
  gamma = A / (A + D.vec)
  gamma.inv = 1/(1 - gamma)
  gamma.inv.w = w * gamma.inv

  ### Create augmented X
  X.a = cbind(X, gamma.inv.w)

  # Variance of y
  V <- D + A * I #%*% t(Z)

  # Calculation for beta.BPE
  B <- A * ginv(V)
  T <- diag(1,nrow(B)) - B
  K <- t(X.a) %*% T^2 %*% X.a
  beta.BPE <- ginv(K) %*% t(X.a) %*% T^2 %*% y

  theta <- y - T %*% (y - X.a %*% beta.BPE)
  return(list(obpAugmented = theta, A.BPE.aug = A, beta.BPE.aug = beta.BPE))
}



