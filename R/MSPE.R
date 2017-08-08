#' Bootstrap MSPE Estimator of OBP
#'
#' This function computes the bootstrap MSPE estimator of OBP.
#' @param theta.OBP A vector of OBPs.
#' @param D Scalar, Vector (or diagonal matrix) of error variances.
#' @param L is the number of bootstrap samples. the default is 200.
#' @return This function will return the Bootstrap MSPE estimator of OBP.
#' @references Jiang J, Nguyen T, and Rao J. S. (2011), "Best Predictive Small Area Estimation", Journal of the American Statistical Association.
#' @import mvtnorm
#' @export
#'
### Bootstrap MSPE Estimator of OBP ###

MSPE_boot <- function(theta.OBP, D, x, L=200){

  if(length(D) == 1){
    D <- D * diag(rep(1,length(theta.OBP)))
    } else if(is.vector(D)){
      D <- diag(D) #If D is a vector, convert it to a diagonal matrix.
    }

  ## Bootstrap Sample
  normal.sample <- mvtnorm::rmvnorm(L, mean = theta.OBP, sigma = D)
  y.boot = t(normal.sample)

  ## Calculate Bootstrap MSPE Estimator
  obpFH.wrapper <- function(y_wrapper){
    for (ii in ls(parent.frame())){
      assign(ii,get(ii,parent.frame()))
    }
    return(obpFH(y_wrapper~x, errorvar=diag(D))$theta.OBP)
  }
  OBP.boot <- apply(y.boot, 2, FUN = obpFH.wrapper)

  theta.OBP.mat <- matrix(rep(theta.OBP, L), ncol = L)

  square.mat <- (OBP.boot - theta.OBP.mat)^2

  MSPE_boot_est <- rowMeans(square.mat)

  return(MSPE_boot_est)

}

#' Bootstrap MSPE estimator of adjusted OBP
#'
#' This function computes the bootstrap MSPE estimator of benchmarked adjusted OBP.
#' @param theta.OBP.adjusted A vector of adjusted OBPs.
#' @param D Scalar, Vector (or diagonal matrix) of error variances.
#' @param x The independant variable values
#' @param weight Weight of each small area for calculating the benchmarked OBP (generally sampling weights of each small area).
#' @param L is the number of bootstrap samples. the default is 200.
#' @return This function will return the Bootstrap MSPE estimator of Adjusted OBP.
#' @references Jiang J, Nguyen T, and Rao J. S. (2011), "Best Predictive Small Area Estimation", Journal of the American Statistical Association.
#' @import mvtnorm
#' @export
#'
### Bootstrap MSPE Estimator of OBP ###

MSPE_boot_adjusted <- function(theta.OBP.adjusted, D, x, weight, L=200){

  if(length(D) == 1){
    D <- D * diag(rep(1,length(theta.OBP.adjusted)))
    } else if(is.vector(D)){
      D <- diag(D) #If D is a vector, convert it to a diagonal matrix.
    }
  normal.sample <- mvtnorm::rmvnorm(L, mean = theta.OBP.adjusted, sigma = D)
  y.boot = t(normal.sample)

  ### Calculate Bootstrap MSPE for Adjusted OBP
  obpFH.adjusted.wrapper <- function(y_wrapper){
    for (ii in ls(parent.frame())){
      assign(ii,get(ii,parent.frame()))
    }
    return(obpFH_adjusted(y_wrapper~x, errorvar=diag(D), weight = weight)$obpAdjusted)
  }
  OBP.boot.adjusted <- apply(y.boot, 2, FUN = obpFH.adjusted.wrapper)

  theta.OBP.adjusted.mat <- matrix(rep(theta.OBP.adjusted, L), ncol = L)

  square.mat.adjusted <- (OBP.boot.adjusted - theta.OBP.adjusted.mat)^2

  MSPE_boot_adjusted_est <- rowMeans(square.mat.adjusted)

  return(MSPE_boot_adjusted_est)

}

#' Bootstrap MSPE estimator of augmented OBP
#'
#' This function computes the bootstrap MSPE estimator of benchmarked augmented OBP.
#' @param theta.OBP.adjusted A vector of augmented OBPs.
#' @param D Scalar, Vector (or diagonal matrix) of error variances.
#' @param x The independant variable values
#' @param weight Weight of each small area for calculating the benchmarked OBP (generally sampling weights of each small area).
#' @param L is the number of bootstrap samples. the default is 200.
#' @return This function will return the Bootstrap MSPE estimator of Augmented OBP.
#' @references Jiang J, Nguyen T, and Rao J. S. (2011), "Best Predictive Small Area Estimation", Journal of the American Statistical Association.
#' @import mvtnorm
#' @export
### Bootstrap MSPE Estimator of OBP ###

MSPE_boot_augmented <- function(theta.OBP.augmented, D, x, weight, L=200){

  if(length(D) == 1){
    D <- D * diag(rep(1,length(theta.OBP.augmented)))
    } else if(is.vector(D)){
      D <- diag(D) #If D is a vector, convert it to a diagonal matrix.
    }
  normal.sample <- mvtnorm::rmvnorm(L, mean = theta.OBP.augmented, sigma = D)
  y.boot = t(normal.sample)

  obpFH.augmented.wrapper <- function(y_wrapper){
    for (ii in ls(parent.frame())){
      assign(ii,get(ii,parent.frame()))
    }
    return(obpFH_augmented(y_wrapper~x, errorvar=diag(D), weight = weight)$obpAugmented)
  }
  OBP.boot.augmented <- apply(y.boot, 2, FUN = obpFH.augmented.wrapper)

  theta.OBP.augmented.mat <- matrix(rep(theta.OBP.augmented, L), ncol = L)

  square.mat.augmented <- (OBP.boot.augmented - theta.OBP.augmented.mat)^2

  MSPE_boot_augmented_est <- rowMeans(square.mat.augmented)

  return(MSPE_boot_augmented_est)

}


#' Naive MSPE Estimator of OBP
#'
#' This function computes the naive MSPE estimator of OBP.
#' @param y A vector of the response variable.
#' @param D Scalar, Vector (or diagonal matrix) of error variances.
#' @param theta.OBP A vector of OBP values
#' @param A.BPE BPE estimate of variance of random effects.
#' @return This function will return the Naive MSPE estimator of the OBP.
#' @references Jiang J, Nguyen T, and Rao J. S. (2011), "Best Predictive Small Area Estimation", Journal of the American Statistical Association.
#' @export
#'



MSPE_naive <- function(y, D, theta.OBP, A.BPE){

  if(is.matrix(D)){
    D = diag(D)
  } else if(length(D)==1){
    D = c(rep(D,length(y)))
  }

  B = A.BPE/(A.BPE+D)

  MSPEn = (theta.OBP - y)^2 + D*(2*B - 1)

  return(MSPEn)
}


#' JNR MSPE Estimator of OBP
#'
#' This function computes the JNR MSPE estimator of OBP.
#' @param y A vector of the response variable.
#' @param D Scalar, Vector (or diagonal matrix) of error variances.
#' @param theta.OBP A vector of OBP values
#' @param A.BPE BPE estimate of variance of random effects.
#' @return This function will return the Naive MSPE estimator of the OBP.
#' @references Jiang J, Nguyen T, and Rao J. S. (2011), "Best Predictive Small Area Estimation", Journal of the American Statistical Association.
#' @import MASS
#' @export


### JNR MSPE Estimator ###
MSPE_JNR <- function(formula, data, errorvar, theta.OBP, A.BPE, beta.BPE){

  if (!missing(data)) {
    formuladata <- model.frame(formula, na.action = na.omit, data)
    X <- model.matrix(formula, data)
  }
  else {
    formuladata <- model.frame(formula, na.action = na.omit)
    X <- model.matrix(formula)
  }
  y <- formuladata[,1]
  if (attr(attributes(formuladata)$terms, "response") == 1) textformula <- paste(formula[2], formula[1], formula[3])else textformula <- paste(formula[1], formula[2])
  if (length(na.action(formuladata)) > 0) stop("Argument formula=", textformula, " contains NA values.")
  if (any(is.na(errorvar))) stop("Argument errorvar=", namevar, " contains NA values.")

  if (missing(theta.OBP) | missing(A.BPE) | missing(beta.BPE)){
    temp = obpFH(formula, data, errorvar)
    theta.OBP = temp$theta.OBP
    A.BPE = temp$A.BPE
    beta.BPE = temp$beta.BPE
  }

  D = errorvar
  if(is.matrix(D))D = diag(D)

  MSPE_JNR_i <- function(i){

#     for (ii in ls(parent.frame())){
#       assign(ii,get(ii,parent.frame()))
#     }
    yi <- y[i]
    xi <- X[i,]

    Di <- D[i]
    Bi = A.BPE/(A.BPE+Di)
    Mi = matrix(0,(length(xi)+1),2)
    Mi[1:length(xi),1] = xi
    Mi[(length(xi)+1),2] = 1/(A.BPE + Di)

    Ui = c(yi - xi%*%beta.BPE, (yi - xi%*%beta.BPE)^2 - (A.BPE + Di))

    w11 = xi%*%t(xi)
    w12 = (2/(A.BPE+Di))*(yi - xi%*%beta.BPE)*xi
    w1 = cbind(w11, w12)
    w21 = t(w12)
    w22 = (3*(yi - xi%*%beta.BPE)^2 - (A.BPE + Di))/(A.BPE + Di)^2
    w2 = c(w21, w22)
    Wi = rbind(w1, w2)

    fi = (-2)*(1-Bi)^2*Mi%*%Ui

    deli = diag(c(rep(0, length(beta.BPE)), 1/(A.BPE + Di)))
    G2i = 2*(1-Bi)^2*(Wi - deli)

    return(list(Bi = Bi, fi = fi, Wi = Wi, G2i = G2i))
  }

  m = length(y)
  t = c(1:m)

  MJNR <- lapply(t, MSPE_JNR_i)

  G2 <- MJNR[[1]]$G2i
  for(i in 2:m){
    G2 = G2 + MJNR[[i]]$G2i
  }
  G2.inv = ginv(G2)
  h2 = G2.inv[nrow(G2.inv),]

  MSPE_i <- function(i){
    MSPE <- (theta.OBP[i] - y[i])^2 + D[i]*(2*MJNR[[i]]$Bi - 1) + 2*(1 - MJNR[[i]]$Bi)^2*(t(h2)%*%MJNR[[i]]$fi) + 4*D[i]*(1 - MJNR[[i]]$Bi)^3*sum(diag(G2.inv %*% MJNR[[i]]$Wi))

    return(MSPE)
  }

  MSPE.JNR.est <- unlist(lapply(t, MSPE_i))

  return(MSPE.JNR.est)

}
