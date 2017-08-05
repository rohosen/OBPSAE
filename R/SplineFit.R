#' Fitting the optimal spline by cross validation method.
#' 
#' This function fits an optimal spline function to the given data set after choosing optimal degree of the polynomial and optimal number of knots by cross-validation. 
#' @param x independent variable values.
#' @param y dependent variable values
#' @return A list containing fitted values, optimal degree and number of knots, coefficients.
#' @import sae
#' @export

fitSpline.cv <- function(x,y){
  optimum.pq <- OptimizeSpline(x,y)
  temp <- fitSpline(x,y, splineDegree = optimum.pq[1], numKnots = optimum.pq[2])
  f.hat <- temp$yhat
  p <- optimum.pq[1]
  q <- optimum.pq[2]
  beta <- temp$beta
  gamma <- temp$gamma
  knots <- temp$knots
  return(list(f.hat=f.hat, p=p, q=q, beta=beta, gamma=gamma, knots=knots))
}

OptimizeSpline <- function(x,y,DegreeRange=c(1:3),numKnotsRange=c(1:5)){
  CVmatrix = matrix(0,length(DegreeRange),length(numKnotsRange))
  for(i in 1:nrow(CVmatrix)){
    for(j in 1:ncol(CVmatrix)){
      temp = fitSpline(x,y,DegreeRange[i],numKnotsRange[j])
      CVmatrix[i,j] = CV(y,temp$yhat,temp$H)
    }
  }
  min.cv <- which(CVmatrix == min(CVmatrix), arr.ind = TRUE)
  return(min.cv)
}


#' Fitting the optimal spline by cross validation method.
#' 
#' This function fits an optimal spline function to the given data set after choosing optimal degree of the polynomial and optimal number of knots by cross-validation. 
#' @param x independent variable values.
#' @param y dependent variable values
#' @return A list containing fitted values, optimal degree and number of knots, coefficients.
#' @import sae
#' @export

fitSpline <- function(x,y,splineDegree,numKnots){
  if(numKnots<=0){
    stop("Number of knots should be a positive integer.")
  }
  knots = seq(min(x),max(x),len=(numKnots+2))[2:(numKnots+1)]
  X = matrix(1,length(x),1)
  if(splineDegree > 0){
    for(i in 1:splineDegree){
      X = cbind(X,x^i)
    }
  }
  
  for(i in 1:length(knots)){
    X = cbind(X, sapply(x-knots[i], function(u)max(u,0))^splineDegree)
  }
  paramEst = ginv(t(X) %*% X) %*% t(X) %*% y
  H = X %*% ginv(t(X) %*% X) %*% t(X)
  yhat = H %*% y
  beta = paramEst[1:(splineDegree+1)]
  gamma = paramEst[(splineDegree+2):(splineDegree+numKnots+1)]
  
  
  return(list(yhat = yhat, H = H, beta = beta, gamma = gamma, knots = knots))
}

CV <- function(y,yhat,H){
  return(sum(((y-yhat)/(1-diag(H)))^2))
}

GCV <- function(y,yhat,H){
  return(sum(((y-yhat)/(1-mean(diag(H))))^2))
}

