######################################################
###  Functions to be called in the main algorithm  ###
######################################################

#' The quantile check function
#'
#' @param r A vector or scaler
#' @param tau Percentile
#'
#' @return Quantile loss
#' @noRd
rq.loss<- function(r, tau){
  (abs(r)+(2*tau-1)*r)/2
}

# tanh loss (not used for current work)
tanh.loss<- function(r, gamma) gamma*log(cosh(r/gamma))




#' The Huber loss function
#'
#' @param r A vector or scaler
#' @param gamma Huber parameter
#'
#' @return Huber loss
#' @noRd
huber.loss <- function(r,gamma){
  le.ind<- which(abs(r) <= gamma)
  if(length(le.ind)==0){
    return( abs(r)-gamma/2)
  } else{
    r[le.ind]<- (r[le.ind])^2/2/gamma
    r[-le.ind]<- abs(r[-le.ind])-gamma/2
    return(r)
  }
}
##############################


## Huberized quantile loss (not used in this implementation)
rq.huber<- function(r, tau, gamma){
  r<- as.vector(r)
  (huber.loss(r,gamma)+(2*tau-1)*r)/2
}
##############################


#' Scale back coefficients
#'
#' @param coefs Coefficient vector
#' @param mu_x 
#' @param sigma_x 
#' @param intercept 
#' @noRd
transform_coefs <- function(coefs,mu_x,sigma_x,intercept=TRUE){
  if(intercept){
    new_coefs<- coefs[-1]/sigma_x
    intercept<- coefs[1] - sum(coefs[-1]*mu_x/sigma_x)
    new_coefs<- c(intercept, new_coefs)
  } else{
    new_coefs<- coefs/sigma_x
  }
  new_coefs
}
############################


#  
#' First order derivative w.r.t. residual
#'
#' @param r residual
#' @param tau Percentile
#' @param gamma Huber parameter
#'
#' @return First order derivative, a vector of n
#' @noRd
#'
rq.huber.deriv<- function(r, tau, gamma){
  r<- as.vector(r)
  le.ind<- which(abs(r) <= gamma)
  if(length(le.ind)!= 0){
    l.vec<- r
    l.vec[le.ind]<- (r[le.ind]/gamma+(2*tau-1))/2
    l.vec[-le.ind]<- (sign(r[-le.ind])+(2*tau-1))/2
  } else{
    l.vec <- (sign(r)+(2*tau-1))/2
  }
  return(l.vec)
} # end of function

## hypertangent approach, currently not used
rq.tanh.deriv<- function(r, tau, gamma){
  (tanh(r/gamma)+2*tau-1)/2
} # end of function
############################


#' Negative gradient of huberized quantile loss (w.r.t. beta)
#'
#' @param r Residual
#' @param weights Observation weight
#' @param tau Huber parameter
#' @param gamma Percentile
#' @param x Design matrix
#' @param n Sample size
#' @param apprx Approximation method, huber or tanh
#'
#' @return Gradient vector
#' @noRd
neg.gradient <- function(r,weights,tau,gamma,x,n,apprx){
  if(apprx=="huber"){
    wt_deriv <- as.vector(weights*rq.huber.deriv(r, tau, gamma))
  }else{
    wt_deriv <- as.vector(weights*rq.tanh.deriv(r, tau, gamma))
  }
  
  if(is.null(dim(x))){
    mean(x*wt_deriv)
  } else{
    apply(x*wt_deriv,2,mean)
  }
} # end of function

## l2 norm
l2norm<- function(x){
  sqrt(sum(x^2))
} # end of function


#' A function used for kkt condition check in order to verify strong rule
#'
#' @param r Residual
#' @param weights Observation weights
#' @param w Weight of shrinkage parameter. Default is square root of each group size.
#' @param gamma Huber parameter
#' @param tau Percentile
#' @param group.index A vector of group index, e.g., (1,1,1,2,2,2,3,3)
#' @param inactive.ind A vector of index for group with zero coefficients according to strong rule, e.g., (1,2,4) 
#' @param lambda Shrinkage parameter
#' @param x Reduced design matrix
#' @param n Sample size
#' @param apprx Approximation method
#'
#' @return True and False to indicate if KKT condition is met.
#' @noRd
kkt_check<- function(r, weights, w, gamma, tau, group.index, inactive.ind, lambda, x, n, apprx){
  grad<- -neg.gradient(r, weights, tau, gamma, x, n, apprx)
  grad.norm<- tapply(grad, group.index, l2norm)/w
  bad_spots <- grad.norm > lambda
  if(sum(bad_spots[inactive.ind])==0){
    list(kkt=TRUE,new_groups=NULL)
  } else{
    list(kkt=FALSE,new_groups=inactive.ind[which(bad_spots[inactive.ind])])
  }
} # end of function

