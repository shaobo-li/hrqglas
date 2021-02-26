
## coefficient
#' Extract coefficients from cv.hrq_glasso object
#'
#' @param object The model object \code{cv.hrq_glasso object}.
#' @param s Value of lambda. 
#' @param ... other input parameters.
#'
#' @return The function extract estimated coefficients from \code{cv.hrq_glasso} object.
#' @export
#'
coef.cv.hrq_glasso<- function(object, s, ...){
  lambda<- object$lambda
  if(missing(s)) s<- "lambda.min"
  if(s == "lambda.min") s<- object$lambda.min
  if(s == "lambda.1se") s<- object$lambda.1se
  ind<- which.min(abs(s-lambda))[1]
  beta<- as.vector(object$beta[,ind])
  names(beta)<- c("Intercept", paste("V", 1:(length(beta)-1), sep = ""))
  return(beta)
}

#' Extract coefficients from hrq_glasso object
#'
#' @param object The model object \code{hrq_glasso object}.
#' @param s Value of lambda.
#' @param ... other input parameters.
#'
#' @return The function extract estimated coefficients from \code{hrq_glasso} object.
#' @export
#'
coef.hrq_glasso<- function(object, s, ...){
  if(missing(s)){
    return(object$beta)
  }else{
    lambda<- object$lambda
    ind<- which.min(abs(s-lambda))[1]
    beta<- as.vector(object$beta[,ind])
    names(beta)<- c("Intercept", paste("V", 1:(length(beta)-1), sep = ""))
    return(beta)
  }
}



