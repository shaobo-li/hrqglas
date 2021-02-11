
############################################################
#' Prediction for the hrq_glasso object
#' 
#' This function provides the prediction of the \code{hrq_glasso} object.
#'
#' @param object The model object of \code{hrq_glasso}.
#' @param newX New design matrix.
#' @param s Value of lambda. Default is NULL, so that the function provides prediction at all lambdas used in \code{hrq_glasso}.
#' @param ... other input parameters.
#'
#' @return The function returns predicted values based on the fitted model from \code{hrq_glasso}.
#' @export
#'
#' @examples
#' n<- 100
#' p<- 10
#' x0<- matrix(rnorm(n*p),n,p)
#' X<- cbind(x0, x0^2, x0^3)[,order(rep(1:p,3))]
#' y<- -2+X[,1]+0.5*X[,2]-X[,3]-0.5*X[,7]+X[,8]-0.2*X[,9]+rt(n,2)
#' group<- rep(1:p, each=3)
#' fit<- hrq_glasso(X, y, group)
#' pred<- predict(fit, newX=X, s=0.3)
predict.hrq_glasso<- function(object, newX, s=NULL, ...){
  lambda<- object$lambda
  if(length(lambda)==1){
    beta<- as.vector(object$beta)
    if(length(beta)==ncol(newX)+1){
      pred<- as.matrix(cbind(1,newX))%*%beta
    }else{
      stop("newX has wrong dimension!")
    }
    
  }else{
    if(is.null(s)){
      beta<- as.matrix(object$beta)
      if(nrow(beta)==ncol(newX)+1){
        pred<- as.matrix(cbind(1,newX))%*%beta
      }else{
        stop("new X has wrong dimension!")
      }
    }else{
      ind<- which.min((lambda-s)^2)
      beta<- as.vector(object$beta[,ind])
      if(length(beta)==ncol(newX)+1){
        pred<- as.matrix(cbind(1,newX))%*%beta
      }else{
        stop("new X has wrong dimension!")
      }
    }
  }
  
  return(pred)
}


