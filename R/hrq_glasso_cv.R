
#########################
### Cross validation  ###
#########################

#' Cross-validation for quantile regression with group lasso
#' 
#' K fold cross-validation is conducted. Four types of loss (squared error (\code{se}), absolute error (\code{ae})
#' quantile check loss (\code{check}) and huber loss (\code{he})) 
#' can be specified as the CV criterion. 
#'
#' @param x Design matrix
#' @param y Response variable
#' @param tau Percentage
#' @param group.index A vector of group index, e.g., (1,1,1,2,2,2,3,3) 
#' @param k Number of folders.
#' @param loss The loss function used for computing the cross-validation error. Supported losses include squared error (\code{se}), absolute error (\code{ae}), quantile check loss (\code{check}) and huber loss (\code{he}).
#' @param method Choice for mean or quantile regression. Default is \code{quantile}.
#' @param folds A vector of folder index for all observations. The procedure random splits if this argument is not specified.
#' @param gamma Huber parameter 
#' @param ... Other inputs of function \code{hrq_glasso()}.
#'
#' @return The full solution path is returned. It also returns the vector of CV score 
#' as well as the optimal values in terms of min and 1se. Corresponding lambda values are also returned. 
#' \item{beta}{The estimated coefficients for all lambdas, stored in sparse matrix format, where each column corresponds to a lambda.}
#' \item{lambda}{The sequence of lambdas.}
#' \item{lambda.min}{The optimal lambda that minimizes the CV error}
#' \item{lambda.1se}{The largest lambda such that CV error is within 1 standard error of the minimum CV error.}
#' \item{cv.all}{The vector of all values of CV error for all lambdas.}
#' \item{cv.min}{The value of CV error corresponding to \code{lambda.min}.}
#' \item{cv.1se}{The value of CV error corresponding to \code{lambda.1se}.}
#' \item{folds}{The vector of indices for k folds split.}
#' \item{cvup}{CV error + 1 standard error}
#' \item{cvlo}{CV error + 1 standard error}
#' \item{n.grp}{The number of selected groups for each lambda.}
#' @export
#'
#' @examples
#' n<- 100
#' p<- 10
#' x0<- matrix(rnorm(n*p),n,p)
#' X<- cbind(x0, x0^2, x0^3)[,order(rep(1:p,3))]
#' y<- -2+X[,1]+0.5*X[,2]-X[,3]-0.5*X[,7]+X[,8]-0.2*X[,9]+rt(n,2)
#' group<- rep(1:p, each=3)
#' fitcv<- cv.hrq_glasso(x=X, y=y, group.index=group, method="quantile")
#' plot(fitcv)
cv.hrq_glasso<- function(x, y, tau=0.5, group.index=NULL, k=5, loss="check", method="quantile", folds=NULL, gamma=0.2, ...){
  
  fullmodel<- hrq_glasso(x=x, y=y, group.index=group.index, ...)
  lambda<- fullmodel$lambda
  
  n <- length(y)
  if(is.null(folds)){
    folds <- cut(1:n, k, labels = F)[sample(n)]
  } else{
    k <- max(folds)
  }
  
  nlambda<- length(lambda)
  
  mse<- mae<- mqe<- mhe <- matrix(0, nlambda, k);
  for(i in 1:k){
    ind<- which(folds==i)
    train_x<- x[-ind,]
    train_y<- y[-ind]
    test_x<- x[ind,]
    test_y<- y[ind]
    
    train_model<- hrq_glasso(x=train_x, y=train_y, group.index=group.index, lambda=lambda, lambda.discard=FALSE,...)
    pred<- predict.hrq_glasso(train_model, newX = test_x)
    if(loss == "se"){
      se<- (test_y-pred)^2
      mse[,i]<- apply(se,2,mean)[1:nlambda]
    } 
    if(loss == "ae"){
      ae<- abs(test_y-pred)
      mae[,i]<- apply(ae,2,mean)[1:nlambda]
    } 
    if(loss == "check"){
      eq<- sapply(1:nlambda, function(xx) rq.loss((test_y-pred[,xx]), tau))
      mqe[,i]<- apply(eq,2,mean)[1:nlambda]
    } 
    if(loss == "he"){
      eh<- sapply(1:nlambda, function(xx) huber.loss((test_y-pred[,xx]),gamma))  
      mhe[,i]<- apply(eh,2,mean)[1:nlambda]
    } 
  }
  cv.mse<- apply(mse, 1, mean)
  cv.mae<- apply(mae, 1, mean)
  cv.mqe<- apply(mqe, 1, mean)
  cv.mhe<- apply(mhe, 1, mean)
  cv.mse.1se<- apply(mse, 1, sd)
  cv.mae.1se<- apply(mae, 1, sd)
  cv.mqe.1se<- apply(mqe, 1, sd)
  cv.mhe.1se<- apply(mhe, 1, sd)
  
  if(loss == "check"){
    cv.ind.min<- which(cv.mqe==min(cv.mqe))
    cv.ind.1se<- which(cv.mqe<=min(cv.mqe)+cv.mqe.1se[cv.ind.min])
    lambda.1se<- max(lambda[cv.ind.1se])
    output<- list(beta=fullmodel$beta, lambda=lambda, lambda.min=lambda[cv.ind.min], lambda.1se=lambda.1se, 
                  cv.all= cv.mqe, folds=folds, cv.min=cv.mqe[cv.ind.min], cv.1se=cv.mqe[which(lambda==lambda.1se)], 
                  cvup=cv.mqe+cv.mqe.1se, cvlo=cv.mqe-cv.mqe.1se, n.grp=fullmodel$n.grp)
    class(output) <- "cv.hrq_glasso"
    return(output)
  }else{
    if(loss=="se"){
      cv.ind.min<- which(cv.mse==min(cv.mse))
      cv.ind.1se<- which(cv.mse<=min(cv.mse)+cv.mse.1se[cv.ind.min])
      lambda.1se<- max(lambda[cv.ind.1se])
      output<- list(beta=fullmodel$beta, lambda=lambda, lambda.min=lambda[cv.ind.min], lambda.1se=lambda.1se, 
                    cv.all= cv.mse, folds=folds, cv.min=cv.mse[cv.ind.min], cv.1se=cv.mse[which(lambda==lambda.1se)], 
                    cvup=cv.mse+cv.mse.1se, cvlo=cv.mse-cv.mse.1se, n.grp=fullmodel$n.grp)
      class(output) <- "cv.hrq_glasso"
      return(output)
    }else{
      if(loss=="ae"){
        cv.ind.min<- which(cv.mae==min(cv.mae))
        cv.ind.1se<- which(cv.mae<=min(cv.mae)+cv.mae.1se[cv.ind.min])
        lambda.1se<- max(lambda[cv.ind.1se])
        output<- list(beta=fullmodel$beta, lambda=lambda, lambda.min=lambda[cv.ind.min], lambda.1se=lambda.1se, 
                      cv.all= cv.mae, folds=folds, cv.min=cv.mae[cv.ind.min], cv.1se=cv.mae[which(lambda==lambda.1se)], 
                      cvup=cv.mae+cv.mae.1se, cvlo=cv.mae-cv.mae.1se, n.grp=fullmodel$n.grp)
        class(output) <- "cv.hrq_glasso"
        return(output)
      }else{
        if(loss == "he"){
          cv.ind.min<- which(cv.mhe==min(cv.mhe))
          cv.ind.1se<- which(cv.mhe<=min(cv.mhe)+cv.mhe.1se[cv.ind.min])
          lambda.1se<- max(lambda[cv.ind.1se])
          output<- list(beta=fullmodel$beta, lambda=lambda, lambda.min=lambda[cv.ind.min], lambda.1se=lambda.1se, 
                        cv.all=cv.mhe, folds=folds, cv.min= cv.mhe[cv.ind.min], cv.1se=cv.mhe[which(lambda==lambda.1se)], 
                        cvup=cv.mhe+cv.mhe.1se, cvlo=cv.mhe-cv.mhe.1se, n.grp=fullmodel$n.grp)
          class(output) <- "cv.hrq_glasso"
          return(output)
        } else{
          stop("'loss' is not supported! Supported loss functions are 'check', 'ae', 'se' and 'he'.")
        }
      }
    }
  }
  
}# end of function
############################



