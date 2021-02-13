
################################
## plot

#' Generating plots for cross-validation
#' @param x The object of function \code{cv.hrq_glasso}.
#' @param ... other input parameters for the generic function \code{plot}.
#' @importFrom graphics axis arrows abline
#' @return Cross-validation plot for the entire solution path.
#' @export
#'
plot.cv.hrq_glasso<- function(x, ...){
  cv.fit<- x
  lambda<- cv.fit$lambda
  lambda.min<- cv.fit$lambda.min
  lambda.1se<- cv.fit$lambda.1se
  loglambda<- log(lambda)
  y<- cv.fit$cv.all
  plot(loglambda, y, xlab="log(lambda)", ylab="CV Error", ylim = range(cv.fit$cvup, cv.fit$cvlo), pch=16, col="red", cex = .6)
  #glmnet::error.bars(log(cv.fit$lambda), cv.fit$cvup, cv.fit$cvlo, width = 0.01, col = "darkgrey")
  arrows(loglambda,cv.fit$cvup,loglambda,cv.fit$cvlo,length=0.01, angle=90, code=3, col = "darkgrey")
  axis(side=3, at=loglambda, labels = cv.fit$n.grp, tick = FALSE, line = 0)
  abline(v = log(cv.fit$lambda.min), lty = 3)
  abline(v = log(cv.fit$lambda.1se), lty = 3)
  invisible()
}






