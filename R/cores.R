
#' Core function for updating beta
#'
#' @param x Design matrix
#' @param y Response variable
#' @param tau Percentile
#' @param gamma Huber parameter
#' @param weights Observation weights
#' @param group.index A vector of group index, e.g., (1,1,1,2,2,2,3,3)
#' @param lambdaj Shrinkage parameter
#' @param w.lambda Weights for Shrinkage parameter of each group
#' @param eigenval Largest eigenvalue of the submatrix of H corresponding to each group
#' @param beta.ini Initial or previously updated beta
#' @param max_iter Maximum number of iteration
#' @param apprx Approximation method
#' @param epsilon Convergence criterion
#'
#' @return
#' \item{beta.update}{the updated coefficient estimates}
#' \item{converge}{indicator of convergence}
#' \item{r}{updated residual}
#' \item{iter}{number of iteration}
#' \item{grad}{gradient vector}
#' 
#' @export
solvebeta<- function(x, y, tau, gamma, weights, group.index, lambdaj, w.lambda, eigenval, beta.ini, max_iter, apprx, epsilon=1e-4){
  
  x<- as.matrix(x)
  np<- dim(x); n<- np[1]; p<- np[2]
  ng<- length(unique(group.index))
  nng<- table(group.index)
  
  int0<- beta.ini[1]
  beta0<- beta1<- beta.ini[-1]
  r0<- y-int0-x%*%beta0
  
  delta<- 2; iter<- 0; #delta.seq<- 2
  while (delta>epsilon & iter<max_iter) {
    
    iter<- iter+1
    
    if(!(apprx %in% c("huber", "tanh"))){
      stop("must use huber or tanh for apprx")
    }
    int1 <- int0 + neg.gradient(r0, weights, tau, gamma, rep(1,n), n, apprx)*gamma
    
    #r0<- as.vector(y-int1-as.matrix(x)%*%beta0)
    r0<- r0-int1+int0
    
    # update for each group
    for(k in unique(group.index)){
      
      ind<- which(group.index==k)
      u<- neg.gradient(r0, weights, tau, gamma, x[,ind], n, apprx)
      temp1<- (u+eigenval[k]*beta0[ind])
      temp2<- l2norm(temp1)
      if(temp2<=lambdaj*w.lambda[k]){
        beta1_k<- 0
      }else{
        beta1_k<- temp1*(1-lambdaj*w.lambda[k]/temp2)/eigenval[k]
      }
      
      beta1[ind]<- beta1_k
      #r0<- as.vector(y-int1-as.matrix(x)%*%beta1)
      r0<- as.vector(r0-as.matrix(x[,ind])%*%(beta1_k-beta0[ind]))
      
    }
    
    delta<- max(abs(c(int1, beta1)- c(int0, beta0)))
    #delta.seq<- c(delta.seq, delta)
    
    beta0<- beta1
    int0<- int1
    
  }
  
  if(iter< max_iter){
    converge<- TRUE
  }else{
    converge<- FALSE
  }
  
  return(list(beta.update=c(int0, beta0), converge=converge, r=r0, iter=iter, grad=-u))
  
} # end of function
########################################
