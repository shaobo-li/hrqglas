
#####################
### main function ###
#####################
### for Huber regression, gamma=IQR(y)/10
### for Quantile regression, gamma is adaptively chosen for each lambda: gamma=min(0.8, max(Q10(|r|), gamma_min)),
###   where gamma_min can be set from (0.01, 0.5) with faster convergence at larger value, default is 0.2 

#' Robust group variable selection for quantile and mean regression
#' 
#' This function conducts group-wise (with known groups) variable selection for quantile and robust mean regression with the group lasso penalty. 
#' The Huber loss is used for both types of regression model, where the quantile check function is approximated by Huber loss. A full solution path
#' is generated unless a single value of the shrinkage parameter is specified.
#' @importFrom Matrix Matrix
#' @importFrom stats sd quantile coefficients
#' @importFrom MASS rlm 
#' @useDynLib hrqglas solve_beta
#' 
#' @param x Design matrix (in matrix format)
#' @param y Response variable
#' @param group.index A vector of group index, e.g., (1,1,1,2,2,2,3,3) 
#' @param tau Percentile
#' @param lambda Shrinkage parameter, default is NULL so that the algorithm chooses a sequence.
#' @param weights Observation weights, default is NULL
#' @param w.lambda Weights for Shrinkage parameter of each group, default is NULL
#' @param gamma Huber parameter. An initial value is 0.2, while the algorithm adaptively tunes the value in each iteration.
#' @param max_iter Maximum number of iteration
#' @param apprx Approximation method. Default is \code{huber}. The other option is \code{tanh} which uses the hypertangent function to approximate the first order derivative of absolute loss. 
#' @param lambda.discard Default is TRUE, meaning that the solution path stops if the relative deviance changes sufficiently small. It usually happens near the end of solution path. However, the program returns at least 70 models along the solution path. 
#' @param method Choice for mean or quantile regression. Default is \code{quantile}.
#' @param scalex Standardize design matrix. Default is TRUE.
#' @param epsilon The epsilon level convergence. Default is 1e-4.
#' @param beta0 Initial estimates. Default is NULL.
#'
#' @return It returns a sequence of estimated coefficients for quantile regression with group feature selection corresponding to a sequence of lambda. 
#' The estimated coefficients are in the sparse matrix format. Returned values also include the sequence of lambda, the null deviance, 
#' values of penalized loss, and unpenalized loss across the sequence of lambda. 
#' \item{beta}{The estimated coefficients for all lambdas, stored in sparse matrix format, where each column corresponds to a lambda.}
#' \item{lambda}{The sequence of lambdas.}
#' \item{null.dev}{The null deviance.}
#' \item{pen.loss}{The value of penalized loss for each lambda.}
#' \item{loss}{The value of unpenalized loss for each lambda.}
#' \item{index.grp}{Group indices that correspond to the estimated coefficient matrix \code{beta}.}
#' \item{n.grp}{The number of selected groups for each lambda.}
#' 
#' @references
#' Sherwood, B., and Li, S. (2021) An Efficient Approach to Feature Selection and Estimation for Quantile Regression with Grouped Variables. \emph{Working paper}.
#' 
#' Yang, Y., and Zou, H., (2015) A Fast Unified Algorithm for Solving Group-lasso Penalize Learning Problems, \emph{Statistics and Computing}, 25 1129-1141.
#' \doi{10.1007/s11222-014-9498-5}.
#' 
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
#' fit$beta[,8]
#' 
hrq_glasso<- function(x, y, group.index, tau=0.5, lambda=NULL, weights=NULL, w.lambda=NULL, gamma=0.2, max_iter=200, apprx="huber", 
                      lambda.discard=TRUE, method="quantile", scalex=TRUE, epsilon=1e-4, beta0=NULL){
  
  if(scalex){
    x <- scale(x)
    mu_x <- attributes(x)$`scaled:center`
    sigma_x <- attributes(x)$`scaled:scale`
  }
  
  np<- dim(x)
  n<- np[1]; p<- np[2]
  
  if(!(apprx %in% c("huber", "tanh"))){
    stop("must use 'huber' or 'tanh' for apprx")
  }
  
  if(min(group.index)!=1 | max(group.index)!=length(unique(group.index)) | !is.numeric(group.index)){
    stop("'group.index' must be numeric sequence starting from 1.")
  }
  if(is.null(weights)) weights<- rep(1, n)
  
  ng<- max(group.index)
  nng<- table(group.index)
  if(is.null(w.lambda)) w.lambda<- sqrt(nng)

  # initial value
  gamma0<- gamma
  if(is.null(beta0)){
    if(method=="quantile"){
      # intercept 
      b.int<- quantile(y, probs = tau)
      r<- y-b.int
      # null devience
      dev0<- sum(rq.loss(r,tau))
      gamma.max<- 4; gamma<- min(gamma.max, max(gamma0, quantile(abs(r), probs = 0.1)))
    } else if(method == "mean"){
      b.int <- coefficients(rlm(y ~ 1, k2=gamma))[1] # gamma stays fixed for huber regression
      r <- y - b.int
      dev0 <- sum(huber.loss(r,gamma))
    } else{
      stop("method must be 'quantile' or 'mean'. ")
    }
    beta0<- c(b.int, rep(0,p)) 
  } else{
    r <- y - beta0[1] - x%*%beta0[-1]
    if(method=="quantile"){
      dev0<- sum(rq.loss(r,tau))
      gamma.max<- 4; gamma<- min(gamma.max, max(gamma0, quantile(abs(r), probs = 0.1)))
    }else if(method == "mean"){
      dev0 <- sum(huber.loss(r,gamma))
    }else{
      stop("method must be 'quantile' or 'mean'.")
    }
  }
  r0<- r
  
  ## required by C code ##
  group_order <- order(group.index)
  group.index <- sort(group.index)
  x <- x[,group_order]
  beta_order <- c(1,group_order+1)
  
  
  ## get sequence of lambda if not supplied
  # l2norm of gradient for each group
  grad_k<- -neg.gradient(r0, weights, tau, gamma=gamma, x, n, apprx)
  grad_k.norm<- tapply(grad_k, group.index, l2norm)
  
  lambda.max<- max(grad_k.norm/w.lambda)
  lambda.flag<- 0
  if(is.null(lambda)){
    lambda.min<- ifelse(n>p, lambda.max*0.001, lambda.max*0.01)
    #lambda<- seq(lambda.max, lambda.min, length.out = 100)
    lambda<- exp(seq(log(lambda.max), log(lambda.min), length.out = 101))
  }else{
    # user supplied lambda
    lambda.discard<- FALSE
    if(lambda.max> max(lambda)){
      lambda<- exp(c(log(lambda.max), log(sort(lambda, decreasing = TRUE))))
    }else{
      if(length(lambda)>1 & min(lambda)<lambda.max){
        lambda.flag<- 1
        lambda.user<- lambda
        lambda<- exp(c(log(lambda.max), log(sort(lambda[lambda<lambda.max], decreasing = TRUE))))
      }else{
        #warning("lambda is too large, all coefficients are shrunk to 0!")
        return(list(beta=matrix(c(b.int, rep(0,p)), nrow = p+1, ncol = length(lambda) )))
      }
    }
  }
  
  ## QM condition in Yang and Zou, Lemma 1 (2) -- PD matrix H
  H<- 2*t(x)%*%diag(weights)%*%x/(n*gamma)	
  # get eigen values of sub matrices for each group
  eigen.sub.H<- rep(0, ng)
  for(k in 1:ng){
    ind<- which(group.index==k)
    sub.H<- H[ind, ind]
    eigen.sub.H[k]<- max(eigen(sub.H)$value)+1e-6
  }
  
  ### for loop over lambda's
  if(length(lambda)==2){
    lambda<- lambda[2]
    # update beta and residuals
    update <- .C("solve_beta", as.double(y), as.double(cbind(1,x)), as.double(tau), as.double(gamma), 
                 as.double(weights), as.double(lambda), as.double(w.lambda), 
                 as.double(eigen.sub.H), as.double(beta0), as.integer(max_iter), as.double(epsilon), as.integer(apprx=="huber"), 
                 as.integer(n), as.integer(p), as.integer(ng), as.integer(nng), as.integer(0), as.double(r0), as.integer(0))
    beta0<- update[[9]]
    iter_num <- update[[19]]
    converge_status <- update[[17]]==1
    update.r <- update[[18]]
    
    if(scalex){
      beta0 <- transform_coefs(beta0, mu_x, sigma_x)
    }
    beta0[beta_order] <- beta0
    
    if(method=="quantile"){
      dev1<- sum(weights*rq.loss(update.r, tau))
    } else if( method == "mean"){
      dev1<- sum(weights*huber.loss(update.r,gamma))	
    }
    
    pen.loss<- dev1/n+lambda*sum(eigen.sub.H*sapply(1:ng, function(xx) l2norm(beta0[-1][group.index==xx])))
    group.index.out<- unique(group.index[beta0[-1]!=0])
    output<- list(beta=beta0, lambda=lambda, null.dev=dev0, pen.loss=pen.loss, loss=dev1/n, tau=tau, 
                  apprx=apprx, n.grp=length(group.index.out), index.grp=group.index.out, x=x, y=y)
    output.hide<- list(converge=update$converge, iter=update$iter, rel_dev=dev1/dev0)
    class(output)<- "hrq_glasso"
    return(output)
    
  }else{
    
    group.index.out<- matrix(0, ng, length(lambda))
    n.grp<- rep(0, length(lambda)); n.grp[1]<- 0
    beta.all<- matrix(0, p+1, length(lambda))
    beta.all[,1]<- beta0
    kkt_seq<- rep(NA, length(lambda)); kkt_seq[1]<- TRUE  # for internal use
    converge<- rep(0, length(lambda)); converge[1]<- TRUE # for internal use
    iter<- rep(0, length(lambda)); iter[1]<- 1
    outer_iter_count <- rep(0,length(lambda)); outer_iter_count[1] <- 0 # for internal use
    rel_dev<- rep(0, length(lambda)); rel_dev[1]<- 1
    loss<- rep(0, length(lambda)); loss[1]<- dev0/n
    pen.loss<- rep(0, length(lambda)); pen.loss[1]<- dev0/n
    gamma.seq<- rep(0, length(lambda)); gamma.seq[1]<- gamma
    active.ind<- NULL
    for(j in 2:length(lambda)){ #
      
      if(length(active.ind)<ng){
        ## use strong rule to determine active group at (i+1) (pre-screening)
        grad_k<- -neg.gradient(r0, weights, tau, gamma=gamma, x, n, apprx)
        grad_k.norm<- tapply(grad_k, group.index, l2norm)
        active.ind<- which(grad_k.norm>=w.lambda*(2*lambda[j]-lambda[j-1])) 
        n.active.ind<- length(active.ind)
        
        if(length(active.ind)==0){
          inactive.ind<- 1:ng
          outer_iter<- 0
          kkt_met<- NA
          update.iter<- 0
          update.converge<- NA
          update.r<- r0
          
        }else{
          if(length(active.ind)>ng/2) max_iter<- 50
          
          inactive.ind<- (1:ng)[-active.ind]
          
          ## outer loop to update beta and check KKT after each update
          kkt_met <- FALSE
          outer_iter <- 0 
          
          while(!kkt_met & length(inactive.ind)>0){
            outer_iter<- outer_iter+1
            
            # reduced data
            reduced.ind<- which(group.index %in% active.ind)
            reduced.group.index<- group.index[reduced.ind]
            u.reduced.group.index <- unique(reduced.group.index)
            reduced.ng<- length(u.reduced.group.index)
            x.sub<- x[,reduced.ind]
            reduced.eigen <- eigen.sub.H[u.reduced.group.index]
            
            
            # update beta and residuals
            x.sub.c <- cbind(1,x.sub)
            residuals.c <- y - x.sub.c %*% beta0[c(1, reduced.ind+1)]
            reduced.w.lambda <- w.lambda[u.reduced.group.index]
            reduced.nng <- table(reduced.group.index)
            update <- .C("solve_beta", as.double(y), as.double(x.sub.c), as.double(tau), as.double(gamma), 
                         as.double(weights), as.double(lambda[j]), as.double(reduced.w.lambda), 
                         as.double(reduced.eigen), as.double(beta0[c(1, reduced.ind+1)]), as.integer(max_iter), as.double(epsilon), as.integer(apprx=="huber"), 
                         as.integer(n), as.integer(p), as.integer(reduced.ng), as.integer(reduced.nng), as.integer(0), as.double(residuals.c), as.integer(0))
            
            beta.update <- update[[9]]
            update.r <- update[[18]]
            update.converge <- update[[17]]
            update.iter <- update[[19]]
            
            beta0[c(1, reduced.ind+1)]<- beta.update
            
            ## check inactive set by KKT condition
            kkt_results <- kkt_check(r=update.r,weights=weights,w=w.lambda,gamma=gamma,tau=tau,group.index=group.index,
                                     inactive.ind=inactive.ind,lambda=lambda[j],x=x,n=n,apprx=apprx)
            if(kkt_results$kkt==TRUE){
              kkt_met <- TRUE
            } else{
              active.ind <- c(active.ind, kkt_results$new_groups)
              inactive.ind<- (1:ng)[-active.ind]
              #break
            }
          }
        }
      }
      else{ # nonsparse estimates, length(active_ind)=ng
        
        outer_iter<- 0
        kkt_met<- NA
        max_iter<- 50
        
        # update beta and residuals
        update <- .C("solve_beta", as.double(y), as.double(cbind(1,x)), as.double(tau), as.double(gamma), 
                     as.double(weights), as.double(lambda[j]), as.double(w.lambda), 
                     as.double(eigen.sub.H), as.double(beta0), as.integer(max_iter), as.double(epsilon), as.integer(apprx=="huber"), 
                     as.integer(n), as.integer(p), as.integer(ng), as.integer(nng), as.integer(0), as.double(r0), as.integer(0))
        beta.update <- update[[9]]
        update.r <- update[[18]]
        update.converge <- update[[17]]
        update.iter <- update[[19]]
        
        beta0<- beta.update
        
      }
      
      outer_iter_count[j] <- outer_iter
      kkt_seq[j]<- kkt_met
      converge[j]<- update.converge
      iter[j]<- update.iter
      gamma.seq[j]<- gamma
      
      beta0[beta_order] <- beta0  # required for C
      
      beta.all[,j]<- beta0
      r0<- update.r 
      if(method=="quantile"){
        gamma.max<- 4; gamma<- min(gamma.max, max(gamma0, quantile(abs(r0), probs = 0.1)))
      } 
      
      # group index and number of groups
      grp.ind<- unique(group.index[beta0[-1]!=0])
      group.index.out[grp.ind,j]<- grp.ind
      n.grp[j]<- length(grp.ind)
      
      # alternative deviance
      if(method=="quantile"){
        dev1<- sum(weights*rq.loss(r0, tau))
      } else if( method == "mean"){
        dev1<- sum(weights*huber.loss(r0,gamma))	
      }
      loss[j]<- dev1/n
      pen.loss[j]<- dev1/n+lambda[j]*sum(eigen.sub.H*sapply(1:ng, function(xx) l2norm(beta0[-1][group.index==xx])))
      rel_dev[j]<- dev1/dev0
      rel_dev_change<- rel_dev[j]-rel_dev[j-1]
      if(abs(rel_dev_change)<1e-3 & j>70 & lambda.discard) break

    } # end of for loop of lambda
    
    
    stop.ind<- which(rel_dev!=0)
    if(length(stop.ind)==0) stop.ind<- 1:length(lambda)
    
    if(lambda.flag==0){
      stop.ind<- stop.ind[-1]
      if(scalex){
        beta.final<- apply(beta.all[,stop.ind], 2, transform_coefs, mu_x,sigma_x)
      } else{
        beta.final <- beta.all[,stop.ind]
      }
      rownames(beta.final)<- c("Intercept", paste("V", 1:p, sep = ""))
      
      output<- list(beta=Matrix(beta.final), lambda=lambda[stop.ind], null.dev=dev0, pen.loss=pen.loss[stop.ind], 
                    loss=loss[stop.ind], n.grp=n.grp[stop.ind], index.grp=Matrix(group.index.out[,stop.ind]),
                    tau=tau, apprx=apprx, group.index=group.index, method=method, x=x, y=y) 
     
      output_hide<- list(iter=iter[stop.ind], rel_dev=rel_dev[stop.ind], outer.iter=outer_iter_count[stop.ind], 
                         kkt=kkt_seq[stop.ind], gamma=gamma.seq[stop.ind])
      #output1<- c(output, output_hide)
      class(output) <- "hrq_glasso"
      return(output)
    }
    
    ####################
    if(lambda.flag==1){
      length.diff<- length(lambda.user)-length(lambda)+1
      beta.all<- cbind(matrix(beta.all[,1], nrow = p+1, ncol = length.diff),beta.all[,-1])
      if(scalex){
        beta.final<- apply(beta.all, 2, transform_coefs, mu_x,sigma_x)
      } else{
        beta.final <- beta.all
      }
      rownames(beta.final)<- c("Intercept", paste("V", 1:p, sep = ""))
      
      output<- list(beta=Matrix(beta.final), lambda=lambda.user, null.dev=dev0, pen.loss=c(rep(pen.loss[1], length.diff), pen.loss[-1]), 
                    loss=c(rep(loss[1], length.diff), loss[-1]), n.grp=c(rep(n.grp[1], length.diff), n.grp[-1]), 
                    index.grp=Matrix(cbind(matrix(group.index.out[,1], nrow = nrow(group.index.out), ncol = length.diff),group.index.out[,-1])), 
                    tau=tau, apprx=apprx, group.index=group.index, method=method, x=x, y=y)
      output_hide<- list(iter=c(rep(iter[1], length.diff), iter), rel_dev=c(rep(rel_dev[1], length.diff), rel_dev), outer.iter=outer_iter_count[stop.ind], 
                         kkt=c(rep(kkt_seq[1], length.diff), kkt_seq[-1]), gamma=c(rep(gamma.seq[1], length.diff), gamma.seq[-1]))
      
      class(output) <- "hrq_glasso"
      warning(paste("first ", length.diff, " lambdas results in pure sparse estimates!", sep = ""))
      return(output)
      
    }
    
  } # end of else condition
  
} # end of function
########################################





