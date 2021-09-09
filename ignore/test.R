rm(list=ls(all=TRUE))
library(devtools)
devtools::unload("rqPen")
devtools::unload("hrqglas")
install_github("bssherwood/hrqglas")
3
#install_github("shaobo-li/rqPen")
install_github("bssherwood/rqPen")
3
library(rqPen)
library(hrqglas)
library(quantreg)
data(barro)
y <- barro$y.net
x <- as.matrix(barro[,-1])

groups <- c(rep(1,4),rep(2,3),rep(3,3),rep(4,3))


#r1 <- rq.pen(x,y)
lamMax <- rqPen:::getLamMaxGroup(x, y, groups, .5, rep(1,4), 
                         penalty = "gLASSO", scalex=TRUE,tau.penalty.factor=1)
eps <- ifelse(nrow(x) < ncol(x), 0.01, 1e-04)
nlambda <- 100
lambda <- exp(seq(log(lamMax), log(eps * lamMax), length.out = nlambda))

h0 <- hrq_glasso(x,y,groups,tau=.5, w.lambda=rep(1,4))#weird cutting off first lambda and also doesn't quite match with these next two
h1 <- hrq_glasso(x,y,groups,tau=.5,lambda=lambda, w.lambda=rep(1,4))
r1 <- rq.group.pen(x,y,groups=groups)


getLamMaxGroup <- function(x,y,group.index,tau=.5,group.pen.factor,gamma=.2,gamma.max=4,gamma.q=.1,penalty="gLASSO",scalex=TRUE){
  # code improvement: Hacky approach to the group.pen.factor issue. 
  returnVal <- 0
  n <- length(y)
  if(scalex){
    x <- scale(x)
  }
  validSpots <- which(group.pen.factor!=0)
  for(tau_val in tau){
    r <- y - quantile(y,tau_val)
    gamma0<- min(gamma.max, max(gamma, quantile(abs(r), probs = gamma.q)))
    
    grad_k<- -neg.gradient(r, rep(1,n), tau_val, gamma=gamma0, x, apprx="huber")
    grad_k.norm<- tapply(grad_k, group.index, l2norm)
    
    lambda.max<- max(c(returnVal,grad_k.norm[validSpots]/group.pen.factor[validSpots]))
  }
  lambda.max
}

n <- length(y)
tau <- .5
gamma.max <- 4
gamma0 <- .2
weights <- rep(1,n)
group.index <- groups
w.lambda <- rep(1,4)
x <- scale(x)

b.int<- quantile(y, probs = tau)
r<- y-b.int
gamma.max<- 4; 
gamma<- min(gamma.max, max(gamma0, quantile(abs(r), probs = 0.1)))

r0 <- r
apprx <- "huber"
grad_k2 <- -hrqglas:::neg.gradient(r0, weights, tau, gamma=gamma, x, n, apprx)
grad_k.norm<- tapply(grad_k2, group.index, hrqglas:::l2norm)

lambda.max<- max(grad_k.norm/w.lambda)


#figure out the difference in the lambda maximum values. 
fit <- hrq_glasso(x,y,groups,tau=.5)

#fit <- hrq_glasso(x,y,groups,tau=.5,lambda=r1$lambda)
x <- as.matrix(x)
r2 <- rq.group.pen(x,y,groups=groups)

weights <- rep(1,length(y))
w.lambda <- rep(1,4)

b.int <- quantile(y, probs = tau)
grad_k <- -neg.gradient(r0, weights, tau, gamma = gamma, 
                        x, n, apprx)
grad_k.norm <- tapply(grad_k, group.index, l2norm)
lambda.max <- max(grad_k.norm/w.lambda)
r <- y - b.int