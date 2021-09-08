rm(list=ls(all=TRUE))
library(devtools)
devtools::unload("rqPen")
devtools::unload("hrqglas")
install_github("shaobo-li/hrqglas")
install_github("bssherwood/rqpen")
library(rqPen)
library(hrqglas)
library(quantreg)
data(barro)
y <- barro$y.net
x <- barro[,-1]

groups <- c(rep(1,4),rep(2,3),rep(3,3),rep(4,3))


#r1 <- rq.pen(x,y)
lamMax <- rqPen:::getLamMaxGroup(x, y, groups, .5, rep(1,4), 
                         penalty = "gLASSO", TRUE)
eps <- ifelse(nrow(x) < ncol(x), 0.01, 1e-04)
nlambda <- 100
lambda <- exp(seq(log(lamMax), log(eps * lamMax), length.out = nlambda))
hrq_glasso(x,y,groups,tau=.5,lambda=lambda)
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