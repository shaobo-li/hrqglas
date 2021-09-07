library(rqPen)
library(hrqglas)
library(quantreg)
data(barro)
y <- barro$y.net
x <- barro[,-1]

groups <- c(rep(1,4),rep(2,3),rep(3,3),rep(4,3))


r1 <- rq.pen(x,y)


fit <- hrq_glasso(x,y,groups,tau=.5,lambda=r1$lambda)
r2 <- rq.group.pen(x,y,groups=groups)