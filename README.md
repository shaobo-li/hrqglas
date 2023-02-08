# hrqglas: An R Package for Group Variable Selection for Quantile and Robust Mean Regression

## Overview

This R package provides a program to conduct group-wise variable selection and estimation for quantile and robust mean regression. The group lasso penalty $p_{\lambda_j}(\beta_j)=\lambda_j\|\beta_j\|_2$, where $j$ is the group index, is used for group-wise variable selection (Yuan and Lin, 2006). For quantile regression, the check loss is approximated by the Huber loss for the median and the tilted version of Huber loss at other quantiles. This approximation overcomes the nondifferentiability at the origin of check loss, which may otherwise cause instable estimation. Statistical consistency has been shown for this approximated quantile regression estimates (Sherwood and Li, 2022). The estimation algorithm follows Yang and Zou (2015), and it is computationally efficient and stable. A robust estimation of mean regression is a byproduct of this implementation as Huber loss, with appropriate choices of the tuning parameter, is intrinsically a robust loss function that is insensitive to outliers. 

## Installation

``` r
devtools::install_github("shaobo-li/hrqglas")
```

## Example

``` r
library(hrqglas)
n<- 200
p<- 30
x0<- matrix(rnorm(n*p),n,p)
X<- cbind(x0, x0^2, x0^3)[,order(rep(1:p,3))]
y<- -2+X[,1]+0.5*X[,2]-X[,3]-0.5*X[,7]+X[,8]-0.2*X[,9]+rt(n,2)
group<- rep(1:p, each=3)

# quantile regression
fit<- hrq_glasso(x=X, y=y, group.index=group, method="quantile", tau=0.3)
fit.cv<- cv.hrq_glasso(x=X, y=y, group.index=group, method="quantile", tau=0.3, loss="check")
plot(fit.cv)

# mean regression
fit1<- hrq_glasso(x=X, y=y, group.index=group, method="mean")
fit.cv1<- cv.hrq_glasso(x=X, y=y, group.index=group, method="mean", loss="se")
plot(fit.cv1)
```



## References

Sherwood, B., & Li, S. (2022). Quantile regression feature selection and estimation with grouped variables using Huber approximation. *Statistics and Computing*, 32(5), 75. [https://doi.org/10.1007/s11222-022-10135-w](https://doi.org/10.1007/s11222-022-10135-w)
 
Yang, Y. and Zou, H., (2015) A Fast Unified Algorithm for Solving Group-lasso Penalize Learning Problems, *Statistics and Computing*, 25 1129-1141. [https://doi.org/10.1007/s11222-014-9498-5](https://doi.org/10.1007/s11222-014-9498-5).

Yuan, M. and Lin, Y., (2005) Model Selection and Estimation in Regression with Grouped Variables, \emph{Journal of the Royal Statistical Society: Series B}, 68 49-67.  [https://doi.org/10.1111/j.1467-9868.2005.00532.x](https://doi.org/10.1111/j.1467-9868.2005.00532.x).
