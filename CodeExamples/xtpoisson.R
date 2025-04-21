library(numDeriv)
library(lme4)
library(dplyr)
source("panelFunctions.r")

logsumexp <- function(x){
  xmax <- max(x)
  x2 <- x-max(x)
  return(xmax + log(sum(exp(x2))))
}

xtpoisson <- function(beta, y, X, id){
  XB <- drop(X %*% beta)
  ll <- y*(XB-ave(XB,factor(id), FUN=logsumexp))
  ll <- split(ll, factor(id))
  ll <- sapply(ll, sum)
  return(-sum(ll))
}

D.xtpoisson <- function(beta, y, X, id){
  XB <- drop(X %*% beta)
  eXB <- exp(XB)
  denom <- ave(eXB,factor(id), FUN=sum)
  numer <- do.call(rbind, 
                   lapply(split.matrix(X*eXB,factor(id)),
                          \(x){t(ave(x, FUN=colSums))}))
  gr <- y*(X -numer/denom)
  return(-colSums(gr))
}




## log likelihood for RE poisson with gamma heterogeneity 
repoisson <- function(theta, X, y, id, byid=FALSE){
  s <- exp(theta[length(theta)])
  beta <- theta[-length(theta)]
  XB <- drop(X%*%beta)
  lam <-  exp(XB)
  yXB <- y*XB
  
  mat <- cbind(y, lam, lgamma(y+1), yXB)
  MAT <- sapply(split.matrix(mat, factor(id)), colSums)
  sumY <- MAT[1,]
  sumLam <- MAT[2,]
  sumlgY <- MAT[3,]
  sumXBy <- MAT[4,]
  
  nu <- s/(s+sumLam)
  
  
  ll <- lgamma(sumY+s) - 
     sumlgY-
    lgamma(s) + 
    s*log(nu) + 
    log(1-nu)*sumY + 
    sumXBy - 
    sumY*log(sumLam)
  if(byid){
    return(-ll)
  }else{
    return(-sum(ll))
  }
}
# 
# ### test data 
# N <- 100
# T <- 2
# 
# 
# alpha <- floor(runif(N, -5, 5))
# tau <- floor(runif(T, -1, 3))
# 
# a <- rep(alpha, each=T) 
# t <- rep(tau, N)
# 
# X <- rnorm(N*T, mean=2, sd=.5) + .25*a + .15*t
# lambda <- exp(X* (.5) + a+t)
# 
# beta <- .5
# 
# dat <- data.frame(id= rep(1:N, each=T),
#                   time= rep(1:T, N),
#                   X=X,
#                   y= rpois(N*T, lambda=lambda))%>%
#   mutate(Xbar=mean(X), .by=id)
# 
# 
# 
# 
# ## for comparison
# mundlak <- glmer(y~X+Xbar+I(time==2)+(1|id), data=dat, family=poisson, nAGQ = 6)
# 
# theta0 <- runif(5)
# names(theta0) <- c("Intercept", "X", "Xbar", "tau2", "scale")
# cre2 <- optim(theta0,
#               repoisson,
#               y=dat$y,
#               X=cbind(1,dat$X, dat$Xbar, 1*(dat$time==2)), 
#               id=dat$id, 
#               method="BFGS",
#               hessian=TRUE)
# c(fixef(mundlak), attributes(VarCorr(mundlak)$id)$stddev)
# c(cre2$par[1:4], 1/exp(cre2$par[5]))
# 
# 
# ## Clustered standard errors ##
# Ji <- jacobian(f=repoisson, x=cre2$par,  
#                y=dat$y, X=cbind(1,dat$X, dat$Xbar, 1*(dat$time==2)), 
#                id=dat$id, 
#                byid=TRUE)
# Meat <- crossprod(Ji)         
# bread <- solve(cre2$hessian)
# Vcl <- bread %*% Meat %*% bread 
# sqrt(diag(Vcl))
# 
# 
