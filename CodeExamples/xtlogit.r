library(matrixStats)
library(fastGHQuad)
library(numDeriv)
library(lme4)
library(dplyr)
## log likelihood for RE logit 
relogit <- function(theta, X, y, id, GH, byid=FALSE){
  s <- exp(theta[length(theta)])
  beta <- theta[-length(theta)]
  XB <- drop( X%*%beta )
  sepXB <- split(XB, factor(id))
  sepY <- split(y, factor(id))
  
  int <- mapply(x=sepXB, y=sepY,
                \(x,y){
                  f <- \(u){return(logit.f(u=u, XB=x, y=y, s=s))}
                  sum(GH$w * f(GH$x))
                })
  if(byid){
    return(-log(int))
  }else{
    return(-sum(log(int)))
  }
}

## Auxillary function for the likelihood to integrate
logit.f <- function(u, XB, y, s){
  u <- matrix(rep(u, each=length(y)), nrow=length(y))
  return(1/sqrt(pi) *colProds(plogis((2*y-1)*(XB+u*s*sqrt(2)))))
}

felogit <- function(beta, X, y, id, byid=FALSE){
  XB <- drop( X%*%beta )
  eXB <- exp(XB)
  k <- lapply(split(y, factor(id)), sum)
  sepXB <- split(eXB, factor(id))
  sepT <-  lapply(sepXB, length)
  lf <- mapply(\(t, k, e){log(fi(t,k,e))},
    t=sepT,
    k=k,
    e=sepXB)
  s1 <- sapply(split(y*XB, factor(id)), sum)
  ll <- s1-lf
  if(byid){
    return(-ll)
  }else{
    return(-sum(ll))
  }
}
## Auxillary function for clogit 
fi <- function(T,k, eXB){
  if(T<k){
    return(0)
  }else{if(k==0){
    return(1)
  }else{
    return(fi(T-1, k, eXB) + fi(T-1, k-1, eXB)*eXB[T] )     
  }
  }
}




## Clustered standard errors ##
vcovCL.relogit <- function(mod){
  par <- c(fixef(mod), log(attributes(VarCorr(mod)$id)$stddev))
  df <- lme4::glFormula(mod@call$formula,data = mod@frame)
  X <- df$X
  id <- df$reTrms$flist$id
  y <- df$fr[,1]
  GH <- gaussHermiteData(cre@call$nAGQ)
  Ji <- jacobian(f=relogit, x=par,  
                 y=y, X=X, 
                 id=id,
                 GH=GH,
                 byid=TRUE)
  Meat <- crossprod(Ji)         
  bread <- solve(mod@optinfo$derivs$Hessian)
  V <- bread %*% Meat %*% bread 
  dBeta <- diag(c(rep(1, length(par)-1), exp(par[length(par)])))
  Vcl <- dBeta %*% V %*% dBeta
  rownames(Vcl) <- colnames(Vcl) <- c(names(fixef(mod)), "sigma")
  return(Vcl)  
}

vcovCL.felogit <- function(mod, y, X, id){
  par <- mod$coefficients
  Ji <- jacobian(f=felogit, x=par,  
                 y=y, X=X, 
                 id=id,
                 byid=TRUE)
  Meat <- crossprod(Ji)         
  bread <- mod$var
  V <- bread %*% Meat %*% bread 
  rownames(V) <- colnames(V) <- names(mod$coefficients)
  return(V)  
}



# Uncomment this to test
# ### test data 
# N <- 100
# T <- 2
# 
# a <- runif(N, -2, 2)
# X <- runif(N*T) + .2*rep(a,each=T)
# beta <- .5
# 
# true.ame <- mean(dlogis(X*beta + rep(a,each=T))*beta)
# dat <- data.frame(id= rep(1:N, each=T),
#                   T= rep(1:T, N),
#                   X=X)%>%
#   mutate(Xbar=mean(X), .by=id)
# dat$y <- 1*(X*beta + rep(a,each=T) + rlogis(N*T)>0)
# 
# 
# 
# 
# ## for comparison
# mundlak <- glmer(y~X+Xbar+(1|id), data=dat, family=binomial(), nAGQ = 6)
# 
# 
# GH <- gaussHermiteData(6)
# theta0 <- rep(0,4)
# names(theta0) <- c("Intercept", "X", "Xbar", "sigma")
# cre2 <- optim(theta0,
#               relogit,
#               y=dat$y, X=cbind(1,dat$X, dat$Xbar), 
#               id=dat$id, GH=GH,
#               hessian=TRUE)
# c(fixef(mundlak), attributes(VarCorr(mundlak)$id)$stddev)
# c(cre2$par[1:3], exp(cre2$par[4]))
# 
# 
# ## Clustered standard errors ##
# Ji <- jacobian(f=relogit, x=cre2$par,  
#          y=dat$y, X=cbind(1,dat$X, dat$Xbar), 
#          id=dat$id, GH=GH,
#          byid=TRUE)
# Meat <- crossprod(Ji)         
# bread <- solve(cre2$hessian)
# Vcl <- bread %*% Meat %*% bread 
# c(sqrt(diag(Vcl))[1:3],
#   exp(cre2$par[4])*sqrt(diag(Vcl))[4])





