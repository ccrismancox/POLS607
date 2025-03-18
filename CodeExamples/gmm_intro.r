library(MASS)
library(fixest)
library(Matrix)
library(car)
library(ivreg)
library(sandwich)
library(gmm)
library(dplyr)
rm(list=ls())
N <- 25
T <- 10
set.seed(1)
errors <- mvrnorm(N*T, c(0,0), 
                  Sigma=matrix(c(4, 2*4*.6,
                                 2*4*.6, 16),
                               nrow=2, ncol=2))

## unobserved heterogeneity
U <- mvrnorm(N, c(-1,1), 
        Sigma=matrix(c(1, .9,
                       .9, 1),
                     nrow=2, ncol=2))
kappa <- rep(U[,1], each=T)
alpha <- rep(U[,2], each=T)

## instruments (exogenous only conditional on the unobserved heterogeneity)
z1 <- rnorm(N*T) + kappa
z2 <- rnorm(N*T) + kappa

Z <- cbind(z1, z2)
gamma <- c(1, 3)

## endogenous x
x1 <- drop(Z %*% gamma + kappa + errors[,1]) 

## exogenous x
x2 <- rnorm(N*T)

## betas
beta <- c(-1, 1)
X <- cbind(x1, x2)
y <- drop(X %*% beta + alpha + errors[,2])

cor(cbind(X,  alpha,errors[,2]))

df <- data.frame(state=rep(1:N, each=T),
                 year=1:T,
                 y=y,
                 x1=x1, x2=x2, z1=z1, z2=z2)

within2sls <- feols(y~x2|state|x1~z1+z2, data=df)
summary(within2sls)



df$u <- within2sls$residuals
aux <- feols(u~z1+z2+x2|state, data=df)

## overidentified Sargan test
S <- N*T*(1- crossprod(aux$residuals)/crossprod(within2sls$residuals))
linearHypothesis(lm(x1~z1+z2+x2+factor(state)-1, data=df),
                 c("z1=0", "z2=0"),
                 vcov=\(x){vcovCL(x,df$state)})

var.names <- c("y", "x1", "x2", "z1", "z2")
df <- df %>% 
  mutate(across(all_of(var.names), 
                \(x){x-mean(x)}, 
                .names = "{col}.within"),
         .by=state)
    
Zdot <- with(df, cbind(z1.within, z2.within,x2.within))
Xdot <- with(df, cbind(x1.within, x2.within))
ydot <- df$y.within
eps.hat <- within2sls$residuals



## 2step gmm 
### 1. start with W from another model (e.g. 2sls or a different gmm)
W <- solve(Reduce( `+`, by(df$state, data=Zdot*eps.hat,
                           \(x){ 
                             tcrossprod(colSums(x))/T
                           }))/N)
### 2. Fit the GMM with this W
beta.hat.gmm <- solve(t(Xdot) %*% Zdot %*% W %*%t(Zdot) %*% Xdot) %*% 
  (t(Xdot) %*% Zdot %*% W %*% t(Zdot) %*% ydot)
### 3. Compute a new W with the GMM output
eps.hat1 <- drop(ydot - Xdot %*% beta.hat.gmm)
W1 <- solve(Reduce( `+`, by(df$state, data=Zdot*eps.hat1,
                            \(x){ 
                              tcrossprod(colSums(x))/T
                            }))/N)
## 4. Wuse the updated W for the standard errors
se.gmm <- sqrt(diag(N*T*solve(t(Xdot) %*% Zdot %*% (W1) %*% t(Zdot) %*% Xdot)))


cbind(beta.hat.gmm, se.gmm, beta.hat.gmm/se.gmm)


## Canned version. Need to supply our own weights 
within.gmm2 <- gmm(y.within~x1.within+x2.within-1, 
                   ~z1.within+z2.within+x2.within-1,
                   data=df, 
                   weightsMatrix = W)
summary(within.gmm2)
eps.hat1 <- drop(ydot - Xdot %*% within.gmm2$coefficients)
W1 <- solve(Reduce( `+`, by(df$state, data=Zdot*eps.hat1,
                            \(x){ 
                              tcrossprod(colSums(x))/T
                            }))/N)





## Sargan
G <- colMeans(Zdot*eps.hat1)
J <-  N*T*(G %*% W1 %*% G)
J
pchisq(J, lower=FALSE, df=ncol(Zdot)-ncol(Xdot))


## What is being minimized? This objective
g <- function(b, y=y, X=X, Z=Z, W=W){
  ehat <- (y-X%*%b)
  return(t(ehat) %*% Z  %*% W %*% t(Z) %*% ehat/(nrow(Z)^2))
}

Dg <- function(b, y=y, X=X, Z=Z, W=W){
  ehat <- (y-X%*%b)
  return( -2* t(ehat) %*% Z  %*% W  %*% t(Z) %*% X/(nrow(Z)^2))
} 

gmm.brute <- optim(within.gmm2$coef, g, gr=Dg,
                   y=ydot,X=Xdot, Z=Zdot, W=W,
                   method="BFGS")
#Sargan
c(J, gmm.brute$value *(N*T))

