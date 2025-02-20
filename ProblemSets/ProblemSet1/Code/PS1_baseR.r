rm(list=ls())

set.seed(1)
source("panelFunctions.r")

## set up simulation conditions and parameters 
conds <- expand.grid(Units=c(50, 200),
                     het=c(1,0))
burnin <- 25
rho.x <- .25
rho.e <- .8
theta <- c(1, 4, -3, 5)

sd.u <- sqrt(7)
sd.x <- 1/2
sd.z2 <- 1/2
a.high <- 1


## How different would we expect the FE and RE results to be at an average T
sigma2.e.true <- sd.u^2/(1-rho.e^2)
sigma2.a.true <- (2*a.high)^2/12+theta[4]^2 * sd.z2^2
omega.T25 <- 1 - sqrt(sigma2.e.true/(25*sigma2.a.true + sigma2.e.true))

## initialize storage for the results
Results <- list()
conds$obs.cor  <-  0
conds$pooled.OVB  <-  0
for(i in 1:nrow(conds)){
  
  
  ## Set the conditions for this experiment (NT)
  N <- conds$Units[i]
  Ti <- rpois(N, 24)+1
  T <- max(Ti) #We'll create it balanced and then trim
  
  
  ## draw fixed effects and exogenous data
  a <- runif(N, -1, 1)
  z2 <- rnorm(N, sd=.5)
  z1 <- rnorm(N, sd=2)
  
  x2 <- rpois(N*T, 5)
  
  if(conds$het[i]==1){
    ## If we're in the main set of simulations
    ## we'll create correlation matrices for our stationary 
    ## x1 and epsilon so we don't have to have an inner sequential loop
    ## over time periods (nothing wrong with looping over time periods, 
    ## but looping over **both** N and T gets slow)
    T2 <- T+ burnin
    Rx <- sapply(T2:1, \(x){c(rep(0,T2-x),rho.x^(0:(x-1)))})
    R.e <-sapply(T2:1, \(x){c(rep(0,T2-x),rho.e^(0:(x-1)))})
    
    ## static component of x1
    mu.x <- matrix(15 -3*rep(a,each=T2)- 2*rep(z2,each=T2), ncol=N)
    ## initial condition plus, discounted by how long ago it was
    x0 <- mu.x * rho.x^(1:T2)
    ## x.star shocks 
    x.star <- matrix(rnorm(N*T2,sd=sd.x), nrow=T2)
    ## the correlation matrix times static term plus shocks 
    ## plus the initial condition discounted by time
    x1 <- Rx %*% (mu.x  + x.star) + x0
    
    ## discard the burn-in
    x1 <- c(x1[(burnin+1):T2,])
  }else{
    ## or exogenous for the other simulations
    x1 <- rnorm(N*T, 15, sd=sd.x)
  }
  
  ### fix the data for each simulation
  dat <- data.frame(unit=rep(1:N, each=T),
                    year=1:T,
                    z1 = rep(z1, each=T),
                    z2 = rep(z2, each=T),
                    x1 = x1,
                    x2 = x2,
                    a = rep(a, each=T))

  ## discard unused years for the unbalanced panel
  dat <- subset(dat, year<=rep(Ti,each=T))
  dat$x1.bar <- with(dat, ave(x1, unit))
  dat$x2.bar <- with(dat, ave(x2, unit))
    

  ## full observed X
  bigX <- with(dat, cbind(1, x1, x2, z1))
  ## full observed X group means
  Xbar <- with(dat, cbind(1, x1.bar, x2.bar, z1))
  
  ## The true omega weights under these conditions
  omega.true <- 1 - sqrt(sigma2.e.true/(rep(Ti,Ti)*sigma2.a.true + sigma2.e.true))
  
  
  ## Omitted variable bias calculations
  Omitted <- with(dat, cbind(a,z2))
  OV.bias <- solve(crossprod(bigX)) %*% t(bigX) %*% Omitted %*% c(1,theta[4])
  OV.bias <- OV.bias[2]

  conds$obs.cor[i] <- with(dat, cor(a+theta[4]*z2, x1))
  conds$pooled.OVB[i] <- OV.bias

  
  ## Start the simulation
  out <- matrix(0, nrow=1000, ncol=14)
  for(b in 1:1000){
    
    
    
    
    
    ## Create error terms based on simulation
    if(conds$het[i]==1){
      u <- matrix(rnorm(N*T2, mean=0, sd=sd.u), nrow=T2)
      e <- R.e %*% u
      dat$e <- c(e[(burnin+1):T2,])[rep(1:T,N) <= rep(Ti,each=T)]
    }else{
      dat$e <- rnorm(sum(Ti), 0, sd=sd.u)
    }
    
    ## Generate y and bar{y}
    dat$y <- with(dat, a +  cbind(x1, x2, z1, z2) %*% theta + e)
    dat$y.bar <- ave(dat$y, dat$unit)
    

    ## pooled model
    pooled <- lm(y~x1+x2+z1,data=dat)
    pooled.se <- c(pooled.se0=sqrt(diag(vcov(pooled)))[2],
                   pooled.se1=sqrt(diag(clusterVCOV(bigX, pooled$residuals, dat$unit)))[2])
    
    ## within transformation
    X.within <- with(dat, cbind(x1-x1.bar, x2-x2.bar))
    y.within <- dat$y-dat$y.bar
    theta.within <- solve(crossprod(X.within)) %*% t(X.within) %*% y.within
    e.within <- y.within - X.within %*% theta.within
    
    s2.e.hat <- drop(crossprod(e.within)/length(e.within))
    V.within <- s2.e.hat * (solve(crossprod(X.within)))
    within.se <- c(within.se0=sqrt(diag(V.within))[1],
                   within.se1=sqrt(diag(clusterVCOV(X.within, e.within, dat$unit)))[1])
    
    
    ## GLS: step1 variance
    ehat <- pooled$residuals
    combo.var <- crossprod(ehat)/length(ehat)
    s2.a.hat <- drop(combo.var - s2.e.hat)
    ## if you wanted to do the covariances approach
    # s2.a.hat <- mean((by(ehat, INDICES = dat$unit,
    #                      \(x){ 
    #                        M<-tcrossprod(x);
    #                        return( mean(M[lower.tri(M)]))})), na.rm=TRUE)
    # s2.e.hat <- drop(combo.var - s2.a.hat)
    
    ## weights and transformations
    w.hat <-  1 - sqrt(s2.e.hat)/sqrt(Ti * s2.a.hat + s2.e.hat)
    w.hat <- rep(w.hat, Ti)
    y.tilde <- (dat$y-w.hat*dat$y.bar)/sqrt(s2.e.hat)
    X.tilde <- (bigX - w.hat *Xbar)/sqrt(s2.e.hat)
    
    ## GLS 
    theta.gls <- solve(crossprod(X.tilde)) %*% t(X.tilde) %*% y.tilde
    e.gls <- y.tilde - X.tilde %*% theta.gls
    s2.gls <- drop(crossprod(e.gls)/length(e.gls))
    V.gls <- solve(crossprod(X.tilde))
    gls.se <- c(gls.se0=sqrt(diag(V.gls))[2],
                gls.se1=sqrt(diag(clusterVCOV(X.tilde, e.gls, dat$unit)))[2])
    
    
 
    ## RE-MLE
    mle.re <- optim(c(theta.gls, log(c(s2.a.hat, s2.e.hat))), #start values
                    fn=relik,  ##loglikelihood
                    gr=re_mle_gr, ##gradient
                    method="BFGS", ##optimization method
                    X=bigX, ## data
                    Ti=Ti,  Xbar=Xbar, ybar=dat$y.bar,
                    y=dat$y)
    
    ## untransform the variance estimates 
    s2a.mle <- exp(mle.re$par[5])
    s2e.mle <- exp(mle.re$par[6])
    resid.re <- dat$y - bigX %*% mle.re$par[1:4]
    omega.mle <- rep(1- sqrt(s2e.mle/(Ti*s2a.mle + s2e.mle)),Ti)
    resid.bar <- rep(tapply(resid.re, dat$unit, mean),Ti)
    X.mle <- (bigX- Xbar*omega.mle)/sqrt(s2e.mle)
    resid.mle <- (resid.re-omega.mle*resid.bar)/sqrt(s2e.mle)
    
    mle.re.se <- c(mle.se0=sqrt(diag(solve(crossprod(X.mle))))[2],
                   mle.se1=sqrt(diag(clusterVCOV(X.mle,resid.mle,dat$unit)))[2])
    
    
    ## Hausman test
    Hstat <- t(theta.gls[2:3]-theta.within) %*% 
      solve(V.within -V.gls[2:3, 2:3]) %*% 
      (theta.gls[2:3]-theta.within)
    Htest <- c(pchisq(Hstat, lower=F, df=2)  < 0.05)
    
    ## Mundlak
    mundlak <- lm(y~x1+x2+ z1+x1.bar + x2.bar, data=dat, x=TRUE)
    mundlak.V <- clusterVCOV(mundlak$x,mundlak$residuals,dat$unit)
    
    A <- cbind(matrix(0, 2,4), diag(2))
    wald.stats <- c(t(A %*% mundlak$coefficients) %*% solve(A %*% mundlak.V %*% t(A)) %*%  (A %*% mundlak$coefficients))
    
    Wald <- c(pchisq(wald.stats, lower=FALSE, df=2)  < 0.05)
    
    
    out[b, ]<- c(bhat.pooled=pooled$coef[2], 
      bhat.within=theta.within[1],
      bhat.mle=mle.re$par[2],
      bhat.gls=theta.gls[2], 
      pooled.se,
      within.se,
      gls.se,
      mle.re.se,
      Hausman=Htest,
      Wald=Wald)
  }
  Results[[i]] <- out
}

save(list=c("conds", "Results", "theta"), file="PS1.rdata")

