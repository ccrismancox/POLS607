## data manipulation packages
library(readstata13)
library(data.table) 

## econometrics packages
library(lmtest)
library(car)
library(sandwich)
library(fixest) 
library(lme4) 
library(clubSandwich)

## tables and figures
library(modelsummary)

## checking out the data
protests <- read.dta13("datasets/Replication_secpol_protestComplete.dta")
protests <- data.table(protests)
protests <- protests[order(ccode, year),]

colnames(protests)

## panel dimenions
length(unique(protests$ccode))
summary(protests[, length(year), by = ccode])

## adjust the variables based on their replication file

## Normalize the latent varaible to be  mean 0, var 1 
protests[, Protest := scale(mean5)] 
protests[, nbr_protest := scale(nbr_mean5)] 

## create the controls: lag(log(pop)), lag(log(gdp_pc), lag(excluded population))
protests[, `:=` (l.ln_pop = shift(log(pop+1)), 
                 l.ln_gdppc = shift(log(gdp_pc)),
                 l.lexclpop = shift(lexclpop)),
         by=ccode]

## model formula
f1 <- Protest~ secretpol_revised + l.ln_pop + l.ln_gdppc  + l12gr+  l.lexclpop+
  nbr_protest+intrastate+attempt


## Fitting with the pooled esetimator
pooled  <- lm(f1, data=protests, x=TRUE)
summary(pooled)

## To make life easy
## We're going to restrict ourselves to just the used sample
protests <- protests[as.numeric(row.names(pooled$model)), ]

## Let's consider the residual autocorrelation
## in choosing standard errors
protests[, e.hat := pooled$residuals]
protests[, L.e.hat := shift(e.hat), by =ccode]
summary(lm(e.hat~L.e.hat, data=protests)) #that's pretty high!


## Clustering the standard errors 
Vcl.pooled <- vcovCL(pooled, cluster=protests$ccode)
round(coeftest(pooled, Vcl.pooled), 5)

## Suppose we wanted to bootstrap we have the clustered bootstrap
pooled.boot <- t(replicate(50, {
  idx <- sample(unique(protests$ccode), 
                size=length(unique(protests$ccode)),
                replace=TRUE)
  d <- copy(protests)
  d <- d[unlist(sapply(idx, \(x){which(d$ccode==x)}))]
  pooled.bs <- lm(f1, dat=d)
  pooled.bs$coef
}))
round(coeftest(pooled, var(pooled.boot)), 5)

## And the clustered bayesian bootstrap
pooled.bayes.boot <- t(replicate(50, {
  Ti <- table(protests$ccode)
  d <- copy(protests)
  weight <- rexp(length(unique(protests$ccode)))
  weight <- weight/sum(weight)
  d$weight <- rep(weight*length(unique(protests$ccode)), Ti)
  lm(f1, dat=d, weights=weight)$coef
}))
round(coeftest(pooled, var(pooled.bayes.boot)), 5)


## We can consider fixed effects estimators too. Starting with the LSDV
lsdv  <- lm(update(f1, .~. -1 + factor(ccode)), data=protests)
Vcl.lsdv <- vcovCL(lsdv, cluster=protests$ccode)
round(coeftest(lsdv, Vcl.lsdv)[1:8,], 4)


## Within transformation
var.names <- colnames(pooled$model)
protests[,paste0(var.names, ".within"):=lapply(.SD, \(x){x- mean(x)}), 
         by=ccode, .SDcols=var.names ]
fwithin <- paste0(var.names[1], ".within ~ -1 + ", 
                  paste0(var.names[-1], ".within", collapse=" + "))
within1 <- lm(fwithin, data=protests)
Vcl.within1 <- vcovCL(within1, cluster=protests$ccode)
round(coeftest(within1, Vcl.within1), 4)


## The fixest is the better way to go here.  It takes 
## a formula of the form y~x|heterogeneity. And automatically
## clusters the variance
within2 <- feols(Protest~ secretpol_revised + l.ln_pop + l.ln_gdppc  + l12gr+  l.lexclpop+
                   nbr_protest+intrastate+attempt|ccode, data=protests)
summary(within2)

## truely the same 
max(abs(lsdv$residuals-within1$residuals))
max(abs(lsdv$residuals-within2$residuals))


## here we can see the difference between the 
## total and within r-squared
c(summary(lsdv)$r.sq, summary(within1)$r.sq)


## why are these different?
## which of these are unbiased estimates? Which are consistent?
c(summary(lsdv)$sigma, summary(within1)$sigma, sqrt(summary(within2)$sigma2))



## build the weights for RE=GLS
Ti <- table(protests$ccode) #unbalanced panel so each unit has different weight
Ti <- rep(Ti, Ti)
sigma2.eps <- within2$sigma2 #unbiased and consistent 
sigma2.a <- mean(pooled$residuals^2) -sigma2.eps
protests$omega.hat <- 1- sqrt(sigma2.eps/(Ti*sigma2.a+sigma2.eps) )
mean(protests$omega.hat) ## fairly similar on this measure
protests[,paste0(var.names, ".gls"):=lapply(.SD, \(x){x-omega.hat*mean(x)}), 
         by=ccode, .SDcols=var.names ]
protests[,const.gls:=1-omega.hat]
fgls <- paste0(var.names[1], ".gls ~ -1 + const.gls + ",
               paste0(var.names[-1], ".gls", collapse=" + "))

gls <- lm(fgls, data=protests)
summary(gls)
Vcl.gls <- vcovCL(gls, cluster=protests$ccode)
round(coeftest(gls, Vcl.gls), 4)



## The lme4 package is the more common way to go here.  It takes 
## a formula of the form y~x+(1|heterogeneity). 
## However, it does not work with the sandwich package, so 
## we move the clubSandwich package for clustering.
## It also doesn't like the lmtest package that much
re  <- lmer(update(f1, . ~ . + (1|ccode)), data=protests)
Vcl.re <- vcovCR(re, cluster=protests$ccode, type="CR1")
coef_test(re, Vcl.re)


## hausman  (with iid)
Htest <- c(within2$coefficients - re@beta[-1]) %*% 
  solve(within2$cov.iid - vcov(re)[-1,-1]) %*%
  c(within2$coefficients - re@beta[-1])
pchisq(drop(Htest), df=length(within2$coefficients), lower=FALSE)


##hausman  (with clustering) but this version is sus
Htest.cl <- c(within2$coefficients - re@beta[-1]) %*% 
  solve(vcov(within2) - Vcl.re[-1,-1]) %*%
  c(within2$coefficients - re@beta[-1])
pchisq(drop(Htest.cl), df=length(within2$coefficients), lower=FALSE)



### Mundlak--pooled
Xnames <- colnames(pooled$model)[-1]
protests[,paste0(var.names, ".bar"):=lapply(.SD, \(x){mean(x)}), 
         by=ccode, .SDcols=var.names ]
mundlak.add <- paste(".~.",paste0("+", Xnames, ".bar",collapse = "" ))
mundlak.formula <- update(f1, mundlak.add)
mundlak <- lm(mundlak.formula, data=protests)
Vcl.m <- vcovCL(mundlak, cluster=protests$ccode)
round(coeftest(mundlak, Vcl.m), 4)
linearHypothesis(mundlak, paste0(Xnames, ".bar=0"), vcov=Vcl.m)


### Mundlak--CRE
cre  <- lmer(update(mundlak.formula, . ~ . + (1|ccode)), data=protests)
Vcl.cre <- vcovCR(cre, cluster=protests$ccode, type="CR1")
coef_test(cre, Vcl.cre)


### between
protests[, Protest.bar := mean(Protest), by=ccode]
f.btwn <- as.formula(paste("Protest.bar ~", 
                           paste0(Xnames, ".bar",collapse = " + " ) ))
btwn <- lm(f.btwn, data=protests)

cbind(within2$coefficients, mundlak$coef[2:9])
cbind(BtwnDiff=btwn$coef[-1]-within2$coefficients, 
      mundlak$coef[10:17]) 


modelsummary(list("Pooled"=pooled,  
                  "RE-GLS"=gls,
                  "RE-MLE"=re, 
                  "LSDV"=lsdv, 
                  "Within-lm"=within1, 
                  "Within-feols"=within2,
                  "Mundlak"=mundlak),
             vcov=list(Vcl.pooled, Vcl.gls, Vcl.re,
                       Vcl.lsdv, Vcl.within1, vcov(within2),
                       Vcl.m),
             fmt=2,
             coef_map=c("secretpol_revised"="Secret police",
                        "secretpol_revised.gls"="Secret police",
                        "secretpol_revised.within"="Secret police"),
             gof_map=c("nobs", "r.squared", "r2.within"))

