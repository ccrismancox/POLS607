library(car)
library(dplyr)
library(fixest)
library(ivreg)
library(sandwich)
library(readstata13)
library(lmtest)
library(Matrix)
rm(list=ls())


data.use <- read.dta13( "datasets/subsample.dta")
data.use <- data.use %>% 
  mutate(nattack = log(sinh(nattack)+1)) %>% 
  mutate(L.nattack = lag(nattack,1), .by=id)




within.ldv <- feols(nattack ~ v2x_corr +sp_pop_totl
                    + ny_gdp_pcap_kd
                    +kg_democracy +statefailure 
                    + L.nattack|id+year,
                    data=data.use)
summary(within.ldv)


##effects? Let's focus on the average within-unit standard deviation 
## increase in corruption. What is that?
data.use %>% 
  summarize(s=sd(v2x_corr, na.rm=TRUE), .by=id) %>% 
  summarize(mean(s))
## about 0.07

## Short-run effect of corruption on terrorism? 
## Approximation: A 1.7% increase  in terrorist attacks 
## for every within standard deviation  (all else equal)
deltaMethod(within.ldv, "0.07*v2x_corr")

## 1 period later  (about a 3% increase)
deltaMethod(within.ldv, "0.07*(v2x_corr*L.nattack+v2x_corr)")

## permanent increase? (about a 6% increase)
deltaMethod(within.ldv, "0.07*(v2x_corr/(1-L.nattack))")


### Anderson  Hsiao
col.names <- c("id", "year", "nattack",
               "v2x_corr", "sp_pop_totl",
               "ny_gdp_pcap_kd", "kg_democracy",
               "statefailure", "L.nattack")

## we lose 2 time dummies to differencing and colinearity
time.dummies <- model.matrix(~ factor(year)-1, 
                             data=data.use)[,-c(1:2)] 
colnames(time.dummies) <- paste0("year", min(data.use$year):max(data.use$year))[-c(1:2)]
data.use <- cbind(data.use, time.dummies)

## Difference the variables 
var.names <- c(col.names[-c(1:2)], colnames(time.dummies))
data.use <- data.use %>% 
  group_by(id) %>%
  mutate(across(all_of(var.names), \(x){x-lag(x)}, .names="D.{col}"),
         L2.nattack = lag(L.nattack)) %>%
  ungroup()
AH <- ivreg(D.nattack ~
              D.v2x_corr + D.sp_pop_totl + D.ny_gdp_pcap_kd + 
              D.kg_democracy + D.statefailure +
              D.L.nattack+
              D.year1973 + D.year1974 + D.year1975 + 
              D.year1976 + D.year1977 + D.year1978 + 
              D.year1979 + D.year1980 + D.year1981 + 
              D.year1982 + D.year1983 + D.year1984 + 
              D.year1985 + D.year1986 + D.year1987 + 
              D.year1988 + D.year1989 + D.year1990 + 
              D.year1991 + D.year1992 + D.year1993 + 
              D.year1994 + D.year1995 + D.year1996 + 
              D.year1997 + D.year1998 + D.year1999 + 
              D.year2000 + D.year2001 + D.year2002 + 
              D.year2003 + D.year2004 + D.year2005 + 
              D.year2006 + D.year2007 + D.year2008 + 
              D.year2009 + D.year2010 + D.year2011 + 
              D.year2012 + D.year2013 + D.year2014 + 
              D.year2015 + D.year2016 + D.year2017 + 
              D.year2018-1 |
              D.v2x_corr + D.sp_pop_totl + D.ny_gdp_pcap_kd + 
              D.kg_democracy + D.statefailure +
              L2.nattack+
              D.year1973 + D.year1974 + D.year1975 + 
              D.year1976 + D.year1977 + D.year1978 + 
              D.year1979 + D.year1980 + D.year1981 + 
              D.year1982 + D.year1983 + D.year1984 + 
              D.year1985 + D.year1986 + D.year1987 + 
              D.year1988 + D.year1989 + D.year1990 + 
              D.year1991 + D.year1992 + D.year1993 + 
              D.year1994 + D.year1995 + D.year1996 + 
              D.year1997 + D.year1998 + D.year1999 + 
              D.year2000 + D.year2001 + D.year2002 + 
              D.year2003 + D.year2004 + D.year2005 + 
              D.year2006 + D.year2007 + D.year2008 + 
              D.year2009 + D.year2010 + D.year2011 + 
              D.year2012 + D.year2013 + D.year2014 + 
              D.year2015 + D.year2016 + D.year2017 + 
              D.year2018-1,
            data=data.use,
            x=TRUE)
V.AH <- vcovCL(AH, cluster=data.use$id)
coeftest(AH, V.AH)[1:6,]


resid.dat <- data.use %>%
  select(id,year,D.nattack) %>%
  mutate(AH.resid = predict(AH, newdata=data.use)- D.nattack) %>% 
  mutate(Lah.resid = lag(AH.resid), 
         L2.ah.resid = lag(AH.resid,2), 
         .by=id)


## first one should be strong and negative (check)
summary(lm(AH.resid~Lah.resid-1, data=resid.dat))
## second should be very zero (could be better but definitely small)
summary(lm(AH.resid~L2.ah.resid-1,  data=resid.dat))


## Short-run effect of corruption on terrorism? 
## Roughly a 1.4% increase  in terrorist attacks 
## for every within standard deviation 
deltaMethod(AH, "0.07*D.v2x_corr", vcov=V.AH)

## 1 period later  (about a 2% increase)
deltaMethod(AH, "0.07*(D.v2x_corr*D.L.nattack+D.v2x_corr)", vcov=V.AH)

## permanent increase? (about a 2% increase)
deltaMethod(AH, "0.07*(D.v2x_corr/(1-D.L.nattack))", vcov=V.AH)


#### Arellano and Bond ####

## step 1 create a balanced the panel using fake data
data.use$real.data <- 1
pseudo.data <- data.frame(id=rep(unique(data.use$id), each=48),
                          year=rep(1971:2018))
data.ab <- merge(data.use, pseudo.data, by=c("id",  "year"), all=TRUE)
data.ab$real.data[is.na(data.ab$real.data)] <- 0 
dropU <- data.ab %>% 
  summarize(across(starts_with("D."),
                   \(x){all(is.na(x))}, 
                   .names="{.col}"),
            .by=id) %>% 
  mutate(id=NULL) %>%
  apply(1, any) %>% 
  which()

data.ab <- data.ab %>% 
  filter(!id %in% unique(id)[dropU])

## step 2 pull out the data of interest
X <- data.ab %>% 
  filter(year >= 1973) %>% 
  select(starts_with("D."))

keep <- (1-apply(X, 1, anyNA)) ## create a real measure of sample size
X[apply(X, 1, anyNA),] <- 0 ## make the dropped rows all 0
y <- X$D.nattack
Ly <- X$D.L.nattack
X$D.nattack <- NULL
X$D.L.nattack <- NULL
X <- as.matrix(X)


N <- length(unique(data.ab[data.ab$year>=1973 ,][keep==1,]$id))
print(N)


## useful function for creating multiple lags 
## with tidy
lag_multiple <- function(x, n_vec){
  Mat <- sapply(n_vec, lag, x = x)
  colnames(Mat) <- paste0("lag", n_vec)
  as_tibble(Mat)
}

## fill in NA lagged attacks with 0s to create
## the instrument matrix 
data.ab$nattack[data.ab$real.data==0] <- 0
lags <- subset(data.ab, select=c("id", "year", "nattack"))
## generate all the lags
lagMat <- lags %>% mutate(lag_multiple(nattack,47:2), .by=id)
lagMat <- lagMat[data.ab$year >= 1973,]
lagMat <- cbind(lagMat, X, keep) #bind in the things that instrument for themselves and keep
lagMat <- lagMat %>% select(!c(year, nattack))
X <- cbind(Ly, X) ## put the lag back into X

## I suspect there's a better way, 
## but this is what I have for now
Z.list <- tapply(lagMat,
                 INDEX=lagMat$id,
                 \(x){
                   zout <- cbind(do.call(bdiag, ## parts of lagMat that are lags of y
                                         apply(x[,2:47],1, 
                                               \(y){
                                                 y <- matrix(na.omit(y), nrow=1)
                                                 return(y)
                                               }
                                         )
                   ),
                   as.matrix(x[,48:(ncol(x)-1)])) ## parts that are self-instruments
                   ## Do not ask me why that needs to be a matrix, but it does 
                   zout[as.matrix(x[,ncol(x)]==0),] <- 0  ##make sure our dropped rows are all 0
                   
                   return(zout)
                 })
Z <- do.call(rbind, Z.list)


Omega1.hat.inv <- solve(Reduce(`+`,
                               lapply(Z.list,
                                      \(z){
                                        if(nrow(z)>0){
                                          Di <- -cbind(Diagonal(nrow(z)),0) +
                                            cbind(0,Diagonal(nrow(z)))
                                          H <- Di %*% t(Di)
                                          t(z) %*% H %*% z
                                        }else{
                                          0}
                                      })))
ab1.hat <- solve(t(X) %*% Z %*%  Omega1.hat.inv %*% t(Z) %*% X) %*%
  (t(X) %*% Z %*%  Omega1.hat.inv %*% t(Z) %*% y)
ab1.hat[1:6,]


e.ab1.hat <- y - X %*% ab1.hat
sigma2.hat <- crossprod(e.ab1.hat)/ (sum(keep)) ## consistent 
V.ab0 <- solve(t(X) %*% Z %*% Omega1.hat.inv %*% t(Z) %*% X)*drop(sigma2.hat) 
se.ab0 <- sqrt(diag(V.ab0))
cbind(ab1.hat, se.ab0, ab1.hat/se.ab0)[1:6,]


## Robust standard errors
Ze <- Z*drop(e.ab1.hat)
Omega2.hat <- Reduce(`+`,
                     lapply(unique(data.ab$id),
                            \(i){
                              Zi <- Ze[data.ab[data.ab$year>=1973,]$id==i,]
                              return(tcrossprod(colSums(Zi)))
                            }
                     )
)

V.ab1 <- solve(t(X) %*% Z %*% Omega1.hat.inv %*% t(Z) %*% X) %*%
  (t(X) %*% Z %*% Omega1.hat.inv %*%Omega2.hat 
   %*% Omega1.hat.inv %*% t(Z) %*% X) %*% 
  solve(t(X) %*% Z %*% Omega1.hat.inv %*% t(Z) %*% X)
se.ab1 <- sqrt(diag(V.ab1))
cbind(ab1.hat, se.ab0, se.ab1)[1:6,]

## AR tests 
ar.dat <- data.ab %>% 
  filter(year>=1973) %>% 
  select(c(id, year)) %>%
  mutate(ehat =drop(e.ab1.hat)) %>% 
  group_by(id) %>% 
  mutate(L.ehat = lag(ehat),
         L2.ehat = lag(ehat,2)) %>% 
  ungroup()
summary(lm(ehat~L.ehat-1, data=ar.dat))
summary(lm(ehat~L2.ehat-1, data=ar.dat))

## Sargan on the first step
sarganJ <- (t(e.ab1.hat) %*% Z %*% Omega1.hat.inv %*% 
              t(Z) %*% e.ab1.hat)/sigma2.hat
sarganJ
pchisq(drop(sarganJ), df=ncol(Z)-ncol(X), lower=FALSE) ## good

##### strength#####
## We won't be able to invert A V_Cluster A' with this many instruments, so 
## we'll use regular robust SEs
g.hat <- solve(crossprod(Z)) %*% t(Z) %*% X[,"Ly"]
Ze <- Z*drop(e.ab1.hat)
bread <- solve(crossprod(Z))
meat <- crossprod(Ze)
V1.white <- bread %*% meat %*% bread

A <- cbind(Diagonal(ncol(Z)-ncol(X)+1), 
           Matrix(0, nrow=ncol(Z)-ncol(X)+1, ncol=ncol(X)-1))
Fstat.white <- t(A %*% g.hat ) %*% 
  solve(A %*% V1.white %*% t(A)) %*% 
  (A %*% g.hat ) / nrow(A)
Fstat.white

## Not great


#### effects ####

## short run: about a 3.5% increase 
c(0.07*ab1.hat["D.v2x_corr",], 0.07*sqrt(se.ab1["D.v2x_corr"]^2))


## One year later about a 5.6% increase
D1 <- c(1+ab1.hat["Ly",], ab1.hat["D.v2x_corr",])
c(0.07*(ab1.hat["D.v2x_corr",]+
          ab1.hat["D.v2x_corr",]*ab1.hat["Ly",]), 
  0.07*drop(sqrt(D1 %*% V.ab1[2:1, 2:1] %*% D1)))



## Permanent increase in corruption: about a 7% increase in terrorism
Dlr <- c(1/(1-ab1.hat["Ly",]), ab1.hat["D.v2x_corr",]/(1-ab1.hat["Ly",])^2)
c(0.07*(ab1.hat["D.v2x_corr",]/(1-ab1.hat["D.v2x_corr",])), 
  0.07*drop(sqrt(Dlr %*% V.ab1[2:1, 2:1] %*% Dlr)))

















######### Using fewer lags ##################
#### Fewer lags but with the two-step ####
Z.list <- tapply(lagMat,
                 INDEX=lagMat$id,
                 \(x){
                   zout <- cbind(
                     do.call(bdiag,
                             apply(x[,2:47],1,
                                   \(y){
                                     y <- matrix(na.omit(y), nrow=1)
                                     if(ncol(y)>2){ ## number of lags
                                       y <- y[,(ncol(y)-1):ncol(y),drop=FALSE]
                                     }
                                     return(y)
                                   }
                             )
                     ),
                     as.matrix(x[,48:(ncol(x)-1)]))
                   zout[as.matrix(x[,ncol(x)]==0),] <- 0 
                   ## Do not ask me why that needs to be a matrix, but it does 
                   return(zout)
                 })
Z <- do.call(rbind, Z.list)

Omega1.hat.inv <- solve(Reduce(`+`,
                               lapply(Z.list,
                                      \(z){
                                        Di <- -cbind(Diagonal(nrow(z)),0) +
                                          cbind(0,Diagonal(nrow(z)))
                                        H <- Di %*% t(Di)
                                        t(z) %*% H %*% z
                                      })))
ab1.hat <- solve(t(X) %*% Z %*%  Omega1.hat.inv %*% t(Z) %*% X) %*%
  (t(X) %*% Z %*%  Omega1.hat.inv %*% t(Z) %*% y)
ab1.hat[1:6,]


e.ab1.hat <- y - X %*% ab1.hat
sigma2.hat <- drop(crossprod(e.ab1.hat))/ (sum(keep)) ## consistent 

V.ab0 <- solve(t(X) %*% Z %*% (Omega1.hat.inv/sigma2.hat) %*% t(Z) %*% X)
se.ab0 <- sqrt(diag(V.ab0))
cbind(ab1.hat, se.ab0, ab1.hat/se.ab0)[1:6,]

Ze <- Z*drop(e.ab1.hat)
Omega2.hat <- Reduce(`+`,
                     lapply(unique(data.ab$id),
                            \(i){
                              Zi <- Ze[data.ab[data.ab$year>=1973,]$id==i,]
                              return(tcrossprod(colSums(Zi)))
                            }
                     )
)


bread <- solve(t(X) %*% Z %*% Omega1.hat.inv %*% t(Z) %*% X)
meat <-  (t(X) %*% Z %*% Omega1.hat.inv %*%Omega2.hat %*% 
            Omega1.hat.inv %*% t(Z) %*% X)
V.ab1 <- bread  %*% meat %*% bread
se.ab1 <- sqrt(diag(V.ab1))
cbind(ab1.hat, se.ab0, se.ab1)[1:6,]


## AR tests 
ar.dat <- data.ab %>% 
  filter(year>=1973) %>% 
  select(c(id, year)) %>%
  mutate(ehat =drop(e.ab1.hat)) %>% 
  group_by(id) %>% 
  mutate(L.ehat = lag(ehat),
         L2.ehat = lag(ehat,2)) %>% 
  ungroup()
summary(lm(ehat~L.ehat-1, data=ar.dat))
summary(lm(ehat~L2.ehat-1, data=ar.dat))

## Sargan on the first step -- Not as good as before
sarganJ <- t(e.ab1.hat) %*% Z %*% 
  (Omega1.hat.inv/sigma2.hat) %*%
  t(Z) %*% e.ab1.hat 
sarganJ
pchisq(drop(sarganJ), df=ncol(Z)-ncol(X), lower=FALSE)

## The two step
Omega2.hat.inv <- solve(Omega2.hat)
ab2.hat <- solve(t(X) %*% Z %*%  Omega2.hat.inv %*% t(Z) %*% X) %*%
  (t(X) %*% Z %*%  Omega2.hat.inv %*% t(Z) %*% y)
ab2.hat[1:6]
e.ab2.hat <- y - X %*% ab2.hat

V.ab2 <- solve(t(X) %*% Z %*% Omega2.hat.inv %*% t(Z) %*% X)
se.ab2 <- sqrt(diag(V.ab2))
cbind(ab2.hat, se.ab2, ab2.hat/se.ab2)[1:6,]

## two step with correction
Dgi.dtheta <- -t(Z) %*% X 
Gn <- Dgi.dtheta
gn <- (t(Z) %*% e.ab2.hat)

D <- matrix(0, ncol(X),ncol(X))
for(j in 1:ncol(X)){
  GAMMA.j <- (t(Z)%*%drop(e.ab2.hat) %*% t(Dgi.dtheta[,j,drop=FALSE]))
  dOmega.j <- (GAMMA.j + t(GAMMA.j))
  D[,j] <- drop(solve(t(Gn) %*% Omega2.hat.inv %*% Gn) %*% 
                  t(Gn) %*% 
                  (Omega2.hat.inv  %*% dOmega.j %*% Omega2.hat.inv) 
                %*% gn )/N
}

V2.corrected <- (V.ab2) + 
  D %*% (V.ab2) + 
  (V.ab2) %*% t(D) + 
  D %*% (V.ab1) %*% t(D)

se.ab2.c <- sqrt(diag(V2.corrected))
cbind(ab2.hat, se.ab2, se.ab2.c)[1:6,]

## AR tests 
ar.dat <- data.ab %>% 
  filter(year>=1973) %>% 
  select(id, year) %>% 
  mutate(ehat = drop(e.ab2.hat)) %>% 
  mutate(L.ehat = lag(ehat),
         L2.ehat = lag(ehat,2),
         .by=id)
summary(lm(ehat~L.ehat-1, data=ar.dat))
summary(lm(ehat~L2.ehat-1, data=ar.dat))


## Sargan on the second step
sarganJ <- t(e.ab2.hat) %*% Z %*% Omega2.hat.inv %*% t(Z) %*% e.ab2.hat 
sarganJ
pchisq(drop(sarganJ), df=ncol(Z)-ncol(X), lower=FALSE)
## better




## strength?
g.hat <- solve(crossprod(Z)) %*% t(Z) %*% X[,"Ly"]
Ze <- Z*drop(e.ab2.hat)
bread <- solve(crossprod(Z))
meat <- Reduce(`+`,
               lapply(unique(data.ab$id),
                      \(i){
                        Zi <- Ze[data.ab[data.ab$year>=1973,]$id==i,]
                        return(tcrossprod(colSums(Zi)))
                      }
               )
)
Vcl <- bread %*% meat %*% bread


A <- cbind(Diagonal(ncol(Z)-ncol(X)+1), 
           Matrix(0, nrow=ncol(Z)-ncol(X)+1, ncol=ncol(X)-1))
Fstat.cl<- (t(A %*% g.hat ) %*% 
              solve(A %*% Vcl %*% t(A)) %*% 
              (A %*% g.hat )) / nrow(A)
Fstat.cl
## Not bad



##### Unit root tests #####
library(readstata13)
library(dplyr)
library(Matrix)
library(fUnitRoots)
source("panelFunctions.r")

terror <- read.dta13("Rcode/corruption_terrorism/subsample.dta")

#### LLC test ####
data.llc <- terror %>% 
  select(id, year, nattack) %>%
  mutate(nattack= nattack-mean(nattack,na.rm=TRUE), .by=year) %>%
  mutate(y=nattack, 
         l.y = lag(y),
         D.y = y-l.y,
         .by=id) %>%
  mutate(Ti = length(na.omit(D.y)), .by=id)


N <- length(unique(data.llc$id))
data.llc %>% summarize(mean(Ti), .by=id) %>% summary()



e.tilde <- v.tilde <- list()
s<- rep(0, N)
corrs <- rep(0, N)
for(i in 1:N){
  data.i <- data.llc %>% filter(id==unique(data.llc$id)[i])
  fit0 <- lm(D.y~l.y, data=data.i)
  corrs[i] <- lm(fit0$residuals[-length(fit0$residuals)]~
                   fit0$residuals[-1]-1)$coef
  
  fit1 <- lm(D.y~1, data=data.i)
  e.i <- predict(fit1, data.i)- data.i$D.y
  fit2 <- lm(l.y~ 1, data=data.i)
  v.i <- predict(fit2, data.i)- data.i$l.y
  fit3 <- lm(e.i~v.i-1)
  s.eps.i <- summary(fit3)$sigma
  e.tilde[[i]] <- e.i/s.eps.i
  v.tilde[[i]] <- v.i/s.eps.i
  lag.coef <- fit0$coefficients[grep(names(fit0$coefficients), pattern="L[1-4].D.y")]
  if(is.null(lag.coef)){lag.coef <- 0}
  s[i] <- abs(1-sum(lag.coef))
  
}
S <- mean(s)
e.tilde <- unlist(e.tilde)
v.tilde <- unlist(v.tilde)
mean(corrs)

fit4 <- lm(e.tilde~v.tilde-1)
NT <-  length(fit4$residuals)
t.stat <- summary(fit4)$coef["v.tilde", "t value"]
se <-  summary(fit4)$coef["v.tilde", "Std. Error"]
s2.eps.tilde <- summary(fit4)$sigma^2

### values for mu and sigma star from LLC table 2
### For "model 2" and average T of about 45

t.adj <-( t.stat -  (NT*S*se*(-0.533))/s2.eps.tilde)/0.837
t.adj
pnorm(t.adj)
## reject the null of unit root



### IPS test ####
data.ips <- terror %>% 
  select(id, year, nattack) %>% 
  mutate(Ti = length(na.omit(nattack)), .by=id) %>%
  filter(Ti > 10) %>% 
  mutate(y=nattack, 
         y= y- mean(y,na.rm=TRUE), 
         .by=year)%>%
  mutate(l.y = lag(y),
         D.y = y-l.y,
         .by=id) 
data.ips.clean <- na.omit(data.ips) 
N <- length(unique(data.ips$id))

## System estimator
DeltaN <- sparse.model.matrix(~factor(id) -1 ,data=data.ips.clean)
y <- data.ips.clean$D.y
X <- cbind(data.ips.clean$l.y* DeltaN, DeltaN)
colnames(X)[1:N] <- paste0("phi",1:N)
XX <- solve(crossprod(X))
ests <- XX %*% t(X) %*% y


## t-normal
si <- lapply(split.matrix(y-X%*%ests, data.ips.clean$id),
             \(x){drop(crossprod(x)/(length(x)-2))})
se <- mapply(\(x,y){sqrt(diag(y*solve(crossprod(x))))[2]},
             x=split.matrix(cbind(1,data.ips.clean$l.y), data.ips.clean$id),
             y=si)
phi <- ests[1:N]
t1 <- phi/se
mean(t1) ## critical value for tbar (based on Table 2) -1.67 and the average Ti
## t-tilde
st <- mapply(\(x,y){sqrt(diag(var(y)*solve(crossprod(x))))[2]},
             x=split.matrix(cbind(1,data.ips.clean$l.y), data.ips.clean$id),
             y=split(data.ips.clean$D.y, data.ips.clean$id))
t.tilde <- phi/st
mean(t.tilde)

## Approximating using table 1 in IPS
E.ti.tilde <- -1.47
V.ti.tilde <- 0.6475
z.stat <- (mean(t.tilde) -  E.ti.tilde)/sqrt(V.ti.tilde/N)
pnorm(z.stat)
## Reject the null



ttp <- matrix(0, nrow=N, ncol=3)
## Or regression by regression approach
for(i in 1:N){
  fit.i <- lm(D.y~l.y, data=data.ips, 
              subset=id==unique(data.ips$id)[i],
              x=TRUE, y=TRUE)
  t.stat <- summary(fit.i)$coeff[,"t value"]["l.y"]
  V.tilde <- var(fit.i$y) * solve(crossprod(fit.i$x))
  t.til <- fit.i$coefficients["l.y"]/sqrt(diag(V.tilde)["l.y"])
  
  ## Let's do the Fisher test while we're in here
  pi <- suppressWarnings(adfTest(data.ips[data.ips$id==unique(data.ips$id)[i], ]$y,
                                 lags=0, type="c"))
  pi <- pi@test$p.value
  ttp[i,] <- c(t.stat, t.til, pi)
}
summary(ttp)
FisherP <- -2*sum(log(ttp[,3]))
FisherZ <-sum(qnorm(ttp[,3]))/sqrt(N)
pchisq(FisherP, lower.tail = FALSE, df=2*N)
pnorm(FisherZ)
