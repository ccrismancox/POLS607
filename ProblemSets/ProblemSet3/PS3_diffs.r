### Casey Crisman-Cox
### Code for if you decided that 
### log(GDP pc) * 100 was not stationary



## packages 
library(readstata13)
library(dplyr)
library(fixest)
library(car)
library(Matrix)
library(xtable)
library(ggplot2)
library(fUnitRoots)

rm(list=ls())

## additional functions and data
source("panelFunctions.r")
source("effects.r")
ddcg <- read.dta13("DDCGdata_final.dta")


ddcg <- ddcg %>% 
  mutate(y0 = y, #save the original measure just in case
         y = y-lag(y), ## replace y with diff(y)
         l.y = lag(y),
         l2.y = lag(y,2),
         l3.y = lag(y,3),
         l4.y = lag(y,4),
         .by=wbcode2)


## within estimator AR(1)
within.l1 <- feols(y~dem+l.y|wbcode2+year, data=ddcg)


within.sign.check <- feols(l.y~dem|wbcode2+year, data=ddcg)$coefficient
sign(within.sign.check) 

## as true rho increases the bias increases in magnitude
plot(-within.sign.check* (-(seq(0,.99, length=100)+1)/(39)) ~ 
       seq(0,.99, length=100),
     ylab=expression("Bias in "~hat(beta)[1]),
     xlab=expression(rho))

## AR(2) and AR(4) with the within
within.l2 <- feols(y~dem+l.y+l2.y|wbcode2+year, data=ddcg)
within.l4 <- feols(y~dem+l.y+l2.y+l3.y+l4.y|wbcode2+year, data=ddcg)

## create a list with the quantities of interest 
eff.within <- list(effects(within.l1$coefficients, vcov(within.l1)),
                   effects(within.l2$coefficients, vcov(within.l2)),
                   effects(within.l4$coefficients, vcov(within.l4)))



#### GMM lag 1 ####
data.ab <- ddcg %>% # subset to relevant columns 
  select(wbcode2, year, y, dem, l.y,l2.y,l3.y,l4.y, starts_with("yy")) %>%
  group_by(wbcode2) %>%
  filter(!all(is.na(y))) %>% ## and remove any all NA units
  ungroup()  %>% 
  filter(year > 1960) # lose 1960 due to differencing in y



data1 <- data.ab %>% 
  arrange(wbcode2, year) %>% 
  group_by(wbcode2) %>%
  mutate(across(c(y, dem, l.y, starts_with("yy")),
                .fn=\(x){x-lag(x)}, 
                .names = "D.{.col}")) %>%
  ungroup() 

## Extract relevant data
X <- data1 %>% 
  filter(year >= 1963) %>% #lose 2 more years to lags and differences
  select( starts_with("D.")) 
keep <- (1-apply(X, 1, anyNA))
X[keep==0,] <- 0
y <- X$D.y #pull the DV
Xendo <- with(X, cbind(D.l.y)) #pull endogenous variables
## delete the non-exogenous variables
X$D.y <- NULL
X$D.l.y <- NULL
X$D.yy1  <- NULL 
X$D.yy2  <- NULL
X$D.yy3  <- NULL # lose first three year dummies

X <- as.matrix(X)
qr(crossprod(X))$rank == ncol(X) #check rank



N <- length(unique(data1[data1$year >= 1963,][keep==1,]$wbcode2))

lag_multiple <- function(x, n_vec){
  Mat <- sapply(n_vec, lag, x = x)
  colnames(Mat) <- paste0("lag", n_vec)
  as_tibble(Mat)
}



lags <- subset(data1, select=c("wbcode2", "year", "y"))
lags[is.na(lags$y),]$y <- 0
lagMat <- lags %>% mutate(lag_multiple(y,49:2), .by=wbcode2) #up to 49 lags
lagMat <- lagMat[data1$year >= 1963,]
lagMat <- cbind(lagMat, X, keep)
lagMat <- lagMat %>% select(!c(year, y))

## instruments 
Z.list <- tapply(lagMat,
                 INDEX=lagMat$wbcode2,
                 \(x){
                   zout <- cbind(do.call(bdiag,
                                         apply(x[,2:49],1,## 2-49 are lags
                                               \(y){
                                                 y <- matrix(na.omit(y), nrow=1)
                                                 return(y)
                                               }
                                         )
                   ),
                   as.matrix(x[,c(50:(ncol(x)-1))]))##50+ are exogenous X and keep
                   zout[as.matrix(x[,ncol(x)]==0),] <- 0 ## Do not ask me why that needs to be a matrix, but it does 
                   return(zout)
                 })
Z <- do.call(rbind, Z.list)
dim(Z)


Di <- -cbind(Diagonal(nrow(Z.list[[1]])),0) +
  cbind(0,Diagonal(nrow(Z.list[[1]])))
H <- Di %*% t(Di)
Omega1.hat.inv <- solve(Reduce(`+`,
                               lapply(Z.list,
                                      \(z){
                                        if(nrow(z)>0){
                                          t(z) %*% H %*% z
                                        }else{
                                          0
                                        }
                                      })))
X <- cbind(X[,1], Xendo, X[,-1])
ab1.hat <- solve(t(X) %*% Z %*%  Omega1.hat.inv %*% t(Z) %*% X) %*%
  (t(X) %*% Z %*%  Omega1.hat.inv %*% t(Z) %*% y)


e.ab1.hat <- y - X %*% ab1.hat
sigma2.hat <- crossprod(e.ab1.hat)/ (sum(keep)) ## consistent 
V.ab0 <- solve(t(X) %*% Z %*% Omega1.hat.inv %*% t(Z) %*% X)*drop(sigma2.hat)
se.ab0 <- sqrt(diag(V.ab0))
cbind(ab1.hat, se.ab0, ab1.hat/se.ab0)[1:2,]

Ze <- Z*drop(e.ab1.hat)
Omega2.hat <- Reduce(`+`,
                     lapply(unique(data1$wbcode2),
                            \(i){
                              Zi <- Ze[data1[data1$year>=1963,]$wbcode2==i,]
                              return(tcrossprod(colSums(Zi)))
                            }
                     )
)


V.ab1 <- solve(t(X) %*% Z %*% Omega1.hat.inv %*% t(Z) %*% X) %*%
  (t(X) %*% Z %*% Omega1.hat.inv %*%Omega2.hat %*% Omega1.hat.inv %*% t(Z) %*% X) %*% 
  solve(t(X) %*% Z %*% Omega1.hat.inv %*% t(Z) %*% X)
se.ab1 <- sqrt(diag(V.ab1))

## store the results
gmm1.out <- list(est=ab1.hat, V=V.ab1,
                 effects=effects(est=ab1.hat[1:2],
                                 V= V.ab1[1:2, 1:2]),
                 dim=c(sum(keep), ncol(Z)))





#### GMM lag 2 #####
## just update the above in key spots
data2 <- data.ab %>% 
  arrange(wbcode2, year) %>% 
  group_by(wbcode2) %>%
  mutate(across(c(y, dem, l.y, l2.y, starts_with("yy")), # 2 lags
                .fn=\(x){x-lag(x)}, 
                .names = "D.{.col}")) %>%
  ungroup() 

X <- data2 %>% 
  filter(year >= 1964) %>% #lose another year
  select(starts_with("D.")) 
keep <- (1-apply(X, 1, anyNA))
X[keep==0,] <- 0
y <- X$D.y
Xendo <- with(X, cbind(D.l.y, D.l2.y))
X$D.y <- NULL
X$D.l.y <- NULL
X$D.l2.y <- NULL
X$D.yy1  <- NULL 
X$D.yy2  <- NULL
X$D.yy3  <- NULL
X$D.yy4  <- NULL #lose another year dummy

X <- as.matrix(X)
qr(crossprod(X))$rank == ncol(X)



N <- length(unique(data1[data1$year >= 1964,][keep==1,]$wbcode2))

lags <- subset(data1, select=c("wbcode2", "year", "y"))
lags[is.na(lags$y),]$y <- 0
lagMat <- lags %>% mutate(lag_multiple(y,49:2), .by=wbcode2)
lagMat <- lagMat[data1$year >= 1964,]
lagMat <- cbind(lagMat, X, keep)
lagMat <- lagMat %>% select(!c(year, y))


## Instruments
Z.list <- tapply(lagMat,
                 INDEX=lagMat$wbcode2,
                 \(x){
                   zout <- cbind(do.call(bdiag,
                                         apply(x[,2:48],1,
                                               \(y){
                                                 y <- matrix(na.omit(y), nrow=1)
                                                 return(y)
                                               }
                                         )
                   ),
                   as.matrix(x[,c(50:(ncol(x)-1))]))
                   zout[as.matrix(x[,ncol(x)]==0),] <- 0 ## Do not ask me why that needs to be a matrix, but it does 
                   return(zout)
                 })
Z <- do.call(rbind, Z.list)
dim(Z)
  

Di <- -cbind(Diagonal(nrow(Z.list[[1]])),0) +
  cbind(0,Diagonal(nrow(Z.list[[1]])))
H <- Di %*% t(Di)
Omega1.hat.inv <- solve(Reduce(`+`,
                               lapply(Z.list,
                                      \(z){
                                        if(nrow(z)>0){
                                          t(z) %*% H %*% z
                                        }else{
                                          0
                                        }
                                      })))
X <- cbind(X[,1], Xendo, X[,-1])

ab1.hat <- solve(t(X) %*% Z %*%  Omega1.hat.inv %*% t(Z) %*% X) %*%
  (t(X) %*% Z %*%  Omega1.hat.inv %*% t(Z) %*% y)

e.ab1.hat <- y - X %*% ab1.hat
sigma2.hat <- crossprod(e.ab1.hat)/ (sum(keep)) ## consistent 
V.ab0 <- solve(t(X) %*% Z %*% Omega1.hat.inv %*% t(Z) %*% X)*drop(sigma2.hat)
se.ab0 <- sqrt(diag(V.ab0))
cbind(ab1.hat, se.ab0, ab1.hat/se.ab0)[1:3,]

Ze <- Z*drop(e.ab1.hat)
Omega2.hat <- Reduce(`+`,
                     lapply(unique(data1$wbcode2),
                            \(i){
                              Zi <- Ze[data1[data1$year>=1964,]$wbcode2==i,]
                              return(tcrossprod(colSums(Zi)))
                            }
                     )
)




V.ab1 <- solve(t(X) %*% Z %*% Omega1.hat.inv %*% t(Z) %*% X) %*%
  (t(X) %*% Z %*% Omega1.hat.inv %*%Omega2.hat %*% Omega1.hat.inv %*% t(Z) %*% X) %*% 
  solve(t(X) %*% Z %*% Omega1.hat.inv %*% t(Z) %*% X)
se.ab1 <- sqrt(diag(V.ab1))
gmm2.out <- list(est=ab1.hat, V=V.ab1,
                 effects=effects(est=ab1.hat[1:3],
                                 V= V.ab1[1:3, 1:3]),
                 dim=c(sum(keep), ncol(Z)))








#### GMM lag 4#####
data4 <- data.ab %>% 
  arrange(wbcode2, year) %>% 
  group_by(wbcode2) %>%
  mutate(across(c(y, dem, l.y, l2.y, l3.y, l4.y, starts_with("yy")),
                .fn=\(x){x-lag(x)}, 
                .names = "D.{.col}")) %>%
  ungroup() 

X <- data4 %>% 
  filter(year >= 1966) %>%
  select(starts_with("D.")) 
keep <- (1-apply(X, 1, anyNA))
X[keep==0,] <- 0
y <- X$D.y
Xendo <- with(X, cbind(D.l.y, D.l2.y,D.l3.y, D.l4.y))
X$D.y <- NULL
X$D.l.y <- NULL
X$D.l2.y <- NULL
X$D.l3.y <- NULL
X$D.l4.y <- NULL
X$D.yy1  <- NULL 
X$D.yy2  <- NULL
X$D.yy3  <- NULL
X$D.yy4  <- NULL
X$D.yy5  <- NULL
X$D.yy6  <- NULL

X <- as.matrix(X)
qr(crossprod(X))$rank == ncol(X)



N <- length(unique(data1[data1$year >= 1966,][keep==1,]$wbcode2))


lags <- subset(data1, select=c("wbcode2", "year", "y"))
lags[is.na(lags$y),]$y <- 0
lagMat <- lags %>% mutate(lag_multiple(y,49:2), .by=wbcode2)
lagMat <- lagMat[data1$year >= 1966,]
lagMat <- cbind(lagMat, X, keep)
lagMat <- lagMat %>% select(!c(year, y))





## Instruments
Z.list <- tapply(lagMat,
                 INDEX=lagMat$wbcode2,
                 \(x){
                   zout <- cbind(do.call(bdiag,
                                         apply(x[,2:46],1,
                                               \(y){
                                                 y <- matrix(na.omit(y), nrow=1)
                                                 return(y)
                                               }
                                         )
                   ),
                   as.matrix(x[,c(50:(ncol(x)-1))]))
                   zout[as.matrix(x[,ncol(x)]==0),] <- 0 ## Do not ask me why that needs to be a matrix, but it does 
                   return(zout)
                 })
Z <- do.call(rbind, Z.list)
dim(Z)

Di <- -cbind(Diagonal(nrow(Z.list[[1]])),0) +
  cbind(0,Diagonal(nrow(Z.list[[1]])))
H <- Di %*% t(Di)
Omega1.hat.inv <- solve(Reduce(`+`,
                               lapply(Z.list,
                                      \(z){
                                        if(nrow(z)>0){
                                          t(z) %*% H %*% z
                                        }else{
                                          0
                                        }
                                      })))
X <- cbind(X[,1], Xendo, X[,-1])

ab1.hat <- solve(t(X) %*% Z %*%  Omega1.hat.inv %*% t(Z) %*% X) %*%
  (t(X) %*% Z %*%  Omega1.hat.inv %*% t(Z) %*% y)






e.ab1.hat <- y - X %*% ab1.hat
sigma2.hat <- crossprod(e.ab1.hat)/ (sum(keep)) ## consistent 
V.ab0 <- solve(t(X) %*% Z %*% Omega1.hat.inv %*% t(Z) %*% X)*drop(sigma2.hat)
se.ab0 <- sqrt(diag(V.ab0))
cbind(ab1.hat, se.ab0, ab1.hat/se.ab0)[1:5,]

Ze <- Z*drop(e.ab1.hat)
Omega2.hat <- Reduce(`+`,
                     lapply(unique(data1$wbcode2),
                            \(i){
                              Zi <- Ze[data1[data1$year>=1966,]$wbcode2==i,]
                              return(tcrossprod(colSums(Zi)))
                            }
                     )
)




V.ab1 <- solve(t(X) %*% Z %*% Omega1.hat.inv %*% t(Z) %*% X) %*%
  (t(X) %*% Z %*% Omega1.hat.inv %*%Omega2.hat %*% Omega1.hat.inv %*% t(Z) %*% X) %*% 
  solve(t(X) %*% Z %*% Omega1.hat.inv %*% t(Z) %*% X)
se.ab1 <- sqrt(diag(V.ab1))
gmm4.out <- list(est=ab1.hat, V=V.ab1,
                 effects=effects(est=ab1.hat[1:5],
                                 V= V.ab1[1:5, 1:5]),
                 dim=c(sum(keep), ncol(Z)))





###### AR test for gmm4######
ar.dat <- data.ab %>% 
  filter(year>=1966) %>% 
  select(c(wbcode2, year)) %>%
  mutate(ehat =drop(e.ab1.hat)) %>% 
  mutate(ehat = ifelse(ehat==0, NA, ehat)) %>% 
  group_by(wbcode2) %>% 
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















###### fewer Instruments #######
## Instruments
Z.list <- tapply(lagMat,
                 INDEX=lagMat$wbcode2,
                 \(x){
                   zout <- cbind(do.call(bdiag,
                                         apply(x[,2:46],1,
                                               \(y){
                                                 y <- matrix(na.omit(y), nrow=1)
                                                 if(ncol(y)>2){
                                                   y <- y[,(ncol(y)-1):ncol(y),drop=FALSE]
                                                 }
                                                 return(y)
                                               }
                                         )
                   ),
                   as.matrix(x[,c(50:(ncol(x)-1))]))
                   zout[as.matrix(x[,ncol(x)]==0),] <- 0 ## Do not ask me why that needs to be a matrix, but it does 
                   return(zout)
                 })
Z <- do.call(rbind, Z.list)
dim(Z)

Di <- -cbind(Diagonal(nrow(Z.list[[1]])),0) +
  cbind(0,Diagonal(nrow(Z.list[[1]])))
H <- Di %*% t(Di)
Omega1.hat.inv <- solve(Reduce(`+`,
                               lapply(Z.list,
                                      \(z){
                                        if(nrow(z)>0){
                                          t(z) %*% H %*% z
                                        }else{
                                          0
                                        }
                                      })))

ab1.hat <- solve(t(X) %*% Z %*%  Omega1.hat.inv %*% t(Z) %*% X) %*%
  (t(X) %*% Z %*%  Omega1.hat.inv %*% t(Z) %*% y)



e.ab1.hat <- y - X %*% ab1.hat
sigma2.hat <- crossprod(e.ab1.hat)/ (sum(keep)) ## consistent 
V.ab0 <- solve(t(X) %*% Z %*% Omega1.hat.inv %*% t(Z) %*% X)*drop(sigma2.hat)
se.ab0 <- sqrt(diag(V.ab0))
cbind(ab1.hat, se.ab0, ab1.hat/se.ab0)[1:5,]

Ze <- Z*drop(e.ab1.hat)
Omega2.hat <- Reduce(`+`,
                     lapply(unique(data1$wbcode2),
                            \(i){
                              Zi <- Ze[data1[data1$year>=1966,]$wbcode2==i,]
                              return(tcrossprod(colSums(Zi)))
                            }
                     )
)




V.ab1 <- solve(t(X) %*% Z %*% Omega1.hat.inv %*% t(Z) %*% X) %*%
  (t(X) %*% Z %*% Omega1.hat.inv %*%Omega2.hat %*% Omega1.hat.inv %*% t(Z) %*% X) %*% 
  solve(t(X) %*% Z %*% Omega1.hat.inv %*% t(Z) %*% X)
se.ab1 <- sqrt(diag(V.ab1))
gmmFewer.out <- list(est=ab1.hat, V=V.ab1,
                     effects=effects(est=ab1.hat[1:5],
                                     V= V.ab1[1:5, 1:5]),
                     dim=c(sum(keep), ncol(Z)))




## AR tests 
ar.dat <- data.ab %>% 
  filter(year>=1966) %>% 
  select(c(wbcode2, year)) %>%
  mutate(ehat =drop(e.ab1.hat)) %>% 
  mutate(ehat = ifelse(ehat==0, NA, ehat)) %>% 
  group_by(wbcode2) %>% 
  mutate(L.ehat = lag(ehat),
         L2.ehat = lag(ehat,2)) %>% 
  ungroup()
summary(lm(ehat~L.ehat-1, data=ar.dat))
summary(lm(ehat~L2.ehat-1, data=ar.dat))

## Sargan on the first step
sarganJ <- (t(e.ab1.hat) %*% Z %*% Omega1.hat.inv %*% 
              t(Z) %*% e.ab1.hat)/sigma2.hat
sarganJ
pchisq(drop(sarganJ), df=ncol(Z)-ncol(X), lower=FALSE) 

## Not great












## Figure 2 for this I'm pulling from my effects 
## function but now saving each incremental output
est <- within.l4$coefficients
V <- vcov(within.l4)


rho <- est[2:length(est)]
rho <- c(rho, rep(0, 4-length(rho)))
rho[is.na(rho)] <- 0
beta <- est[1]
eff.out <- rep(0, 30)
D.out <- matrix(0, 30, 5)

eff1 <- beta
eff2 <- rho[1]*eff1 + beta
eff3 <- rho[1]*eff2 + rho[2]*eff1 + beta
eff4 <- rho[1]*eff3 + rho[2]*eff2 + rho[3]*eff1 + beta
eff <-  c(eff4, eff3, eff2, eff1)
eff.out[4:1] <- eff

Deff1 <- c(1, rep(0, 4))
Deff2 <- c(1+rho[1]*Deff1[1], eff1,  rep(0, 3))
Deff3 <- c(1+rho[1]*Deff2[1] + rho[2]*Deff1[1], 
           eff2+ rho[1]*Deff2[2],
           eff1,
           rep(0,2))
Deff4 <- c(1+rho[1]*Deff3[1] + rho[2]*Deff2[1]+rho[3]*Deff1[1], 
           eff3 +rho[1]*Deff3[2] + rho[2]*Deff2[2],
           eff2 + rho[1]*Deff3[3]+rho[2]*Deff2[3],
           eff1, 
           0)

D <- rbind(Deff4, Deff3, Deff2, Deff1)
D.out[4:1,] <- D[1:4, ]

for(j in 5:30){
  eff.j <- rho %*% eff +beta
  
  Dj <- c(1+ rho %*% D[,1],
          rho %*% D[,-1] + eff)
  D <- rbind(Dj,D[-4,])    
  eff <- c(eff.j, eff[-4])
  eff.out[j] <- eff.j
  D.out[j,]  <- Dj
  
}
eff.out <- cumsum(eff.out)
D.out <- apply(D.out, 2, cumsum)
Vout <- D.out %*% V %*% t(D.out)
se.out <- sqrt(diag(Vout))


plot.df <- data.frame(est=eff.out,
                      lo=eff.out-1.96*se.out,
                      hi=eff.out+1.96*se.out,
                      years=0:29)

fig2 <- ggplot(plot.df)+
  geom_ribbon(aes(x=years, ymin=lo, ymax=hi), alpha=.25)+
  geom_line(aes(x=years, y=est)) +
  theme_bw(16)+
  xlab("Years since democratization") +
  ylab("Change in logged GDP per capita (levels)")
ggsave(fig2, file="figure2.pdf", width=3, height=3.5)



## create tables
within.list <- list(within.l1, within.l2, within.l4)
within.out <- do.call(cbind, 
                      mapply(cbind,
                             x=lapply(within.list, 
                                      \(x){cbind(rbind(x$coefficients,
                                                       sqrt(diag(vcov(x)))),
                                                 matrix(NA, ncol=5-length(x$coefficients),
                                                        nrow=2)
                                      )}),
                             y=lapply(eff.within, t), SIMPLIFY = FALSE)
)

within.out <- num2str(within.out,digits=2)
within.out[2,] <- paste0("(", within.out[2,], ")")
within.out[grep(within.out, pattern="NA")] <- NA

within.out <- matrix(within.out, ncol=3)
within.out <- rbind(within.out, 
                    paste0("\\multicolumn{1}{c}{",
                           sapply(within.list, \(x){x$nobs}),
                           "}"),
                    rep(NA,3))

gmm.list <- list(gmm1.out, gmm2.out, gmm4.out, gmmFewer.out)
gmm.out <- do.call(cbind, 
                   lapply(gmm.list, 
                          \(x){
                            idx <- c(1, grep(rownames(x$est), pattern="D.l"))
                            cbind(rbind(x$est[idx],
                                        sqrt(diag(x$V)[idx])),
                                  matrix(NA, ncol=5-length(x$est[idx]),
                                         nrow=2),
                                  t(x$effects))}))

gmm.out <- num2str(gmm.out,digits=2)
gmm.out[2,] <- paste0("(", gmm.out[2,], ")")
gmm.out[grep(gmm.out, pattern="NA")] <- NA
gmm.out <- matrix(gmm.out, ncol=4)
gmm.out <- rbind(gmm.out,
                 matrix(paste0("\\multicolumn{1}{c}{",
                        sapply(gmm.list, \(x){x$dim}),
                        "}"), nrow=2))
names <- c("Democracy","", 
           "log GDP pc, lag 1", "",
           "log GDP pc, lag 2", "",
           "log GDP pc, lag 3", "",
           "log GDP pc, lag 4", "",
           "Long-run effect of democracy", "",
           "Effect of democracy after 25 years", "",
           "Persistence", "",
           "Sample size", "Number of instruments")

out <- cbind(names, within.out, gmm.out)
colnames(out) <- paste0("\\multicolumn{1}{c}{",
                        c(" ",
                          "Within (1)", "Within (2)", "Within (3)",
                          "GMM (1)", "GMM (2)","GMM (3)", "GMM (4)"),
                        "}")
print(xtable(out,
             align="rrddddddd",
             caption="Regression results",
             label="tab:regTab"),
      sanitize.text.function=\(x){x},
      booktabs=TRUE,
      include.rownames=FALSE,
      caption.placement="top")



