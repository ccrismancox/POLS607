rm(list=ls())
library(ggplot2)
## A short note on bootstraping 


## The basic bootstrap is justified based on
## the ecdf and how well it matches the CDF
set.seed(1)
y <- rpois(1000, lambda=3)
plot(ecdf(y))


## Because the ECDF is a type of CDF, we can use 
## inverse uniform sampling 
U <- runif(1000)
y2 <- quantile(y, probs=U, type=1)
barplot(table(y2)/1000, col="blue")
points(dpois(seq(0, 20,by=1), lambda=3), col="red", pch=16, cex=2)
lines(dpois(seq(0, 20,by=1), lambda=3), col="red", lwd=2)

## Which is just a convoluted way of 
## using sample with replacement
y3 <- sample(y, size=1000, replace=TRUE)
barplot(table(y3)/1000, col="blue")
points(dpois(seq(0, 20,by=1), lambda=3), col="red", pch=16, cex=2)
lines(dpois(seq(0, 20,by=1), lambda=3), col="red", lwd=2)


####### Bayesian bootstrap ##### 
##  First I want to convince you that the 
## Dirchelet is a reasonable, continuous approx
## for the multinomial distribution
N <- 50
B <- 10000

## draw a bunch of samples. Each row is a sample 
multinomial <- matrix(sample(1:N, size=N*B, replace=TRUE), nrow=B)
m1 <- rowSums(multinomial==1) 

## The continous version of the same thing
gamma11 <-  matrix(rexp(N*B), nrow=B)
dircihlet <- gamma11/rowSums(gamma11)
d1 <- dircihlet[1,]*N 

ggplot()+
  geom_bar(aes(x=m1, y=after_stat(prop)), fill="navyblue")+
  geom_density(aes(x=d1), color="orangered")+
  ggtitle("Distribution of weights on unit 1")+
  xlab("weights")+
  theme_bw(12)+
  ylab("Density/frequency")




## And we can see that either method gives us
## a correct approximation of the variance of the 
## sample mean for y (3/1000) = 0.003

boot <- replicate(50, {
  mean(y[sample(1:1000, size=1000, replace=TRUE)])
})
var(boot)

bboot <- replicate(50, {
  w <- rexp(1000)
  w <- w/sum(w)
  weighted.mean(y, w=w*1000)
})
var(bboot)

