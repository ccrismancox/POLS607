library(dplyr)
library(tidyr)
library(sandwich)
library(lmtest)
library(fixest)
rm(list=ls())
ck <- read.csv("Rcode/datasets/card_krueger_full.csv")


ck <- ck %>% 
  mutate(fte=empft+nmgrs+(0.5*emppt),
         fte2=empft2+nmgrs2+(0.5*emppt2),
         pmeal = psoda+pfry+pentree,
         pmeal2 = psoda2+pfry2+pentree2,
         Dy=fte2-fte,
         restID=1:nrow(.)) %>% 
  filter(!(is.na(fte) | is.na(fte2) | is.na(wage_st) | is.na(wage_st2)
           is.na(pmeal) | is.na(pmeal2)))
dim(ck)



fd <- lm(Dy~state, data=ck)
summary(fd)
t.test(Dy~factor(1-state), data=ck)
diff(t.test(Dy~factor(1-state), data=ck)$estimate)


did <- mean(ck$Dy[ck$state==1]) -mean(ck$Dy[ck$state==0]) 
v0 <- var(ck$Dy[ck$state==0])/sum(ck$state==0)
v1 <- var(ck$Dy[ck$state==1])/sum(ck$state==1)

did/sqrt(v1+v0)
coeftest(fd, vcovCL(fd, cluster=~restID))



## What do the chain coefficients represent here?
fd2<- lm(Dy~state+factor(chain)-1, data=ck)
coeftest(fd2, vcovCL(fd2, cluster=~restID))




## what if we did it as twfe
ck.twfe <- ck %>% 
  select(c(state, restID, 
           co_owned, chain, 
           fte, wage_st, nmgrs,hrsopen,pmeal,
           fte2, wage_st2,nmgrs2,hrsopen2,pmeal2))%>%
  pivot_longer(cols=c(fte:pmeal2),
               names_to = c(".value"),
               names_pattern = "([^0-9]*)",
               cols_vary = "slowest") %>%
  mutate(wave=rep(0:1, each=nrow(ck)),
         d = as.numeric(wave==1 & state==1)) %>% 
  arrange(restID, wave)

twfe1 <- lm(fte~d + factor(state) + factor(wave)-1, data=ck.twfe)
twfe2 <- feols(fte~d|state+wave,data=ck.twfe, vcov=~restID)

coeftest(twfe1, vcovCL(twfe1, cluster=~restID))
summary(twfe2)


twfe3 <- feols(fte~d|restID+wave,
               data=ck.twfe)
summary(twfe3)


## placebos?
twfe4<- feols(nmgrs~d|state+wave,
              data=ck.twfe, vcov=~restID)

twfe5<- feols(hrsopen~d|state+wave,
              data=ck.twfe, vcov=~restID)
summary(twfe4)
summary(twfe5)