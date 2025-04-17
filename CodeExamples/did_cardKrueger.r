library(dplyr)
library(tidyr)
library(sandwich)
library(lmtest)
library(fixest)
rm(list=ls())


#### base line 2x2###
ck <- read.csv("datasets/card_krueger_full.csv")


ck <- ck %>% 
  mutate(fte=empft+nmgrs+(0.5*emppt),
         fte2=empft2+nmgrs2+(0.5*emppt2),
         pmeal = psoda+pfry+pentree,
         pmeal2 = psoda2+pfry2+pentree2,
         Dy=fte2-fte,
         restID=1:nrow(.)) %>% 
  filter(!(is.na(fte) | is.na(fte2) | is.na(wage_st) | is.na(wage_st2)
           | is.na(pmeal) | is.na(pmeal2)))
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



##### DiDiD (triple diff)####
ck.twfe <- ck.twfe %>% 
  mutate(h = max(wage_st >= 5 & wave == 0), .by=restID) %>%
  mutate(d2 = 1*(state==1 & wave==1 & h==0),
         NJ.lo = 1*(state==1 &  h==0),
         NJ.hi = 1*(state==1 & h==1),
         PA.lo = 1*(state==0 & h==0),
         PA.hi = 1*(state==0 & h==1),
         NJt = 1*(state==1 & wave==1),
         ht =  1*(h==1 & wave==1))



tripleD <- lm(fte~ d2+ 
                NJ.lo + NJ.hi+
                PA.lo + PA.hi+
                wave+
                ht + 
                NJt-1,
              data=ck.twfe)


NJ.l <- mean(ck.twfe[ck.twfe$state==1 & ck.twfe$wave==0 & ck.twfe$h==0, ]$fte)
## \alpha_\text{NJ,low} =  19.6269  


NJ.h <- mean(ck.twfe[ck.twfe$state==1 & ck.twfe$wave==0 & ck.twfe$h==1, ]$fte)
## \alpha_\text{NJ,high} =  21.443  



NJ.T.B.l <- mean(ck.twfe[ck.twfe$state==1 & ck.twfe$wave==1 & ck.twfe$h==0, ]$fte)
## \beta+ \alpha_\text{NJ,low}+ \alpha'_{NJ,2}+ \tau'
## 0.2541 + 19.6269 +1.9279 + -1.0526

NJ.T.h <- mean(ck.twfe[ck.twfe$state==1 & ck.twfe$wave==1 & ck.twfe$h==1, ]$fte)
## \alpha_\text{NJ,high}+ \alpha'_{NJ,2}+\tau'+\gamma
# 21.4430+ 1.9279 -1.0526 + -2.6648  


PA.l <- mean(ck.twfe[ck.twfe$state==0 & ck.twfe$wave==0 & ck.twfe$h==0, ]$fte)
## \alpha_\text{PA,low} = 23.3289

PA.h <- mean(ck.twfe[ck.twfe$state==0 & ck.twfe$wave==0 & ck.twfe$h==1, ]$fte)
## \alpha_\text{PA,high} =  24.0435

PA.T.l <- mean(ck.twfe[ck.twfe$state==0 & ck.twfe$wave==1 & ck.twfe$h==0, ]$fte) 
## \alpha_\text{PA,low} + \tau' + 
# 23.329+ -1.0526

PA.T.h <- mean(ck.twfe[ck.twfe$state==0 & ck.twfe$wave==1 & ck.twfe$h==1, ]$fte)
## \alpha_\text{PA,high} + \gamma + \tau'
## 24.0435+  -2.6648 -1.0526 


## D1 removes \alpha_{gi,h} but not tau', gamma, or alpha_{NJ,2}
NJ2.T.h <- NJ.T.h -  NJ.h  #-1.0526 +1.9279  -2.6648  
NJ2.T.B.l <- NJ.T.B.l-  NJ.l # -1.0526 +1.9279+0.2541
PA2.T.h <- PA.T.h -  PA.h #  -1.0526   -2.6648
PA2.T.l <- PA.T.l-  PA.l #-1.0526

## D2 removes tau and alpha_{NJ,2}, but not gamma
B.h.l <- NJ2.T.B.l-NJ2.T.h # 0.2541 - -2.6648  
h.l <- PA2.T.l-PA2.T.h #2.6648 

## D3 removes gamma
D3 <- B.h.l - h.l
D3





#### 2x2 with covariates ####
library(did)


ck <- ck %>% 
  mutate(h = ifelse(wage_st>=5, 1, 0))
ck.twfe$h <- rep(ck$h, each=2)

## Recheck the base line##
summary(fd)
twfe <- feols(fte~d|state+wave,data=ck.twfe, vcov=~restID)
summary(twfe)

### the questionable approach
twfe2 <- feols(fte~d+co_owned+h|state+wave,data=ck.twfe, vcov=~restID)

###### A better approach (Regression adjustment)#######
att_gt(yname="fte",
       tname="wave",
       idname="restID",
       gname="state",
       data=ck.twfe,
       est_method="reg",
       xformla = ~co_owned+h)

twfe2


## RA by hand
ra <- lm(Dy~co_owned+h, data=ck, subset=state==0)
summary(ra)
ck$Dra <- predict(ra, newdata=ck)

ck %>% 
  filter(state==1) %>%
  summarize(mean(Dy-Dra))


###### A better approach (IPW)#######

att_gt(yname="fte",
       tname="wave",
       idname="restID",
       gname="state",
       data=ck.twfe,
       est_method="ipw",
       xformla = ~co_owned+h)

twfe2


## IPW by hand


ip <- glm(state~co_owned+h, 
          data=ck,family=binomial("logit"))
summary(ip)
ck$phat <- predict(ip, newdata=ck, type="response")

ck %>% 
  mutate(w1 = state/mean(state),
         w0 = ((1-state)*phat/(1-phat))/mean((1-state)*phat/(1-phat)) ) %>% 
  summarize(mean((w1-w0)*(Dy)))

###### A better approach (Double robust)#######


att_gt(yname="fte",
       tname="wave",
       idname="restID",
       gname="state",
       data=ck.twfe,
       est_method="dr",
       xformla = ~co_owned+h)

twfe2

## DR by hand
ck %>% 
  mutate(w1 = state/mean(state),
         w0 = ((1-state)*phat/(1-phat))/mean((1-state)*phat/(1-phat)) ) %>% 
  summarize(mean((w1-w0)*(Dy-Dra)))
