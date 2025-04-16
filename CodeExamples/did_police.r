library(readstata13)
library(did)
library(bacondecomp)
library(fixest)
library(dplyr)
library(tidyr)
library(gridExtra)
library(ggplot2)

rm(list=ls())

policing <- read.dta13("datasets/underReport.dta")
p2017 <- read.csv("datasets/statePop2017.csv")
policing[policing$year==2017,]$pop_100000 <- p2017$pop2017/100000
policing$ln_police_homicide <- with(policing,log(police_homicide_fatal/pop_100000+1))

## First any issue with just using the full data staggered?
decomp <- bacon(ln_police_homicide ~ invLaw,
                data = policing[!is.na(policing$ln_police_homicide),],
                id_var = "state",
                time_var = "year")
print(decomp)

policing <- policing %>%
  mutate(G= min(ifelse(invLaw==0, Inf, pmin(year*invLaw))), .by=state) 

policing %>% 
  with(., table(G, year))

## Let's consider the ATT(t) for the 2015 group with controls for
policing15 <- policing %>% 
  filter( G >= 2015) %>% #last pre treated and no already treated
  mutate(treated = 1*(G==2015))
table(policing15$treated)



## let's start by looking at balance among observables
policing15  %>% 
  filter(year==2014) %>%
  select(treated, 
         std_squire_capacity, ## state "capacity"
         termEnact, ## term limits
         std_med_income, ## income
         std_employees, ## employees
         std_police_ideal, ## police ideal point
         std_black, ## Black population
         demGov, #Dem gov
         demLeg, #unified dem leg
         repLeg) %>% #unified gop leg
  reframe(across(everything(),
                 .fns=list(mean, sd)),
          .by=c(treated)) %>%
  pivot_longer(!c(treated)) %>% 
  pivot_wider(names_from=c(treated),
              names_prefix = "Group") %>% 
  mutate(Stat=rep(c("mean", "sd"), 9),
         var=rep(1:9, each=2),
         D=ifelse(Stat=="mean", 
                  Group1-Group0, 
                  sqrt((Group1^2+Group0^2)/2))) %>% 
  mutate(Norm = D[1]/D[2], .by=var) %>% 
  filter(Stat=="mean") %>% 
  select(!c(var, D, Stat))


### ATT(t) 
UNatt <-RAatt <- IPWatt <- DRatt <-rep(0, 5)

i <- 0
for(t in 2013:2017){
  i <- i+1
  means <- policing15 %>% 
    filter(year %in% c(2014, t)) %>% 
    summarize(ybar=mean(ln_police_homicide, na.rm=TRUE),
              .by=c(treated,year)) 
  UNatt[i] <- (means[means$treated==1 & means$year==t,]$ybar-
                 means[means$treated==1 & means$year==2014,]$ybar)-
    (means[means$treated==0 & means$year==t,]$ybar-
       means[means$treated==0 & means$year==2014,]$ybar)
  
  
  Z <- policing15 %>% 
    filter(year ==2014) %>%
    select(treated, 
           std_squire_capacity, ## state "capacity"
           termEnact, ## term limits
           std_med_income, ## income
           std_employees, ## employees
           std_police_ideal, ## police ideal point
           std_black, ## Black population
           demGov, #Dem gov
           demLeg, #unified dem leg
           repLeg)
  Z$Dy <- policing15$ln_police_homicide[policing15$year==t]-
    policing15$ln_police_homicide[policing15$year==2014]
  
  ra <- lm(Dy~std_squire_capacity+termEnact+std_med_income+
             std_employees+std_police_ideal+std_black+
             demGov+demLeg+repLeg, data=Z, subset=treated==0)
  Z$Dra <- predict(ra, newdata=Z)
  
  RAatt[i]  <- Z %>% 
    filter(treated==1) %>%
    summarize(mean(Dy-Dra)) %>% 
    as.numeric()
  
  
  
  
  ip <- glm(treated~std_squire_capacity+termEnact+std_med_income+
              std_employees+std_police_ideal+std_black+
              demGov+demLeg+repLeg, data=Z,family=binomial("logit"))
  summary(ip)
  Z$phat <- predict(ip, 
                    newdata=Z, 
                    type="response")
  
  IPWatt[i] <- Z %>% 
    mutate(w1 = treated/mean(treated),
           w0 = ((1-treated)*phat/(1-phat))/mean((1-treated)*phat/(1-phat)) ) %>% 
    summarize(mean((w1-w0)*(Dy))) %>% 
    as.numeric()
  
  DRatt[i] <- Z %>% 
    mutate(w1 = treated/mean(treated),
           w0 = ((1-treated)*phat/(1-phat))/mean((1-treated)*phat/(1-phat)) ) %>% 
    summarize(mean((w1-w0)*(Dy-Dra))) %>% 
    as.numeric()
  
}


cbind(2013:2017,
      UNatt, RAatt, IPWatt, DRatt)

policing15 <- policing15 %>% 
  mutate(time.m2 = ifelse(treated==1&year==2013, 1,0),,
         time0 = ifelse(treated==1&year==2015, 1,0),
         time1 = ifelse(treated==1&year==2016, 1,0),
         time2 = ifelse(treated==1&year==2017, 1,0))
event.fit <-feols(ln_police_homicide~time.m2+time0+time1+time2
                  |state+year,
                  data=policing15)
event.fit ## match UNatt

event.fit2 <-feols(ln_police_homicide~time.m2+time0+time1+time2+
                     std_squire_capacity+termEnact+
                     std_med_income+
                     std_employees+
                     std_police_ideal+std_black+
                     demGov+demLeg+repLeg
                   |state+year,
                   data=policing15)
event.fit2 ## matches nothing!


policing15$stateCode <- as.numeric(as.factor(policing15$state))
policing15$treated <- 2015*(policing15$G==2015)
unAtt <- att_gt(yname = "ln_police_homicide",
                tname="year",
                idname="stateCode",
                gname="treated",
                base_period = "universal",
                data=policing15)
attRA <- att_gt(yname = "ln_police_homicide",
                tname="year",
                idname="stateCode",
                base_period = "universal",
                gname="treated",
                est_method="reg",
                xformla=~std_squire_capacity+termEnact+
                  std_med_income+
                  std_employees+
                  std_police_ideal+std_black+
                  demGov+demLeg+repLeg,
                data=policing15)

attIPW <- att_gt(yname = "ln_police_homicide",
                 tname="year",
                 idname="stateCode",
                 gname="treated",
                 est_method="ipw",
                 base_period = "universal",
                 xformla=~std_squire_capacity+termEnact+
                   std_med_income+
                   std_employees+
                   std_police_ideal+std_black+
                   demGov+demLeg+repLeg,
                 data=policing15)

attDR <- att_gt(yname = "ln_police_homicide",
                tname="year",
                idname="stateCode",
                gname="treated",
                est_method="dr",
                base_period = "universal",
                xformla=~std_squire_capacity+termEnact+
                  std_med_income+
                  std_employees+
                  std_police_ideal+std_black+
                  demGov+demLeg+repLeg,
                data=policing15)

grid.arrange(
  ggdid(unAtt)+ggtitle("No controls"),
  ggdid(attRA)+ggtitle("Reg. Adjust."),
  ggdid(attIPW)+ggtitle("IPW"),
  ggdid(attDR)+ggtitle("Double robust"),
  nrow=2)