library(fixest)
library(car)
library(ggplot2)
rm(list=ls())
## Start with just using the 2x2 that is the transition from 2013 to 2014
medicaid <- read.csv("datasets/medicaid_expansion.csv")
medicaid$G <- ifelse(is.na(medicaid$Implemented), 0, medicaid$Implemented)

df.2by2 <- medicaid %>% 
  filter(G <= 2014) %>% 
  filter(year ==2013 | year==2014) %>% 
  mutate(treated = ifelse(is.na(Implemented),  0, 1*(Implemented==2014)),
         post=treated*(year==2014))
did.2by2 <- lm(death.rate~post+factor(treated)+factor(year)-1, data=df.2by2)
summary(did.2by2)



## 2 by T
ggplot(medicaid %>% 
         filter(G<=2014) %>% 
         summarize(death.rate=mean(death.rate,na.rm=TRUE), .by=c(G,year)))+
  geom_line(aes(x=year, y=death.rate, color=factor(G)))+
  theme_bw(14)+
  theme(legend.position = "bottom")+
  ylab("Infant mortality rate") +
  xlab("year")+
  labs(color = "Treatment group")



### ATT(t) 
att <- rep(0, 14)
i <- 0
for(t in 2009:2022){
  i <- i+1
  means <- medicaid %>% 
    filter(G <= 2014) %>% 
    filter(year %in% c(2013, t)) %>% 
    summarize(ybar=mean(death.rate, na.rm=TRUE),
              .by=c(G,year))
  att[i] <- (means[means$G==2014 & means$year==t,]$ybar-
               means[means$G==2014 & means$year==2013,]$ybar)-
    (means[means$G==0 & means$year==t,]$ybar-
       means[means$G==0 & means$year==2013,]$ybar)
}


did.2byT <- medicaid %>% 
  filter(G<=2014) %>% 
  mutate(time2 = ifelse(is.na(Implemented),NA, (year-2013)))
dummies <- model.matrix.lm(~factor(time2)-1, data=did.2byT, na.action=na.pass)
dummies[is.na(dummies)] <- 0
colnames(dummies) <- c(paste0("time.m", 4:1),paste0("time.", 0:9))
did.2byT <- cbind(did.2byT,dummies)
event <- feols(death.rate~time.m4+time.m3+time.m2+time.m1+
                 time.1+time.2+time.3+time.4+time.5+time.6+
                 time.7+time.8+time.9|state+year,data=did.2byT) 



eventPlot <- data.frame(time=c(-4:-1, 1:9,0),
                        effect= c(event$coefficients,0),
                        hi = c(confint(event)[,1],0),
                        lo = c(confint(event)[,2],0))

ggplot(eventPlot)+
  geom_pointrange(aes(x=time, y=effect, ymin=lo, ymax=hi))+
  geom_hline(aes(yintercept = 0)) +
  theme_bw(14)+
  xlab("Time")+
  scale_x_continuous(breaks=-4:9)+
  ylab("Marginal effect")+
  geom_vline(aes(xintercept = 0), linetype="dashed", alpha=.4)


### A packaged version  with standard error correction ## 
library(did)
did.2byT$stateCode <- as.numeric(as.factor(did.2byT$state))
out1 <- att_gt(yname="death.rate",
               tname="year",
               idname="stateCode",
               gname="G",
               base_period="universal",
               data=did.2byT)
summary(out1)



cbind(att, out1$att, c(event$coefficients[1:4],0,event$coefficients[5:13]))




did.2byT <- did.2byT %>% 
  mutate(post=1*(G==2014 & year>=2014))
twfe <- feols(death.rate~post|state+year,data=did.2byT) 
summary(twfe)

##versus 
mean(att[6:14])
## or 
Att.combo <- deltaMethod(event,
                         "(time.1+time.2+time.3+time.4+time.5+time.6+
        time.7+time.8+time.9)/9")
c(Att.combo[1:2])