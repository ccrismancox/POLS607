library(fixest)
library(ggplot2)
library(dplyr)
rm(list=ls())
mayData <- read.csv("did_May.csv")
table(mayData$treatment1, mayData$wave)

mayData$post <- mayData$treatment1 *( mayData$wave==12)
mayData$time <- mayData$wave -11

table(mayData$post, mayData$treatment1)
table(mayData$post, mayData$treatment1, mayData$time)

mayData$treatment.m1 <- mayData$treatment1*(mayData$time==-1)
mayData$treatment.1 <- mayData$treatment1*(mayData$time==1)


groups <- mayData %>% 
  group_by(treatment1, time)  %>% 
  summarize(Approval=mean(likeMayW)) %>% 
  mutate(Group=ifelse(treatment1==1, "Treated", "Untreated"))

ggplot(groups)+
  geom_line(aes(x=time, y=Approval, color=Group))+
  theme_bw(14)+
  xlab("Time")+
  geom_vline(aes(xintercept = 0), linetype="dashed", alpha=.4)
theme(legend.position = "bottom")


event <- feols(likeMayW~ treatment.m1+treatment.1|treatment1+wave, 
               data=mayData, vcov="hetero")
summary(event)

eventPlot <- data.frame(time=c(-1,1,0),
                        effect= c(event$coefficients,0),
                        hi = c(confint(event)[,1],0),
                        lo = c(confint(event)[,2],0))

ggplot(eventPlot)+
  geom_pointrange(aes(x=time, y=effect, ymin=lo, ymax=hi))+
  geom_hline(aes(yintercept = 0)) +
  theme_bw(14)+
  xlab("Time")+
  ylab("Marginal effect")+
  geom_vline(aes(xintercept = 0), linetype="dashed", alpha=.4)


did <- feols(likeMayW~post|treatment1+wave+id, data=mayData, vcov=~id)

summary(did)