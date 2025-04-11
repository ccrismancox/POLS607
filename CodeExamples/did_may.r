library(fixest)
library(ggplot2)
library(dplyr)
rm(list=ls())

mayData <- read.csv("datasets/did_May.csv")
table(mayData$treatment1, mayData$wave) ##three waves

mayData$post <- mayData$treatment1 *( mayData$wave==12)
mayData$time <- mayData$wave -12 ## periods -2, -1, 0

table(mayData$post, mayData$treatment1)
table(mayData$post, mayData$treatment1, mayData$time)

mayData$treatment.m2 <- mayData$treatment1*(mayData$time==-2)
mayData$treatment.0 <- mayData$treatment1*(mayData$time==0)

trendDat <- mayData %>% 
  group_by(time, treatment1) %>%
  summarize(likeMayW=mean(likeMayW)) %>%
  mutate(Treatment=ifelse(treatment1==0, "Untreated", "Treated")) %>%
  ungroup()%>%
  arrange(Treatment, time)

pre.slopes <- trendDat %>% 
  group_by(Treatment) %>% 
  filter(time < 0) %>% 
  summarize(slope=diff(likeMayW))
pre.slopes

## plotting just the trend in outcomes over time
## does it look parallel pre-treatment? Hard to say
ggplot(trendDat, aes(x=time,y=likeMayW, color=Treatment)) +
  geom_point()+
  geom_line()+
  theme_bw(14)+
  theme(legend.position = "bottom")+
  ylab("Approves May") +
  xlab("Time")+
  geom_vline(aes(xintercept = -1), linetype="dashed", alpha=.2)


mayData$female <- 1*(mayData$gender==1)

Xdemo <- model.matrix(~white+female+conservative-1, data=mayData)
Xbig <- model.matrix(~post+factor(treatment1)+factor(wave)-1,
                     data=mayData)
e1 <- Xdemo  - Xbig %*% (solve(t(Xbig) %*% Xbig) %*% t(Xbig) %*% Xdemo)
y <- mayData$likeMayW
e2 <- y  - e1 %*% (solve(t(e1) %*% e1) %*% t(e1) %*% y)

cond.Dat <- data.frame(treatment1=Xbig[,3],
                       time= (Xbig[,4] + 2*Xbig[,5])-2,
                       likeMayW=e2)


trendDat2 <- cond.Dat %>% 
  group_by(time, treatment1) %>%
  summarize(likeMayW=mean(likeMayW)) %>%
  mutate(Treatment=ifelse(treatment1==0, "Untreated", "Treated")) %>%
  ungroup()%>%
  arrange(Treatment, time)

pre.slopes2 <- trendDat2 %>% 
  group_by(Treatment) %>% 
  filter(time < 0) %>% 
  summarize(slope=diff(likeMayW))
pre.slopes2 ## about the same gap, just flipped

## plotting just the trend in outcomes over time
## does it look parallel pre-treatment? Hard to say
ggplot(trendDat2, aes(x=time,y=likeMayW, color=Treatment)) +
  geom_point()+
  geom_line()+
  theme_bw(14)+
  theme(legend.position = "bottom")+
  ylab("Approves May") +
  xlab("Time")+
  geom_vline(aes(xintercept = -1), linetype="dashed", alpha=.2)





event <- feols(likeMayW~ treatment.m1+treatment.1|treatment1+wave, 
               data=mayData, vcov=~id)
summary(event)

eventPlot <- data.frame(time=c(-2,0,-1),
                        effect= c(event$coefficients,0),
                        hi = c(confint(event)[,1],0),
                        lo = c(confint(event)[,2],0))

ggplot(eventPlot)+
  geom_pointrange(aes(x=time, y=effect, ymin=lo, ymax=hi))+
  geom_hline(aes(yintercept = 0))+  
  theme(legend.position = "bottom")+
  ylab("Marginal effect on May's approval") +
  xlab("Time")

## with individual level controls
event2 <- feols(likeMayW~ treatment.m1+treatment.1|id+wave, 
                data=mayData, vcov=~id)
summary(event2)

eventPlot <- data.frame(time=c(-2,0,-1),
                        effect= c(event2$coefficients,0),
                        hi = c(confint(event2)[,1],0),
                        lo = c(confint(event2)[,2],0))

ggplot(eventPlot)+
  geom_pointrange(aes(x=time, y=effect, ymin=lo, ymax=hi))+
  geom_hline(aes(yintercept = 0))+  
  theme(legend.position = "bottom")+
  ylab("Marginal effect on May's approval") +
  xlab("Time")

did <- feols(likeMayW~post|treatment1+wave, data=mayData, vcov=~id)
did2 <- feols(likeMayW~post|id+wave, data=mayData, vcov=~id)

summary(did)
summary(did2)