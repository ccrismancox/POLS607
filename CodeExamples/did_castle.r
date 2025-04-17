# library(remotes)
# install_github("johnson-shuffle/mixtape")
library(mixtape) 
library(ggplot2)
library(did)
library(bacondecomp)
library(fixest)
library(dplyr)

rm(list=ls())

data("castle_doctrine_2000_2010")
castle <- castle_doctrine_2000_2010


castle <- castle %>% 
  mutate(across(all_of(c("homicide", "exp_subsidy", "exp_pubwelfare",
                         "police", "income", "prisoner", "lagprisoner")),
                .fns = log,
                .names="l_{col}"),
         post=1*(cdl==1)) %>%
  mutate(trend = 1:length(state), .by=state)


ggplot(castle %>% 
         summarize(Homicide = mean(l_homicide), 
                   .by=c(effyear, year))) +
  geom_line(aes(x=year, y=Homicide, color=factor(effyear)))+
  scale_x_continuous(breaks=2000:2010)+
  geom_vline(aes(xintercept=effyear, color=factor(effyear)),
             linetype="dotted")+
  theme_bw(14)+
  theme(legend.position = "bottom")+
  labs(color = "Treatment year")

### Parallel trends 1
### ATT(g,t) 
castle$G <- castle$effyear
castle$G[is.na(castle$G)] <- 0

att <- matrix(0, 5, 11)
i <- 0
for(g in sort(unique(castle$effyear))){
  i <- i +1
  j <- 0
  for(t in 2000:2010){
    j <- j+1
    means <- castle %>% 
      filter(G %in% c(0,g) & year %in% c(t, g-1)) %>%
      summarize(ybar=mean(l_homicide, na.rm=TRUE),
                .by=c(G,year))
    att[i,j] <- (means[means$G==g & means$year==t,]$ybar-
                   means[means$G==g & means$year==g-1,]$ybar)-
      (means[means$G==0 & means$year==t,]$ybar-
         means[means$G==0 & means$year==g-1,]$ybar)
  }
}

colnames(att) <- paste0("year", 2000:2010)
rownames(att) <- paste0("treated", 2005:2009)
att


### Parallel trends 2
### ATT(g, t) 
att2 <- matrix(0, 5, 11)
i <- 0
for(g in sort(unique(castle$effyear))){
  i <- i +1
  j <- 0
  for(t in 2000:2010){
    j <- j+1
    means <- castle %>% 
      mutate(G = ifelse(G==0, Inf, G)) %>% 
      filter( (G ==g  | G > pmax(g,year)) &
                (year %in% c(t, g-1))) %>% 
      mutate(G = ifelse(G==g, g, 0))%>% 
      summarize(ybar=mean(l_homicide, na.rm=TRUE),
                .by=c(G,year))
    att2[i,j] <- (means[means$G==g & means$year==t,]$ybar-
                    means[means$G==g & means$year==g-1,]$ybar)-
      (means[means$G==0 & means$year==t,]$ybar-
         means[means$G==0 & means$year==g-1,]$ybar)
  }
}

colnames(att2) <- paste0("year", 2000:2010)
rownames(att2) <- paste0("treated", 2005:2009)
att2



### Parallel trends 3
### ATT(g.t) 
att3 <- matrix(0, 5, 11)
i <- 0
for(g in sort(unique(castle$effyear))){
  i <- i +1
  j <- 0
  for(t in 2000:2010){
    j <- j+1
    means <- castle %>% 
      mutate(pre = 1*(year < g),
             t = 1*(year==t),
             G= 1*(G==g))
    preMeans <- means %>% 
      filter(pre==1) %>%
      summarize(ybar=mean(l_homicide), .by=G)
    tMeans <- means %>% 
      filter(t==1) %>% 
      summarize(ybar=mean(l_homicide), .by=G)
    
    att3[i,j] <- (tMeans[tMeans$G==1,]$ybar-
                    preMeans[preMeans$G==1,]$ybar)-
      (tMeans[tMeans$G==0,]$ybar-
         preMeans[preMeans$G==0,]$ybar)
  }
}

colnames(att3) <- paste0("year", 2000:2010)
rownames(att3) <- paste0("treated", 2005:2009)
att3


plot.df <- data.frame(year=2000:2010,
                      group=rep(2005:2009, each=11),
                      ATTgt = c(t(att), t(att2), t(att3)),
                      Assumptions=rep(1:3, each=55))
ggplot(plot.df)+
  geom_line(aes(x=year, y=ATTgt, color=factor(Assumptions)))+
  geom_point(aes(x=year, y=ATTgt, color=factor(Assumptions), shape=factor(Assumptions)))+
  facet_wrap(~group)+
  geom_vline(aes(xintercept = group), linetype="dotted")+
  scale_x_continuous(breaks=seq(2000, 2010, by=3))+
  ylab("ATT(g,t)")+
  theme_bw(14)+
  theme(legend.position = "bottom")+
  labs(color = "PT Assumption", shape="PT Assumption")




library(did)
castle <- castle %>% 
  mutate(ifelse(G==0, Inf, G))

## matches att_1(g,t)
att_gt("l_homicide",
       tname = "year",
       idname = "sid",
       gname = "G", 
       control_group="nevertreated",
       base_period = "universal",
       data = castle)

## matches att_2(g,t)
att_gt("l_homicide",
       tname = "year",
       idname = "sid",
       gname = "G", 
       control_group="notyettreated",
       base_period = "universal",
       data = castle)


castle$time2 <- castle$year-castle$effyear
gdummies <- model.matrix.lm(~factor(G):factor(time2)-1, 
                            data=castle,
                            na.action="na.pass")
gdummies[is.na(gdummies)] <- 0
gdummies <- gdummies[,grep(pattern="factor\\(G\\)0", 
                           x=colnames(gdummies), invert = TRUE)]
gdummies <- gdummies[,grep(pattern="factor\\(time2\\)-1", 
                           x=colnames(gdummies), invert = TRUE)]
fe <- model.matrix(~factor(state)+factor(year)-1, data=castle)
SAx <- cbind(gdummies, fe)
SAx <- SAx[,colSums(SAx)!=0]
qr(crossprod(SAx))$rank
SA.hat <- solve(crossprod(SAx)) %*% t(SAx) %*% castle$l_homicide
beta.hat <- SA.hat[1:50,]
SA.hat <- rbind(beta.hat[grep(pattern="factor\\(G\\)2005", names(beta.hat))],
                beta.hat[grep(pattern="factor\\(G\\)2006", names(beta.hat))],
                beta.hat[grep(pattern="factor\\(G\\)2007", names(beta.hat))],
                beta.hat[grep(pattern="factor\\(G\\)2008", names(beta.hat))],
                beta.hat[grep(pattern="factor\\(G\\)2009", names(beta.hat))])


## OR with the fixest package
castle$effyear2 <- castle$effyear
castle$effyear2[is.na(castle$effyear2)] <- 10000 #set the treatment date for
## the untreated to something out of sample
sa <- feols(l_homicide~sunab(effyear2, year)|state+year,data=castle)
summary(sa)


sa$coefficients #same as beta.hat



att.sa <- sapply(c(-9:-2, 0:5), 
                 \(ell){beta.hat[grep(pattern=paste0("factor\\(time2\\)", ell),
                                      names(beta.hat))]})
weights <- castle %>% 
  summarize(obs=length(l_homicide), .by=c(G, time2)) %>% 
  mutate(size=sum(obs), .by=time2) %>% 
  mutate(weight=obs/size) %>% 
  arrange(time2, G) %>%
  filter(!(time2==-1 | G==0)) 
weights <- split(weights$weight, weights$time2)
mapply(\(x,w)weighted.mean(x,w), x=att.sa, w=weights)









#Fixed effect regression using post as treatment variable 
dd_reg1 <- feols(l_homicide ~ post|sid + year, data = castle)
summary(dd_reg1)
decomp <- bacon(l_homicide ~ post,
                data = castle,
                id_var = "state",
                time_var = "year")

## Early v. late: Maybe good
## Late v. early: Maybe bad? maybe ok?
## Treated v. untreated: good!


## event study version 
dummies <- model.matrix.lm(~factor(time2)-1, data=castle, na.action="na.pass")
dummies[is.na(dummies)] <- 0 ## this is the interaction for the untreated
colnames(dummies) <- c(paste0("time.m", 9:1), paste0("time.",0:5))
castle <- cbind(castle, dummies)
event.fit <- feols(l_homicide ~ time.m9 + time.m8 + time.m7 + 
                     time.m6 + time.m5 + time.m4 + time.m3 + 
                     time.m2 + time.0 +  time.1 + time.2 + 
                     time.3 + time.4 + time.5
                   | state+year, data=castle)


eventPlot <- data.frame(time=c(-9:-1, 1:5,0),
                        effect= c(event.fit$coefficients[1:14],0),
                        hi = c(confint(event.fit)[1:14,1],0),
                        lo = c(confint(event.fit)[1:14,2],0))

g1 <- ggplot(eventPlot)+
  geom_pointrange(aes(x=time, y=effect, ymin=lo, ymax=hi))+
  geom_hline(aes(yintercept = 0)) +
  theme_bw(14)+
  xlab("Time")+
  scale_x_continuous(breaks=-9:5)+
  ylab("Marginal effect")+
  geom_vline(aes(xintercept = 0), linetype="dashed", alpha=.4)
g1

eventPlot.sa <- data.frame(time=c(-9:-1, 1:5,0),
                           effect= c(sa$coeftable[,1],0),
                           hi = c(confint(sa)[1:14,1],0),
                           lo = c(confint(sa)[1:14,2],0))

g2 <- ggplot(eventPlot.sa)+
  geom_pointrange(aes(x=time, y=effect, ymin=lo, ymax=hi))+
  geom_hline(aes(yintercept = 0)) +
  theme_bw(14)+
  xlab("Time")+
  scale_x_continuous(breaks=-9:5)+
  ylab("Marginal effect")+
  geom_vline(aes(xintercept = 0), linetype="dashed", alpha=.4)
g2







