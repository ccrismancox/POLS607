## Packages 
library(dplyr)
library(fixest)
library(bacondecomp)
library(ggplot2)
library(gridExtra)
library(did)

#Clear workspace
rm(list=ls())
## data
use_df <- read.csv("afghanData.csv")

## Cohort-level trends
ggplot(use_df %>% 
         mutate(cohort = ifelse(cohort==10000, NA, cohort)) %>% 
         summarize(Influence = mean(gov_influence, na.rm=TRUE), 
                   .by=c(cohort, year))) +
  geom_line(aes(x=year, y=Influence, color=factor(cohort)))+
  scale_x_continuous(breaks=2008:2014)+
  geom_vline(aes(xintercept=cohort, color=factor(cohort)),
             linetype="dotted")+
  ylab("Government influence")+
  theme_bw()


## Baseline TWFE
twfe <- feols(gov_influence~treated|DISTID+year ,data=use_df)



## Decompose the TWFE
decomp <- bacon(gov_influence ~ treated,
                data = use_df,
                id_var = "DISTID",
                time_var = "year")

decomp


## Individual cohort 2x2s
df.2by2 <- use_df %>%  
  filter(cohort %in% c(2011, 10000)) %>% 
  filter(year ==2010 | year==2011) 

did.2by2.2011 <- feols(gov_influence~treated|cohort +year,
                       data=df.2by2)

df.2by2 <- use_df %>%  
  filter(cohort %in% c(2012, 10000)) %>% 
  filter(year ==2011 | year==2012) 

did.2by2.2012 <- feols(gov_influence~treated|cohort +year,
                       data=df.2by2)

df.2by2 <- use_df %>%  
  filter(cohort %in% c(2013, 10000)) %>% 
  filter(year ==2012 | year==2013) 

did.2by2.2013 <- feols(gov_influence~treated|cohort +year,
                       data=df.2by2)

round(rbind(did.2by2.2011$coeftable,
      did.2by2.2012$coeftable,
      did.2by2.2013$coeftable),2 )


## ATT(g,t) under different PT assumptions
att <- att2 <- att3 <-  matrix(0, 3,7)
i <- 0
for(g in sort(unique(use_df$cohort))[-4]){
  i <- i +1
  j <- 0
  for(t in 2008:2014){
    j <- j+1
    means <- use_df %>% 
      filter(cohort %in% c(10000,g) & year %in% c(t, g-1)) %>%
      summarize(ybar=mean(gov_influence, na.rm=TRUE),
                .by=c(cohort,year))
    att[i,j] <- (means[means$cohort==g & means$year==t,]$ybar-
                   means[means$cohort==g & means$year==g-1,]$ybar)-
      (means[means$cohort==10000 & means$year==t,]$ybar-
         means[means$cohort==10000 & means$year==g-1,]$ybar)
    
    
    means <- use_df %>% 
      filter( (cohort ==g  | cohort > max(g,t)) &
                (year %in% c(t, g-1))) %>% 
      mutate(cohort = ifelse(cohort==g, g, 0))%>% 
      summarize(ybar=mean(gov_influence, na.rm=TRUE),
                .by=c(cohort,year))
    att2[i,j] <- (means[means$cohort==g & means$year==t,]$ybar-
                    means[means$cohort==g & means$year==g-1,]$ybar)-
      (means[means$cohort==0 & means$year==t,]$ybar-
         means[means$cohort==0 & means$year==g-1,]$ybar)
    
    
    
    att3[i,j] <- use_df %>% 
      arrange(DISTID, year) %>%
      mutate(Ig = 1*(cohort==g) + 2*(cohort>t &cohort!=g),
             const=1,
             pre=ifelse(year < g, 1, NA)) %>% 
      mutate(ybar = mean(gov_influence*pre, na.rm=TRUE), .by=DISTID) %>% 
      filter(year ==t & Ig > 0) %>% 
      summarize(Dt = mean(gov_influence - ybar), .by=Ig) %>%
      summarize(att3 = diff(Dt)) %>%
      as.numeric()
  }
}

plot.df <- data.frame(year=2008:2014,
                      group=rep(2011:2013, each=7),
                      ATTgt = c(t(att), t(att2), t(att3)),
                      Assumptions=rep(1:3, each=21))
ggplot(plot.df)+
  geom_line(aes(x=year, y=ATTgt, color=factor(Assumptions)))+
  geom_point(aes(x=year, y=ATTgt, color=factor(Assumptions), shape=factor(Assumptions)))+
  facet_wrap(~group)+
  theme_bw(14)+
  geom_vline(aes(xintercept = group), linetype="dotted")+
  theme(legend.position = "bottom")+  
  scale_x_continuous(breaks=seq(2000, 2010, by=3))



## Sun and Abrams (same as ATT 1)
sa <- feols(gov_influence~sunab(cohort, year)|DISTID+year,data=use_df)
summary(sa)



use_df$time2 <- ifelse(use_df$cohort==10000, NA, use_df$year-use_df$cohort)
dummies <- model.matrix.lm(~factor(time2)-1, data=use_df, na.action="na.pass")
dummies[is.na(dummies)] <- 0 ## this is the interaction for the untreated
colnames(dummies) <- c(paste0("time.m", 5:1), paste0("time.",0:3))
event.df <- cbind(use_df, dummies)
event.fit <- feols(gov_influence ~ time.m5 + time.m4 + time.m3 + time.m2 + 
                     time.0 +  time.1 + time.2 + 
                     time.3 | DISTID+year, data=event.df)

## Event plot "Standard"
eventPlot <- data.frame(time=c(-5:-2, 0:3,-1),
                        effect= c(event.fit$coefficients,0),
                        hi = c(confint(event.fit)[,1],0),
                        lo = c(confint(event.fit)[,2],0))

g1 <- ggplot(eventPlot)+
  geom_pointrange(aes(x=time, y=effect, ymin=lo, ymax=hi))+
  geom_hline(aes(yintercept = 0)) +
  theme_classic(14)+
  xlab("Time")+
  scale_x_continuous(breaks=-4:3)+
  ylab("Marginal effect")


## Same as ATT1 & SA
did_att_gt <- att_gt(yname="gov_influence",
                     tname="year",
                     idname="DISTID",
                     gname="cohort",
                     data=use_df,
                     base_period = "universal",
                     bstrap=FALSE,
                     cband=FALSE)
summary(did_att_gt)
did_es <- aggte(did_att_gt, type="dynamic")
summary(did_es)



## Same as ATT2
did_att_gt2 <- att_gt(yname="gov_influence",
                     tname="year",
                     idname="DISTID",
                     gname="cohort",
                     data=use_df,
                     base_period = "universal",
                     control_group = "notyettreated",
                     bstrap=FALSE,
                     cband=FALSE)
summary(did_att_gt2)
did_es2 <- aggte(did_att_gt2, type="dynamic")
summary(did_es2)








## Wooldridge set up (akin to ATT3)
wooldridge <- use_df %>% 
  mutate(gs = year*(year>=cohort)*(cohort!=10000),
         gs =  ifelse(cohort==10000 | gs==0, NA, gs),
         G= ifelse(cohort==10000, NA, cohort),
         gs= ifelse(is.na(gs), NA, 
                    paste0(G,".", gs)))
wdummies <- model.matrix.lm(~factor(gs)-1, data=wooldridge, na.action="na.pass")
wdummies[is.na(wdummies)] <-0 

GS <- expand.grid(G=2011:2013, s=2011:2014)
GS <- GS[GS$s >= GS$G,]
GS <- GS[order(GS$G),]
colnames(wdummies) <- paste0("g", GS$G, ".s", GS$s)
wooldridge <- cbind(wooldridge, wdummies)
wool.fit <- feols(gov_influence ~ g2011.s2011 + g2011.s2012 + g2011.s2013 + 
        g2011.s2014 + g2012.s2012 + g2012.s2013 + g2012.s2014 +
        g2013.s2013 + g2013.s2014| DISTID+year, data=wooldridge)





## same kind of aggregation based on relative group sizes
weights <- wooldridge %>% 
  summarize(obs=length(gov_influence), .by=c(G,time2)) %>% 
  mutate(size=sum(obs), .by=time2) %>% 
  mutate(weight=obs/size) %>% 
  arrange(G, time2) %>%
  filter(!(time2 < 0| is.na(G))) 


## I'll make a transformation matrix to make the SEs a little 
## easier to compute
woolA <- t(replicate(4, weights$weight))
woolA[1,weights$time2!=0] <- 0
woolA[2,weights$time2!=1] <- 0
woolA[3,weights$time2!=2] <- 0
woolA[4,weights$time2!=3] <- 0

## now the weighted means are just A\hat{\beta}
# and the variance is AV(\hat{\beta}) A'

wool.plot <- data.frame(ATT=woolA %*% wool.fit$coefficients,
                        Time = 0:3,
                        se= sqrt(diag(woolA %*% vcov(wool.fit) %*% t(woolA)))) %>% 
  mutate(lo = ATT-1.96*se,
         hi = ATT+1.96*se)


g.wool <-ggplot(wool.plot)+
  geom_pointrange(aes(x=Time, y=ATT, ymin=lo, ymax=hi))+
  theme_classic(14)+
  xlab("Time since treatment")+
  scale_x_continuous(breaks=-9:5)+
  ylab("Marginal effect")+
  geom_hline(aes(yintercept = 0), linetype="dashed", alpha=.4)



grid.arrange(
  g1 + ggtitle("Classic event study"),
  ggdid(did_es)+ggtitle("Sun & Abraham"),
  ggdid(did_es2)+ggtitle("C&S (PT2)"),
 g.wool+ggtitle("Wooldridge"),
  nrow=2)







