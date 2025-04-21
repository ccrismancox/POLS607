library(mixtape)
library(ggplot2)
library(did)
library(bacondecomp)
library(fixest)
library(dplyr)
library(gridExtra)
rm(list=ls())

data("castle_doctrine_2000_2010")
castle <- castle_doctrine_2000_2010


castle <- castle %>%
  mutate(across(all_of(c("homicide", 
                         "exp_subsidy",
                         "exp_pubwelfare",
                         "police", "income",
                         "prisoner", "lagprisoner")),
                .fns = log,
                .names="l_{col}"),
         post=1*(cdl==1)) %>%
  mutate(trend = 1:length(state), .by=state)


ctrl <-  ~ poverty+
  unemployrt+
  l_income+
  l_exp_pubwelfare


castle$G <- castle$effyear
castle$G[is.na(castle$G)] <- Inf


UNatt <-RAatt <- IPWatt <- DRatt <-
  NAtracker<- matrix(0, 5, 11)
i <- 0
for(g in sort(unique(castle$effyear))){
  i <- i +1
  j <- 0
  for(t in 2000:2010){
    j <- j+1
    means <- castle %>%
      filter( (G ==g  | G > max(g,t)) &
                (year %in% c(t, g-1))) %>%
      mutate(G = ifelse(G==g, g, 0))%>%
      summarize(ybar=mean(l_homicide, na.rm=TRUE),
                .by=c(G,year))
    UNatt[i,j] <- (means[means$G==g & means$year==t,]$ybar-
                     means[means$G==g & means$year==g-1,]$ybar)-
      (means[means$G==0 & means$year==t,]$ybar-
         means[means$G==0 & means$year==g-1,]$ybar)
    
    
    Z <- castle %>%
      mutate(treated=1*(G==g))%>%
      filter( (G ==g  | G > max(g,t)) &
                (year == min(t,  g-1))) %>%
      select(treated,  state, year, G,
             poverty,
             unemployrt,
             l_income,
             l_exp_pubwelfare)
    Z.Dy <- castle%>%
      mutate(treated=1*(G==g))%>%
      filter( (G ==g  | G > max(g,t)) &
                (year %in% c(t,g-1))) %>%
      select(state,year,l_homicide) 
    Z$Dy <- Z.Dy$l_homicide[Z.Dy$year==t]-
      Z.Dy$l_homicide[Z.Dy$year==g-1]
    
    
    ra <- lm(update(ctrl, Dy~.),
             data=Z, subset=treated==0)
    Z$Dra <- predict(ra, newdata=Z)
    
    RAatt[i,j]  <- Z %>%
      filter(treated==1) %>%
      summarize(mean(Dy-Dra)) %>%
      as.numeric()
    
    
    ip <- glm(update(ctrl, treated~.),
              data=Z,family=binomial("logit"))
    
    
    Z$phat <- predict(ip,
                      newdata=Z,
                      type="response")
    ## note that did returns an NA if we get 
    NAtracker[i,j] <- ifelse(any(Z$phat>0.999),
                             NA, 1)
    ##
    
    IPWatt[i,j] <- Z %>%
      mutate(w1 = treated/mean(treated),
             w0 = ((1-treated)*phat/(1-phat))/mean((1-treated)*phat/(1-phat)) ) %>%
      summarize(mean((w1-w0)*(Dy))) %>%
      as.numeric()
    
    DRatt[i,j] <- Z %>%
      mutate(w1 = treated/mean(treated),
             w0 = ((1-treated)*phat/(1-phat))/mean((1-treated)*phat/(1-phat)) ) %>%
      summarize(mean((w1-w0)*(Dy-Dra))) %>%
      as.numeric()
    
  }
  
  
  
  
}

colnames(UNatt) <- colnames(RAatt) <- 
  colnames(IPWatt) <- colnames(DRatt) <- paste0("year", 2000:2010)
rownames(UNatt) <- rownames(RAatt) <- 
  rownames(IPWatt) <- rownames(DRatt) <- paste0("treated", 2005:2009)


UNatt
RAatt

IPWatt
DRatt



plot.df <- data.frame(year=2000:2010,
                      group=rep(2005:2009, each=11),
                      ATTgt = c(t(UNatt), t(RAatt), t(IPWatt),t(DRatt) ),
                      Strategy=rep(c("No controls", "RA", "IPW", "DR"), each=55))
ggplot(plot.df)+
  geom_line(aes(x=year, y=ATTgt, color=Strategy))+
  geom_point(aes(x=year, y=ATTgt, color=Strategy, shape=Strategy))+
  facet_wrap(~group)+
  geom_vline(aes(xintercept = group), linetype="dotted")+
  scale_x_continuous(breaks=seq(2000, 2010, by=3))+
  ylab("ATT(g,t)")+
  theme_bw(14)+
  theme(legend.position = "bottom")

att.names <- expand.grid(2005:2009,2000:2010)%>%
  arrange(Var1,Var2) %>% 
  mutate(name=paste0(Var1, ".", Var2)) %>% 
  select(name) %>% 
  unlist()

cannedRA <- att_gt("l_homicide",
                   tname = "year",
                   idname = "sid",
                   gname = "G", 
                   xformla = ~ poverty+
                     unemployrt+
                     l_income+
                     l_exp_pubwelfare,
                   control_group="notyettreated",
                   base_period = "universal",
                   est_method = "reg",
                   data = castle)

RA.est <- cbind(cannedRA$att, 
                c(t(RAatt)))
rownames(RA.est) <- att.names
print(RA.est)

cannedIPW <- att_gt("l_homicide",
                    tname = "year",
                    idname = "sid",
                    gname = "G", 
                    xformla = ~ poverty+
                      unemployrt+
                      l_income+
                      l_exp_pubwelfare,
                    control_group="notyettreated",
                    base_period = "universal",
                    est_method = "ipw",
                    data = castle)


IPW.est <- cbind(cannedIPW$att,
                 c(t(IPWatt)),
                 c(t(IPWatt*NAtracker)))
rownames(IPW.est) <- att.names
print(IPW.est)


cannedDR <- att_gt("l_homicide",
                   tname = "year",
                   idname = "sid",
                   gname = "G", 
                   xformla = ~ poverty+
                     unemployrt+
                     l_income+
                     l_exp_pubwelfare,
                   control_group="notyettreated",
                   base_period = "universal",
                   est_method = "dr",
                   data = castle)


DRest <- cbind(cannedDR$att,
               c(t(DRatt)),
               c(t(DRatt*NAtracker)))
row.names(DRest) <- att.names
print(DRest)

g1 <- ggdid(aggte(cannedDR, "dynamic",na.rm=TRUE))+ggtitle("Double robust")
g2 <- ggdid(aggte(cannedRA, "dynamic",na.rm=TRUE))+ggtitle("Regression adjusted")
g3 <- ggdid(aggte(cannedIPW, "dynamic",na.rm=TRUE))+ggtitle("IPW")


grid.arrange(g1,g2,g3, nrow=3)

weights <- castle %>% 
  mutate(time2=year-G) %>% 
  summarize(obs=length(l_homicide), .by=c(G, time2)) %>% 
  mutate(size=sum(obs), .by=time2) %>% 
  mutate(weight=obs/size) %>% 
  arrange(G,time2) %>%
  filter(!is.infinite(G))
weights <- split(weights$weight, weights$time2)


ests <- split(cannedRA$att, factor(cannedRA$t-cannedRA$group))
cbind(mapply(weighted.mean,
             x=ests,
             w=weights),
      aggte(cannedRA, "dynamic",na.rm=TRUE)$att.egt)
