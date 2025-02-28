library(dplyr)
library(readstata13)

library(sandwich)
library(fixest)
library(ivreg)
library(lmtest)

library(modelsummary)
rm(list=ls())
aid <- read.dta13("datasets/finaldata.dta")

### setup in the paper###
aid  <- aid %>% 
  mutate(lcommit3a=log(1000000*commit3a+1),
           lpopulation =log(1000*population),
           listock = log(istock+1),
           lgdpcap = log(gdpcap+1),
           lexports = log(exports+1),
           ldist = log(distance+1),
           lusmil = log(usmil+1),
           ldisaster = log(disaster +1)) %>% 
  mutate(lpopulation_lag= lag(lpopulation), 
           listock_lag=lag(listock),
           lgdpcap_lag=lag(lgdpcap),
           lexports_lag = lag(lexports),
           lusmil_lag = lag(lusmil),
           fh_lag = lag(fh),
           civilwar_lag = lag(civilwar),
           ldisaster_lag = lag(ldisaster),
         .by=dyad) %>%
  filter(year > 1992 & year < 2009)

## Main outcome
# lcommit3a: foreign aid commitments from donor to receipient (USD log)

## Main regressors
# listock_lag: Size of the migrant population from the recipient country 
#               in the donor (log, lag)
# lgdpcap_lag: Recipient GDP per capita (USD/person log, lag)
# lpopulation_lag: Recipient population (log, lag)
# lexports_lag: Exports from donor to the recipient  (USD log, lag)
# ldist: Distance from donor to recipient
# colony: Recipient is a former colony of the donor
# lusmil_lag: US military aid (log lag)
# fh_lag: 1-7 measure of democracy (lag)
# civilwar_lag: Binary, is there civil war (lag)
# ldisaster_lag: Number of people affected by a natural disaster (log, lag)



m1 <- feols(lcommit3a~listock_lag+lgdpcap_lag+lpopulation_lag+ lexports_lag +
              ldist+ colony+ lusmil_lag+ fh_lag+ civilwar_lag+ldisaster_lag|
              donor+year, data=aid)
summary(m1)

m2 <- feols(lcommit3a~listock_lag+lgdpcap_lag+lpopulation_lag+ lexports_lag +
              ldist+ colony+ lusmil_lag+ fh_lag+ civilwar_lag+ldisaster_lag|
              dyad+year, data=aid)
summary(m2)

aid.sam <- aid[m2$obs_selection$obsRemoved,]


Years <- model.matrix(~factor(year)-1, data=aid.sam)[,-1]
colnames(Years) <- paste0("year", 1994:2008)
aid.sam <- cbind(aid.sam, Years)


var.names<- c("listock_lag", "lgdpcap_lag", "lpopulation_lag", 
              "lexports_lag", "lusmil_lag", "fh_lag", "civilwar_lag",
              "ldisaster_lag",paste0("year", 1994:2008) )
aid.sam <- aid.sam %>%
  group_by(dyad) %>%
  mutate(across(all_of(var.names), \(x){x-mean(x)},
                .names= "{col}.within")) %>% 
  ungroup()



fx <- ~listock_lag+lgdpcap_lag+lpopulation_lag+ lexports_lag +
  ldist+ colony+ lusmil_lag+ fh_lag+ civilwar_lag+ldisaster_lag +
  year1994+year1995+year1996+year1997+year1998+year1999+
  year2000+year2001+year2002+year2003 +year2004+
  year2005+year2006+year2007+year2008-1


fz <- ~lgdpcap_lag.within+lgdpcap_lag.within+lpopulation_lag.within+  
  lexports_lag.within +
  lusmil_lag.within+ 
  fh_lag.within+ civilwar_lag.within+ldisaster_lag.within +
  year1994.within+year1995.within+year1996.within+year1997.within+
  year1998.within+year1999.within+
  year2000.within+year2001.within+year2002.within+year2003.within +
  year2004.within+ year2005.within+year2006.within+year2007.within+
  year2008.within +
  listock_lag.within + 
  ldist+ colony-1#z1 
## (no z2 here, so no need to include xbar as additional instruments)

ht <- ivreg(update(fx, lcommit3a ~.), fz, data=aid.sam)
ht.vcl <-vcovCL(ht, aid.sam$dyad)
coeftest(ht, vcov=ht.vcl)[1:10,]




aid.sam <- aid.sam %>%
  group_by(dyad) %>%
  mutate(across(all_of(var.names), \(x){mean(x)},
                .names= "{col}.bar")) %>% 
  ungroup()

mundlak <- lm(lcommit3a~listock_lag+lgdpcap_lag+
                lpopulation_lag+ lexports_lag +
                ldist+ colony+
                lusmil_lag+ fh_lag+ civilwar_lag+ldisaster_lag +
                year1994+year1995+year1996+year1997+year1998+year1999+
                year2000+year2001+year2002+year2003 +year2004+
                year2005+year2006+year2007+year2008+
                listock_lag.bar+
                lgdpcap_lag.bar+
               lgdpcap_lag.bar+lpopulation_lag.bar+  
               lexports_lag.bar +
               lusmil_lag.bar+ 
               fh_lag.bar+ civilwar_lag.bar+ldisaster_lag.bar +
               year1994.bar+year1995.bar+year1996.bar+year1997.bar+
               year1998.bar+year1999.bar+
               year2000.bar+year2001.bar+year2002.bar+year2003.bar +
               year2004.bar+ year2005.bar+year2006.bar+year2007.bar+
               year2008.bar,
              data=aid.sam)
m.vcl <-vcovCL(mundlak, aid.sam$dyad)
coeftest(mundlak, vcov=m.vcl)[2:11,]


modelsummary(list("Donor-FE"=m1,  
                  "Dyad-FE"=m2,
                  "Dyad-FE (HT)"=ht,
                  "Dyad-FE (Mundlak)"=mundlak),
             vcov=list(vcov(m1),vcov(m2),
                       ht.vcl, m.vcl),
             fmt=2,
             coef_map=c("listock_lag"="Migrant population (log)",
                        "colony"="Former colony"),
             gof_map=c("nobs"))
