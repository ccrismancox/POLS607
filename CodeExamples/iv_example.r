library(dplyr)
library(readstata13)
library(fixest)
library(ivreg)
library(sandwich)
library(lmtest)
library(car)
library(Matrix) #sparse matrices

terror <- read.dta13("datasets/subsample.dta")


## lead the outcome by one (lag everything else)
terror <- terror %>% 
  mutate(f.nattack = lead(nattack, 1), .by=id) ## move the DV forward one

baseline.ols <- feols(f.nattack~v2x_corr+sp_pop_totl+ny_gdp_pcap_kd+
                        kg_democracy+statefailure|id+year, data=terror)
summary(baseline.ols)

within.2sls <- feols(f.nattack~sp_pop_totl+ny_gdp_pcap_kd+
                       kg_democracy+statefailure|id+year|
                       v2x_corr~iv_region, data=terror)
summary(within.2sls, stage=1:2)

## subset to used data
terror2 <- terror[as.numeric(within.2sls$obs_selection$obsRemoved),]

## panel dimensions
length(unique(terror2$id))
summary(terror2 %>% summarize(length(year), .by=id))


## Sparse functions
DeltaN <- sparse.model.matrix(~factor(id)-1, data=terror2)
DeltaT <- sparse.model.matrix(~factor(year)-1, data=terror2)[,-1]
M <- Diagonal(nrow(terror2)) - DeltaN %*% solve(crossprod(DeltaN)) %*% t(DeltaN)

## build two-way transformation (often a drag)
M2 <- M - M %*% DeltaT  %*% solve(t(DeltaT) %*% M %*% DeltaT) %*% t(DeltaT) %*% M



## alternative is to include the time dummies with X and just do the regular 
## within transformation
var.names<- c("f.nattack", "v2x_corr", "sp_pop_totl", 
              "ny_gdp_pcap_kd", "kg_democracy", "statefailure", "iv_region")
terror2 <- terror2 %>%
  mutate(across(all_of(var.names), \(x){as.numeric(M2 %*%x)},
         .names= "{col}.within"))
        



within.2sls2 <- ivreg(f.nattack.within~v2x_corr.within+
                        sp_pop_totl.within+ny_gdp_pcap_kd.within+
                        kg_democracy.within+statefailure.within-1|
                        iv_region.within+
                        sp_pop_totl.within+ny_gdp_pcap_kd.within+
                        kg_democracy.within+statefailure.within-1,
                      data=terror2)
summary(within.2sls2, vcov=\(x){vcovCL(x,cluster=terror2$id)})

within.first <- lm(v2x_corr.within~
                     iv_region.within+
                     sp_pop_totl.within+ny_gdp_pcap_kd.within+
                     kg_democracy.within+statefailure.within-1,
                   data=terror2)
coeftest(within.first, vcov=vcovCL(within.first, terror2$id))
linearHypothesis(within.first, "iv_region.within",
                 vcov=vcovCL(within.first, terror2$id))
