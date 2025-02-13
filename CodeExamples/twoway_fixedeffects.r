library(data.table)
library(readstata13)
library(fixest)
library(ivreg)
library(sandwich)
library(lmtest)
library(car)
library(Matrix) #sparse matrices
rm(list=ls())
terror <- read.dta13("datasets/corruption_terrorism_subsample.dta")

terror <- data.table(terror)
#####
within1 <- feols( nattack ~ v2x_corr+sp_pop_totl+ ny_gdp_pcap_kd
                  +kg_democracy +statefailure|id,data=terror)
summary(within1)


within2 <- feols( nattack ~ v2x_corr+sp_pop_totl+ ny_gdp_pcap_kd
                  +kg_democracy +statefailure|id+year,data=terror)
summary(within2)


within2a <- feols( nattack ~ v2x_corr+sp_pop_totl+ ny_gdp_pcap_kd
                   +kg_democracy +statefailure|id+year,
                   cluster=~year,
                   data=terror)
summary(within2a)

within2b <- feols( nattack ~ v2x_corr+sp_pop_totl+ ny_gdp_pcap_kd
                   +kg_democracy +statefailure|id+year,
                   vcov="DK",
                   panel.id=c("id", "year"),
                   data=terror)
summary(within2b)


within2c <- feols( nattack ~ v2x_corr+sp_pop_totl+ ny_gdp_pcap_kd
                   +kg_democracy +statefailure|id+year,
                   vcov="twoway",
                   data=terror)
summary(within2c)



sqrt(diag(vcov(within2) + vcov(within2a) 
          - vcov(update(within2, vcov="hetero"))))








####
baseline.ols <- feols(nattack~v2x_corr+sp_pop_totl+ny_gdp_pcap_kd+
                        kg_democracy+statefailure|id+year, data=terror)

within.2sls <- feols(nattack~sp_pop_totl+ny_gdp_pcap_kd+
                       kg_democracy+statefailure|id+year|
                       v2x_corr~iv_region, data=terror)
summary(within.2sls, stage=1:2)

terror2 <- terror[as.numeric(within.2sls$obs_selection$obsRemoved)]
length(unique(terror2$id))
summary(terror2[,length(year), by=id]$V1)


## Sparse functions
DeltaN <- sparse.model.matrix(~factor(id)-1, data=terror2)
DeltaT <- sparse.model.matrix(~factor(year)-1, data=terror2)[,-1]
M <- Diagonal(nrow(terror2)) - DeltaN %*% solve(crossprod(DeltaN)) %*% t(DeltaN)


## build two-way transformation
M2 <- M - M %*% DeltaT  %*% solve(t(DeltaT) %*% M %*% DeltaT) %*% t(DeltaT) %*% M

var.names<- c("nattack", "v2x_corr", "sp_pop_totl", 
              "ny_gdp_pcap_kd", "kg_democracy", "statefailure", "iv_region")
terror2[ , paste0(var.names, ".within") := lapply(.SD, \(x){as.numeric(M2 %*%x)}), 
         .SDcols=var.names ]


within.2sls2 <- ivreg(nattack.within~v2x_corr.within+
                        sp_pop_totl.within+ny_gdp_pcap_kd.within+
                        kg_democracy.within+statefailure.within-1|
                        iv_region.within+
                        sp_pop_totl.within+ny_gdp_pcap_kd.within+
                        kg_democracy.within+statefailure.within-1,
                      data=terror2)

within.2sls2 <- ivreg(nattack.within~v2x_corr.within+
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
