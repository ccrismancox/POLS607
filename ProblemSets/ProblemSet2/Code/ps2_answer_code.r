library(readstata13)
library(Matrix)
library(fixest)
library(lme4)
library(sandwich)
library(car)
library(modelsummary)
library(xtable)
library(dplyr)


source("panelFunctions.r")


## read in the data
cw <- read.dta13("conflict.dta", convert.dates = FALSE)
cw <- cw %>% 
  arrange(ccode, year) 

## summmary stats
sum.stats <- cw %>%
  select(c(any_prio, gdp_g, gdp_g_l,
         y_0, polity2l, ethfrac, relfrac,
         Oil, lmtnest, lpopl1)) %>% 
  reframe(across(everything(),
                 \(x){
                   c(mean(x, na.rm=TRUE),
                     sd(x, na.rm=TRUE),
                     length(na.omit(x)))
                 }
  )) %>% 
  t() %>% 
  round(digits=3)

rownames(sum.stats) <- c("PRIO conflict",
                         "Growth",
                         "Lagged growth",
                         "Log(GDP p.c.), 1979",
                         "Democracy (lag)",
                         "Eth. Frac.",
                         "Rel. Frac.",
                         "Oil exporter",
                         "Log(Mountainous terrain)",
                         "Log(Population) (lag)")
colnames(sum.stats) <- c("Mean", "St.\ Dev.", "Obs.")
sum.stats[,3] <- as.character(sum.stats[,3])

print(xtable(sum.stats, align="rlll",
             caption="Summary statistics",
             label="tab:sumstats"),
      sanitize.text.function=\(x){x},
      booktabs=TRUE,
      caption.placement="top",
      hline.after=c(-1,0,1,3,10))





m2 <- lm(any_prio~gdp_g + gdp_g_l+y_0+polity2l+ethfrac+relfrac+
           Oil+lmtnest+ lpopl1 +I(year-1978), data=cw)
V2 <- vcovCL(m2, ~ccode)
m3 <- lm(any_prio~gdp_g + gdp_g_l+y_0+polity2l+ethfrac+relfrac+
           Oil+lmtnest+ lpopl1+factor(ccode):I(year-1978), data=cw)
V3 <- vcovCL(m3, ~ccode)
m4 <- feols(any_prio~gdp_g + gdp_g_l+factor(ccode):year|ccode, data=cw)


others <- rbind.data.frame(c("Country fixed effects", "no", "no", "yes"),
                           c("Country-specific time trend", "no", "yes", "yes"))
attr(others,"position") <- 19:20
modelsummary(list(m2, m3, m4),
             fmt=2,
             title="Replicating the OLS models from Table 4 of MSS\\label{tab:ols}",
             escape=FALSE,
             output="latex",
             vcov=list(V2, V3, NULL),
             coef_map=c("gdp_g"="Growth",
                        "gdp_g_l"="Lagged growth",
                        "y_0"="Log(GDP p.c.), 1979",
                        "polity2l"="Democracy ",
                        "ethfrac"= "Eth. Frac.",
                        "relfrac"="Rel. Frac.",
                        "Oil"="Oil exporter",
                        "lmtnest"="Log(Mountainous terrain)",
                        "lpopl1"="Log(Population)"),
             gof_map=c("r.squared", "rmse","nobs"),
             add_rows=others
)




m4 <- feols(any_prio~gdp_g + gdp_g_l|ccode+year, data=cw)
m4a <- update(m4, vcov=~year)
m4b <- update(m4, vcov="DK", panel.id=c("ccode", "year"))
m4c <- update(m4, vcov="twoway", panel.id=c("ccode", "year"))


se.comp <- rbind(m4$coefficients,
                 sqrt(diag(vcov(m4))),
                 sqrt(diag(vcov(m4a))),
                 sqrt(diag(vcov(m4b))),
                 sqrt(diag(vcov(m4c))))

rownames(se.comp) <- c("Estimate",
                       "Clustered by country",
                       "Clustered by year",
                       "Driscoll-Kraay",
                       "2-way clustering")
colnames(se.comp) <- c("Growth", "Lagged growth")
print(xtable(se.comp, align="rll",
             caption="Comparing standard errors",
             label="tab:se_comp"),
      sanitize.text.function=\(x){x},
      booktabs=TRUE,
      caption.placement="top")


me <- feols(any_prio~polity2l+lpopl1|ccode+year|gdp_g~GPCP_g,
            data=cw)


## one approach the system
sys.dat <- cw %>% 
  select(ccode, year, any_prio,lpopl1, 
         polity2l, gdp_g, GPCP_g)
sys.dat <- rbind(sys.dat, sys.dat)
sys.dat <-sys.dat %>%
  mutate(eq2 =c(rep(0, nrow(cw)), rep(1, nrow(cw))),
         eq1 = 1-eq2,
         ccode2 = ccode+eq2*500,
         year2 = year+eq2*500,
         y=c(cw$any_prio, cw$gdp_g))

sys.est <- feols(y~GPCP_g:eq1+polity2l:eq1+lpopl1:eq1+
                   GPCP_g:eq2+ polity2l:eq2+lpopl1:eq2|ccode2+year2,
                 data=sys.dat,
                 cluster=~ccode)
sys1 <- deltaMethod(sys.est, "`GPCP_g:eq1`/`GPCP_g:eq2`")
### pretty good match


## another system approach
y <- sys.dat$y
X <- with(cw, cbind(GPCP_g,polity2l, lpopl1))
Di <- sparse.model.matrix(~factor(ccode) + factor(year)-1, data=cw)
X.sys <- bdiag(X,X)
D.sys <- bdiag(Di, Di)
bigX <- cbind(X.sys,D.sys)
## Note the column ordering here is:
## cols 1:3 are eq1 rain, democracy, population
## cols 4:6 are eq2 rain, democracy, population


rf.hat <- drop(solve(crossprod(bigX)) %*% t(bigX) %*% y)
ehat <- y-bigX %*% rf.hat
V <- clusterVCOV(bigX, ehat, rep(cw$ccode, 2) )
se.delta <- sqrt(c(1/rf.hat[4], -rf.hat[1]/rf.hat[4]^2) %*%  
                   V[c(1,4), c(1,4)] %*% 
                   c(1/rf.hat[4], -rf.hat[1]/rf.hat[4]^2))

## pretty close!
out <- rbind(summary(me)$coeftable[1,1:2],
             c(rf.hat[1]/rf.hat[4], drop(se.delta)))
rownames(out) <- c("2SLS", "System")
print(xtable(out, align="rll",
             caption="Identifying $\\beta_1$",
             label="tab:id"),
      sanitize.text.function=\(x){x},
      booktabs=TRUE,
      caption.placement="top")





m5 <- feols(any_prio ~  y_0 + polity2l + ethfrac + 
              relfrac + Oil + lmtnest + lpopl1 +
              factor(ccode):I(year-1978)|
              gdp_g + gdp_g_l ~ GPCP_g + GPCP_g_l, 
            vcov=~ccode,
            data = cw)
m6 <- feols(any_prio ~  
              factor(ccode):year|ccode|
              gdp_g + gdp_g_l ~ GPCP_g + GPCP_g_l, 
            vcov=~ccode,
            data = cw)

others <- rbind.data.frame(c("Country fixed effects", "no", "yes"),
                           c("Country-specific time trend", "yes", "yes"))
attr(others,"position") <- 19:20
modelsummary(list(m5, m6),
             fmt=2,
             title="Replicating 2SLS results from MSS\\label{tab:2sls}",
             escape=FALSE,
             output="latex",
             coef_map=c("fit_gdp_g"="Growth",
                        "fit_gdp_g_l"="Lagged growth",
                        "y_0"="Log(GDP p.c.), 1979",
                        "polity2l"="Democracy ",
                        "ethfrac"= "Eth. Frac.",
                        "relfrac"="Rel. Frac.",
                        "Oil"="Oil exporter",
                        "lmtnest"="Log(Mountainous terrain)",
                        "lpopl1"="Log(Population)"),
             gof_map=c("r.squared", "rmse","nobs"),
             add_rows=others
)




## first stage F
Fstats <- c(linearHypothesis(m5$iv_first_stage$gdp_g,
                             c("GPCP_g", "GPCP_g_l"), test="F")$F[2],
            linearHypothesis(m5$iv_first_stage$gdp_g_l, 
                             c("GPCP_g", "GPCP_g_l"), test="F")$F[2],
            linearHypothesis(m6$iv_first_stage$gdp_g,
                             c("GPCP_g", "GPCP_g_l"), test="F")$F[2],
            linearHypothesis(m6$iv_first_stage$gdp_g_l, 
                             c("GPCP_g", "GPCP_g_l"), test="F")$F[2])
Fstats <- matrix(Fstats, nrow=2)
colnames(Fstats) <- c("Model 5", "Model 6")
rownames(Fstats) <- c("Growth", "Lagged growth")

print(xtable(Fstats, align="rll",
             caption="First-stage $F$ statistics",
             label="tab:F"),
      sanitize.text.function=\(x){x},
      booktabs=TRUE,
      caption.placement="top")


## over id
mi <- feols(any_prio ~ polity2l +  lpopl1|ccode+year|
              gdp_g ~ GPCP_g + GPCP_g_l, 
            vcov=~ccode,
            data = cw)

cw$ehat <- cw$any_prio- predict(mi, newdata=cw)
aux <- feols(ehat~GPCP_g + GPCP_g_l+polity2l +  lpopl1|ccode+year, data=cw)

## overidentified Sargan test for 2sls
J.2sls <- (1- crossprod(aux$residuals)/crossprod(mi$residuals))*743
J.2sls
pchisq(J.2sls,df=1, lower=FALSE)





## GMM
gmm.dat <- cw %>%
  select(ccode, year, 
         any_prio, gdp_g,
         polity2l, lpopl1,
         GPCP_g, GPCP_g_l)

tmat <- model.matrix(~factor(year)-1, data=cw)[,-1]
colnames(tmat) <- paste0("yy", 1982:1999)

gmm.dat <- cbind(gmm.dat,tmat)
gmm.dat <- gmm.dat %>%
  mutate(across(!year,
                \(x){x-mean(x)}, 
                .names = "{col}.within"),
         .by=ccode) %>% 
  select(ccode, ends_with(".within"))

Zdot <- gmm.dat %>% 
  select(GPCP_g.within, GPCP_g_l.within,
         polity2l.within, lpopl1.within,
         starts_with("yy")) %>% 
  as.matrix()
Xdot <- gmm.dat %>% 
  select(gdp_g.within,
         polity2l.within, lpopl1.within,
         starts_with("yy")) %>% 
  as.matrix()
ydot <- gmm.dat$any_prio.within
eps.hat <- mi$residuals



## 2-step gmm 
### 1. start with W from another model (e.g. 2sls or a different gmm)
W <- solve(Reduce( `+`, by(cw$ccode, data=Zdot*eps.hat,
                           \(x){ 
                             tcrossprod(colSums(x))
                           })))


### 2. Fit the GMM with this W
beta.hat.gmm <- solve(t(Xdot) %*% Zdot %*% W %*%t(Zdot) %*% Xdot) %*% 
  (t(Xdot) %*% Zdot %*% W %*% t(Zdot) %*% ydot)


### 3. Standard errors with the new W
e.gmm.hat <- gmm.dat$any_prio.within- Xdot %*%beta.hat.gmm
Ze <- Zdot*drop(e.gmm.hat)
W1 <- solve(Reduce( `+`, by(gmm.dat$ccode, data=Ze,
                            \(x){ 
                              tcrossprod(colSums(x))
                            })))
V.gmm <- solve(t(Xdot) %*% Zdot %*% (W1) %*% t(Zdot) %*% Xdot)
se.gmm <- sqrt(diag(V.gmm))

cbind(beta.hat.gmm, se.gmm)[1:3,]
## Sargan
G <- colSums(Ze)
J <-  (G %*% W %*% G)
pVal.J <- pchisq(J, lower=FALSE, df=ncol(Zdot)-ncol(Xdot))




out <- cbind(beta.hat.gmm, se.gmm)[1:3,]
colnames(out) <- c("Est", "St. Err.")
rownames(out) <- c("Growth", "Democracy", "Population")

J.test <- paste0("\\midrule\nSargan Test & \\multicolumn{2}{c}{",
                 round(J,2), " ($p \\approx ", round(pVal.J,2),"$)}\\\\")

print(xtable(out, digits=3,
             caption="GMM results",
             label="tab:gmm"),
      sanitize.text.function=\(x){x},
      booktabs=TRUE,
      caption.placement="top",
      hline.after=c(-1,0,3),
      add.to.row=list(pos=list(3), 
                       command=J.test))



X<- cw %>% 
  select(ccode, year, polity2, lpopl1,gdp_g_l) %>%
  na.omit()
tmat <- model.matrix(~factor(year)-1, data=X)[,-1]
colnames(tmat) <- paste0("yy", 1982:1999)
X <- cbind(X, tmat)
X.dot <- X %>% 
  mutate(across(!year,
                \(x){x-mean(x)}, 
                .names = "{col}.within"),
         .by=ccode) %>% 
  select(ends_with(".within"))
X <- X %>% 
  select(!c(ccode, year))
ht.dat <- cw %>% 
  filter(!is.na(polity2)) %>%
  select(ccode, year, any_prio, ethfrac) %>% 
  cbind(., X, X.dot)


ht <- ivreg(any_prio~ethfrac + polity2 + lpopl1 + gdp_g_l + 
              yy1982 + yy1983 + yy1984 + yy1985 + yy1986 + 
              yy1987 + yy1988 + yy1989 + yy1990 + yy1991 + 
              yy1992 + yy1993 + yy1994 + yy1995 + yy1996 + 
              yy1997 + yy1998 + yy1999|
              ethfrac + polity2.within + lpopl1.within + gdp_g_l.within + 
              yy1982.within + yy1983.within + yy1984.within + yy1985.within + 
              yy1986.within + yy1987.within + yy1988.within + yy1989.within + 
              yy1990.within + yy1991.within + yy1992.within + yy1993.within + 
              yy1994.within + yy1995.within + yy1996.within + yy1997.within + 
              yy1998.within + yy1999.within,
            data=ht.dat)
            


m.dat <- ht.dat %>%
  select(!ends_with(".within")) %>%
  mutate(across(!c(year, any_prio),
                \(x){mean(x)}, 
                .names = "{col}.bar"),
         .by=ccode)

m.dat <- m.dat[, !duplicated(round(as.matrix(m.dat),15),MARGIN=2)]

mundlak <- lm(any_prio~ethfrac + polity2 + lpopl1 + gdp_g_l + 
                yy1982 + yy1983 + yy1984 + yy1985 + yy1986 + 
                yy1987 + yy1988 + yy1989 + yy1990 + yy1991 + 
                yy1992 + yy1993 + yy1994 + yy1995 + yy1996 + 
                yy1997 + yy1998 + yy1999+
                polity2.bar + lpopl1.bar + gdp_g_l.bar + 
                yy1982.bar + yy1991.bar + yy1992.bar,
              data=m.dat)
                

mundlak2 <- lmer(any_prio~ethfrac + polity2 + lpopl1 + gdp_g_l + 
  yy1982 + yy1983 + yy1984 + yy1985 + yy1986 + 
  yy1987 + yy1988 + yy1989 + yy1990 + yy1991 + 
  yy1992 + yy1993 + yy1994 + yy1995 + yy1996 + 
  yy1997 + yy1998 + yy1999+
  polity2.bar + lpopl1.bar + gdp_g_l.bar + 
  yy1982.bar + yy1991.bar + yy1992.bar+(1|ccode),
data=m.dat)


Vht <- vcovCL(ht, ~ccode)
Vmundlak <- vcovCL(mundlak,~ccode)



modelsummary(list(HT=ht, Mundlak=mundlak),
             fmt=2,
             vcov=~ccode,
             title="Working with time invariant covariates\\label{tab:invar}",
             escape=FALSE,
             output="latex",
             coef_map=c("polity2"="Democracy ",
                        "ethfrac"= "Eth. Frac.",
                        "gdp_g_l"="Lagged growth",
                        "lpopl1"="Log(Population)"),
             gof_map=c("r.squared","nobs"))







