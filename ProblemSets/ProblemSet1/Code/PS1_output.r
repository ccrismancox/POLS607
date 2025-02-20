library(matrixStats)
library(xtable)
rm(list=ls())
load("PS1.rdata")
source("panelFunctions.r")

## Expected Omitted variable bias in pooled estimator
conds$pooled.OVB[1:2]


## Main simulation
estimates <- t(sapply(Results[1:2], \(x){colMeans(x[,1:4])}))
obs.sd <- t(sapply(Results[1:2], \(x){colSds(x[,1:4])}))
classical.se <- t(sapply(Results[1:2], \(x){colMeans(x[,seq(5,11, by=2)])}))
cluster.se <- t(sapply(Results[1:2], \(x){colMeans(x[,seq(6,12, by=2)])}))
power <- t(sapply(Results[1:2], \(x){colMeans(x[,13:14])}))
size <- t(sapply(Results[3:4], \(x){colMeans(x[,13:14])}))


table1 <- matrix(rbind(c(estimates), 
                       c(obs.sd), 
                       c(classical.se), 
                       c(cluster.se)), ncol=4)
table1 <- num2str(table1)


table1 <- cbind(c("N=50", rep("", 3),
                  "N=200", rep("", 3)),
                c("Average estimate", 
                  "Observed st.~dev.", 
                  "Classical st.~err.",
                  "Clustered st.~err."),
                table1)
colnames(table1) <- c(" ", " ", "Pooled", "Within", "RE-MLE", "RE-GLS")

print(xtable(table1,
             caption="Different approaches for estimating $\\beta_1=1$",
             label="tab:sim1",
             align="lllcccc"),
      include.rownames=FALSE,
      booktabs=TRUE,
      caption.placement = "top",
      sanitize.text.function = \(x){x},
      hline.after=c(-1,0,4,8))





table2 <- matrix(rbind(c(power), c(size)), ncol=2)
table2 <- round(table2, digits=3)
table2 <- cbind(c("N=50", "",
                  "N=200", ""),
                c("Power", 
                  "Size"),
                table2)
colnames(table2) <- c(" ", " ", "Hausman", "Mundlak")


print(xtable(table2,
             caption="Comparing tests for fixed effects versus efficient estimation",
             label="tab:sim2",
             align="lllcc"),
      include.rownames=FALSE,
      booktabs=TRUE,
      caption.placement = "top",
      sanitize.text.function = \(x){x},
      hline.after=c(-1,0,2,4))


