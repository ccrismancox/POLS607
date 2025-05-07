effects <- function(est, V, differenced=TRUE){
  ## Substantive effects for democracy on the outcome
  ## inputs: 
  ##   est: a vector of length 2-5. The first element is the 
  ##        coefficient of democracy. The next 1-4 elements are 
  ##        are the coefficients on lagged outcomes
  ##   V: A matrix of est by est. Variance-covariance of matrix 
  ##      the elements in est
  ## output:
  ##   outMat: A 3 by 2 matrix of substantive effects and standard errors 
  
  
  ## setup coefficients
  rho <- est[2:length(est)]
  rho <- c(rho, rep(0, 4-length(rho)))
  rho[is.na(rho)] <- 0
  beta <- est[1]
  
  ## first four effects
  eff1 <- beta
  eff2 <- rho[1]*eff1 + beta
  eff3 <- rho[1]*eff2 + rho[2]*eff1 + beta
  eff4 <- rho[1]*eff3 + rho[2]*eff2 + rho[3]*eff1 + beta
  eff <- c(eff4, eff3, eff2, eff1)
  
  ## If we're dealing with differences, we'll sum these
  ## for the total change
  level <- sum(eff)
  
  
  ## Derivatives
  Deff1 <- c(1, rep(0, 4))
  Deff2 <- c(1+rho[1]*Deff1[1], eff1,  rep(0, 3))
  Deff3 <- c(1+rho[1]*Deff2[1] + rho[2]*Deff1[1], 
             eff2+ rho[1]*Deff2[2],
             eff1,
             rep(0,2))
  Deff4 <- c(1+rho[1]*Deff3[1] + rho[2]*Deff2[1]+rho[3]*Deff1[1], 
             eff3 +rho[1]*Deff3[2] + rho[2]*Deff2[2],
             eff2 + rho[1]*Deff3[3]+rho[2]*Deff2[3],
             eff1, 
             0)
  
  D <- rbind(Deff4, Deff3, Deff2, Deff1)
  
  ## If we're dealing with differences, we'll sum these
  ## for the total change
  Dlevel <- colSums(D)
  for(j in 5:25){
    ## loop for the remaining periods 
    eff.j <- rho %*% eff +beta
    
    Dj <- c(1+ rho %*% D[,1],
            rho %*% D[,-1] + eff)
    D <- rbind(Dj,D[-4,])    
    eff <- c(eff.j, eff[-4])
    
    
    level <- level + eff.j
    Dlevel <- Dlevel + Dj
  }
  
  if(differenced){
    ## if the model is in differences
    ## we want to return that summed value (level)
    Dlevel <- Dlevel[1:length(est)]
    eff25 <- level
    V25 <- drop(Dlevel %*% V %*% Dlevel)
  }else{
    ## otherwise, it's already in levels
    Dj <- Dj[1:length(est)]
    eff25 <- eff.j
    V25 <- drop(Dj %*% V %*% Dj)
  }
  
  
  ## asked for long run effect on the outcome
  ## so this is right either way
  lr <- beta/(1-sum(rho))
  D.lr <- c(1/(1-sum(rho)),
            rep(beta/((1-sum(rho))^2),4))[1:length(est)]
  
  V.lr <- drop(D.lr %*% V %*% D.lr     )
  
  
  persist <- sum(rho)
  D.per <- c(0, rep(1,4))[1:length(est)]
  V.p <- drop(D.per %*% V %*% D.per)
  
  outMat <- cbind(c(lr,eff25, persist), sqrt(c(V.lr,V25, V.p)))
  rownames(outMat) <- c("Long-run effect of democracy",
                        "Effect of democracy after 25 years",
                        "Persistence of GDP")
  colnames(outMat) <- c("Est", "SE")
  return(outMat)
}




