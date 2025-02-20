
relik <- function(x0, y, X, Ti, Xbar, ybar){
  ## Log-likelihood for the linear random effects model
  ## inputs: 
  ##   x0: length-(k+2) matrix of parameter vector (theta, ln(sigma^2_alpha), ln(sigma^2_epsilon))
  ##   y: length-sum(Ti) vector of the outcome variable
  ##   X: sum(Ti) by k matrix of independent variables
  ##   Ti: length-N vector recording the number of times each unit is observed
  ##   Xbar: sum(Ti) by k matrix of within-unit means of X
  ##   ybar: length-sum(Ti) vector matrix of within-unit means of y
  ## output:
  ##   ll: -1 * log-likelihood at these parameters
  theta <- x0[1:ncol(X)]
  s2.a <- exp(x0[ncol(X)+1])
  s2.e <- exp(x0[ncol(X)+2])
  
  ln2.detSIGMA.i <- 1/2 * log(Ti * s2.a + s2.e) + (Ti-1)/2 * log(s2.e)
  
  w <- 1 - sqrt(s2.e)/sqrt(Ti * s2.a + s2.e)
  w <- rep(w, Ti)
  ehat <- drop(y-X %*%theta)
  ebar <- drop(ybar-Xbar %*%theta)
  e.tilde <- (ehat-w*ebar)/sqrt(s2.e)
  ll <- sum( ln2.detSIGMA.i) + (1/2) * crossprod(e.tilde)
  return(ll)
}





re_mle_gr <- function(x0, y, X, Ti,  Xbar, ybar){
  ## Gradient of the log-likelihood for the linear random effects model
  ## inputs: 
  ##   x0: length-(k+2) matrix of parameter vector (theta, ln(sigma^2_alpha), ln(sigma^2_epsilon))
  ##   y: length-sum(Ti) vector of the outcome variable
  ##   X: sum(Ti) by k matrix of independent variables
  ##   Ti: length-N vector recording the number of times each unit is observed
  ##   Xbar: sum(Ti) by k matrix of within-unit means of X
  ##   ybar: length-sum(Ti) vector matrix of within-unit means of y
  ## output:
  ##   gr: -1 * Derivative of the log-likelihood wrt to x0
  theta <- x0[1:ncol(X)]
  s2.a <- exp(x0[ncol(X)+1])
  s2.e <- exp(x0[ncol(X)+2])
  
  sigma2 = Ti*s2.a+s2.e
  w <- 1 - sqrt(s2.e)/sqrt(Ti * s2.a + s2.e)
  w <- rep(w, Ti)
  ehat  <-  drop(y-X %*%theta)
  ebar <- drop(ybar-Xbar %*%theta)
  
  
  inner <- ((sqrt(s2.e)/(2*rep(sigma2,Ti)^(3/2)))-
              (1/(2*sqrt(s2.e)*rep(sqrt(sigma2),Ti))))
  
  Dtheta  <- -colSums((ehat-w*ebar)*(X - w*Xbar))/s2.e
  Ds2a <- sum(Ti/(2*sigma2)) - 
    sum(rep(Ti,Ti)*ebar*(ehat-w*ebar)/(2*sqrt(s2.e)*rep(sigma2,Ti)^(3/2)))
  Ds2a <- Ds2a*s2.a ## Chain rule for the transformation
  
  Ds2e <- sum((Ti-1)/(2*s2.e)) + sum(1/(2*sigma2)) - 
    sum(((ehat-ebar*w)^2)/(2*s2.e^2)) -
    sum( (ebar*(ehat-ebar*w)*inner))/s2.e
  Ds2e <- Ds2e*s2.e ## Chain rule for the transformation
  gr <- c(Dtheta,Ds2a,Ds2e)
  return(gr)
}

clusterVCOV <- function(X, residuals, idx){
  ## Return clustered standard errors for a linear model
  ## inputs: 
  ##   X: sum(Ti) by k matrix of (transformed if necessary) independent variables. 
  ##   Residuals: length-sum(Ti) vector of (transformed if necessary) residuals from the model 
  ##   idx: Length-sum(Ti) vector that records with unit each observation belongs to
  ## output:
  ##   V: The clustered covariance matrix
  ## NOTE: Requires the split.matrix or split.Matrix function to work properly (see below)
  bread <- solve(crossprod(X))
  meat <- Reduce(`+`, 
                 lapply(split.matrix(X*drop(residuals), idx),
                        \(x){tcrossprod(colSums(x))})
  )
  
  G <- length(unique(idx))
  n <- nrow(X)
  k <- ncol(X)
  V <- (bread %*% meat %*% bread)* (G/(G-1))
  return(V)
}

split.matrix <- function(x, f){
  lapply(split(x = seq_len(nrow(x)), f = f),
         function(ind) x[ind, , drop = FALSE])
  
}


###  a useful way to format numbers for tables particularly when you have 
### leading zeros to worry about 
num2str <- function(x, digits=2){
  digitCounter <- rep(digits, length(x))
  while(any(formatC(abs(x), digits=digits, format='f') == paste0("0.", paste0(rep("0", digits), collapse="")))){
    idx <- which(formatC(abs(x), digits=digits, format='f') == paste0("0.", paste0(rep("0", digits), collapse="")))
    out <- x*0
    out[idx] <-  formatC(x[idx], digits=digits+1, format='f')
    out[-idx] <- formatC(x[-idx], digits=digits, format='f')
    
    digitCounter[idx] <- digitCounter[idx] +1
    digits <- digits +1
    # cat("digits=", digits, "\n")
  }
  out <- x*0
  for(i in 1:length(x)){
    out[i] <- formatC(x[i], digits=digitCounter[i], format='f')
  }
  return(out)
}
