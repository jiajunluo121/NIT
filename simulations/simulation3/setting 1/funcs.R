


# Oracle case: simulation 1

TFOR.simul1.func <- function(y, vary, var,se, p){

  fy <- (1-p)*dnorm(y, 0, sqrt(vary)) + p*dnorm(y, se, sqrt(vary))
  dfy <- (1-p)*(0-y)/vary*dnorm(y, 0, sqrt(vary)) + p*(se-y)/vary*dnorm(y, se, sqrt(vary))
  
  delta <- y + var * dfy/fy
  return(delta)
}

TFSIOR.simul1.func <- function(y, s, vary, vars, var, se, p){
  
  fys <- (1-p)*dnorm(y, 0, sqrt(vary))*dnorm(s[,1], 0, sqrt(vars))*dnorm(s[,2], 0, sqrt(vars))*
    dnorm(s[,3], 0, sqrt(vars))*dnorm(s[,4], 0, sqrt(vars)) + 
    p*dnorm(y, se, sqrt(vary))*dnorm(s[,1], se, sqrt(vars))*dnorm(s[,2], se, sqrt(vars))*
    dnorm(s[,3], se, sqrt(vars))*dnorm(s[,4], se, sqrt(vars))
  
  dyfys <- (1-p)*(0-y)*dnorm(y, 0, sqrt(vary))/vary*dnorm(s[,1], 0, sqrt(vars))*dnorm(s[,2], 0, sqrt(vars))*
    dnorm(s[,3], 0, sqrt(vars))*dnorm(s[,4], 0, sqrt(vars)) + 
    p*(se-y)*dnorm(y, se, sqrt(vary))/vary*dnorm(s[,1], se, sqrt(vars))*dnorm(s[,2], se, sqrt(vars))*
    dnorm(s[,3], se, sqrt(vars))*dnorm(s[,4], se, sqrt(vars)) 
    
  delta <- y + var * dyfys/fys
  return(delta)
}

TFKER.func_2 <- function(y, sigma)
{
  n <- length(y)
  h.grids <- c(0.5, 1, 1.5, 2, 2.5, 3)
  minRisk <- 2000
  minh <- 0
  K <- 5
  l <- as.integer(n/K)
  for(j in length(h.grids))
  {
    h <- h.grids[j]
    curRisk <- 0
    
    for(k in 1:K)
    {
      le <- (k-1)*l+1
      ri <- k*l
      y_test <- y[c(le: ri)]
      y_train <- y[-c(le: ri)]
      n_test <- length(y_test)
      score <- rep(0, n_test)
      for(i in 1:n_test)
      {
        diff <- y_test[i] -y_train
        diff_exp <- exp(-0.5*h^(2)*diff^2)
        top <- sum(-diff*diff_exp*h^2)
        bottom <- sum(diff_exp)
        score[i] <- top/bottom
      }
      
      cur_delta <- y_test + sigma*score
      
      curRisk <- curRisk +  risk.estimator(y_test, cur_delta, sigma)
    }
    
    if(curRisk < minRisk)
    {
      minRisk <- curRisk
      minh <- h
    }
  }
  
  h <- minh
  score <- rep(0, n)
  for(i in 1:n)
  {
    diff <- y[i] -y
    diff_exp <- exp(-0.5*h^(2)*diff^2)
    top <- sum(-diff*diff_exp*h^2)/n
    bottom <- sum(diff_exp)/n
    score[i] <- top/bottom
  }
  delta <- y + sigma*score
  return(delta)
}




