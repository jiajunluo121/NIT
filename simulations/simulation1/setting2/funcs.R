


# Oracle case: simulation 1

TFOR.simul1.func <- function(y, var, sd){
  var_total <- var + 1 + sd^2
  fy <- 0.5*dnorm(y, 0, sqrt(var_total)) + 0.5*dnorm(y, 1, sqrt(var_total))
  
  dfy <- 0.5*(0-y)/var_total*dnorm(y, 0, sqrt(var_total)) + 0.5*(1-y)/var_total*dnorm(y, 1, sqrt(var_total))
  
  delta <- y + var * dfy/fy
  return(delta)
}

TFSIOR.simul1.func <- function(y, s, vary, sd){
  v1 <- 1 + sd^2
  v2 <- 1 + sd^2
  print(v2)
  fys <-0.5/(2*pi)*sqrt(1/(v1*v2+v1+v2))*exp((-(v1*v2*1+v2*y^2+v1*s^2)*(v1*v2+v1+v2)+(v1*v2*1+v2*y+v1*s)^2)/(2*(v1*v2+v1+v2)*v1*v2)) + 
    0.5/(2*pi)*sqrt(1/(v1*v2+v1+v2))*exp((-(v2*y^2+v1*s^2)*(v1*v2+v1+v2)+(v2*y+v1*s)^2)/(2*(v1*v2+v1+v2)*v1*v2))
  
  dyfys <- 0.5/(2*pi)*sqrt(1/(v1*v2+v1+v2))*exp((-(v1*v2*1+v2*y^2+v1*s^2)*(v1*v2+v1+v2)+(v1*v2*1+v2*y+v1*s)^2)/(2*(v1*v2+v1+v2)*v1*v2))*
    (-v2*2*y*(v1*v2+v1+v2)+2*v2*(v1*v2*1+v2*y+v1*s))/(2*(v1*v2+v1+v2)*v1*v2) +
    0.5/(2*pi)*sqrt(1/(v1*v2+v1+v2))*exp((-(v2*y^2+v1*s^2)*(v1*v2+v1+v2)+(v2*y+v1*s)^2)/(2*(v1*v2+v1+v2)*v1*v2))*
    (-v2*2*y*(v1*v2+v1+v2)+2*v2*(v2*y+v1*s))/(2*(v1*v2+v1+v2)*v1*v2)
  delta <- y + vary * dyfys/fys
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




