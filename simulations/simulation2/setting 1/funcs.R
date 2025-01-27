


# Oracle case: simulation 2

TFOR.simul2.func <- function(y, k, n, var){
  p <- k/n
  fy <- p*dnorm(y, 1.5, sqrt(var)) + (1-p)*dnorm(y, 0, sqrt(var))
  
  dfy <- p*(1.5-y)/var*dnorm(y, 1.5, sqrt(var)) + (1-p)*(0-y)/var*dnorm(y, 0, sqrt(var))
  
  delta <- y + var * dfy/fy
  return(delta)
}

TFSIOR.simul2.func <- function(y, s, k, n, var){
  p <- k/n
  fy <- p*dnorm(y, 1.5, sqrt(var))*dnorm(s, 3.5, sqrt(var)) +
    p*dnorm(y, 0, sqrt(var))*dnorm(s, 2, sqrt(var)) +
    (1-2*p)*dnorm(y, 0, sqrt(var))*dnorm(s, 0, sqrt(var)) 
  
  
  dfy <- p*(1.5-y)/var*dnorm(y, 1.5, sqrt(var))*dnorm(s, 3.5, sqrt(var)) +
    p*(0-y)/var*dnorm(y, 0, sqrt(var))*dnorm(s, 2, sqrt(var)) +
    (1-2*p)*(0-y)/var*dnorm(y, 0, sqrt(var))*dnorm(s, 0, sqrt(var)) 
  delta <- y + var * dfy/fy
  return(delta)
}


TFKER.func_2 <- function(y, sigma)
{
  n <- length(y)
  
  
  h <- 1/density(y)$bw
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
# TFKER.func_2 <- function(y, sigma)
# {
#   n <- length(y)
#   h.grids <- c(0.5, 1, 1.5, 2, 2.5, 3)
#   minRisk <- Inf
#   minh <- 0
#   K <- 5
#   l <- as.integer(n/K)
#   for(j in length(h.grids))
#   {
#     h <- h.grids[j]
#     curRisk <- 0
#     
#     for(k in 1:K)
#     {
#       le <- (k-1)*l+1
#       ri <- k*l
#       y_test <- y[c(le: ri)]
#       y_train <- y[-c(le: ri)]
#       n_test <- length(y_test)
#       score <- rep(0, n_test)
#       for(i in 1:n_test)
#       {
#         diff <- y_test[i] -y_train
#         diff_exp <- exp(-0.5*h^(2)*diff^2)
#         top <- sum(-diff*diff_exp*h^2)
#         bottom <- sum(diff_exp)
#         score[i] <- top/bottom
#       }
#       
#       cur_delta <- y_test + sigma*score
#       
#       curRisk <- curRisk +  risk.estimator(y_test, cur_delta, sigma)
#     }
#     
#     if(curRisk < minRisk)
#     {
#       minRisk <- curRisk
#       minh <- h
#     }
#   }
#   
#   h <- minh
#   score <- rep(0, n)
#   for(i in 1:n)
#   {
#     diff <- y[i] -y
#     diff_exp <- exp(-0.5*h^(2)*diff^2)
#     top <- sum(-diff*diff_exp*h^2)/n
#     bottom <- sum(diff_exp)/n
#     score[i] <- top/bottom
#   }
#   delta <- y + sigma*score
#   return(delta)
# }




