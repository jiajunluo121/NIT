library(caret)
library(xgboost)
library(mgcv)

# NIT estimator.
ksd_MahDist_est.1 <- function(y, s, h, var, monotonicity=FALSE)
{
  # Input:
  ## y: primary observations
  ## s: auxiliary sequences
  ## h: bandwidth, note that h = 1/lambda
  ## var: standard deviation of primary observation based on normal means

  n <- length(y)
  K <- matrix(0, n, n)
  dK <- matrix(0, n, n)
  X <- cbind(y, s)
  covM <- cov(X)
  invCov <- solve(covM)
  for(i in 1:n)
  {
    tmat <- matrix(rep(X[i,],each=n),n)-X
    K[i,] <- exp(-0.5*h^(2)*mahalanobis(tmat, center=FALSE, invCov, inverted = TRUE))
    dK[i,] <- h^2*tmat %*% invCov[1,]*K[i,]
  }

  w <- Variable(n)
  constVec<- matrix(rep(1,n),nrow=n,ncol=1)

  # Object function for optimization function
  obj <- quad_form(w, K) +2*t(w) %*% dK %*% constVec

  # Constraints
  ## Unbiaseness
  constraints <- list(t(constVec)%*%w==0)

  ## Monotonicity.
  if (monotonicity) {
    yest <- y + var*w
    idxs <- sort.list(y)
    for(i in 2:n)
    {
      idx2 <- idxs[i]
      idx1 <- idxs[i-1]
      append(constraints, yest[idx2] >=yest[idx1])
    }
  }

  prob <- Problem(Minimize(obj), constraints)
  res<- solve(prob)

  w.hat <- res$getValue(w)
  val<- res$value
  return(list("w.hat"=w.hat,"val"=val, "metric"=res$metrics))

}

ksd_MahDist_est <- function(y, s, h, var, monotonicity=FALSE)
{
  # Input:
  ## y: primary observations
  ## s: auxiliary sequences
  ## h: bandwidth, note that h = 1/lambda
  ## var: standard deviation of primary observation based on normal means
  
  n <- length(y)
  K <- matrix(0, n, n)
  dK <- matrix(0, n, n)
  X <- y
  covM <- var(X)
  invCov <- 1.0/covM
  
  disc.dist <- discrete_mand_dist(s)
  for(i in 1:n)
  {
    tmat <- matrix(rep(X[i],each=n),n)-X
    K[i,] <- exp(-0.5*h^(2)*(mahalanobis(tmat, center=FALSE, invCov, inverted = TRUE)+ disc.dist[i,]))
    dK[i,] <- h^2*tmat %*% invCov*K[i,]
  }
  
  w <- Variable(n)
  constVec<- matrix(rep(1,n),nrow=n,ncol=1)
  
  # Object function for optimization function
  obj <- quad_form(w, K) +2*t(w) %*% dK %*% constVec
  
  # Constraints
  ## Unbiaseness
  constraints <- list(t(constVec)%*%w==0)
  
  ## Monotonicity.
  if (monotonicity) {
    yest <- y + var*w
    idxs <- sort.list(y)
    for(i in 2:n)
    {
      idx2 <- idxs[i]
      idx1 <- idxs[i-1]
      append(constraints, yest[idx2] >=yest[idx1])
    }
  }
  
  prob <- Problem(Minimize(obj), constraints)
  res<- solve(prob)
  
  w.hat <- res$getValue(w)
  val<- res$value
  return(list("w.hat"=w.hat,"val"=val, "metric"=res$metrics))
  
}



# Bandwidth choice based on MCV.
bd.est <-function(y, s, var){
  n <- length(y)
  alpha <- 0.1
  ep <- rnorm(n,mean=0,sd=1)
  U <- y + sqrt(var)*alpha*ep
  V <- y- sqrt(var)*ep/alpha
  h.grids <- c(0.0001, 0.001, 0.01)

  minRisk <- Inf
  minh <- 0
  for (i in 1:length(h.grids)){
    h <- h.grids[i]
    delta <- U + (var+var*alpha^2)*ksd_MahDist_est(U, s, h, var+var*alpha^2)$w.hat
    risk <- mse(delta, V)
    if (risk < minRisk){
      minRisk <- risk
      minh <- h
    }
  }

  return(minh)
}


TFSIKSD_MahDist.MCV.func <- function(y,s, var)
{
  h <- bd.est(y, s, var)

  delta <- y + var * ksd_MahDist_est(y, s, h, var)$w.hat

  return(delta)
}

# Unbasied risk estimator for estimating kernel bandwidth h for univariate Tweedie
risk.estimator <- function(y, delta, var)
{
  n <- length(y)
  risk <- -sum((delta-y)^2)/n + sum(var)/n
  return(risk)
}

JS.func <- function(y, sigma){
  m <- length(y)
  delta <- y - (m-2)*sigma^2/sum(y^2)*y
  return(delta)
}



#
EBCF.func <- function(y, s, var) {
  s <- matrix(s)
  n <- length(y)
  est <- rep(0, n)
  kfolds <- createFolds(y, k=5)

  for (fold_data in kfolds) {
    y_train <- y[-fold_data]
    s_train <- s[-fold_data, ]
    
    y_test <- y[fold_data]
    s_test <- s[fold_data, ]
    
    data <- matrix(s_train)
    bst <- xgboost(data = data, label = y_train, nrounds=30, objective = "reg:squarederror", verbose = 0)
    test_data <- matrix(s_test)
    m_hat <- predict(bst, test_data)
    
    # estinate A
    A_hat <- mse(m_hat, y_test) - var
    A_hat <- pmax(A_hat, 0.0001)
    pr <- A_hat/(var + A_hat )
    t <- pr * y_test + (1 - pr) * m_hat
    est[fold_data] <- t
  }
  
  return(est)
}


discrete_mand_dist <- function(x) {
  
  n <- length(x)
  
  
  
  dists <- matrix(rep(0, n*n), nrow=n, ncol=n)
  
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (x[i] != x[j]) {
        dists[i, j] = 1
      }
      
    }
  }
  return(dists)
}




