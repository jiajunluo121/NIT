# example of Tweedie's formula with side information

source("NIT.R")
suppressMessages(suppressWarnings(library(CVXR)))
library(Rfast)
library(mvtnorm)
library(Metrics)
library(REBayes)
library(reshape)
library(ggplot2)

set.seed(5)
time1 <- Sys.time()


n <- 1000
np <- length(k.vec)

k <- 50
# two sample setup
mu1 <- rep(0, n)
mu2 <- rep(0, n)
mu1[1:k] <- 2.5
be <- k + 1 
end <- 2*k
mu1[be: end] <- 1
mu2[1:k] <- 1
mu2[be:end] <- 1

n1 <- rnorm(n,0,sd=1)
n2 <- rnorm(n,0,sd=1)
U <- mu1 + n1
V <- mu2 + n2
muy <- mu1-mu2
variance <- 2

y <- U-V
s <- U+V


TFSIKSD_MahDist.delta <-TFSIKSD_MahDist.MCV.func(y,s, variance)
# Mean square error
TFSIKSD_MahDist.mse <- mse(TFSIKSD_MahDist.delta, muy)
    
