# Simulation of Tweedie's formula with side information

source("funcs.R")

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


nrep <- 100
k.vec <- seq(from=50, to=450,by=50)
n <- 1000
np <- length(k.vec)
print(n)

# Mean square error of different methods
JS.mse <- rep(0, np)
TFOR.mse <- rep(0, np)
TFKER.mse <- rep(0, np)
TFSIOR.mse <- rep(0, np)
TFSIKSD_MahDist.mse <- rep(0, np)
EBCF.mse <- rep(0, np)
NPMLE.mse <- rep(0, np)

JS.se <- matrix(rep(0, nrep*np), np, nrep)
TFOR.se <- matrix(rep(0, nrep*np), np, nrep)
TFKER.se <- matrix(rep(0, nrep*np), np, nrep)
TFSIOR.se <- matrix(rep(0, nrep*np), np, nrep)
EBCF.se <- matrix(rep(0, nrep*np), np, nrep)
NPMLE.se <- matrix(rep(0, nrep*np), np, nrep)
TFSIKSD_MahDist.se <- matrix(rep(0, nrep*np), np, nrep)

for(i in 1:np)
{
  k <- k.vec[i]
  for (j in 1:nrep)
  { 

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
    y <- U-V
    s <- U+V
    
    
    TFOR.delta <- TFOR.simul2.func(y, k, n, 2)
    
    TFKER.delta <- TFKER.func_2(y, 2)
    
    TFSIOR.delta <- TFSIOR.simul2.func(y, s, k, n, 2)

    EBCF.delta <- EBCF.func(y, s, 2)
    
    NPMLE.density <- GLmix(y, sigma=sqrt(2))
    NPMLE.delta <- predict(NPMLE.density, y)
    
    JS.delta <- JS.func(y, sqrt(2))
    
    
    TFSIKSD_MahDist.delta <-TFSIKSD_MahDist.MCV.func(y,s, 2)
    
    
    JS.se[i,j] <- mse(JS.delta, muy)
    EBCF.se[i,j] <- mse(EBCF.delta, muy)
    NPMLE.se[i,j] <- mse(NPMLE.delta, muy)
    TFOR.se[i,j] <- mse(TFOR.delta, muy)
    TFKER.se[i,j] <- mse(TFKER.delta, muy)
    TFSIOR.se[i,j] <- mse(TFSIOR.delta, muy)
    TFSIKSD_MahDist.se[i,j] <- mse(TFSIKSD_MahDist.delta, muy)
    
  }
  
  JS.mse[i] <- mean(JS.se[i, ])
  EBCF.mse[i] <- mean(EBCF.se[i, ])
  NPMLE.mse[i] <- mean(NPMLE.se[i, ])
  TFOR.mse[i] <- mean(TFOR.se[i, ])
  TFKER.mse[i] <- mean(TFKER.se[i, ])
  TFSIOR.mse[i] <- mean(TFSIOR.se[i, ])
  TFSIKSD_MahDist.mse[i] <- mean(TFSIKSD_MahDist.se[i, ])
}

# (NM.mse, TFOR.mse, TFKER.mse, TFSIOR.mse, TFSIKSD_1dkernel.mse, TFSIKSD_2dkernel.mse, TFSIKSD_MahDist.mse)
mse<-cbind(JS.mse, EBCF.mse, NPMLE.mse, TFKER.mse, TFSIOR.mse, TFSIKSD_MahDist.mse)
time2 <- Sys.time()
print(Sys.time()-time1)

dat = cbind(k.vec, mse)
df <- as.data.frame(dat)
names(df) <- c("sd","JS", "EBCF", "NPMLE", "EBT", "NIT.OR","NIT.DD" )
meltR = melt(df, id = "sd")

plt <- ggplot(meltR, aes(x=sd, y = value, group = variable, colour = variable)) + geom_point(aes(shape=variable))+  geom_line(size=0.4) +theme_bw()
plt <- plt + ylab("MSE") + xlab("k") 
plt + theme(legend.title = element_blank(),legend.position ='top',legend.text=element_text(size=10),legend.margin=margin(0,0,-13,0),legend.background = element_blank(),legend.key = element_blank())

ggsave('figure.png')
