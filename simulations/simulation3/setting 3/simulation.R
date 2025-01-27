# Simulation of Tweedie's formula with side information

source("funcs.R")

source("NIT.R")
suppressMessages(suppressWarnings(library(CVXR)))
library(Rfast)
library(mvtnorm)
library(Metrics)
library(reshape)
library(ggplot2)
time1 <- Sys.time()

set.seed(5)

nrep <- 100
sd.vec <- seq(from=0.1, to=1,by=0.1)
n <- 1000
np <- length(sd.vec)


# Mean square error of different methods
JS.mse <- rep(0, np)
TFOR.mse <- rep(0, np)
TFKER.mse <- rep(0, np)
TFSIOR.mse <- rep(0, np)
TFSIKSD_MahDist.mse <- rep(0, np)
EBCF.mse <- rep(0, np)
NPMLE.mse <- rep(0, np)
TFSIKSD_one.mse <- rep(0, np)

JS.se <- matrix(rep(0, nrep*np), np, nrep)
TFOR.se <- matrix(rep(0, nrep*np), np, nrep)
TFKER.se <- matrix(rep(0, nrep*np), np, nrep)
TFSIOR.se <- matrix(rep(0, nrep*np), np, nrep)
EBCF.se <- matrix(rep(0, nrep*np), np, nrep)
NPMLE.se <- matrix(rep(0, nrep*np), np, nrep)
TFSIKSD_MahDist.se <- matrix(rep(0, nrep*np), np, nrep)
TFSIKSD_one.se <- matrix(rep(0, nrep*np), np, nrep)


for(i in 1:np)
{
  sd1 <- 0.5
  sd2 <- sd.vec[i]
  se <- 2
  for (j in 1:nrep)
  { 
    #print(j)
    eta <- 2*rbinom(n, 1, 0.5)
    
    mu1 <- rep(0, n)
    mu2 <- rep(0, n)
    mu3 <- rep(0, n)
    mu4 <- rep(0, n)
    mu5 <- rep(0, n)
    
    mu1[1:500] <- eta[1:500] 
    mu1 <-  mu1 + rnorm(n,0,sd1)
    
    mu2[501:1000] <- eta[501:1000] 
    mu2 <-  mu2 + rnorm(n,0,sd1)
    
    mu3[1:500] <- eta[1:500] 
    mu3 <-  mu3 + rnorm(n,0,sd1)
    
    mu4[501:1000] <- eta[501:1000]
    mu4<- mu4 + rnorm(n,0,sd1)
    
    mu5[1:500] <-  eta[1:500]
    mu5<- mu5 + rnorm(n,0,sd1)
    
    
    x1 <- mu1 + rnorm(n,0, 1)

    x2 <- mu2 + rnorm(n,0,sd2)

    x3 <- mu3 + rnorm(n,0,sd2)

    x4 <- mu4 + rnorm(n,0,sd2)
    
    x5 <- mu5 + rnorm(n,0,sd2)
    
    
    muy <- mu1
    
    y <- x1
    s1 <- (x2 +x3+x4+x5)/4
    s <- cbind(x2, x3,x4, x5)
    
    
    
    
    NM.delta <- y
    TFOR.delta <- TFOR.simul1.func(y, 1+sd1^2, 1,se, 0.5)
    
    TFKER.delta <- TFKER.func_2(y, 1)
    
    TFSIOR.delta <- TFSIOR.simul1.func(y, s, 1+sd1^2, sd1^2+sd2^2, 1, se, 0.5)
    
    EBCF.delta <- EBCF.func(y, s, 1)
    
    NPMLE.density <- GLmix(y, sigma=1)
    NPMLE.delta <- predict(NPMLE.density, y)
    
    JS.delta <- JS.func(y, 1)
    
    TFSIKSD_one.delta <- TFSIKSD_MahDist.MCV.func(y,s1, 1) 
    TFSIKSD_MahDist.delta <- TFSIKSD_MahDist.MCV.func(y,s, 1) 
    
    
    JS.se[i,j] <- mse(JS.delta, muy)
    EBCF.se[i,j] <- mse(EBCF.delta, muy)
    NPMLE.se[i,j] <- mse(NPMLE.delta, muy)
    TFOR.se[i,j] <- mse(TFOR.delta, muy)
    TFKER.se[i,j] <- mse(TFKER.delta, muy)
    TFSIOR.se[i,j] <- mse(TFSIOR.delta, muy)
    TFSIKSD_MahDist.se[i,j] <- mse(TFSIKSD_MahDist.delta, muy)
    TFSIKSD_one.se[i,j] <- mse(TFSIKSD_one.delta, muy)

    
  }
  
  JS.mse[i] <- mean(JS.se[i, ])
  EBCF.mse[i] <- mean(EBCF.se[i, ])
  NPMLE.mse[i] <- mean(NPMLE.se[i, ])
  TFOR.mse[i] <- mean(TFOR.se[i, ])
  TFKER.mse[i] <- mean(TFKER.se[i, ])
  TFSIOR.mse[i] <- mean(TFSIOR.se[i, ])
  TFSIKSD_one.mse[i] <- mean(TFSIKSD_one.se[i, ])
  TFSIKSD_MahDist.mse[i] <- mean(TFSIKSD_MahDist.se[i, ])
}


mse<-cbind(JS.mse, EBCF.mse, NPMLE.mse, TFKER.mse, TFSIOR.mse, TFSIKSD_one.mse, TFSIKSD_MahDist.mse)
time2 <- Sys.time()
print(Sys.time()-time1)


dat = cbind(sd.vec, mse)
df <- as.data.frame(dat)
names(df) <- c("sd","JS", "EBCF", "NPMLE", "EBT", "NIT.OR","NIT1.DD", "NIT.DD" )
meltR = melt(df, id = "sd")

plt <- ggplot(meltR, aes(x=sd, y = value, group = variable, colour = variable)) + geom_point(aes(shape=variable))+scale_shape_manual(values=c(16,17,15,3,7,4, 8))+  geom_line(size=0.4) +theme_bw()
plt <- plt + ylab("MSE") + xlab(expression(sigma[s])) 
plt + theme(legend.title = element_blank(),legend.position ='top',legend.text=element_text(size=9),legend.margin=margin(0,0,-13,0),legend.background = element_blank(),legend.key = element_blank())


ggsave('figure.png')
