# Simulation of Tweedie's formula with side information

source("funcs.R")

source("NIT.R")
suppressMessages(suppressWarnings(library(CVXR)))
library(Rfast)
library(mvtnorm)
library(Metrics)
library(REBayes)
time1 <- Sys.time()


nrep <- 100
var.vec <- seq(from=0.1, to=1,by=0.1)
n <- 1000
np <- length(var.vec)


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
  sd <- var.vec[i]
  for (j in 1:nrep)
  { 
    
    ep <- rbinom(n, 1, 0.5) + rnorm(n,0,1)
    
    muy <- ep + rnorm(n,0,sd)
    y <- muy + rnorm(n,0,1)
    s <- ep + rnorm(n,0,sd) + rnorm(n,0,1)
    
    
    NM.delta <- y
    TFOR.delta <- TFOR.simul1.func(y, 1, sd)
    
    TFKER.delta <- TFKER.func_2(y, 1)
    
    TFSIOR.delta <- TFSIOR.simul1.func(y, s, 1, sd)
    
    EBCF.delta <- EBCF.func(y, s, 1)
    
    NPMLE.density <- GLmix(y, sigma=1)
    NPMLE.delta <- predict(NPMLE.density, y)
    
    JS.delta <- JS.func(y, 1)
    
    sdy <- sqrt(var(y))
    sds <- sqrt(var(s))
    TFSIKSD_MahDist.delta <- TFSIKSD_MahDist.MCV.func(y,s, 1)
    
    
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

mse<-cbind(JS.mse, EBCF.mse, NPMLE.mse, TFKER.mse, TFSIOR.mse, TFSIKSD_MahDist.mse)
time2 <- Sys.time()
print(Sys.time()-time1)


dat = cbind(var.vec, mse)
df <- as.data.frame(dat)
names(df) <- c("sd","JS", "EBCF", "NPMLE", "EBT", "NIT.OR","NIT.DD" )
meltR = melt(df, id = "sd")

plt <- ggplot(meltR, aes(x=sd, y = value, group = variable, colour = variable)) + geom_point(aes(shape=variable))+  geom_line(size=0.4) +theme_bw()
plt <- plt + ylab("MSE") + xlab(expression(sigma)) 
plt + theme(legend.title = element_blank(),legend.position ='top',legend.text=element_text(size=10),legend.margin=margin(0,0,-13,0),legend.background = element_blank(),legend.key = element_blank())

ggsave('figure.png')

