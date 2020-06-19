rm(list = ls())

Ohio_Data<- read.csv("Ohio_Data.csv")
head(Ohio_Data)
summary(Ohio_Data)
SMRs <- Ohio_Data$Obs/Ohio_Data$Exp
SMRs
Ohio_Data$X
Ohio_Data <- cbind(Ohio_Data,SMRs)
Ohio_Data
hist(Ohio_Data$SMRs,breaks = 40)
plot(Ohio_Data$X,SMRs)

source("OhioMap.R")
library(maps)
map.text("county", "ohio")
SMRs <- runif(88) # need to read in the OhioMap function
SMRs
OhioMap(SMRs,ncol=8,type="e",figmain="Ohio SMRs for each county",lower=0,upper=2)

# Fit the model
library(tidyverse)
library(R2jags)
library(coda)
library(lattice)

# We have observations from 88 counties, hence the for loop.
jags.mod <- function(){
  for (i in 1:88) {
    Obs[i]~dpois(mu[i]) # Poisson likelihood
    log(mu[i]) <- log(Exp[i])+b0+log(theta[i])
    RR[i]<- exp(b0)*theta[i]
  }
  # priors
  b0~dunif(-100,100)
  alpha~ dgamma(1,1)
  for (k in 1:88) {
    theta[k]~dgamma(alpha,alpha)
  }
}

dat <- Ohio_Data
dat
X <- dat$X
Obs <- dat$Obs
Exp <- dat$Exp
SMRs <- dat$SMRs
jags.data <- list("Obs","Exp","SMRs")

# Set the parameters we want to monitor
jags.param <- c("b0","RR")

# Specify initial values
inits1 <- list("b0" = 0.2,"alpha" = 10)
inits2 <- list("b0" = 10,"alpha" = 2000)
jags.inits <- list(inits1,inits2)

# Fit the JAGS model,n.thin is thinning rate, 
# a positive integer, used for generating sims.array
jags.mod.fit <- jags(data = jags.data, inits = jags.inits,
                     parameters.to.save = jags.param, n.chains = 2, n.iter = 10000,
                     n.burnin = 5000,n.thin=1,model.file = jags.mod)

# Look at oucome, summaries, trace plots, densities
# This gives us summary statistics of the posterior of the monitored nodes (mean, sd, quantiles).
print(jags.mod.fit)
plot(jags.mod.fit)
traceplot(jags.mod.fit)
# Check conbergence
jagsfit.mcmc <- as.mcmc(jags.mod.fit)
gelman.diag(jagsfit.mcmc)

# Extract RR and calcluate its means
RR <- as.data.frame(jags.mod.fit$BUGSoutput$sims.list$RR)
head(RR)
RRmeans <- colMeans(RR)
RRmeans
summary(jagsfit.mcmc)

# Mapping RR means
source("OhioMap.R")
library(maps)
map.text("county", "ohio")
RRmeans <- runif(88) # need to read in the OhioMap function
RRmeans
OhioMap(RRmeans,ncol=8,type="e",figmain="RRmeans",lower=0,upper=2)

# Probabilities RR in each area exceed 1.2
# To estimate the probabilities, we have to modify the model definition
jags.mod2 <- function(){
  for (i in 1:88) {
    Obs[i]~dpois(mu[i]) # Poisson likelihood
    log(mu[i]) <- log(Exp[i])+b0+log(theta[i])
    RR[i]<- exp(b0)*theta[i]
  }
  # priors
  b0~dunif(-100,100)
  alpha~ dgamma(1,1)
  for (k in 1:88) {
    theta[k]~dgamma(alpha,alpha)
  }
  #rank
  rankRR <- rank(RR[])
  p1.2 <- ifelse(RR>1.2,1,0)
}

jags.param <- c("rankRR","p1.2")

jags.mod.fit <- jags(data = jags.data, inits = jags.inits,
                     parameters.to.save = jags.param, n.chains = 2, n.iter = 10000,
                     n.burnin = 5000,n.thin=1,model.file = jags.mod2)
# Ordered ranks
r <- as.data.frame(jags.mod.fit$BUGSoutput$sims.list$rankRR)

# Probabilities
p1.2 <- colMeans(jags.mod.fit$BUGSoutput$sims.list$p1.2)
p1.2

# Mapping the Probabilities
source("OhioMap.R")
library(maps)
map.text("county", "ohio")
p1.2 <- runif(88) # need to read in the OhioMap function
OhioMap(p1.2,ncol=8,type="e",figmain="Probabilities for RR exceeds 1.2",lower=0,upper=2)

# Repeat the analysis with different priors for p(b0) and alpha
# I set b0 to be standard uniform distributions and alpha to be poisson distributions
jags.mod3 <- function(){
  for (i in 1:88) {
    Obs[i]~dpois(mu[i]) # Poisson likelihood
    log(mu[i]) <- log(Exp[i])+b0+log(theta[i])
    RR[i]<- exp(b0)*theta[i]
  }
  # priors
  b0~dnorm(0,1)
  alpha~ dpois(2)
  for (k in 1:88) {
    theta[k]~dgamma(alpha,alpha)
  }
}

dat <- Ohio_Data
dat
X <- dat$X
Obs <- dat$Obs
Exp <- dat$Exp
SMRs <- dat$SMRs
jags.data <- list("Obs","Exp","SMRs")

# Set the parameters we want to monitor
jags.param <- c("b0","RR")

# Specify initial values
inits1 <- list("b0" = 0.2,"alpha" = 10)
inits2 <- list("b0" = 10,"alpha" = 2000)
jags.inits <- list(inits1,inits2)

# Fitting the model
jags.mod.fit <- jags(data = jags.data, inits = jags.inits,
                     parameters.to.save = jags.param, n.chains = 2, n.iter = 10000,
                     n.burnin = 5000,n.thin=1,model.file = jags.mod3)

RR <- as.data.frame(jags.mod.fit$BUGSoutput$sims.list$RR)
head(RR)
RRmeans <- colMeans(RR)
# Mapping RR means
source("OhioMap.R")
library(maps)
map.text("county", "ohio")
RRmeans <- runif(88) # need to read in the OhioMap function
OhioMap(RRmeans,ncol=8,type="e",figmain="RRmeans for changed priors model",lower=0,upper=2)

jagsfit.mcmc <- as.mcmc(jags.mod.fit)
summary(jagsfit.mcmc)


