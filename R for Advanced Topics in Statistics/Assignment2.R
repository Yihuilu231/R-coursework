install.packages("chron")
install.packages("date")
install.packages("naniar")
install.packages("padr")
install.packages("zoo")
install.packages("imputeTS")
install.packages("timeSeries")

library(chron)
library(date)
library(lubridate)
library(dplyr)
library(ggplot2)
library(naniar)
library(padr)
library(zoo)
library(imputeTS)
library(tidyverse)
library(timeSeries)
library(spdep)
library(sf)
library(CARBayes)
library(rgdal)
library(rgeos)
library(RColorBrewer)
library(R2jags)
library(coda)
library(lattice)

London_Pollution <- read.csv("London_Pollution.csv")
# Have a look at the dataset
summary(London_Pollution)
sum(is.na(London_Pollution$Bloomsbury))
sum(is.na(London_Pollution$Barking))
London_Pollution$Date
London_Pollution$Bloomsbury[1]
Date <- as.data.frame(London_Pollution$Date)
# Extract year from the dataset
year <- format(as.Date(London_Pollution$Date, format="%d/%m/%Y"),"%Y")
NewLondon_Pollution <- cbind(year,London_Pollution)
head(NewLondon_Pollution)
year[1]

# Bloomsburry area
N2000 <- 0
N2001 <- 0
N2002 <- 0
N2003 <- 0
N2004 <- 0
N2005 <- 0
for (i in 1:1827) {
  if(year[i]==2000&is.na(NewLondon_Pollution$Bloomsbury[i])==TRUE)
    N2000 <- N2000+1
  else if(year[i]==2001&is.na(NewLondon_Pollution$Bloomsbury[i])==TRUE)
    N2001 <- N2001+1
  else if(year[i]==2002&is.na(NewLondon_Pollution$Bloomsbury[i])==TRUE){
    N2002 <- N2002+1
  }
  else if(year[i]==2003&is.na(NewLondon_Pollution$Bloomsbury[i])==TRUE){
    N2003 <- N2003+1
  }
  else if(year[i]==2004&is.na(NewLondon_Pollution$Bloomsbury[i])==TRUE){
    N2004 <- N2004+1
  }
  else if(year[i]==2005&is.na(NewLondon_Pollution$Bloomsbury[i])==TRUE){
    N2005 <- N2005+1
  }
}
# Checking the answer
N2000+N2001+N2002+N2003+N2004+N2005

# Barking area
BarkingN2000 <- 0
BarkingN2001 <- 0
BarkingN2002 <- 0
BarkingN2003 <- 0
BarkingN2004 <- 0
BarkingN2005 <- 0
for (i in 1:1827) {
  if(year[i]==2000&is.na(NewLondon_Pollution$Barking[i])==TRUE)
    BarkingN2000 <- BarkingN2000+1
  else if(year[i]==2001&is.na(NewLondon_Pollution$Barking[i])==TRUE)
    BarkingN2001 <- BarkingN2001+1
  else if(year[i]==2002&is.na(NewLondon_Pollution$Barking[i])==TRUE){
    BarkingN2002 <- BarkingN2002+1
  }
  else if(year[i]==2003&is.na(NewLondon_Pollution$Barking[i])==TRUE){
    BarkingN2003 <- BarkingN2003+1
  }
  else if(year[i]==2004&is.na(NewLondon_Pollution$Barking[i])==TRUE){
    BarkingN2004 <- BarkingN2004+1
  }
  else if(year[i]==2005&is.na(NewLondon_Pollution$Barking[i])==TRUE){
    BarkingN2005 <- BarkingN2005+1
  }
}
# Checking the answer
BarkingN2000+BarkingN2001+BarkingN2002+BarkingN2003+BarkingN2004+BarkingN2005

# Question 2 -- Plot the PM10 measurements against time for the two sites
# exploring missing data
vis_miss(London_Pollution)
# The plot shows that Bloomsbury missed 23.37% data and 
# Barking missed 12.7% data
gg_miss_upset(London_Pollution)
ggplot(aes(x=Date,y=Bloomsbury),data = NewLondon_Pollution) + geom_miss_point()+
  labs(x="Date",y="PM10 measurements",title = "Bloomsbury")

class(London_Pollution$Date)
# Change the time format
Time <- as.Date(parse_date_time(London_Pollution$Date,"dmy"))
Bloomsbury <- na.ma(London_Pollution$Bloomsbury, k = 4, weighting = "exponential")
Precip1 <- data.frame(Time, Bloomsbury, check.rows=TRUE)
ggplot(data = Precip1, mapping= aes(x= Time, y= Bloomsbury)) + geom_line()+
  labs(title = "Bloomsbury")
# Barking
Barking <- na.ma(London_Pollution$Barking, k = 4, weighting = "exponential")
Precip2 <- data.frame(Time, Barking, check.rows=TRUE)
ggplot(data = Precip2, mapping= aes(x= Time, y= Barking)) + geom_line()+
  labs(title = "Barking")

tsAirgap
plotNA.distribution(tsAirgap)
plotNA.distribution(London_Pollution$Bloomsbury)
plotNA.distribution(London_Pollution$Barking)
# Convert the Bloomsbury column into a ts object, 
# specifying on what date the measurements start
London_Pollution$Date <- as.Date(London_Pollution$Date,format="%d/%m/%Y")
Bloomsbury <- ts(London_Pollution$Bloomsbury,start=c(2000,as.numeric(format(London_Pollution$Date[1],"%j"))),frequency=365)
plotNA.distribution(Bloomsbury)
# Barking
Barking <- ts(London_Pollution$Barking,start=c(2000,as.numeric(format(London_Pollution$Date[1],"%j"))),frequency=365)
plotNA.distribution(Barking)

# Question 3
# Reading in Scotland shapefiles
London <- readOGR(dsn = '.',
                    layer = 'London')
# DataFrame
Site <- c("Bloomsburry","Barking")
Easting <- c(530123,548030)
Northing <- c(182014,183363)
London_data <- data_frame(Site,Easting,Northing)
# Convert shapefiles into something ggplot can work with
# i.e. use the st_as_sf function
Ind <- st_as_sf(London)
ggplot(Ind)+geom_sf(colour=NA)+geom_point(data=London_data,aes(x=Easting,y=Northing))

# Question 4
# Prepare the dataset
London_Pollution <- read.csv("London_Pollution.csv")
head(London_Pollution)
Bloomsbury <-London_Pollution$Bloomsbury
Bloomsbury <- Bloomsbury[1:1461]
head(Bloomsbury)

N <- 1461
jags.mod <- function(){
  B.pred[1] ~ dnorm(0, 1.0E-3)
  for (i in 2 : N) {
    Bloomsbury[i] ~ dnorm(B.pred[i],tau.v)
    B.pred[i] ~ dnorm( B.pred[i-1], tau.w)
  }
  # priors
  tau.w ~ dgamma(1,0.01)
  sigma.w2 <- 1/tau.w
  tau.v ~ dgamma(1,0.01)
  sigma.v2 <- 1/tau.v
}
jags.data <- list("Bloomsbury","N")
# parameters we want to monitor
jags.param <- c("Bloomsbury","sigma.w2","tau.w","sigma.v2","tau.v","B.pred")
# Initial values
miss.v <- is.na(Bloomsbury)
b.init1 <- rep(1,times=N)
b.init1[miss.v==FALSE] <- NA
b.init1[miss.v==TRUE] <- 22
b.init2 <- rep(1,times=N)
b.init2[miss.v==FALSE] <- NA
b.init2[miss.v==TRUE] <- 20
B.pred.inits1 = rep(20, N)
B.pred.inits2 = rep(20, N)
inits1 <- list( "tau.w" = 1, "tau.v" = 1,"Bloomsbury" = b.init1, "B.pred" = B.pred.inits1)
inits2 <- list( "tau.w" = 1, "tau.v" = 1,"Bloomsbury"=b.init2 , "B.pred" = B.pred.inits2)
jags.inits <- list(inits1, inits2)
# Fit JAGS
jags.mod.fit <- jags(data = jags.data, inits = jags.inits,
                     parameters.to.save = jags.param, n.chains = 2, n.iter = 10000,
                     n.burnin = 5000,n.thin=1,model.file = jags.mod)
# Look at outcome - point/interval estimates
print(jags.mod.fit)
# Check convergence
jagsfit.mcmc <- as.mcmc(jags.mod.fit)
gelman.diag(jagsfit.mcmc)
gelman.diag(jagsfit.mcmc,multivariate = FALSE)
traceplot(jags.mod.fit)
# Summaries
#Bloomsburys <- jags.mod.fit$BUGSoutput$sims.list$Bloomsbury
#summary(Bloomsburys)
#sigma.w2s <- jags.mod.fit$BUGSoutput$sims.list$sigma.w2
#summary(sigma.w2s)
tau.ws <- jags.mod.fit$BUGSoutput$sims.list$tau.w
#summary(tau.ws)
#sigma.v2s <- jags.mod.fit$BUGSoutput$sims.list$sigma.v2
#summary(sigma.v2s)
#tau.ws <- jags.mod.fit$BUGSoutput$sims.list$tau.w
#summary(tau.ws)
#B.preds <- jags.mod.fit$BUGSoutput$sims.list$B.pred
#summary(B.preds)
# Convert into an MCMC object for more diagnostics
summary(jagsfit.mcmc)
plot(jagsfit.mcmc)
xyplot(jagsfit.mcmc)
densityplot(jagsfit.mcmc)
autocorr.plot(jagsfit.mcmc)
# Extract posteriors
B.preds <- jags.mod.fit$BUGSoutput$sims.list$B.pred
#Predict
N <- 1468
jags.mod <- function(){
  B.pred[1] ~ dnorm(0, 1.0E-3)
  for (i in 2 : N) {
    Bloomsbury[i] ~ dnorm(B.pred[i],tau.v)
    B.pred[i] ~ dnorm( B.pred[i-1], tau.w)
  }
  # priors
  tau.w ~ dgamma(1,0.01)
  sigma.w2 <- 1/tau.w
  tau.v ~ dgamma(1,0.01)
  sigma.v2 <- 1/tau.v
}
jags.data <- list("Bloomsbury","N")
# parameters we want to monitor
jags.param <- c("Bloomsbury","sigma.w2","tau.w","sigma.v2","tau.v","B.pred")
# Initial values
miss.v <- is.na(Bloomsbury)
b.init1 <- rep(1,times=N)
b.init1[miss.v==FALSE] <- NA
b.init1[miss.v==TRUE] <- 22
b.init2 <- rep(1,times=N)
b.init2[miss.v==FALSE] <- NA
b.init2[miss.v==TRUE] <- 20
B.pred.inits1 = rep(20, N)
B.pred.inits2 = rep(20, N)
inits1 <- list( "tau.w" = 1, "tau.v" = 1,"Bloomsbury" = b.init1, "B.pred" = B.pred.inits1)
inits2 <- list( "tau.w" = 1, "tau.v" = 1,"Bloomsbury"=b.init2 , "B.pred" = B.pred.inits2)
jags.inits <- list(inits1, inits2)
# Fit JAGS
jags.mod.fit <- jags(data = jags.data, inits = jags.inits,
                     parameters.to.save = jags.param, n.chains = 2, n.iter = 10000,
                     n.burnin = 5000,n.thin=1,model.file = jags.mod)
b.preds <- jags.mod.fit$BUGSoutput$sims.list$B.pred
b.predsmeans <- colMeans(b.preds)
b.predsmeans[,1461:1468]
# RW(2)
jags.mod <- function(){
  B.pred[1] ~ dnorm(0, 1.0E-3)
  B.pred[2] ~ dnorm(0, 1.0E-3)
  for (i in 3 : N) {
    Bloomsbury[i] ~ dnorm(B.pred[i],tau.v)
    B.pred[i] ~ dnorm(2 * B.pred[i-1] - B.pred[i-2], tau.w)
  }
  # priors
  tau.w ~ dgamma(1,0.01)
  sigma.w2 <- 1/tau.w
  tau.v ~ dgamma(1,0.01)
  sigma.v2 <- 1/tau.v
}

# data
jags.data <- list("Bloomsbury","N")
# parameters we want to monitor
jags.param <- c("Bloomsbury","sigma.w2","tau.w","sigma.v2","tau.v","B.pred")
# Initial values
miss.v <- is.na(Bloomsbury)
b.init1 <- rep(1,times=N)
b.init1[miss.v==FALSE] <- NA
b.init1[miss.v==TRUE] <- 22
b.init2 <- rep(1,times=N)
b.init2[miss.v==FALSE] <- NA
b.init2[miss.v==TRUE] <- 20
B.pred.inits1 = rep(20, N)
B.pred.inits2 = rep(20, N)
inits1 <- list( "tau.w" = 1, "tau.v" = 1,"Bloomsbury" = b.init1, "B.pred" = B.pred.inits1)
inits2 <- list( "tau.w" = 1, "tau.v" = 1,"Bloomsbury"=b.init2 , "B.pred" = B.pred.inits2)
jags.inits <- list(inits1, inits2)
# Fit JAGS
jags.mod.fit <- jags(data = jags.data, inits = jags.inits,
                     parameters.to.save = jags.param, n.chains = 2, n.iter = 10000,
                     n.burnin = 5000,n.thin=1,model.file = jags.mod)
# Look at outcome - point/interval estimates
print(jags.mod.fit)
# Check convergence
jagsfit.mcmc <- as.mcmc(jags.mod.fit)
gelman.diag(jagsfit.mcmc,multivariate = FALSE)
traceplot(jags.mod.fit)
# Summaries
Bloomsburys <- jags.mod.fit$BUGSoutput$sims.list$Bloomsbury
summary(Bloomsburys)
sigma.w2s <- jags.mod.fit$BUGSoutput$sims.list$sigma.w2
summary(sigma.w2s)
tau.ws <- jags.mod.fit$BUGSoutput$sims.list$tau.w
summary(tau.ws)
sigma.v2s <- jags.mod.fit$BUGSoutput$sims.list$sigma.v2
summary(sigma.v2s)
tau.ws <- jags.mod.fit$BUGSoutput$sims.list$tau.w
summary(tau.ws)
B.preds <- jags.mod.fit$BUGSoutput$sims.list$B.pred
summary(B.preds)

# Convert into an MCMC object for more diagnostics
jagsfit.mcmc <- as.mcmc(jags.mod.fit)
summary(jagsfit.mcmc)
plot(jagsfit.mcmc)
xyplot(jagsfit.mcmc)
densityplot(jagsfit.mcmc)
autocorr.plot(jagsfit.mcmc)

# Predict the first week in 2004 for both models





