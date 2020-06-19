require(mnormt)
require(ggplot2)

#  Setup mean function
GPmean <- function(x,a=0,b=0,c=0) a+b*x+c*x*x

#  and a covariance function

GPexppcov <- function(x1,x2,sigma=1,delta=.1,p=2) sigma*(exp((-abs(x1-x2)^p)/delta))
GPMatern <- function(x1,x2,sigma=1,delta=1,nu=1.5) sigma*(2^(1-nu)/gamma(nu))*(sqrt(2*nu)*abs(x1-x2)/delta)^nu*besselK(sqrt(2*nu)*abs(x1-x2)/delta,nu)

# sample from a GP

# produce x the values where we are going to evaluate the GP

y=c(1:100)/10

condx <- c(20,40,60,80)

# Calculate the prior covariance of these points


covmaty <- matrix(nrow=length(y),ncol=length(y))
for (i in 1:length(y)){
  for (j in i:length(y)){
    #    covmaty[i,j] <- GPexppcov(y[i],y[j],delta=0.15 ,p=1.5)
    covmaty[i,j] <- GPMatern(y[i],y[j],delta=1.5,nu=2.5)
    if (is.na(covmaty[i,i]) == T) covmaty[i,i] <- 1
    covmaty[j,i] <- covmaty[i,j]
  }
}

# add nugget

# nugget variance
sigman=0.01

covmaty <- covmaty + diag(length(y))*sigman

#  extract the data

sigzz <- matrix(nrow=length(y)-4,ncol=length(y)-4)

sigyy <- matrix(nrow=4,ncol=4)
sigzy <- matrix(nrow=length(y)-4, ncol=4)

sigzz <- covmaty[-condx,-condx]
sigyy <- covmaty[condx,condx]
sigzy <- covmaty[-condx,condx]

sigyyinv <- chol2inv(chol(sigyy))

# generate some data

muy <- GPmean(y[condx],b=0)
p <- 1
data <- rmnorm(p,muy,sigyy)

#  the number of realisations

p <- 4

mu <- GPmean(y[-condx])
mu <- mu + t((sigzy) %*% sigyyinv %*% as.matrix(data-muy))
mu <- as.vector(mu)
condvar <-matrix(nrow=length(y)-4,ncol=length(y)-4)

condvar <-sigzz-(sigzy %*% sigyyinv %*% t(sigzy))

# Recreate the mean and variance in a dataframe

meanvar <- as.data.frame(matrix(nrow=100,ncol=2))
names(meanvar) <- c('mean','var')

meanvar$mean[-condx] <- mu
meanvar$mean[condx] <- data

meanvar$var[-condx] <- diag(condvar)
meanvar$var[condx] <- c(rep(sigman,4))

#  plot the mean and +/- 2 sd

ggplot(data=meanvar) + geom_line(aes(x=seq(1,100,1),y=mean)) +
  geom_line(aes(x=seq(1,100,1),y=mean + 2*sqrt(var)),colour='red',linetype='dotted') +
  geom_line(aes(x=seq(1,100,1),y=mean - 2*sqrt(var)),colour='red',linetype='dotted') +
  xlab('x')

#  Sample from the conditional GP

real <- rmnorm(p,mu,condvar)

#  Put these into a dataframe with the conditioning values

realised <- as.data.frame(matrix(nrow=100,ncol=p))
realised[-condx,]<-t(real)
realised[condx,] <- data
names(realised) <- c('sample1','sample2','sample3','sample4')
#  Now plot the realisations 

ggplot(data=realised) + 
  geom_line(aes(x=seq(1,100,1),y=sample1)) + 
  geom_line(aes(x=seq(1,100,1),y=sample2),colour='red') + 
  geom_line(aes(x=seq(1,100,1),y=sample3),colour='green') +
  geom_line(aes(x=seq(1,100,1),y=sample4),colour='blue') +
  xlab('x') + ylab('y')

