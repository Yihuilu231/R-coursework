require(mnormt)
require(ggplot2)

#  Setup mean function
GPmean <- function(x,a=0,b=0,c=0) a+b*x+c*x*x

#  and a covariance function

GPexppcov <- function(x1,x2,sigma=1,delta=.1,p=2) sigma*(exp((-abs(x1-x2)^p)/delta))
GPMatern <- function(x1,x2,sigma=1,delta=.1,nu=1.5) sigma*(2^(1-nu)/gamma(nu))*(sqrt(2*nu)*abs(x1-x2)/delta)^nu*besselK(sqrt(2*nu)*abs(x1-x2)/delta,nu)

# sample from a GP

# produce x the values where we are going to evaluate the GP

x <- c(1:100)/10

mu <- GPmean(x,b=1,c=.5)
covmat <- matrix(nrow=length(x),ncol=length(x))
for (i in 1:length(x)){
  for (j in i:length(x)){
    covmat[i,j] <- GPexppcov(x[i],x[j],delta=0.5,p=1.9)
    #   covmat[i,j] <- GPMatern(x[i],x[j],delta=0.25,nu=3.5)
    if (is.na(covmat[i,i]) == T) covmat[i,i] <- 1
    covmat[j,i] <- covmat[i,j]
  }
}

#  Generate p realisations
p=5
real <- rmnorm(p,mu,covmat)

#  Now plot the realisations 

# first create a data frame

df <- data.frame(cbind(as.matrix(x),t(real)))
names(df) <- c('x','realisation_1','realisation_2','realisation_3','realisation_4','realisation_5')
ggplot(df) + geom_line(aes(x=x,y=realisation_1)) + geom_line(aes(x=x,y=realisation_2),colour='red')+ geom_line(aes(x=x,y=realisation_3),colour='green')+ geom_line(aes(x=x,y=realisation_4),colour='blue') + geom_line(aes(x=x,y=realisation_5),colour='orange')


#plt1 <- ggplot(df) + geom_line(aes(x=x,y=X2)) 

