########################################################
### Illusrate scatter plot smoothers (moving average) ##
########################################################

## load the library containing function to simulate appropriate data
library(mgcv) 
## this ensures the same data are simulated every time
set.seed(2) 
## simulate some data, a respose y and some predictors x that are non-lineraly related to y
dat <- gamSim(1,n=400,dist="normal",scale=2)
## Set up plotting parameters
#x11(width=8,height=6)
par(mar = c(4, 4, 1, 1),cex=1.2,lwd=2)
## Plot y against x
plot(dat$x2,dat$y,pch=20,xlab="x",ylab="y",lwd=1)

## Estimate the relationship using kernel smoothing with bandwidth=1
KS <- ksmooth(dat$x2,dat$y,bandwidth=1,n.points=500,kernel="normal")
## superimpose the kernel smoother estimate on the existing plot
lines(KS$x,KS$y,col="red",lwd=3)
## try smaller bandwidth to increase flexibility
KS <- ksmooth(dat$x2,dat$y,bandwidth=0.01,n.points=500,kernel="normal")
lines(KS$x,KS$y,col="blue",lwd=3)
## try bigger bandwidth as line is too wiggly
KS <- ksmooth(dat$x2,dat$y,bandwidth=0.1,n.points=500,kernel="normal")
lines(KS$x,KS$y,col="purple",lwd=3)
## So bandwidth choice is arbitrary! Would be nice to have a way to estimate that from the data.

## Estimate the relationship using LOESS (localised regression)
plot(dat$x2,dat$y,pch=4,xlab="x",ylab="y",lwd=1)
## The loess fit with span 0.1
L <- loess(y~x2,dat,span=0.1)
## Create a sequence of x values to predict from the LOESS fit
xx <- seq(0,1,len=500)
## Add the LOESS estimate (predictions)
lines(xx,predict(L,newdata=data.frame(x2=xx)),col="blue",lwd=3)
## Increase span as the line seems too wiggly
L <- loess(y~x2,dat,span=1)
lines(xx,predict(L,newdata=data.frame(x2=xx)),col="red",lwd=3)
## Decrease span a bit to get the "bump"
L <- loess(y~x2,dat,span=0.5)
lines(xx,predict(L,newdata=data.frame(x2=xx)),col="darkgreen",lwd=3)






##############################################
### Illusrate the polynomial basis function ##
##############################################

## We will simulate data to illustrate this so we use set.seed to ensure simulated data are same evry time
set.seed(1)
## Simulate an Uniform x on [0,1] and sort the values in ascending order (to create a "relationship").
## (Don't worry about the details of simulating this data, it's just for illustration purposes)
x <- sort(runif(40)*10)^0.5
## Simulate an Uniform x on [0,1] and sort the values in ascending order
y <- sort(runif(40))^0.1
# Plot the data
x11(width=8,height=6)
par(mar = c(4, 4, 1, 1),cex=1.2)
plot(x,y,pch=4,ylim=c(0.65,1.1))

## Fit an order 5 polynomial using linear regression
model <- lm(y~poly(x,5))
## Create a sequence of x values for which to predict from the fitted model
xp <- seq(min(x),max(x),len=500)
## predict from model
preds <- predict(model,newdata=data.frame(x=xp))
## Add the predictions to the plot
lines(xp,preds,lwd=3)

## Now try to increase the order of the polynomial to 10 and see what happens
# Fit the model
model2 <- lm(y~poly(x,10))
# redict from it
preds2 <- predict(model2,newdata=data.frame(x=xp))
# and plot the predicted line
lines(xp,preds2,col="blue",lwd=3)
# which illustrates that polynomial splines can be quite unstable!




###################################################
### Illusrate the cubic spline basis function    ##
### and how this can be fitted as a linear model ##
###################################################
# Simulate some data:
library(mgcv) 
set.seed(2) 
dat <- gamSim(1,n=400,dist="normal",scale=2)
# The covariate x:
x <- dat$x2 
# The response y:
y <- dat$y
# Write two functions, to compute cubic splines as in lecture notes, and create the associated model matrix. 
# I will not expect you to understand what this are doing and I only use them here for illustration. We will 
# be using the function gam() to actually fit models.
rk <- function(x,z){ 
  abs((x-z)^3) 
}
spl.X <- function(x,xk){
  q <- length(xk)+2
  n <- length(x)
  X <- matrix(1,n,q)
  X[,2] <- x
  X[,3:q] <- outer(x,xk,FUN=rk)
  X
}
# Additive model with a cubic spline with 3 knots placed at 0.25, 0.5, 0.75
xk <- c(0.25,0.5,0.75)
X <- spl.X(dat$x2,xk) # Calculate the model matrix
model <- glm(y~X-1,data=dat,family=gaussian) # Fit the model (basically a linear model)
summary(model)
# Now estimate the mean of the model at a fine grid between 0 and 1
xp <- seq(0,1,length=500)
mu_hat <- spl.X(xp,xk)%*%coef(model) ## This is X'beta
par(mfrow=c(2,2))
plot(dat$x2,dat$y,pch=4,lwd=1,main="3 knots",xlab="x",ylab=expression(f(x)))
lines(xp,mu_hat,col="red",lwd=2)
points(xk,spl.X(xk,xk)%*%coef(model),col="blue",pch=19) # indicate where the knots are
# Increase the number of knots to 5
xk <- 1*((1:5)/6) ## 5 knots
X <- spl.X(dat$x2,xk) # model matrix
model <- glm(y~X-1,data=dat,family=gaussian)
mu_hat <- spl.X(xp,xk)%*%coef(model) ## This is X'beta
plot(dat$x2,dat$y,pch=4,lwd=1,main="5 knots",xlab="x",ylab=expression(f(x)))
lines(xp,mu_hat,col="red",lwd=2)
points(xk,spl.X(xk,xk)%*%coef(model),col="blue",pch=19)
# Increase the number of knots to 9
xk <- 1*((1:9)/10) ## 9 knots
X <- spl.X(dat$x2,xk)
model <- lm(y~X-1,data=dat)
mu_hat <- spl.X(xp,xk)%*%coef(model) ## This is X'beta
plot(dat$x2,dat$y,pch=4,lwd=1,main="9 knots",xlab="x",ylab=expression(f(x)))
lines(xp,mu_hat,col="red",lwd=2)
points(xk,spl.X(xk,xk)%*%coef(model),col="blue",pch=19)
xk <- 1*((1:11)/12) ## 11 knots
X <- spl.X(dat$x2,xk)
model <- lm(y~X-1,data=dat)
mu_hat <- spl.X(xp,xk)%*%coef(model) ## This is X'beta
plot(dat$x2,dat$y,pch=4,lwd=1,main="11 knots",xlab="x",ylab=expression(f(x)))
lines(xp,mu_hat,col="red",lwd=2)
points(xk,spl.X(xk,xk)%*%coef(model),col="blue",pch=19)





##############################################################
### Revisit question 2 on sheet 1 and fit as Additive Model ##
##############################################################

# View the data
nlmodel[1:20,]

# Recall, this data are 100 obs from the model y~N(mu,sig^2)
# where mu = (theta1*x)/(theta2 + x)

# Set up plot parameters for nicer looking plot
x11(width=8,height=6)
par(mar = c(4, 4, 1, 1),cex=1.2,lwd=2)
plot(y~x,data=nlmodel,pch=20)

# Estimates from fitting model using nlm where
theta1 <- 214.6501
theta2 <- 0.0635

# Add estimated line
xx <- seq(0,1,len=200) ## a sequence of x values on which to produce fitted values
yfit <- (theta1*xx)/(theta2 + xx) ## compute fitted values
lines(xx, yfit,lwd=2,col="red") ## add the estimated line

## Now fit using additive model

# load mgcv library
library(mgcv)
# Fit the model y~N(mu,sig^2) where mu = beta0 + f(x).
# Assume a cubic spline basis representation for f(x).
# Set the rank (basis dimension) of mu to be 4 (1 parameter for the intercept and 3 for f(x)).
# We can think of the rank as the total number of "beta" parameters in the model. I.e. this model
# has four parameters to estimate how the mean depeds on x using a smooth line.
Amodel <- gam(y~s(x,k=4,bs="cs"),data=nlmodel,family=gaussian(link="identity"))
summary(Amodel)

## Let's look at what's inside the fitted model
# Estimate of sig^2
Amodel$sig2
# The rank of f(x)
Amodel$rank
# The beta estimates
Amodel$coefficients
# Estimate of the smoothing parameter
Amodel$sp

## Produce fitted line from this model using perdict
yfitAM <- predict(Amodel,newdata=data.frame(x=xx))
lines(xx, yfitAM,lwd=2,col="blue") ## add the estimated line
## Looks like we could do with a bit more flexibility

## What about the rank of the model being large enough? 
## We can use gam.check() to check this
## and also look at residual plots
x11()
par(mfrow=c(2,2))
gam.check(Amodel,pch=20)
# k' is very close to edf so might need to increase the rank (equiv, k-index is less than 1)

## So, increase the rank to 10. This is equivalent to having 9 parameters in the smooth function. Remember that
## the idea is to have a number that is larger than we think reasonable and the smoothing parameter will ensure
## this is appropriately penalised.
Amodel2 <- gam(y~s(x,k=10,bs="cs"),data=nlmodel,family=gaussian(link="identity"))
# Produce fitted line from this model using perdict
yfitAM2 <- predict(Amodel2,newdata=data.frame(x=xx))
x11()
par(mar = c(4, 4, 1, 1),cex=1.2,lwd=2)
plot(y~x,data=nlmodel,pch=20)
lines(xx, yfit,lwd=2,col="red") ## add the estimated line from original model
lines(xx, yfitAM,lwd=2,col="blue") ## add the estimated line from Amodel
lines(xx, yfitAM2,lwd=2,col="purple") ## add the estimated line from Amodel2
legend("bottomright",c("Original","Amodel","Amodel2"),lty=c(1,1,1),col=c("red","blue","purple"))
# Looks like a better fit, would still need to use gam.check to ensure we have a big enough rank.

### Illustrate inference and effective degrees of freedom

## Go back to Amodel and check what the effective degrees of freedom are for each coefficient
Amodel$edf
# All of these look very close to 1 so maybe we need a few more knots (increase the rank)
# The model degrees of freedom (analogous to n-p-1 in GLMs) is
Amodel$df.residual
# which is the same as
100 - sum(Amodel$edf) ## (where n=100)

## So conditional on the smoothing parametere estimates, we could now do a deviance goodness of fit test
## exactly as in GLMs (although don't need to for a Gaussian of course). The scaled deviance should be chi-sq.
sc.deviance <- Amodel$deviance/Amodel$sig2
sc.deviance
Amodel$df.residual
pchisq(sc.deviance,Amodel$df.residual,lower.tail=F)
## So model is good fit CONDITIONAL ON PENALTY PARAMETERS. Of course we didn't need to do this since the
## Gaussian has a dispersion parameter so it will always fit! Would still need to do this test if the model
## was Poisson or Binomial.






###########################################################
### Revisit question 5 on sheet 1 and fit as Poisson GAM ##
###########################################################

# Plot the data
x11(width=8,height=6)
par(mar = c(4, 4, 1, 1),cex=1.2,lwd=2)
plot(cases~date,data=aids,pch=20)

## Fit Poisson GAM with a cubic spline for time with rank 4
model <- gam(cases~s(date,k=4,bs="cs"),data=aids,family=poisson(link="log"))
x11()
par(mfrow=c(2,2))
gam.check(model,pch=20)
## The edf is very close to k' so probably best to icrease the rank.

model2 <- gam(cases~s(date,k=10,bs="cs"),data=aids,family=poisson(link="log"))
x11()
par(mfrow=c(2,2))
gam.check(model2,pch=20)
## Not such a good fit, but at least we have enough degrees of freedom

## Produce estimates of the fitted line along with confidence intervals (CIs)
## sequence of time values to predict on
pred_x <- seq(min(aids$date),max(aids$date),len=100)
# Now predict the LINEAR PREDICTOR and get its standard errors. The linear predictor
# is of course Gaussian so we will calculate CIs and then exponentiate those
preds <- predict(model,newdata=data.frame(date=pred_x),
                 se.fit=T,type="link")
x11(width=8,height=6)
par(mfrow=c(1,1),mar = c(4, 4, 1, 1),cex=1.2,lwd=2)
plot(cases ~ date,data=aids,pch=20)
lines(pred_x,exp(preds$fit),col="blue",lwd=2)
lines(pred_x,exp(preds$fit+1.96*preds$se.fit),col="blue",lwd=2,lty=2)
lines(pred_x,exp(preds$fit-1.96*preds$se.fit),col="blue",lwd=2,lty=2)

## However, we forgot something impotant! It's a Poisson GAM so does it fit (conditional on penalty estimates)?
model$deviance
model$df.residual
pchisq(model$deviance,model$df.residual,lower.tail=F)
# It does not! So need alternative formulation to make it fit. (see question 5 on sheet 1)






####################################
### Play with rent data in Munich ##
####################################
munichrent[1:20,]

# Fit a linear model
model1 <- lm(rent~yearc+location+yearc:location,data=munichrent)
summary(model1)
## Check residuals
x11(width=16,height=6)
par(mfrow=c(1,2),mar = c(4, 4, 1, 1),cex=1.2,lwd=2)
plot(model1,1,pch=20)
plot(model1,2,pch=20)
# Heavily skewed standardised deviance residuals. Also, rent takes only positive values so fit a 
# Gamma GLM with identity link

model2 <- glm(rent~yearc+location+yearc:location,data=munichrent,
              family=Gamma(link="identity"))
summary(model2)
## Check residuals
x11(width=16,height=6)
par(mfrow=c(1,2),mar = c(4, 4, 1, 1),cex=1.2,lwd=2)
plot(model2,1,pch=20)
plot(model2,2,pch=20)
## Look at predicted values for each location
years <- 1918:1997
loc1 <- predict(model2,newdata=data.frame(yearc=years,location=as.factor(1)),se.fit=T)
loc2 <- predict(model2,newdata=data.frame(yearc=years,location=as.factor(2)),se.fit=T)
loc3 <- predict(model2,newdata=data.frame(yearc=years,location=as.factor(3)),se.fit=T)
x11(width=8,height=6)
par(mar = c(4, 4, 1, 1),cex=1.2,lwd=2)
plot(years,loc1$fit,type="l",ylim=c(320,900),xlab="time (years)",ylab="mean rent")
lines(years,loc1$fit+1.96*loc1$se.fit,lty=2,lwd=1)
lines(years,loc1$fit-1.96*loc1$se.fit,lty=2,lwd=1)
lines(years,loc2$fit,col="red")
lines(years,loc2$fit+1.96*loc2$se.fit,lty=2,lwd=1,col="red")
lines(years,loc2$fit-1.96*loc2$se.fit,lty=2,lwd=1,col="red")
lines(years,loc3$fit,col="blue")
lines(years,loc3$fit+1.96*loc3$se.fit,lty=2,lwd=1,col="blue")
lines(years,loc3$fit-1.96*loc3$se.fit,lty=2,lwd=1,col="blue")
legend("topleft",c("Location 1","Location 2","Location 3"),col=c("black","red","blue"),lty=c(1,1,1))
## So average rent seems to go up for all locations, but faster for location 1 and 3. However, why would we think that 
## the increase is linear? Maybe best to try a GAMMA GAM and see if there were any fluctuations in average rent

library(mgcv)
## Interactions with factors in gam() are handled with the "by" argument
model3 <- gam(rent~s(yearc, by=location, k=10),data=munichrent, 
              family=Gamma(link="identity"))
x11()
par(mfrow=c(2,2))
gam.check(model3,pch=20)
## So need a few more degrees of freedom. After trial and error, settle to k=30
model3 <- gam(rent~s(yearc, by=location, k=30),data=munichrent, 
              family=Gamma(link="identity"))
x11()
par(mfrow=c(2,2))
gam.check(model3,pch=20)
summary(model3)

## Now plot the fitted values
pred1 <- predict(model3,newdata=data.frame(yearc=years,location=1),se.fit=T)
pred2 <- predict(model3,newdata=data.frame(yearc=years,location=2),se.fit=T)
pred3 <- predict(model3,newdata=data.frame(yearc=years,location=3),se.fit=T)
x11(width=8,height=6)
par(mar = c(4, 4, 1, 1),cex=1.2,lwd=2)
plot(years,pred1$fit,type="l",ylim=c(0,1100))
lines(years,pred1$fit+2*pred1$se.fit,lty=2,lwd=1) ## approx. conf. int.
lines(years,pred1$fit-2*pred1$se.fit,lty=2,lwd=1) ## approx. conf. int.
lines(years,pred2$fit,col="red")
lines(years,pred2$fit+2*pred1$se.fit,lty=2,col="red",lwd=1) ## approx. conf. int.
lines(years,pred2$fit-2*pred1$se.fit,lty=2,col="red",lwd=1) ## approx. conf. int.
lines(years,pred3$fit,col="blue")
lines(years,pred3$fit+2*pred1$se.fit,lty=2,col="blue",lwd=1) ## approx. conf. int.
lines(years,pred3$fit-2*pred1$se.fit,lty=2,col="blue",lwd=1) ## approx. conf. int.
legend("topleft",c("Location 1","Location 2","Location 3"),lty=c(1,1,1),col=c("black","red","blue"))

## Clearly the link function needs to be changed as we are getting nonsensical values
model3 <- gam(rent~s(yearc, by=location, k=20),data=munichrent, 
              family=Gamma(link="log"))
x11()
par(mfrow=c(2,2))
gam.check(model3,pch=20)
pred1 <- predict(model3,newdata=data.frame(yearc=years,location=1),se.fit=T)
pred2 <- predict(model3,newdata=data.frame(yearc=years,location=2),se.fit=T)
pred3 <- predict(model3,newdata=data.frame(yearc=years,location=3),se.fit=T)
x11(width=8,height=6)
par(mar = c(4, 4, 1, 1),cex=1.2,lwd=2)
plot(years,exp(pred1$fit),type="l",ylim=c(0,1100))
lines(years,exp(pred1$fit+2*pred1$se.fit),lty=2,lwd=1) ## approx. conf. int.
lines(years,exp(pred1$fit-2*pred1$se.fit),lty=2,lwd=1) ## approx. conf. int.
lines(years,exp(pred2$fit),col="red")
lines(years,exp(pred2$fit+2*pred1$se.fit),lty=2,col="red",lwd=1) ## approx. conf. int.
lines(years,exp(pred2$fit-2*pred1$se.fit),lty=2,col="red",lwd=1) ## approx. conf. int.
lines(years,exp(pred3$fit),col="blue")
lines(years,exp(pred3$fit+2*pred1$se.fit),lty=2,col="blue",lwd=1) ## approx. conf. int.
lines(years,exp(pred3$fit-2*pred1$se.fit),lty=2,col="blue",lwd=1) ## approx. conf. int.
legend("topleft",c("Location 1","Location 2","Location 3"),lty=c(1,1,1),col=c("black","red","blue"))


model3 <- gam(rent~s(yearc, by=location, k=5),data=munichrent, 
              family=Gamma(link="log"))
x11()
par(mfrow=c(2,2))
gam.check(model3,pch=20)
pred1 <- predict(model3,newdata=data.frame(yearc=years,location=1),se.fit=T)
pred2 <- predict(model3,newdata=data.frame(yearc=years,location=2),se.fit=T)
pred3 <- predict(model3,newdata=data.frame(yearc=years,location=3),se.fit=T)
x11(width=8,height=6)
par(mar = c(4, 4, 1, 1),cex=1.2,lwd=2)
plot(years,exp(pred1$fit),type="l",ylim=c(0,1100))
lines(years,exp(pred1$fit+2*pred1$se.fit),lty=2,lwd=1) ## approx. conf. int.
lines(years,exp(pred1$fit-2*pred1$se.fit),lty=2,lwd=1) ## approx. conf. int.
lines(years,exp(pred2$fit),col="red")
lines(years,exp(pred2$fit+2*pred1$se.fit),lty=2,col="red",lwd=1) ## approx. conf. int.
lines(years,exp(pred2$fit-2*pred1$se.fit),lty=2,col="red",lwd=1) ## approx. conf. int.
lines(years,exp(pred3$fit),col="blue")
lines(years,exp(pred3$fit+2*pred1$se.fit),lty=2,col="blue",lwd=1) ## approx. conf. int.
lines(years,exp(pred3$fit-2*pred1$se.fit),lty=2,col="blue",lwd=1) ## approx. conf. int.
legend("topleft",c("Location 1","Location 2","Location 3"),lty=c(1,1,1),col=c("black","red","blue"))