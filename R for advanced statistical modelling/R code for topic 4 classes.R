###################################################################
################### Play with aircraft airconditioning data #######
###################################################################

## Input the data in R and create indicator variable for 
## observations that are censored
times <- c(143,164,188,188,190,192,206,209,213,216,220,227,
           230,234,246,265,304,216,244)
cens <- rep(1,19)
cens[18:19] <- 0

# Estimate Exponential distribution using first principles.
log_likelihood <- function(theta,times){
  ell <- sum(log(theta*exp(-theta*times[1:17]))) + 
    sum(log(exp(-theta*times[18:19])))
  return(-ell)
}
## Or use the R built-in functions to calculate the pdf and tail area
## probabilities of the Exponential distr.:
#log_likelihood <- function(theta,times,cens){
#  -sum( dexp(times[1:17],theta,log=T) ) - 
#   sum( pexp(times[18:19],theta,log=T,lower.tail=F) )
#}
out <- nlm(log_likelihood,1,times=times,hessian = T)
theta <- out$estimate
1/theta # the mean
theta
sqrt(1/(out$hessian)) # std error
# What if we simpl excluded the censored obs as missing?
library(MASS)
# Use function fitdistr to fit simple distributions using 
# maximum likelihood
fitdistr(times[1:17],densfun = "exponential")
# Note the bigger standard error (less information)


## Repeat using survreg()
## Load the survival library to do the modelling
library(survival)
## Create "Surv" object to be used with the survreg function
Times <- Surv(times,cens,type="right")
Times
##### Now use the survreg function in the survival package, 
##### noting that it defines theta=exp(-beta0)
##### where beta0 is "the intercept"
## Fit an Exponential model
Emodel <- survreg(Times~1,dist="exponential")
summary(Emodel)
## and a Weibull model 
Wmodel <- survreg(Times~1,dist="weibull")
summary(Wmodel) # scale is 1/gamma

## Plot the failure rates. For the Exponential is easy, we just 
## need the estimate of theta
## which is just exp(-beta0)
thetaE <- exp(-Emodel$coefficients)
# For the Weibull, we also need the gamma parameter. 
## What R calls the scale is actually 1/gamma
thetaW <- exp(-Wmodel$coefficients)
gamma <- 1/Wmodel$scale
## create a time sequence
tt <- seq(1,350,len=500)
x11(width=8,height=6)
plot(tt,thetaW*gamma*(thetaW*tt)^(gamma-1),type="l",lwd=3,
     ylab="failure rate",xlab="time (hours)") ## Weibull failure rate
## Exponential has constant failure rate
lines(tt,rep(thetaE,500),lwd=3,col="blue") 
legend("topleft",c("Exponential","Weibull"),col=c("blue","black"),
       lwd=c(3,3),lty=c(1,1))

## Plot the Survivor functions
x11(width=8,height=6)
plot(tt,exp(-(thetaW*tt)^gamma),type="l",lwd=3,
     ylab="Survivor function",xlab="time (hours)") ## Weibull survivor
lines(tt,exp(-thetaE*tt),lwd=3,col="blue") ## Exponential survivor
legend("topright",c("Exponential","Weibull"),col=c("blue","black"),
       lwd=c(3,3),lty=c(1,1))

## Plot the pdfs
x11(width=8,height=6)
hist(times,freq=F,xlim=c(0,350),ylab="Density",xlab="time (hours)",lwd=2) # data histogram
lines(tt,dweibull(tt,shape=gamma,scale=1/thetaW),type="l",lwd=3) ## Weibull pdf
lines(tt,dexp(tt,thetaE),lwd=3,col="blue") ## Exponential pdf
legend("topright",c("Exponential","Weibull"),col=c("blue","black"),lwd=c(3,3),lty=c(1,1))

## What about the proabability that a airconditioning 
## unit will last for more than 250 hours?
## Well just get estimates of theta and gamma, and 
## plug in formula for survivor function
# For the Exponential
exp(-thetaE*250)
# Can also use 
1-pexp(250,rate=thetaE)
# For the Weibull
exp(-(thetaW*250)^gamma) # as per slide 17
# Can also use 
1-pweibull(250,shape=gamma,scale=1/thetaW)
# but note that pweibull confusingly calls gamma 
# the shape and 1/theta the scale


lambda <- function(t,x,beta0,beta1){
  exp(-beta0-beta1*x)*exp(-(t*exp(-beta0-beta1*x)))
}
tt <- seq(0,20,len=200)
plot(tt,lambda(tt,x=0,beta0=0.1,beta1=0.1),type="l")
lines(tt,lambda(tt,x=10,beta0=0.1,beta1=0.1),type="l",col="red")


## But which model is better? The Exponential is a special case 
## of the Weibull, i.e. where gamma = 1.
## Function survreg() produces an estimate of log(1/gamma), 
## so testing for gamma=1 is same as testing 
## for log(1/gamma)=0 so we just look at p-value of log(scale) 
## in the summary of Wmodel, which is <0.05 
## and thus reject the test that it is zero and therefore stick 
## with the Weibull.

## Of course we can also use likelihood ratio test as the two models are nested
loglikE <- logLik(Emodel)
loglikW <- logLik(Wmodel)
## get rid of "labels" introduced by logLik function
loglikE <- as.numeric(logLik(Emodel))
loglikW <- as.numeric(logLik(Wmodel))
## the LRT test:
LRT <- -2*(loglikE-loglikW)
## The Weibull is only larger by one parameter so 
## we test using chi-sq with 1 d.o.f.
pchisq(LRT,1,lower.tail=F)
## Extremely low p-value so we Exponential and 
## Weibull DO NOT fit the same, so we choose 
## Weibull as it is bigger.

## Of course, we can always use the AIC:
aicE <- -2*loglikE + 2*1 ## Exp. has two parameters
aicW <- -2*loglikW + 2*2 ## Weibull has two parameters
aicE
aicW
## So Weibull better. Note the AIC() 
## function returns the same
AIC(Emodel)
AIC(Wmodel)

## So we choose Weibull model, but how does the QQ plot of the 
## residuals look?
x11()
qqnorm( residuals(Wmodel,type="deviance"),pch=20 )
qqline( residuals(Wmodel,type="deviance"),pch=20 )
## Not so convincing at the upper tails. Maybe try a 
## different distribution? See Q2 on sheet 2.








######################################################
############### Play with Lung cancer data to  #######
############### illustrate the Weibull         #######
############### accelerated failure time model #######
######################################################
# Look at the data (survival time is in days)
Lung[1:20,]

# Fit a Weibull AFT to Lung data, with Karno and gender
library(survival)
model <- survreg(Surv(time, status) ~ Karno + gender, data=Lung, 
                 dist="weibull")
summary(model)

# Plot estimates of the failure rate with conf. int.
gamma <- 1/model$scale ## gamma estimate
## Predict the linear predictor, X*beta, with 
## the associated standard error.
pred <- predict(model,newdata=data.frame(gender="Female",Karno=50),
                se.fit=T,type="lp")
## Use it to calculate theta = exp(-X*beta)
theta <- exp(-pred$fit)
## and it's confidence interval
thetaL <- exp(-(pred$fit+1.96*pred$se.fit))
thetaU <- exp(-(pred$fit-1.96*pred$se.fit))
## Create a sequence ot time values in the range of the data
tt <- seq(5,1022,len=500) 
## Use theta and its CI along with the estimate 
## of gamma to compute the failure rate
lambda <- theta*gamma*(tt*theta)^(gamma-1)
lambdaU <- thetaU*gamma*(tt*thetaU)^(gamma-1)
lambdaL <- thetaL*gamma*(tt*thetaL)^(gamma-1)
## Plot the failure rate
x11()
par(lwd=2)
plot(tt,lambda,lwd=3,type="l",ylim=c(0,0.01),
     xlab="time (days)", ylab="failure rate") 
lines(tt,lambdaU,lty=2)
lines(tt,lambdaL,lty=2)

## Create a function to automatically plot the 
## failurerate for various covariate combinations
plot.weib <- function(model,gender="Male",Karno=100,
                      add=F,col="black",ylim=0.04){
  gamma <- 1/model$scale ## gamma estimate
  ## predict the linear predictor
  pred <- predict(model,newdata=data.frame(gender=gender,Karno=Karno),
                  se.fit=T,type="lp")
  ## and use it to calculate theta
  theta <- exp(-pred$fit)
  ## and it's confidence interval
  thetaL <- exp(-(pred$fit+1.96*pred$se.fit))
  thetaU <- exp(-(pred$fit-1.96*pred$se.fit))
  tt <- seq(5,1022,len=500)
  lambda <- theta*gamma*(tt*theta)^(gamma-1)
  lambdaU <- thetaU*gamma*(tt*thetaU)^(gamma-1)
  lambdaL <- thetaL*gamma*(tt*thetaL)^(gamma-1)
  if(add==F){
    plot(tt,lambda,lwd=3,col=col,type="l",ylim=c(0,ylim),xlab="time (days)",ylab="failure rate") 
    lines(tt,lambdaU,col=col,lty=2,lwd=1)
    lines(tt,lambdaL,col=col,lty=2,lwd=1)
    text(900,max(lambda)+sd(lambda)/3,paste(gender,", K=",Karno,sep=""),col=col)
  }
  else{
    lines(tt,lambda,lwd=3,col=col)
    lines(tt,lambdaU,col=col,lty=2,lwd=1)
    lines(tt,lambdaL,col=col,lty=2,lwd=1)
    text(900,max(lambda)+sd(lambda)/3,paste(gender,", K=",Karno,sep=""),col=col)
  }
}

x11()
plot.weib(model,Karno=0) # baseline (Male with Karno=0)
plot.weib(model,Karno=0,gender="Female",col="Red",add=T)
plot.weib(model,Karno=20,gender="Male",col="blue",add=T)
#
x11()
plot.weib(model,gender="Male",Karno=0,ylim=0.04)
plot.weib(model,gender="Male",Karno=40,add=T,col="red")
plot.weib(model,gender="Male",Karno=80,add=T,col="blue")
plot.weib(model,gender="Male",Karno=100,add=T,col="green")
#
x11()
plot.weib(model,gender="Male",Karno=100,ylim=0.005)
plot.weib(model,gender="Female",Karno=100,ylim=0.01,add=T,col="red")

## function to calculate the survivor function for 
## different predictor values
calc.survivor <- function(model,gender="Male",Karno=100,Time=365){
  gamma <- 1/model$scale ## gamma estimate
  pred <- predict(model,newdata=data.frame(gender=gender,Karno=Karno),
                  se.fit=T,type="lp")
  theta <- exp(-pred$fit)
  lambda <- theta*gamma*(tt*theta)^(gamma-1)
  ## the survivor function (as in lecture slides)
  S <- exp(-(theta*Time)^gamma) 
  print(paste("Survival prob. over ",Time," days, for a ",gender," with Karno=",Karno,":",sep=""))
  print(S)
}

calc.survivor(model,gender="Male",Karno=0,Time=365)
calc.survivor(model,gender="Female",Karno=0,Time=365)
#
calc.survivor(model,gender="Female",Karno=0,Time=365)
calc.survivor(model,gender="Female",Karno=25,Time=365)
calc.survivor(model,gender="Female",Karno=75,Time=365)
calc.survivor(model,gender="Female",Karno=100,Time=365)
#
calc.survivor(model,gender="Male",Karno=100,Time=365)
calc.survivor(model,gender="Male",Karno=100,Time=730)
calc.survivor(model,gender="Male",Karno=100,Time=1095)

## Note that we should really be providing confidence intervals 
## on these probabilities as they are effectively
## functions of the unknown parameters beta and gamma. 

## Compare nested models (see if gender is significant)
model1 <- survreg(Surv(time, status) ~ Karno, data=Lung)
model2 <- survreg(Surv(time, status) ~ Karno + gender, data=Lung)
loglik1 <- as.numeric(logLik(model1)) 
loglik2 <- as.numeric(logLik(model2))
LRT <- -2*(loglik1-loglik2)
LRT
## model2 is only larger by one parameter so we test using chi-sq with 1 d.o.f.
pchisq(LRT,1,lower.tail=F)
# p-value smaller than 0.05 so gender is significant

## Residual plots
x11()
par(mfrow=c(1,2)) ## set the plot window to two panels
fitted.values <- predict(model2,type="response")
std.deviance.redids <- residuals(model2,type="deviance")
## Resids vs fitted values
plot(fitted.values,std.deviance.redids,pch=20)
abline(h=0)
## and QQ plot
qqnorm( std.deviance.redids, pch=20 )
qqline( std.deviance.redids )

## Fitted values. Note these are defined as exp(X*beta)
## where X*beta is the mean of the log time-to-failure.
## Note however exp(X*beta) is NOT the mean of the 
## time-to-failure (but it is proportional to it).
## To get predictions of time-to-failure, we need to 
## get the mean of the Weibull which is (slide 17)
## mu = (1/theta)Gamma(1+1/gamma)
?predict.survreg
lin_pred <- predict(model2,type="lp")
theta <- exp(-lin_pred)
Mean <- (1/theta)*gamma(1+1/gamma)
## The mean for observarions that were censored, 
## are actually estimates of the true (unobserved) value
corrected <- Mean[Lung$status==0]
x11()
plot(corrected,Lung$time[Lung$status==0],ylab="Obs",xlab="corrected values from Weibull model",
     ylim=c(100,1100),xlim=c(100,1100),pch=20)
abline(0,1)
## Some of the corrected values seem too low. This is an aspect 
## of the model that we may think implies an unconvincing fit.


## We can also think of a model where gender specifies a different 
## baseline failure rate.
## We can do this by making the gamma parameter depend on 
## gender, using the function strata()
model3 <- survreg(Surv(time, status) ~ Karno + strata(gender), data=Lung)
summary(model3)
## Note that assessing the gender effect is not so straightforward 
## from the output. If we want to know the significance 
## of the effect, we need to use the LRT
model4 <- survreg(Surv(time, status) ~ Karno, data=Lung)
loglik3 <- as.numeric(logLik(model3)) 
loglik4 <- as.numeric(logLik(model4))
LRT <- -2*(loglik4-loglik3)
## model3 is only larger by one parameter so we test 
## using chi-sq with 1 d.o.f.
pchisq(LRT,1,lower.tail=F)
# p-value smaller than 0.05 so gender is significant 
# (but just about)

## Use similar code to above to calculate failure rates for men and women. 
### Note we now have one gamma for each
pred <- predict(model3,newdata=data.frame(Karno=100,gender="Male"),
                se.fit=T,type="lp",xlab="time (days)",ylab="failure rate")
theta <- exp(-pred$fit)
thetaU <- exp(-(pred$fit+1.96*pred$se.fit))
thetaL <- exp(-(pred$fit-1.96*pred$se.fit))
gamma1 <- 1/model3$scale[1]
gamma2 <- 1/model3$scale[2]
lambda.male <- theta*gamma1*(tt*theta)^(gamma1-1)
lambda.maleU <- thetaU*gamma1*(tt*thetaU)^(gamma1-1)
lambda.maleL <- thetaL*gamma1*(tt*thetaL)^(gamma1-1)
lambda.female <- theta*gamma2*(tt*theta)^(gamma2-1)
lambda.femaleU <- thetaU*gamma2*(tt*thetaU)^(gamma2-1)
lambda.femaleL <- thetaL*gamma2*(tt*thetaL)^(gamma2-1)
x11() ## plot from original model
plot.weib(model2,Karno=100,ylim=0.007) # baseline (Male with Karno=0)
plot.weib(model2,Karno=100,gender="Female",col="Red",add=T,ylim=0.007)
x11() ## new plot window
tt <- seq(5,1022,len=500)
plot(tt,lambda.male,lwd=3,type="l",ylim=c(0,0.007),xlab="time (days)")
lines(tt,lambda.maleU,lty=2,lwd=1)
lines(tt,lambda.maleL,lty=2,lwd=1)
lines(tt,lambda.female,lwd=3,col="red")
lines(tt,lambda.femaleU,col="red",lty=2,lwd=1)
lines(tt,lambda.femaleL,col="red",lty=2,lwd=1)
legend("topright",c("male","female"),col=c("black","red"),lty=c(1,1),lwd=3)
## Interesting how the curves look rather different. Does this model
## fit better? 
AIC(model3)
AIC(model2)

# Are the corrected values any better?
lin_pred <- predict(model3,type="lp")
theta <- exp(-lin_pred)
gammas <- rep(0.831, nrow(Lung)) ## males gamma
gammas[Lung$gender=="Female"] <- 0.630
Mean2 <- (1/theta)*gamma(1+1/gammas)
corrected2 <- Mean2[Lung$status==0]
x11()
par(mfrow=c(1,2))
plot(corrected,Lung$time[Lung$status==0],ylab="Obs",
     xlab="corrected values from Weibull model",
     ylim=c(100,1100),xlim=c(100,1100),pch=20)
abline(0,1)
plot(corrected2,Lung$time[Lung$status==0],ylab="Obs",
     xlab="corrected values from Weibull model",
     ylim=c(100,1100),xlim=c(100,1100),pch=20)
abline(0,1)
# Looks a little better actually. 

# Also illustrate log-Logistic
model5 <- survreg(Surv(time, status) ~ Karno + gender, data=Lung, 
                  dist="loglogistic")
summary(model5)
# residuals
x11()
par(mfrow=c(1,2)) ## set the plot window to two panels
fitted.values <- predict(model5,type="response") 
std.deviance.redids <- residuals(model5,type="deviance")
plot(fitted.values,std.deviance.redids,pch=20)
abline(h=0)
qqnorm( std.deviance.redids, pch=20 )
qqline( std.deviance.redids )
# Not particularly better than Weibull.
###########################################################################################################################

