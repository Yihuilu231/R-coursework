install.packages("mgcv")
install.packages("survival")
library(survival)
library(mgcv)
library(lme4)
rm(list=ls())
load("ECM3712.RData")
times <- c(143,164,188,188,190,192,206,209,213,216,220,227,230,234,
           246,265,304,216,244)
cens <- rep(1,19)
cens[18:19] <- 0

log_likelihood <- function(theta,times){
  ell <- sum(log(theta*exp(-theta*times[1:17]))) + 
    sum(log(exp(-theta*times[18:19])))
  return(-ell)
}

out <- nlm(log_likelihood,1,times=times,hessian = T)
theta <- out$estimate
1/theta # the mean
theta
sqrt(1/(out$hessian)) # std error

library(survival)
Times <- Surv(times,cens,type="right")
model <- survreg(Times~1,dist="exponential")
summary(model)
theta <- exp(-model$coefficients)
theta
exp(-theta*300)

head(Lung)
model1 <- survreg(Surv(time,status)~Karno+gender+age,data = Lung,dist="weibull")
summary(model1)

model2 <- survreg(Surv(time,status)~Karno+gender,data = Lung,dist="weibull")
summary(model2)

ll_model1 <- logLik(model1)
ll_model2 <- logLik(model2)
ll_model1 <- as.numeric(logLik(model1))
ll_model2 <- as.numeric(logLik(model2))
teststat <- -2*(ll_model2-ll_model1)
teststat
pchisq(teststat,df=1,lower.tail=FALSE)

model3 <- survreg(Surv(time,status)~Karno+gender,data = Lung,dist="exponential")
summary(model3)
gamma <- 1/model3$scale
gamma

#Check for “Normality” of the standardised deviance residuals.
qqnorm(residuals(model2,type = "deviance"))
qqline(residuals(model2,type = "deviance"))

# Question 3d
pred.model<- predict(model2,newdata = data.frame(gender="Male",Karno=100),se.fit=T,type = "lp")
## Use it to calculate theta = exp(-X*beta)
theta <- exp(-pred.model$fit)
## and it's confidence interval
thetaU <- exp(-(pred.model$fit+1.96*pred.model$se.fit))
thetaL <- exp(-(pred.model$fit-1.96*pred.model$se.fit))
gamma1 <- 1/model3$scale[1]
gamma2 <- 1/model3$scale[2]
## Create a sequence ot time values in the range of the data
tt <- seq(5,1022,len=500)
## Use theta and its CI along with the estimate 
## of gamma to compute the failure rate
lambda1 <- theta*gamma*(tt*theta)^(gamma-1)
lambda1
lambdaU <- thetaU*gamma*(tt*thetaU)^(gamma-1)
lambdaL <- thetaL*gamma*(tt*thetaL)^(gamma-1)

## Plot the failure rate
x11()
par(lwd=2)
plot(tt,lambda1,lwd=3,type="l",ylim=c(0,0.01),
     xlab="time (days)", ylab="failure rate") 
lines(tt,lambdaU,lty=2)
lines(tt,lambdaL,lty=2)

#Female
pred.model<- predict(model2,newdata = data.frame(gender="Female",Karno=100),se.fit=T,type = "lp")
## Use it to calculate theta = exp(-X*beta)
theta <- exp(-pred.model$fit)
## and it's confidence interval
thetaU <- exp(-(pred.model$fit+1.96*pred.model$se.fit))
thetaL <- exp(-(pred.model$fit-1.96*pred.model$se.fit))
gamma1 <- 1/model3$scale[1]
gamma2 <- 1/model3$scale[2]
## Create a sequence ot time values in the range of the data
tt <- seq(5,1022,len=500)
## Use theta and its CI along with the estimate 
## of gamma to compute the failure rate
lambda2 <- theta*gamma*(tt*theta)^(gamma-1)
lambda2
lambdaU <- thetaU*gamma*(tt*thetaU)^(gamma-1)
lambdaL <- thetaL*gamma*(tt*thetaL)^(gamma-1)

## Plot the failure rate
x11()
par(lwd=2)
plot(tt,lambda2,lwd=3,type="l",ylim=c(0,0.01),
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
plot.weib(model2,Karno=0) # baseline (Male with Karno=0)
plot.weib(model2,Karno=100,gender="Female",col="Red",add=T)
plot.weib(model2,Karno=100,gender="Male",col="blue",add=T)

#Question3 e
model4 <- survreg(Surv(time,status)~Karno+strata(gender),data = Lung,dist="weibull")
summary(model4)
## Note that assessing the gender effect is not so straightforward 
## from the output. If we want to know the significance 
## of the effect, we need to use the LRT
model5 <- survreg(Surv(time, status) ~ Karno, data=Lung)
loglik4 <- as.numeric(logLik(model4)) 
loglik5 <- as.numeric(logLik(model5))
LRT <- -2*(loglik5-loglik4)
## model3 is only larger by one parameter so we test 
## using chi-sq with 1 d.o.f.
pchisq(LRT,1,lower.tail=F)
# p-value smaller than 0.05 so gender is significant 
# (but just about)

x11()
plot.weib(model4,Karno=0) # baseline (Male with Karno=0)
plot.weib(model4,Karno=100,gender="Female",col="Red",add=T)
plot.weib(model4,Karno=100,gender="Male",col="blue",add=T)

AIC(model4)
AIC(model2)

#Question 4
head(nlmodel)
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
xx <- seq(0,1,len=100) ## a sequence of x values on which to produce fitted values
yfit <- (theta1*xx)/(theta2 + xx) ## compute fitted values
lines(xx, yfit,lwd=2,col="red") ## add the estimated line

# Now fitting the model
model <- gam(y~s(x,k=9,bs="cr"),data = nlmodel,family = gaussian(link="identity"))
gam.check(model)

## Let's look at what's inside the fitted model
# Estimate of sig^2
model$sig2
# The rank of f(x)
model$rank
# The beta estimates
model$coefficients
# Estimate of the smoothing parameter
model$sp

## Looks like we could do with a bit more flexibility
## What about the rank of the model being large enough? 
## We can use gam.check() to check this
## and also look at residual plots
x11()
par(mfrow=c(2,2))
gam.check(model,pch=20)
# k' is close to edf so might need to increase the rank (equiv, k-index is less than 1)

model$edf
model$df.residual
# which is the same as
100 - sum(model$edf) ## (where n=100)


## Produce fitted line from this model using perdict
yfitAM <- predict(model,newdata=data.frame(x=xx),se.fit=T)
x11()
par(mar = c(4, 4, 1, 1),cex=1.2,lwd=2)
plot(y~x,data=nlmodel,pch=20)
lines(xx,yfitAM$fit,col="blue",lwd=2)
lines(xx,yfitAM$fit+1.96*yfitAM$se.fit,col="blue",lwd=2,lty=2)
lines(xx,yfitAM$fit-1.96*yfitAM$se.fit,col="blue",lwd=2,lty=2)

# Question 5
head(aids)
model <- gam(cases~s(date,k=10,bs="cs"),data=aids,family=poisson(link = "log"))
x11()
par(mfrow=c(2,2))
gam.check(model,pch=20)

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

model$deviance
model$df.residual
pchisq(model$deviance,model$df.residual,lower.tail = F)

# Improve the model
head(aids)
Ipmodel1 <- gam(cases~s(date,k=10,bs="cs")+quarter,data=aids,family=gaussian(link = "log"))
AIC(model)
AIC(Ipmodel1)
# Same model checking as the native one
x11()
par(mfrow=c(2,2))
gam.check(Ipmodel1,pch=20)
# Now predict the LINEAR PREDICTOR and get its standard errors. The linear predictor
# is of course Gaussian so we will calculate CIs and then exponentiate those
# Quarter 1
preds1 <- predict(Ipmodel1,newdata=data.frame(date=pred_x,quarter=as.factor(1)),
                 se.fit=T,type="link")
x11(width=8,height=6)
par(mfrow=c(1,1),mar = c(4, 4, 1, 1),cex=1.2,lwd=2)
plot(cases ~ date,data=aids,pch=20)
lines(pred_x,exp(preds1$fit),col="blue",lwd=2)
lines(pred_x,exp(preds1$fit+1.96*preds1$se.fit),col="blue",lwd=2,lty=2)
lines(pred_x,exp(preds1$fit-1.96*preds1$se.fit),col="blue",lwd=2,lty=2)
# Quarter 2
preds2 <- predict(Ipmodel1,newdata=data.frame(date=pred_x,quarter=as.factor(2)),
                  se.fit=T,type="link")
x11(width=8,height=6)
par(mfrow=c(1,1),mar = c(4, 4, 1, 1),cex=1.2,lwd=2)
plot(cases ~ date,data=aids,pch=20)
lines(pred_x,exp(preds2$fit),col="blue",lwd=2)
lines(pred_x,exp(preds2$fit+1.96*preds2$se.fit),col="blue",lwd=2,lty=2)
lines(pred_x,exp(preds2$fit-1.96*preds2$se.fit),col="blue",lwd=2,lty=2)
# Quarter 3
preds3 <- predict(Ipmodel1,newdata=data.frame(date=pred_x,quarter=as.factor(3)),
                  se.fit=T,type="link")
x11(width=8,height=6)
par(mfrow=c(1,1),mar = c(4, 4, 1, 1),cex=1.2,lwd=2)
plot(cases ~ date,data=aids,pch=20)
lines(pred_x,exp(preds3$fit),col="blue",lwd=2)
lines(pred_x,exp(preds3$fit+1.96*preds3$se.fit),col="blue",lwd=2,lty=2)
lines(pred_x,exp(preds3$fit-1.96*preds3$se.fit),col="blue",lwd=2,lty=2)
# Quarter 4
preds4 <- predict(Ipmodel1,newdata=data.frame(date=pred_x,quarter=as.factor(4)),
                  se.fit=T,type="link")
x11(width=8,height=6)
par(mfrow=c(1,1),mar = c(4, 4, 1, 1),cex=1.2,lwd=2)
plot(cases ~ date,data=aids,pch=20)
lines(pred_x,exp(preds4$fit),col="blue",lwd=2)
lines(pred_x,exp(preds4$fit+1.96*preds4$se.fit),col="blue",lwd=2,lty=2)
lines(pred_x,exp(preds4$fit-1.96*preds4$se.fit),col="blue",lwd=2,lty=2)

Ipmodel2 <- gam(cases~s(date,k=10,bs="cs")+quarter,data=aids,family=poisson(link = "log"))
AIC(Ipmodel2)
x11()
par(mfrow=c(2,2))
gam.check(Ipmodel2,pch=20)

# Question 8
head(penicillin)
model1 <- lmer(yield~treat+(1|blend),data=penicillin,REML=F)
model2 <- lmer(yield~(1|blend),data=penicillin,REML=F)
#Model 3 is a linear model
model3 <- lm(yield~treat,data=penicillin)
LRT <- -2*(logLik(model2)-logLik(model1))
LRT
1-pchisq(LRT,1)
anova(model1,model2)

LRT <- -2*(logLik(model3)-logLik(model1))
LRT
pchisq(LRT,1,lower.tail = F)

# 8(b)--Alternative way
sim_LRT <- 1:1000 # vector to store values of the LRT for each iteration
Dat <- simulate(model3,1000) ### Simulate 1000 data sets from model1 (the 
## smaller model, i.e. the NULL hypothesis)
for(i in 1:1000){ ### start loop
  Mod1 <- lmer(Dat[,i]~(1|blend),data=penicillin,REML=F) ### Fit model1 to simulated data
  Mod2 <- lmer(Dat[,i]~treat+(1|blend),data=penicillin,REML=F) ### Fit model2 simulated data
  sim_LRT[i] <- -2*(logLik(Mod1)-logLik(Mod2)) ### Calculate and store LRT
} ## end loop
mean(sim_LRT>LRT) ### Calculate p-value (i.e. prop of sim_LRT values greater than LRT)
plot(density(sim_LRT),xlim=c(0,40),main="",xlab="LRT")
abline(v=LRT,col="red")


#Question 9
head(pupils)
model1 <- glm(test~IQ+ses+Class,data=pupils,family=gaussian(link="identity"))
summary(model1)
# both IQ and sess are highly significant (both have +ve effect on test)

# Do a LRT to see is class is significant
model0 <- glm(test~IQ+ses,data=pupils,family=gaussian(link="identity"))
anova(model0,model1,test="F")

# F-test of difference between a model with class in and class out gives
# value of 4.783 to be tested against an F_(130,2154) distribution
# gives a very small p-value (<.00001) so class is highly significant

# Plot the class fixed effects with error bars
# just to see roughly how classes (i.e. schools) compare
class <- predict(model1,newdata=data.frame(Class=as.character(1:131),IQ=10,ses=23),
                 terms="Class",type="terms",se.fit=T)
class_eff <- as.vector(class$fit)
class_se <- as.vector(class$se.fit)
X <- cbind(class_eff-1.96*class_se,class_eff,class_eff+1.96*class_se)
plot(1:131,X[,2],ylim=c(min(X[,1]),max(X[,3])),pch=20,xlab="class",ylab="fixed class effects",main="" )
for(i in 1:131){ lines(rep(i,3),X[i,],lwd=2) }
abline(h=0)

#9b (iii)
model2 <- lmer(test~IQ+ses+(1|Class),data=pupils)
summary(model2)
# The conditional varince of the response 40.049 is roughly the same 
# as the residual variance of model1:
summary(model1)$dispersion
# Variance of the random efects is 9.212 and we can test if this 
# is significant using the LRT:
model2 <- lmer(test~IQ+ses+(1|Class),data=pupils,REML=F) # refit the model using max. likelihood
model0 <- glm(test~IQ+ses,data=pupils,family=gaussian(link="identity"))
LRT <- -2*(logLik(model0)-logLik(model2))
LRT <- as.vector(LRT)
LRT
pchisq(LRT,1,lower.tail=F)

# Extract the random effects from model 2 and draw a density plot of them
# superimpose on it the theoretical density of a Gaussian distribution 
# with mean zero that
# depends the class random effect standard deviation estimated from model2
class_ran <- ranef(model2)$Class[,1]
x11()
plot(density(class_ran),lwd=2,xlab="Class random effect predictions",
     ylab="density")
lines(density(rnorm(100000,0,3.035)),col="red",lwd=2) ## add theoretical Gaussian distr. with sd 3.035
## looks good enough. Better though to look at a Q-Q plot with qqnrom() 
## and qqline().

# Produce a 'caterpillar' plot of predicted class 
# random effects with error bars
# just to see roughly how classes (i.e. schools) 
# compare. Luckily dotplot() knows how to do this!
library(lattice)
dotplot(ranef(model2,condVar=TRUE))

# Residual plots of the model:
x11(width=16,height=6)
par(mfrow=c(1,2),mar = c(4, 4, 1, 1),cex=1.2,lwd=2)
fits <- fitted(model2)
resids <- resid(model2)
plot(fits,resids,pch=20)
qqnorm(resids,pch=20)
qqline(resids)
# Not so so good, especially the one on the left...

#9(c)
## Now we want to fit a model where the IQ effect is different per 
## class. In other words we want
## an interaction between pupil IQ and class. 
model3 <- lmer(test~IQ+ses+(1+IQ|Class),data=pupils)
summary(model3)
# The fixed effect of IQ is 2.29407 which is very similar to the 
# fixed effect of 2.25325 from model2.
# Significance is roughly assessed by conditioning of the 
# random effect parameters so we look to t-tests
# with n-p-1 where p=2.

# Are the IQ random slopes significant? Use the LRT:
model2 <- lmer(test~IQ+ses+(1|Class),data=pupils,REML=F) # model2 using max likelihood
model3 <- lmer(test~IQ+ses+(1+IQ|Class),data=pupils,REML=F) # model3 using max likelihood
LRT <- -2*(logLik(model2)-logLik(model3))
LRT <- as.vector(LRT)
LRT
pchisq(LRT,1,lower.tail=F)
# So there is a significant class-IQ interaction (although 
# same comments apply as above for using the LRT for 
# testing random effects)

# Note if this was a fixed effects model, we would have 
# 263 parameters as opoosed to 5.
fitted(model2)

