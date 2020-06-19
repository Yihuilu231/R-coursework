install.packages("gridExtra")
library(mgcv)
library(gridExtra)
rm(list = ls())

load("ECM3712.RData")

#Question 5
aids
head(aids)
plot(aids$cases,aids$date,pch=20,xlab="date",ylab = "cases")

model1 <- glm(cases~date,data=aids,family=poisson(link="log"))
summary(model1)

model2 <- glm(cases~date,data=aids,family=gaussian(link="log"))
summary(model2)
model2$deviance
#add trends of each model
aids
plot(aids$date,aids$cases,xlab="date",ylab="cases",
     ylim=c(0,500),xlim=c(80,100),pch=20)
# Create a sequence of 200 temperature values between 30 and 82:
temps <- seq(80,100,length=200)
# Use the predict command to predict the probability of failure for each of these temps. The 
# argument type="response" implies we are getting predictions of the mean and not the linear predictor.
preds1 <- predict(model1,newdata=data.frame(date=temps),type="response")
lines(temps,preds1,lwd=2,col="blue") 

preds2 <- predict(model2,newdata=data.frame(date=temps),type="response")
lines(temps,preds2,lwd=2,col="blue") 
lines(temps,preds1$fit + 1.96*preds1$se.fit,lty=2,col="red",lwd=3)
lines(temps,preds1$fit - 1.96*preds1$se.fit,lty=2,col="red",lwd=3)

# Now to add approx. 95% confidence envelope around this line
# Predict again but at the linear predictor level along with standard errors
lin_preds <- predict(model1,newdata=data.frame(date=temps),se.fit=T,type="link")
str(lin_preds)

with(preds1, lines(0:1000, exp(fit+1.96*se.fit)/(1+exp(fit+1.96*se.fit)), lty=2))
with(preds2, lines(0:1000, exp(fit-1.96*se.fit)/(1+exp(fit-1.96*se.fit)), lty=2))
# Calculate upper confidence interval limit since linear predictor is approx. Gaussian
upper <- lin_preds$fit+1.96*lin_preds$se.fit
# Transform the upper confidence interval limit to get one at the level of the mean
upper <- exp(upper)/(1+exp(upper))
# Add it to the plot
lines(temps,upper,lty=2,lwd=2,col="blue")
# Do the same with the lower interval limit
lower <- lin_preds$fit-1.96*lin_preds$se.fit
lower <- exp(lower)/(1+exp(lower))
lines(temps,lower,lty=2,lwd=2,col="blue")

##########################AIC###################
AIC(model1)   #1153.873
AIC(model2)   #446.3331
######the second model's aic is much lower than the first one,so it is preferred

#######deviance residuals#######################
deviance(model1)
deviance(model2)
model1$fitted.values
plot(model1,which=1,pch=20)
plot(model1,which=2,pch=20)
plot(model1,which=5,pch=20)

plot(model2,which=1,pch=20)
plot(model2,which=2,pch=20)
plot(model2,which=5,pch=20)

####################proposed extensions##############
model3 <- glm(cases~date+quarter,data=aids,family=gaussian(link="log"))
summary(model3)
drop1(model1)
drop1(model2)
drop1(model3)
