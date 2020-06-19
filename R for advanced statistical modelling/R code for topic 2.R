load("ECM3712.RData")

###############################################################################################################
################### Look at very simple linear model data 'chol' when viewed as a glm rather than an lm #######
###############################################################################################################
head(chol)
names(chol)
plot(chol$age,chol$chol,pch=20)

# first get the linear model ouput
model1 <- lm(chol~age,data=chol)
summary(model1) # Residual standard error is the estimate of sigma from the model.

# now fit the same model as a glm
model2 <- glm(chol~age,data=chol,family=gaussian(link="identity"))
summary(model2)

# p is number of covariance, n-p-1 is degrees of freedom
n <- nrow(chol)
p <- 1
plot(seq(-4,4,len=200),dt(seq(-4,4,len=200),n-p-1),type = "l",lwd=3)
# quantile
qt(0.975,n-p-1)
qt(0.025,n-p-1) # The same as qt(0.975,n-p-1),but negative
abline(v=qt(0.975,n-p-1))
abline(v=qt(0.025,n-p-1))
# pt calcluates the aera of the curve to the left of the line
pt(10.136,n-p-1) #10.136 is the t-value of model2
1-pt(10.136,n-p-1) # want the area to the right of my line
2*(1-pt(10.136,n-p-1)) # both tails(is the value of pr of model2)

# Now match up some quantities in the output in order to reinforce 
# concepts in the glm theory when a two parameter member of the exponential 
# family is involved.

# First note that the 'beta' estimates and standard errors and t-values from output of model1
# and model2 are the same.

# Next note the the 'dispersion parameter' reported from model2
# is just the deviance/(n-p-1) from model 2 - this is an estimate of phi for this model
# and in the case of the gaussian this is just an estimate of sigma^2
# so it should be pretty much same as the square of the 'residual standard error'
# (i.e. the 'least squares estimate' of sigma^2) coming from model 1 
model2$deviance # Dm (ppt-19)
model2$deviance/22 # estimate of faer ,for nomal and poisson distribution,faer
                   # is just variance,22 is n-p-1 in this model

model2$deviance/model2$df.residual 
0.334^2  # sigama squared of linear model

# Note that what R calls residual deviance is what is returned by model$deviance. This is
# what we call D_M in the lecture slides. To get the scaled deviance D_M/phi in R, we would
# have to do it manually. First get the phi estimate from 
phi <- summary(model2)$dispersion
phi
# and then 
model2$deviance/phi
# to get the scaled deviance.

# No need to check that the model fits by testing whether the scaled deviance
# comes from a chi^2 with n-p-1 degrees of freedom since phi was estimated
# in order to ensure just that!


# What about the F-statisic? Well this is a likelihood ratio test of the fitted model
# and a model with just the intercept in it (the "NULL" model). So fit this NULL model first
model_null <- glm(chol~1,data=chol,family=gaussian(link="identity"))
# and then calculate the LRT as in slide 19 of topic 2:
n <- nrow(chol)
p1 <- 0 # NULL model has no predictors
p2 <- 1 # Only one predictor in model2
DM1 <- model_null$deviance
DM2 <- model2$deviance
Fstat <- ( (DM1-DM2)/(p2-p1) )/( DM2/(n-p2-1) )
Fstat
plot(seq(0,110,len=200),dt(seq(110,len=200),p2-p1,n-p2-1),type = "l",lwd=3) #right side,the left side I don't care because it is close to zero.
# 22 degrees of freedom, do not use dt(seq(0,110,len=200))
# because I don't care the left side.Get the F distr of 22 degrees of freedom.

# Is there a difference between the fit of the two models? 
1-pf(Fstat,p2-p1,n-p2-1)  # p for probability, f for distribution, extremly small
                          # just is the p-value
# Yes, model2 is a much better fit.
# Note we can use R to do these tests with the anova command:
anova(model_null,model2,test="F")
# p-value is exactly the same

# Next note that the AIC from model 1 is the same as that reported for model 2.
# There are 3 parameters in both models (beta_0, beta_1 and sigma^2) so to get
# the AIC manually we do
as.numeric(-2*logLik(model2)+2*3)
model2$aic
model_null$aic



############################################################
################# Play with challenger data - Binomial GLM #
############################################################
orings

# Naice estimate of failure probability
p <- mean(orings$Failures)/6
p

# Plot observed prob of failure against temperature
plot(orings$Temperature,orings$Failures/orings$Total,xlab="Temperature",
     ylab="Proportion Failures",ylim=c(0,1),xlim=c(30,82),pch=20)

# Try fitting a binomial GLM with temperature 
model1 <- glm(cbind(Failures,Total-Failures)~Temperature,data=orings,family=binomial(link="logit"))
summary(model1) #residual deviance is Dm
# Deviance goodness of fit of model 1 looks OK:
pchisq(model1$deviance,model1$df.residual,lower.tail=F)  #faer is one, model1$df.residual is n-p-1
# This is p-value of chi_sq distr with model1$df.residual degrees of freedom. The argument lower.tail=F means
# that we are calculating the area under the upper tail of the chi-squared. Note, equivalent to doing 
1-pchisq(model1$deviance,model1$df.residual) 

# plot fitted logistic line,pi(i) is the function of temperatures, use the
# below functions to predict
plot(orings$Temperature,orings$Failures/orings$Total,xlab="Temperature",ylab="Proportion Failures",
     ylim=c(0,1),xlim=c(30,82),pch=20)
# Create a sequence of 200 temperature values between 30 and 82:
temps <- seq(30,82,length=200)
# Use the predict command to predict the probability of failure for each of these temps. The 
# argument type="response" implies we are getting predictions of the mean and not the linear predictor.
preds <- predict(model1,newdata=data.frame(Temperature=temps),type="response")
lines(temps,preds,lwd=2,col="blue") #"S" shape

# Now to add approx. 95% confidence envelope around this line
# Predict again but at the linear predictor level along with standard errors
lin_preds <- predict(model1,newdata=data.frame(Temperature=temps),se.fit=T,type="link")
str(lin_preds)
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

###################################################
############## A little bit do not understand######
###################################################
## Alternatively we can use the standard error for the mean mu_hat (as in slides) and 
## assume that it is approx. Gaussian
prob_preds <- predict(model1,newdata=data.frame(Temperature=temps),se.fit=T,type="response")
str(prob_preds)
upper <- prob_preds$fit+1.96*prob_preds$se.fit
lower <- prob_preds$fit-1.96*prob_preds$se.fit
lines(temps,upper,lty=2,lwd=2,col="red")
lines(temps,lower,lty=2,lwd=2,col="red")
## The intervals are quite different from before and in fact exceed the value of 1 in some cases.
## This is probably due to the Gaussian approximation being poor as the number of trials (n_i=6) is low.

# add a line at the actual 1986 launch temperature (31 degrees)
abline(v=31,col="red",lwd=2)
# make a prediction and derive 95% confidence interval at 31 degrees using the same idea as above
lin_preds <- predict(model1,newdata=data.frame(Temperature=32),se.fit=T)
fit <- lin_preds$fit
upper <- fit+1.96*lin_preds$se.fit
lower <- fit-1.96*lin_preds$se.fit
fit <- exp(fit)/(1+exp(fit))
upper <- exp(upper)/(1+exp(upper))
lower <- exp(lower)/(1+exp(lower))
print(paste(round(fit,2)," (",round(lower,2),",",round(upper,2),")"))

## Check the fit also with the Pearson goodness of fit test since n_i is small here
chi2 <- sum(residuals(model1,type="pearson")^2)
chi2
model1$df.residual ## degrees of freedom (n-p-1)
pchisq(chi2,model1$df.residual,lower.tail=F)
## p-value smaller than the one for the Deviance test but still indicates good fit at 5% level

# Compare with smaller model
modeltheta <- glm(cbind(Failures,Total-Failures)~1,data = orings,
                  family = binomial(link = "logit"))
test<- modeltheta$deviance -model1$deviance
test
1-pchisq(test,1)
# model1 fits better
anova(modeltheta,model1,test="Chisq")
#Big model is always better

## Let's also check residual plots
par(mfrow=c(2,2),pch=20) ## Note pch=20 gives solid circles to make plot look nicer
plot(model1,1)
plot(model1,2)
plot(model1,3)
plot(model1,5)
## Top left plot and bottom right plots are essentially useless. Top right (QQ) plot is worrying as
## the points deviate quite a lot from the line. Bottom right one indicates observation 1 is an 
## influential outlier with high leverage

# just for interest repeat all above excluding observation 1
orings1 <- orings[-1,]
model2 <- glm(cbind(Failures,Total-Failures)~Temperature,data=orings1,family=binomial(link="logit"))
summary(model2)
# This has made a big difference in the results, with Temperature now being insignificant at the 5% level.
# check Deviance goodness of fit 
pchisq(model2$deviance,model2$df.residual,lower.tail=F)
# this is OK at the 5% level

# residual plots are not really any better 
par(mfrow=c(2,2),pch=20)
plot(model2,1)
plot(model2,2)
plot(model2,3)
plot(model2,5)
# But no real justification for leaving observation 1 so need to think of other ways of improving the model... 
# Are we justified in assuming all the 23 shuttles were different only due to temperature during launch?






##################################################
###### Titanic data - question 6 on sheet 1 ######
##################################################
head(titanic)  #real data 
head(titanic[,-3])
head(titanic2)

## Start by fitting model with all two-way interactions and main effects
model1 <- glm(survived~(pclass+gender+age)^2,data=titanic,family=binomial(link="logit"))
summary(model1)
## The main effect for age is not significant but some interactions with age are. To see if we
## can reduce this model, we first try to remove each interaction term and compare the AIC each time
## to decide if the term is to be dropped. Function drop1() will do that for us.
drop1(model1)
# Up to you decide if it should be reduced or not


## Now look at the same data, after grouping them into a binomial structure
head(titanic2)

## Fit binomial model equivalent to model1
model2 <- glm(cbind(survived,total-survived)~(age_group+pclass+gender)^2,data=titanic2,family=binomial(link="logit"))
summary(model2)
# Look at the value of the deviance D_M
deviance(model2) 
# Effectively zero. What does this imply about the model? Go back to the general expression of the Deviance D_M 
# for a GLM on top of slide 18 in Topic 2 and think about what values of theta_hat would be needed to make this zero.

## Can we reduce this model based on the AIC?
drop1(model2)
# Well, no according to the table, however we should do so anyway so reduce to a model that will procuce the least worst fit
# in terms of the AIC. Up to you to do the rest...



#Lecture 7

##############################################################
###### Poisson data with offset - question 8 on sheet 1 ######
##############################################################

# Need a model which has cromosonal abnormalities as the response (Poisson)
# but which allows for numbers of cells exposed for each observation (so 
# this is treated as an offset in the model i.e. we are effectively modelling the
# abnormality rate per hundred cells.
help(dicentric,package="faraway")
head(dicentric)

# This model should do the initial job (note the I() function here protects the model 
# formula from any confusion about mathematical transformations within it).
pois.model <- glm(ca~offset(I(log(cells)))+doseamt+doserate, data=dicentric, family=poisson(link="log"))
summary(pois.model)

# Note that this models the rate of abnormalities per 100 cells. To model the rate per cell:
dicentric$log.cells <- log(dicentric$cells*100)
pois.model2 <- glm(ca~offset(log.cells)+doseamt+doserate,family=poisson,data=dicentric)
summary(pois.model2)
# Notice only the intercept changes. The coefficient effects are the same as they are 
# act multiplicatively on the "base rate" which is exp(beta_0).

# This should be similar to a binomial model of the probability of abnormality given the 
# abnormality events are relatively rare (e.g. rate of 1 per 100, i.e. 0.01, will be similar
# to the odds 0.01/0.99. This approximation breaks down for non-rare events, i.e. high rates).
response.matrix <- cbind(dicentric$ca,dicentric$cells*100-dicentric$ca)
bin.model <- glm(response.matrix~doseamt+doserate, data=dicentric, family=binomial(link="logit"))
summary(bin.model)

# However, the residual deviance (of the Poisson) looks very high relative to the degrees of freedom
1 - pchisq(pois.model$deviance,pois.model$df.residual)

# Sometimes, we can give the Poisson much needed flexibility by trying interactions between the predictors.
# We want to make sure, however, that the interaction term makes sense w.r.t. the application. Here, the
# interaction term is sensible as it allows for the effect of the dose rate and amount to not just 
# be multiplicative
pois.model3 <- glm(ca~offset(log.cells)+doseamt*doserate, data=dicentric, family=poisson(link="log"))
summary(pois.model3)
# well the interaction is significant but it doesn't really solve the lack of fit problem
1-pchisq(pois.model3$deviance,pois.model3$df.residual)
# What about residuals? Not very convincing, residuals vs fitted values may indicate non-linear predictor effects.
par(mfrow=c(2,2))
plot(pois.model3,1)
plot(pois.model3,2)
plot(pois.model3,5)
# So try a model with squared terms in it.
pois.model4 <- glm(ca~offset(log.cells)+doseamt*doserate+I(doserate^2)+I(doseamt^2), data=dicentric, family=poisson(link="log"))
summary(pois.model4)
# All significant, but does it fit?
1 - pchisq(pois.model4$deviance,pois.model$df.residual)
# Well, not at the 5% level but better than before. What about residuals?
par(mfrow=c(2,2))
plot(pois.model4,1)
plot(pois.model4,2)
plot(pois.model4,5)
# QQ plot not great, hard to tell with residuals vs fitted values but some signs of funelling.
# Let's try a model with cubic terms in it
pois.model5 <- glm(ca~offset(log.cells)+doseamt*doserate+I(doserate^2)+I(doseamt^2)+I(doserate^3)+I(doseamt^3), data=dicentric, family=poisson(link="log"))
summary(pois.model5)
# Well, not enough unique data points for doesamt (might have even been sensible to treat this as a factor), 
# but cubic term for doserate is significant.  Does it fit?
1 - pchisq(pois.model5$deviance,pois.model$df.residual)
# Yes it does. What about residuals?
par(mfrow=c(2,2))
plot(pois.model5,1)
plot(pois.model5,2)
plot(pois.model5,5)
# QQ plot looks OK. Again hard to tell with residuals vs fitted. Observation 27 may be an outlier but we don't 
# really have any justification for excluding it. So model with cubic term for doserate and quadratic term for 
# doseamt fits the data and also points to non-linear relationships.

# Let's try to undestand the inference by predicting from the model. Fit model5 again without the cubic doseamt
pois.model6 <- glm(ca~offset(log.cells)+doseamt*doserate+I(doserate^2)+I(doseamt^2)+I(doserate^3), data=dicentric, family=poisson(link="log"))
Drate <- seq(0.1,4,len=200)
Damt <- rep(1,200)
Offset <- rep(0,200) # so we are predicting the rate abnormal cells per (exp(0)=1) cell
preds <- predict(pois.model6,newdata=data.frame(doserate=Drate,doseamt=Damt,log.cells=Offset),type="link",se.fit=T)
plot(Drate,exp(preds$fit),type="l",lwd=2,ylim=c(0,0.017))
lines(Drate,exp(preds$fit+preds$se.fit),lwd=2,lty=2)
lines(Drate,exp(preds$fit-preds$se.fit),lwd=2,lty=2)
Damt <- rep(2.5,200)
preds <- predict(pois.model6,newdata=data.frame(doserate=Drate,doseamt=Damt,log.cells=Offset),type="link",se.fit=T)
lines(Drate,exp(preds$fit),lwd=2,col="red")
lines(Drate,exp(preds$fit+preds$se.fit),lwd=2,lty=2,col="red")
lines(Drate,exp(preds$fit-preds$se.fit),lwd=2,lty=2,col="red")
Damt <- rep(5,200)
preds <- predict(pois.model6,newdata=data.frame(doserate=Drate,doseamt=Damt,log.cells=Offset),type="link",se.fit=T)
lines(Drate,exp(preds$fit),lwd=2,col="blue",type="l")
lines(Drate,exp(preds$fit+preds$se.fit),lwd=2,lty=2,col="blue")
lines(Drate,exp(preds$fit-preds$se.fit),lwd=2,lty=2,col="blue")
# So quite complicted non-linear relationship between the rate of abnormalities per cell, and dose rate and amount. As 
# expected, both the abn. rate significantly increases with dose rate and amount (all p-values are <0.05).

# Note we are using the covariates to expand the model. In some cases we may not be able to do this especially if it does
# not make sense w.r.t. the application. In such cases we may need to resort to a quasi-Poisson to "patch up" the standard 
# errors
qpois.model <- glm(ca~offset(log.cells)+doseamt*doserate, data=dicentric, family=quasipoisson(link="log"))
summary(qpois.model)
# now appreciate that the interaction is not really doing much so try to leave it out
drop1(qpois.model,test="Chi")
# Indeed then we can leave the interaction term out
qpois.model2 <- glm(ca~offset(log.cells)+doseamt+doserate, data=dicentric, family=quasipoisson(link="log"))
summary(qpois.model2)
# Looks like that's about the best we can do. Let's look at the residual plots and then interpret the model
par(mfrow=c(2,2))
plot(qpois.model2,1)
plot(qpois.model2,2)
plot(qpois.model2,5)
## Residuals vs fitted not very convincing, may indicate non-linear predictor effects, however rest look acceptable.
# Alternatively, we can try fitting a Negative Binomial GLM. (Better than quasi-Poisson as we know what the distr. 
# of the response is.)
library(MASS) # contains the function glm.nb
nb.model1 <- glm.nb(ca~offset(log.cells)+doseamt*doserate, data=dicentric)
summary(nb.model1)
# As with quasi-Poisson, interaction not significant. Can we reduce?
drop1(qpois.model,test="F")
# Yes, so let's do this:
nb.model2 <- glm.nb(ca~offset(log.cells)+doseamt+doserate, data=dicentric)
summary(nb.model2)
# Let's look at the residual plots and then interpret the model
par(mfrow=c(2,2),pch=20)
plot(nb.model2,1)
plot(nb.model2,2)
plot(nb.model2,5)
## Residuals vs fitted not very convincing, may indicate non-linear predictor effects, however rest look acceptable.
## Summarising, we have three models that fit well, however we would
## go for the Poisson one with the interactions as it has the residual
## plots looking right and also because the interactions make sense
## in terms of the problem at hand.


#####################################################################
###### Poisson log-linear model. Data from R package "faraway" ######
#####################################################################
# In 1972-74, a survey of one in six residents of Whickham, near Newcastle, England was made. 
# Twenty years later, this data recorded in a follow-up study. Only women who are current smokers 
# or who have never smoked are included. Resulting data set comprises 28 obs on the following
# 4 variables: y=observed count for given combination, smoker=a factor with levels yes no,
# dead=a factor with levels yes no, age=a factor with age-group levels 18-24 25-34 35-44 45-54 55-64 65-74 75+
# Interest here lies on the effects of age and smoking on the probability of death. Obviously we could
# fit a binomial model since death is binary, however we fit this here as a Poisson log-linear model.
# Assuming all margins are random, implies no particular terms need to be included by design. Start
# with biggest model possible and reduce if necessary
library(faraway)
head(femsmoke)
# Fit biggest model possible, i.e. one with the 3-way interaction, all 2-way intercations and all main effects:
model1 <- glm(y~smoker*dead*age,data=femsmoke,family=poisson)
# Note that we can also be explicit about what interactions we include. E.g., we can also fit this model as:
# model1 <- glm(y~smoker+dead+age+smoker:age+smoker:dead+age:smoker+smoker:age:dead,data=femsmoke,family=poisson)
summary(model1)
deviance (model1)
# This is obviously the saturated model so it's of no use. Let's reduce to second biggest model, i.e. one 
# with all 2-way intercations and all main effects:
model2 <- glm(y~(smoker+dead+age)^2,data=femsmoke,family=poisson)
summary(model2)
# Model fits?
1-pchisq(deviance(model2),model2$df.residual)
# Yes it does. Could we reduce further? Can use drop1() function to try dropping any one of the 2-way interactions
# performing an LRT each time
drop1(model2,test="Chi")
# So should not really reduce further. Now to interpret the model. Well, the model output shows interactions and main effect
# for dead=no, so we are going to perform inference on the probability of not dying. We always start from the highest order effects
# (in this case the two-way effects) of "dead". Start with "smokerno:deadno". The effect is positive and significant, so not smoking 
# has a positive effect on the prob. of not dying (which makes sense). Note that the model is additive, so that we can make this
# interpretation without caring what the other variables are doing (e.g. deadno:age or smokeno:age). Now, for deadno:age. No difference
# in the prob. of not dying between age group 25-34 and 18-24 (the baseline). For the age groups 35-44, 45-54, 55-64 and 65-74, the 
# probability of not dying gets progressively smaller compared to the baseline 18-24 (which also makes sense). For the age group 75+, 
# there is probably a lack of data (notice the huge standard error on deadno:age75+). The interaction smoker:age doesn't tell us 
# anything about dying, however it gives information on the association between the propensity to smoke and age at the time of the
# study. E.g. there are significantly more smokers in age groups 65-74 and 75+ compared to 18-24.




###########################################################################
###### Handling censored survival data (question 10 Problems 1) ###########
###########################################################################

head(gehan)

# Extract the uncensored data and fit an exponential model in R this is done by fitting a Gamma model and then
# setting the dispersion parameter to 1 note the default link for the Gamma 
# is inverse so we are fitting a model where time_i~exp(lambda_i) and 
# lambda_i=1/mu_i=b_0+b1*treat_i (the control effect is the intercept b_0 
# the treatment effect is b_0+b1)

# Create new data frame with only uncensored data in it
gehan$treat <- relevel(gehan$treat, "control")
gehan2 <- gehan[gehan$cens==1,]
exp.model <- glm(time~treat,family=Gamma(link="inverse"),data=gehan2)
summary(exp.model,dispersion=1)

# The b1 is not significantly different from zero so
# suggestion is no diffence between treatment and control
# but only if model fits.

# Now fit a model to the full censored data set using the Poisson 'trick'
# i.e. fit the model cens_i~Pois(lambda_i * time_i) using an appropriate model
# formulat in R such that lambda_i = b_0 + b1*treat_i. 
# Remember it's a Poisson so might want to check the deviance before 
# interpreting.

# Hint: use 
gehan$treat <- relevel(gehan$treat, "6-MP")
# to ensure that the Poisson model 'buries' the control in the intercept 
# rather than the treatment.


