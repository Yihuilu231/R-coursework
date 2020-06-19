## Illustrate effect of making a parameter of a distribution random
## Say, Y|psi ~ N(mu,sigma^2_y) where mu=beta_0+psi and psi ~ N(0,sigma^2_psi)
sig_psi <- 0.4
beta0 <- 1
sig_y <- 1
psi <- rnorm(10000,0,sig_psi)
mu <- beta0 + psi
y <- rnorm(10000,mu,sig_y)
var(y) # we know this is sigma^2_psi+sigma^2_y from the notes
sig_psi^2 + sig_y^2
# But what is the distribution of y?
x11()
plot(density(y),lwd=3)
# In general we would have to integrate the psi out but fortunately for a
# Gaussian distribution with a Gaussian mean, we get another Gaussian.
lines(seq(-5,5,len=200),dnorm(seq(-5,5,len=200),beta0,sqrt(sig_psi^2 + sig_y^2)),
      col="red",lwd=3)

################################################################
###### Play with the tyres data (Question 7 on sheet 2) ########
################################################################

tyres[1:20,]
levels(tyres$wheel)
levels(tyres$car_type)

## Fit a linear model with wheel and car_type as factors
model1 <- lm(wear~car_type+wheel,data=tyres) 
summary(model1)

## Now fit an equivalent random effects model where car_type 
## is random (which makes sense as
## the 3 cars we have a a random sample of all cars). 
## Note that by default, lmer will fit using
## restricted maximum likelihood (REML).
library(lme4)
model2 <- lmer(wear~wheel+(1|car_type),data=tyres)
summary(model2)
# Note that this model has only 6 parameters: beta0, beta1, 
# beta2, beta3, sigma^2_y and sigma^2_gamma. The estimate #
# of sigma^2_y (the CONDITIONAL variance of the y) 
# is what R calls "Residual". The estimate of sigma^2_gamma 
# (the variance of the car random effects) is what R 
# calls "car_type".

# Note also that the fixed effects estimates (of wheel) are not 
# the same since the models have different intercepts. If we want 
# to compare the coefficient estimates for the two models we need 
# to be careful because the intercept in model1 is front nearside 
# in car type A whereas the intercept in model2 is average front 
# nearside over all three car types. To get a proper comparison, 
# we can use the predict() function for model1 to produce the 
# estimates of the wheel effects for the average car:
pred_wheels <- levels(tyres$wheel) # a character vector of all wheel types
pred_wheels
# now predict each wheel with type="terms" argument to
predict(model1,newdata=data.frame(wheel=pred_wheels,car_type="A"),
        type="terms",terms="wheel")
## The attr(,"constant") of 24.85861 is the overall mean wear of all cars, 
## so to get the intercept of model2 (average front nearside over all 
## three car types), we need the overall mean plus the 
## wheel effects of front nearside:
24.85861 - 7.339722 ## front nearside
## and now results are comparable to the fixed effects of the LMM
fixef(model2)


### Now do some inference

## Test for the significance of the fixed effects (overall significance of 
## wheel position) in the LMM. Is wheel position significant overall? 
## Fit a model with and withouy the wheel fixed effects, but using maximum
## likelihood rather than REML and thus use the LRT
model1 <- lmer(wear~(1|car_type),data=tyres,REML=F)
model2 <- lmer(wear~wheel+(1|car_type),data=tyres,REML=F)
# get the respective log-likelihoods
ll1 <- logLik(model1)
ll2 <- logLik(model2)
LRT <- -2*(ll1-ll2)
LRT <- as.numeric(LRT)
LRT
## 4 wheels so the models differ by 3 parameters (one wheel is buried in the intercept)
1 - pchisq(LRT,3)
# extremely small p-value so wheel position is overall significant. 

## Mow let's test the fixed effects of wheel using parameteric bootstrapping
sim_LRT <- 1:1000 # vector to store values of the LRT for each iteration
Dat <- simulate(model1,1000) ### Simulate 1000 data sets from model1 (the 
## smaller model, i.e. the NULL hypothesis)
for(i in 1:1000){ ### start loop
	Mod1 <- lmer(Dat[,i]~(1|car_type),data=tyres,REML=F) ### Fit model1 to simulated data
	Mod2 <- lmer(Dat[,i]~wheel+(1|car_type),data=tyres,REML=F) ### Fit model2 simulated data
	sim_LRT[i] <- -2*(logLik(Mod1)-logLik(Mod2)) ### Calculate and store LRT
} ## end loop
mean(sim_LRT>LRT) ### Calculate p-value (i.e. prop of sim_LRT values greater than LRT)
plot(density(sim_LRT),xlim=c(0,40),main="",xlab="LRT")
abline(v=LRT,col="red")
## which confirms that the wheel fixed effects are significant, however we 
## would trust these results more since the 
## approximations behind the LRT can be poor 
## sometimes, and this approach does not rely on them.

## More crudely, we can also just do a summary of model2 and treat 
## the random effects variance as fixed (known) and so 
## peform partial t-tests and F-tests (using analogous 
## degrees of freedom)
summary(model2)
# If random effects parameters are known and n=36, then n-p-1 is 36-3-1=32 so 
## we can do t-tests on 32 d.o.f.
qt(0.975,32) ## upper 97.5% quantile of the t-distribution
# Can also do an F-test to compare model1 and model2 under the assumption 
# of known random effects. Not going to do it here
# since it involves manual calculation of Deviances 
# (residual sum of squares). Just note that it can be done, but also that
# it is crude as it ingores the fact that random effects had to be estimated.

# So how about testing for "significance of the random effects"?  
# I.e. sigma^2_gamma=0.
# We can start by simply using the LRT 
# so fit a model without 
# random effects (i.e. a linear model)
model3 <- lm(wear~wheel,data=tyres)
## and compare this to model2 
model2 <- lmer(wear~wheel+(1|car_type),data=tyres,REML=F)
# using the LRT:
ll3 <- logLik(model3)
ll2 <- logLik(model2)
LRT <- -2*(ll3-ll2)
LRT <- as.numeric(LRT)
LRT
# There is only one parameter to separate these models (the 
# variance of the random effects) so
pchisq(LRT,1,lower.tail=F)
# implying that the Null hypothesis is not rejected (i.e. sigma^2_gamma is zero). 
# In other words, car type is not significant.
# However, we know that we shouldn't be using the LRT for testing parameters 
# on the boundary of the parameter space. So let's do a parametric bootstrap LRT:
sim_LRT <- 1:1000 # vector to store values of the LRT for each iteration
Dat <- simulate(model3,1000) ### Simulate 1000 data sets from model3 (the smaller model)
for(i in 1:1000){ ### start loop
	Mod3 <- lm(Dat[,i]~wheel,data=tyres) ### Fit model3 to simulated data
	Mod2 <- lmer(Dat[,i]~wheel+(1|car_type),data=tyres,REML=F) ### Fit model2 simulated data
	sim_LRT[i] <- -2*(logLik(Mod3)-logLik(Mod2)) ### Calculate and store LRT
} ## end loop
mean(sim_LRT>LRT) ### Calculate p-value (i.e. prop of sim_LRT values greater than LRT)
plot(density(sim_LRT),xlim=c(0,10),main="",xlab="LRT")
abline(v=LRT,col="red")
## The p-value is still bigger than 0.05 but only just, illustrating 
## that the LRT is indeed conservative. So the evidence 
## against the random effects being "significant" is not as overwhelming 
## as the LRT might have us think.

### So overall conclusion is that both fixed and random effects 
### are significant, so we go for model2
### but remember to fit this using REML=T to get unbiased 
### estimates for the variances.
model2 <- lmer(wear~wheel+(1|car_type),data=tyres,REML=T)

## OK but what about the equivalent of a partial t-test so that 
## we can say something about differences
## between wheels? 
summary(model2)
## First, we could assume the parameters of the random effects are 
## fixed (i.e. ignore the fact that we had
## to estimate them), and just "read-off" the values of the 
## t-tests and their implied significance. As we saw earlier,
## the upper 97.5% quantile of the t-distribution with 32 
## degrees of freedom is
qt(0.975,32) 
## so we look for t-test values whose absolute value is above this, 
## to judge significance. The "front offside" is not
## significant, so its effect is not different from "front nearside" 
## which is buried in the intercept. Both rear wheels 
## have a positive significant effect implying rear tyres wear 
## less than front ones (response is reduction in depth of thread).
## "Rear offside" tyres wear less than "rear nearside" (which makes sense).

## For more accurate inference (i.e. one that allows for the estimation 
## uncertainty of the random effect parameters), we can
## use parameteric bootstrap again. Fortunately, lme4 has the built-in 
## function confint() to do just that:
confint(model2,method="boot")
## These are 95% confidence intervals, so we look for zero being in the 
## intervals to judge whether a parameter 
## is not significantly different from zero. Inference is the same as 
## before, however with this we can be more
## confident of the accuracy of the significance.

## And what about predicted values for the random effects? 
# First, simple linear model for wear with car_type as a 
# categorical predictor (factor).
model1 <- lm(wear~car_type,data=tyres) 
summary(model1)
## Now fit this as an LMM where car_type is random effect
model2 <- lmer(wear~(1|car_type),data=tyres)
summary(model2)
## Get the estimates of the mean for each car from model1
predict(model1,newdata=data.frame(car_type=c("A","B","C")))
## We can get predictions of the random effects (as per lecture notes).
ranef(model2)
## which we can add to the estimate of the intercept beta0 to 
## get equivalent estimates of the mean
fixef(model2)
fixef(model2) + ranef(model2)$car_type
# Qualitatively (in terms of magnitude) there is agreement, 
# however the random effects model 
# "shrinks" the effects towards the mean

# What about some uncertainty about the random effects? 
# Well these are predictions (not estimates) so no point 
# in talking about partial t-tests. But we can still get 
# 95% confidence intervals and visualise those (note 
# that this is beyond the lecture material, just 
# showing it for illustration):
library(lattice)
model2 <- lmer(wear~wheel+(1|car_type),data=tyres)
dotplot(ranef(model2,condVar=T))
## So car A has least amount of wear and car C is worst. 
## Note that to get a "heuristic" equivalent of "significane",
## we look at whether the predition of a random effect 
## lies in the intervals of another. So prediction of car B
## is just outside the CI of car A, while car C is well outside 
## of the CI of car A. Compare with:
lmodel <- lm(wear~wheel+car_type,data=tyres)
summary(lmodel)
## Note however that it now much easier to compare cars B and C, 
## i.e. they are not significantly different.

## What about predictions from the model? Well, let's look 
## at the fixed and random effects:
model2 <- lmer(wear~wheel+(1|car_type),data=tyres,REML=T)
fixef(model2)
ranef(model2)
## Let's predict mean wear of front offside wheel for car A:
y_hat <- 17.518889 + 4.102222 - 1.9320844
y_hat
## which is the same as 
predict(model2,newdata=data.frame(wheel="front offside",car_type="A"))
## Now predict the wear of front offside wheel on a new car:
y_hat <- 17.518889 + 4.102222
y_hat
## which is the same as 
predict(model2,newdata=data.frame(wheel="front offside"),re.form=~0)
## However we do not get any standard errors with these predictions 
## as they are not easy obtain in general. We can use bootstrapping 
## again but we don't do it here.


## What about residuals from this model? Function residuals() will provide 
## deviance residuals by default, which should be Normally distributed:
resids <- residuals(model2)
## Note, that the fitted values (y_hat) used for getting the residuals 
## are based on the conditional mean,
## i.e. estimates of fixed effects plus predictions of random effects.
## So to get these fitted values, let's use the predict() function, 
## which by default will produce y_hat as above.
fitted <- predict(model2)
# Now produce redidual plots:
x11(width=16,height=6)
par(mfrow=c(1,2),mar = c(4, 4, 1, 1),cex=1.2,lwd=2)
plot(fitted,resids,ylab="Deviance residuals",xlab="fitted values",
     main="Residuals vs fitted values",pch=20)
abline(h=0)
qqnorm(resids,pch=20)
qqline(resids)
## Resids vs fits looks OK in terms even scatter about the x=0 line. 
## However, points seem to cluster in groups of 3
## so maybe we need an interaction between wheel and car type. 
## QQ plot not looking very good.

###########################################################################
## Note the interaction model below is NOT EXAMINABLE and I will
## not go over this code in the lectures, it's just
## here for our reference.
# So let's try a model with interaction, i.e. a model 
# where the wheel effect is random:
model3 <- lmer(wear~wheel+(1+wheel|car_type),data=tyres,REML=T)
summary(model3)
## Let's look at the random effect predictions
ranef(model3)
## which shows what the model is doing: three car_type effects for each wheel position.
dotplot(ranef(model3,condVar=T))

# What do the residuals look like?
resids <- residuals(model3)
fitted <- predict(model3)
x11(width=16,height=6)
par(mfrow=c(1,2),mar = c(4, 4, 1, 1),cex=1.2,lwd=2)
plot(fitted,resids,ylab="Deviance residuals",xlab="fitted values",
     main="Residuals vs fitted values",pch=20)
abline(h=0)
qqnorm(resids,pch=20)
qqline(resids)
## Much better than before so model3 looks like the better option. 
## We should of course formally test
## whether model3 is better than model2 using parametric bootstrapping.
## END OF NON-EXAMINABLE CODE FOR INTERACTIONS
###########################################################################





#############################################################################
########## Play with the pupils data (Question 9 on sheet 2 #################
#############################################################################
library(lme4)
head(pupils)

# These data involve language scores in Dutch schools. 
# (example of a multilevel model). Specifically, the data considers
# 131 schools (1 class per school) a total of 2287 students.
# Interest lies in assessing the impact on language scores (test) 
# of a mixture of individual pupil factors such as IQ (IQ) and social
# status (ses) as well as class effects.

# First fit a simple Gaussian (fixed effects) GLM to individual pupil 
# variables and class effects. 
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


# Now switch to a model with class as a random effect
model2 <- lmer(test~IQ+ses+(1|Class),data=pupils)
summary(model2)

# Assuming the random effects variance is equal to its estimate, both IQ 
# and ses (fixed effects) are significant with 
# similar estimates (and standard errors) to those obtained from 
# model1. Significance is assessed by assuming the quantiles
# 2.5% and 97.5% quantiles of the t-distr. are roughly -2 and +2. 
# More formally need to find these using n-p-1 degrees of freedom 
# with n=2287 and p=2.

# The conditional varince of the response 40.049 is roughly the same 
# as the residual variance of model1:
summary(model1)$dispersion
# (they are different due to numerical approximations). 

# Variance of the random efects is 9.212 and we can test if this 
# is significant using the LRT:
model2 <- lmer(test~IQ+ses+(1|Class),data=pupils,REML=F) # refit the model using max. likelihood
model0 <- glm(test~IQ+ses,data=pupils,family=gaussian(link="identity"))
LRT <- -2*(logLik(model0)-logLik(model2))
LRT <- as.vector(LRT)
LRT
pchisq(LRT,1,lower.tail=F)
# so the variance of the random effects is not zero and hence 
# Class is significant. Strictly the test is not valid 
# but also conservative so we can be fairly certain that 
# the true p-value is even smaller than the one above.


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



##################################################################
## Play with hip fracture data from Portugal.
## 
head(hip)
## At municipality we have knowledge only of the socio-economic
## status and the population. Allowing for a municipality-level
## random effect will capture the aggregate effect of other 
## covariates we might be missing. In addition, the random effects
## will allow for correlation of values withing each municipality.

## So fit Poisson model with log population as offset and municipality
## random effect.
rmodel <- glmer(Nfract~offset(log(Npop))+sex+ses+(1|municipality),
                family=poisson,data=hip)
summary(rmodel)
## Parameters are estimated using maximum likelihood so unlike 
## LMMs, we get z-tests and associated p-values as in a Poisson
## GLM. The rate of incidence of hip fractures is significantly larger 
## for females. Also, the incidence rate in “average” ses (=2) is 
## significantly smaller than the incidence rate in poor (ses=1). 
## Surprisingly, the incidence rate for affluent municipalities (ses=3)
## is the same as for poor, since the effect of ses3 is not significant.

## Check residuals
par(mfrow=c(1,2))
resids <- resid(rmodel)
qqnorm(resids,pch=20,main="Residuals")
qqline(resids)
r.effects <- ranef(rmodel)$municipality[,1]
qqnorm(r.effects,pch=20,main="Random effects")
qqline(r.effects)
## Unfortunately QQ plot not great at lower tails but not much we can do 
## about it here.

## Now let's test the significance of the random effects using the LRT
## First fit Poisson GLM with no random effects:
model <- glm(Nfract~offset(log(Npop))+sex+ses,family=poisson,data=hip)
## Now we use a chi^2 test with 1 degree of freedom as the only parameter
## to separate the two models is the variance of the random effects:
LRT <- 2*( logLik(rmodel)-logLik(model) )
LRT
1-pchisq(as.numeric(LRT) ,1)
## So random effects are significant. But note we should be strictly speaking
## using a parametric bootstrap. We also know however that the test is 
## conervative, so we can be sure that the p-value is even smaller.