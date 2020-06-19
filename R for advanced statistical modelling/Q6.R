rm(list = ls())

load("ECM3712.RData")

titanic
## Start by fitting model with all two-way interactions and main effects
model1 <- glm(survived~(pclass+gender+age)^2,data=titanic,family=binomial(link="logit"))
summary(model1)
## The main effect for age is not significant but some interactions with age are. To see if we
## can reduce this model, we first try to remove each interaction term and compare the AIC each time
## to decide if the term is to be dropped. Function drop1() will do that for us.
AIC(model1) #1137.395
deviance(model1) #1117.395
drop1(model1)
# Up to you decide if it should be reduced or not
# 
## Now look at the same data, after grouping them into a binomial structure
head(titanic2)
head(titanic[,-3])
head(titanic2)
## Fit binomial model equivalent to model1
model2 <- glm(cbind(survived,total-survived)~(age_group+pclass+gender)^2,data=titanic2,family=binomial(link="logit"))
summary(model2)
# Look at the value of the deviance D_M
deviance(model2) #3.508092e-10
AIC(model2) # 56.43098
# Effectively zero. What does this imply about the model? Go back to the general expression of the Deviance D_M 
# for a GLM on top of slide 18 in Topic 2 and think about what values of theta_hat would be needed to make this zero.

## Can we reduce this model based on the AIC?
drop1(model2)
# Well, no according to the table, however we should do so anyway so reduce to a model that will procuce the least worst fit
# in terms of the AIC. Up to you to do the rest...

head(titanic)
head(titanic[,-3])
head(titanic2)
head(titanic3)
model3 <- glm(survived~(age+pclass+sex)^2,data=titanic3,family=binomial(link="logit"))
AIC(model3) #53.27106 i.e. model3 is better than model2
summary(model3)
deviance(model3)
model3$residuals
drop1(model3)

model4 <- glm(cbind(survived,total-survived)~(pclass+gender)^2,data=titanic2,family=binomial(link="logit"))
summary(model4)
AIC(model4)
drop1(model4)
plot(model3,which=1,pch=20)
plot(model3,which=2,pch=20)
plot(model3,which=5,pch=20)


