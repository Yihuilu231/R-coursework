rm(list = ls())

load("ECM3712.RData")

contraceptive
contraceptive$age

model1 <- glm(count~(age+education+wantsMore)^2,data=contraceptive,family=poisson(link="log"))
summary(model1)
deviance(model1) #380.2587
drop1(model1)

model2 <- glm(count~(age*education*wantsMore),data=contraceptive,family=poisson(link="log"))
summary(model2)
deviance(model2) #378.6216
drop1(model2)

head(contraceptive2)
model3 <- glm(using/notUsing~(age*education*wantsMore),data=contraceptive2,family=poisson(link="log"))
summary(model3)
deviance(model3) #-3.486719e-16
drop1(model3)

model4 <- glm(using/notUsing~(age+education+wantsMore)^2,data=contraceptive2,family=poisson(link="log"))
summary(model4)
deviance(model4) #0.07261349
drop1(model4)
