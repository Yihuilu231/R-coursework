require(MASS)
rm(list = ls())

load("ECM3712.RData")

head(gehan)

# Extract the uncensored data and fit an exponential model in R this is done by fitting a Gamma model and then
# setting the dispersion parameter to 1 note the default link for the Gamma 
# is inverse so we are fitting a model where time_i~exp(lambda_i) and 
# lambda_i=1/mu_i=b_0+b1*treat_i (the control effect is the intercept b_0 
# the treatment effect is b_0+b1)

# Create new data frame with only uncensored data in it
gehan$treat <- relevel(gehan$treat, "control")
gehan
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
model1 <- glm(cens~(-1+time+time:treat),family=poisson(link="identity"),data=gehan2)
summary(model1,dispersion=1)


