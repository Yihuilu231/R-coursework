library(geoR)

data(parana)

# plot the data

points(parana)

# plot the variogram

plot(variog(parana,option='bin'))

# Limit the distance
plot(variog(parana,option='bin',max.dist=400))

# plot the cloud

plot(variog(parana,option='cloud',max.dist=400))

#  Try Hawkins and Cressie

plot(variog(parana,option='cloud',max.dist=400,estimator.type='modulus'))

plot(variog(parana,option='bin',max.dist=400,estimator.type='modulus'))

vari.para <- variog(parana,option='bin',max.dist=400,estimator.type='modulus',bin.cloud=T)
plot(vari.para,bin.cloud=T)

vari.default <- variofit(vari.para)
summary(vari.default)

vari.mat1.5 <- variofit(vari.para,kappa=1.5,fix.kappa=T)
summary(vari.mat1.5)

vari.mat1.5 <- variofit(vari.para,kappa=1.5,fix.kappa=T,weights='cressie')
summary(vari.mat1.5)

