### R code from vignette source 'dlm.Rnw'

#  Comments modified by PC 27/2/14

###################################################
### code chunk number 1: dlm.Rnw:26-29
###################################################
library(dlm)
options(digits=3)
library(dlm)
set.seed(1963)


# defining models

dlm(FF = 1, V = 0.8, GG = 1, W = 0.1, m0 = 0, C0 = 100)


dlmModPoly(order = 1, dV = 0.8, dW = 0.1, C0 = 100)


# Lots of defaults

myMod <- dlmModPoly()


#  see what's in the model and modify
FF(myMod)
W(myMod)
m0(myMod)
V(myMod) <- 0.8


#Adding models

myMod <- dlmModPoly() + dlmModSeas(4)


#  An alternative way of adding models

dlmModPoly(dV = 0.2, dW = c(0, 0.5)) %+% 
  (dlmModSeas(4, dV = 0, dW = c(0, 0, 0.35)) + 
     dlmModPoly(1, dV = 0.1, dW = 0.03))

#  Dynamic linear regression 

u <- rnorm(25)
myMod <- dlmModReg(u, dV = 14.5)
myMod$JFF
head(myMod$X)


#  Estimation for the Nile River data

data(Nile)
plot(Nile)

#  Build model (random walk + noise; unknown variances - log scale)

buildFun <- function(x) {
  dlmModPoly(1, dV = exp(x[1]), dW = exp(x[2]))
}


#  fit the data

fit <- dlmMLE(Nile, parm = c(0,0), build = buildFun)
fit$conv
dlmNile <- buildFun(fit$par)
V(dlmNile)
W(dlmNile)


# Atternative code

StructTS(Nile, "level")


#  Build a new model with jump (change of variance) at 1898

buildFun <- function(x) {
  m <- dlmModPoly(1, dV = exp(x[1]))
  m$JW <- matrix(1)
  m$X <- matrix(exp(x[2]), nc = 1, nr = length(Nile))
  j <- which(time(Nile) == 1899)
  m$X[j,1] <- m$X[j,1] * (1 + exp(x[3]))
  return(m)
}
fit <- dlmMLE(Nile, parm = c(0,0,0), build = buildFun)
fit$conv # See whether is converged or not
dlmNileJump <- buildFun(fit$par)
V(dlmNileJump) # Variance of dlmNileJump What can we use it for?
dlmNileJump$X[c(1, which(time(Nile) == 1899)), 1]

# Filtering and smoothing

# assume we know the parameters and filter

nileJumpFilt <- dlmFilter(Nile, dlmNileJump)
plot(Nile, type = 'o', col = "seagreen")
lines(dropFirst(nileJumpFilt$m), type = 'o', 
      pch = 20, col = "brown")


#  stores variance matrices as SVD 

attach(nileJumpFilt)
v <- unlist(dlmSvd2var(U.C, D.C))
pl <- dropFirst(m) + qnorm(0.05, sd = sqrt(v[-1]))
pu <- dropFirst(m) + qnorm(0.95, sd = sqrt(v[-1]))
detach()
lines(pl, lty = 2, col = "brown") 
lines(pu, lty = 2, col = "brown") 




#  means and variances of the smoothing distributions


nileJumpSmooth <- dlmSmooth(nileJumpFilt)
plot(Nile, type = 'o', col = "seagreen")
attach(nileJumpSmooth)
lines(dropFirst(s), type = 'o', pch = 20, col = "brown")
v <- unlist(dlmSvd2var(U.S, D.S))
pl <- dropFirst(s) + qnorm(0.05, sd = sqrt(v[-1]))
pu <- dropFirst(s) + qnorm(0.95, sd = sqrt(v[-1]))
detach()
lines(pl, lty = 2, col = "brown") 
lines(pu, lty = 2, col = "brown") 




# Second example UK Gas consumption
# linear trend + quarterly seasonal term

plot(UKgas)
lGas <- log(UKgas)
plot(lGas)
dlmGas <- dlmModPoly() + dlmModSeas(4)
buildFun <- function(x) {
  diag(W(dlmGas))[2:3] <- exp(x[1:2])
  V(dlmGas) <- exp(x[3])
  return(dlmGas)
}

# fit by MLE

(fit <- dlmMLE(lGas, parm = rep(0, 3), build = buildFun))$conv
dlmGas <- buildFun(fit$par)
drop(V(dlmGas))
diag(W(dlmGas))[2:3]


# smoothing estimates of states
gasSmooth <- dlmSmooth(lGas, mod = dlmGas)
x <- cbind(lGas, dropFirst(gasSmooth$s[,c(1,3)]))
colnames(x) <- c("Gas", "Trend", "Seasonal")
plot(x, type = 'o', main = "UK Gas Consumption",cex=5)


#  Forecasting
plot(x, type = 'o', main = "UK Gas Consumption",cex=5)

gasFilt <- dlmFilter(lGas, mod = dlmGas)
gasFore <- dlmForecast(gasFilt, nAhead = 20)
gasFore
sqrtR <- sapply(gasFore$R, function(x) sqrt(x[1,1]))
pl <- gasFore$a[,1] + qnorm(0.05, sd = sqrtR)
pu <- gasFore$a[,1] + qnorm(0.95, sd = sqrtR)
x <- ts.union(window(lGas, start = c(1982, 1)), 
              window(gasSmooth$s[,1], start = c(1982, 1)),
              gasFore$a[,1], pl, pu) 
plot(x, plot.type = "single", type = 'o', pch = c(1, 0, 20, 3, 3),
     col = c("darkgrey", "darkgrey", "brown", "yellow", "yellow"),
     ylab = "Log gas consumption")
legend("bottomright", legend = c("Observed", 
                                 "Smoothed (deseasonalized)", 
                                 "Forecasted level", "90% probability limit"),
       bty = 'n', pch = c(1, 0, 20, 3, 3), lty = 1,
       col = c("darkgrey", "darkgrey", "brown", "yellow", "yellow"))




# Bayesian Estimation

#  Simulate from the dlm with MLEs
#  model without jump

plot(Nile, type = 'o', col = "seagreen",cex=5)
nileFilt <- dlmFilter(Nile, dlmNile)
for (i in 1:10) # 10 simulated "true" levels 
  lines(dropFirst(dlmBSample(nileFilt)), col = "brown") 




###################################################
### code chunk number 25: dlm.Rnw:658-662
###################################################
lmixnorm <- function(x, weights, means, sds) {
  log(crossprod(weights, exp(-0.5 * ((x - means) / sds)^2 
                             - log(sds))))
}


###################################################
### code chunk number 26: dlm.Rnw:669-681
###################################################
y <- arms(0, myldens = lmixnorm, 
          indFunc = function(x,...) (x > (-100)) * (x < 100), 
          n = 5000, weights = c(1, 3, 2), 
          means = c(-10, 0, 10), sds = c(7, 5, 2))
summary(y)
library(MASS)
truehist(y, prob = TRUE, ylim = c(0, 0.08), bty = 'o')
curve(colSums(c(1, 3, 2) / 6 *
                dnorm(matrix(x, 3, length(x), TRUE), 
                      mean = c(-10, 0, 10), sd = c(7, 5, 2))), 
      add = TRUE)
legend(-25, 0.07, "True density", lty = 1, bty = 'n')


###################################################
### code chunk number 27: dlm.Rnw:685-690
###################################################
truehist(y, prob = TRUE, ylim = c(0, 0.08), bty = 'o')
curve(colSums(c(1, 3, 2) / 6 *
                dnorm(matrix(x, 3, length(x), TRUE), 
                      mean = c(-10, 0, 10), sd = c(7, 5, 2))), add = TRUE)
legend(-25, 0.07, "True density", lty = 1, bty = 'n')


# Gibbs sampler
outGibbs <- dlmGibbsDIG(lGas, dlmModPoly(2) + dlmModSeas(4),
                        a.y = 1, b.y = 1000, a.theta = 1, 
                        b.theta = 1000,
                        n.sample = 1100, ind = c(2, 3),
                        save.states = FALSE)


# diagnostics

burn <- 100
attach(outGibbs)
par(mfrow=c(2,3), mar=c(3.1,2.1,2.1,1.1))
plot(dV[-(1:burn)], type = 'l', xlab="", ylab="", main=expression(sigma^2))
plot(dW[-(1:burn),1], type = 'l', xlab="", ylab="", main=expression(sigma[beta]^2))
plot(dW[-(1:burn),2], type = 'l', xlab="", ylab="", main=expression(sigma[s]^2))
use <- length(dV) - burn
from <- 0.05 * use
at <- pretty(c(0,use),n=3); at <- at[at>=from]
plot(ergMean(dV[-(1:burn)], from), type="l", xaxt="n",xlab="", ylab="")
axis(1, at=at-from, labels=format(at))
plot(ergMean(dW[-(1:burn),1], from), type="l", xaxt="n",xlab="", ylab="")
axis(1, at=at-from, labels=format(at))
plot(ergMean(dW[-(1:burn),2], from), type="l", xaxt="n",xlab="", ylab="")
axis(1, at=at-from, labels=format(at))
detach()


###################################################
### code chunk number 30: dlm.Rnw:805-829
###################################################
burn <- 100
attach(outGibbs)
dV <- dV[-(1:burn)]
dW <- dW[-(1:burn), ]
detach()
par(mfrow=c(2,3), mar=c(3.1,2.1,2.1,1.1))
plot(dV, type = 'l', xlab = "", ylab = "", 
     main = expression(sigma^2))
plot(dW[ , 1], type = 'l', xlab = "", ylab = "", 
     main = expression(sigma[beta]^2))
plot(dW[ , 2], type = 'l', xlab = "", ylab = "", 
     main = expression(sigma[s]^2))
use <- length(dV) - burn
from <- 0.05 * use
at <- pretty(c(0, use), n = 3); at <- at[at >= from]
plot(ergMean(dV, from), type = 'l', xaxt = 'n',
     xlab = "", ylab = "")
axis(1, at = at - from, labels = format(at))
plot(ergMean(dW[ , 1], from), type = 'l', xaxt = 'n',
     xlab = "", ylab = "")
axis(1, at = at - from, labels = format(at))
plot(ergMean(dW[ , 2], from), type = 'l', xaxt = 'n',
     xlab = "", ylab = "")
axis(1, at =  at - from, labels = format(at))


###################################################
### code chunk number 31: dlm.Rnw:836-837
###################################################
mcmcMean(cbind(dV[-(1:burn)], dW[-(1:burn), ]))


###################################################
### code chunk number 32: dlm.Rnw:839-840
###################################################
rm(dV, dW)
