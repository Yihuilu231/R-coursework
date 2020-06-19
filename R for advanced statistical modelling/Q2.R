load("ECM3712.RData")

#(b)
nlmodel
mylike <- function(theta,x,y){ 
  sum((y - (theta[1]*x)/(theta[2] + x))^2) 
}
# Now need to explore and find sensible initial estimates of the parameters 
# by plotting the data, guessing some initial parameter values, and 
# superimposing the model curve using those values
plot(nlmodel$x, nlmodel$y,xlab="x",ylab="y",pch=19) #non-linear
xfit <- seq(.2, 1.1)
theta.guess <- c(100,0.1)
yfit <- (theta.guess[1] * xfit) / (theta.guess[2] + xfit)
lines(xfit, yfit, lwd=3,col="red",lty=3)
# Not very close, could do better, may need some trial and error.
# Then need to fit using nlm(), find standard errors etc. from Hessian of the 
# log-likelihood evaluated at it's maximum correcting for the sigma^2 
# term - not done here - up to you.


out <- nlm(mylike, p = theta.guess, hessian = T, gradtol = 1e-10, 
           iterlim = 1000)

