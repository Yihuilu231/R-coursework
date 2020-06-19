install.packages("TSA")
install.packages("astsa")
library(TSA)
library(astsa)
# Simulate from the AR(1)
# ARMA and ARIMA models are suitable for stationary dataset
plot(arima.sim(n=128,list(ar=0.5,sd=1)))
plot(arima.sim(n=128,list(ar=-0.5,sd=1)))
plot(arima.sim(n=128,list(ar=0.95,sd=1)))
plot(arima.sim(n=128,list(ar=-0.95,sd=1)))
plot(arima.sim(n=128,list(ar=1.05,sd=1)))

par(mfrow = c(2,1))
plot(arima.sim(list(order=c(0,0,1), ma=.5), n=100), ylab="x",
     main=(expression(MA(1)~~~theta==+.5)))
plot(arima.sim(list(order=c(0,0,1), ma=-.5), n=100), ylab="x",
     main=(expression(MA(1)~~~theta==-.5)))

ACF = ARMAacf(ar=c(0.5,0.1), ma=c(1.0,0.05), 24)[-1]
PACF = ARMAacf(ar=c(0.5,0.1), ma=c(1.0,0.05), 24, pacf=TRUE) 
par(mfrow=c(1,2))
plot(ACF, type="h", xlab="lag", ylim=c(-.8,1));  abline(h=0)
plot(PACF, type="h", xlab="lag", ylim=c(-.8,1));  abline(h=0)

require(TSA)
ARMAspec(model=list(ma=-0.9))

require(astsa)
arma.spec(ar=0,ma=0.9)

