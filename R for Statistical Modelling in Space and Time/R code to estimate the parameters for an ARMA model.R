test1 <- arima.sim(n=1000,list(ar=c(.5,-.1),sd =1))
par(mfrow=c(1,2))
acf(test1,lag.max=20)
pacf(test1,lag.max=20)

require(astsa)
data(rec)
plot(rec)
par(mfrow=c(1,2))
acf(rec,lag.max=30)
pacf(rec,lag.max=30)

rec.yw <- ar.yw(rec, order=2)
rec.yw$x.mean
rec.yw$ar
sqrt(diag(rec.yw$asy.var.coef))
rec.yw$var.pred

rec.mle = ar.mle(rec, order=2) 
rec.mle$x.mean
rec.mle$ar 
sqrt(diag(rec.mle$asy.var.coef))  
rec.mle$var.pred 

rec.mle2 <- arima(rec, order=c(2,0,0))

data(soi)
plot(soi)
par(mfrow=c(1,2))
acf(soi,lag.max=20)
pacf(soi,lag.max=20)

require(tseries)
rec.arma <- arma(rec)
rec.arma2 <- arma(rec,order=c(2,0))
soi.arma <- arma(soi)

