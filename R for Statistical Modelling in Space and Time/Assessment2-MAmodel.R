#an MA model is an ARIMA(0, 0, 1) model

#Fitting the MA model 
MA <- arima(rapidmean.ts, order = c(0,0,1))
print(MA)

#plotting the series along with the MA fitted values
ts.plot(rapidmean.ts)
MA_fit <- rapidmean.ts - resid(MA)
points(MA_fit, type = "l", col = 2, lty = 2)

# Forecasting using MA model
#Making a 1-step forecast based on MA
predict_MA <- predict(MA)

#Obtaining the 1-step forecast using $pred[1]
predict_MA$pred[1]
#Alternately Making a 1-step through 6-step forecast based on MA
predict(MA,n.ahead=6)

#Plotting the series plus the forecast and 95% prediction intervals
ts.plot(rapidmean.ts, xlim = c(2004, 2015))
MA_forecasts <- predict(MA, n.ahead = 6)$pred
MA_forecast_se <- predict(MA, n.ahead = 6)$se
points(MA_forecasts, type = "l", col = 2)
points(MA_forecasts - 2*MA_forecast_se, type = "l", col = 2, lty = 2)
points(MA_forecasts + 2*MA_forecast_se, type = "l", col = 2, lty = 2)

#Goodness of fit
# (AIC) and Bayesian information criterion (BIC) are used for Time series Models.
# Find correlation between AR_fit and MA_fit
cor(AR_fit, MA_fit)
AIC(AR)
AIC(MA)
BIC(AR)
BIC(MA) #CHoosing the lower one

# MA model is better
#Fitting ARMA( p,q) models with arima()

# Question 2
# linear trend + a seasonal value
dlmrapid <- dlmModPoly(order = 3) + dlmModSeas(4)
buildFun <- function(x) {
  diag(W(dlmrapid))[2:4] <- exp(x[1:4])
  V(dlmrapid) <- exp(x[4])
  return(dlmrapid)
}

# fit by MLE
(fit <- dlmMLE(rapidmean.ts, parm = rep(0, 4), build = buildFun))$conv
dlmrapid <- buildFun(fit$par)
drop(V(dlmrapid))
diag(W(dlmrapid))[2:3]


# smoothing estimates of states
rapidSmooth <- dlmSmooth(rapidmean.ts, mod = dlmrapid)
xs <- cbind(rapidmean.ts, dropFirst(rapidSmooth$s[,c(1,4)]))
colnames(xs) <- c("Transports", "Trend", "Seasonal")
plot(xs, type = 'o', main = "moc_Transports")

rapidfilt <- dlmFilter(rapidmean.ts,mod=dlmrapid)
rapidfore <- dlmForecast(rapidfilt,nAhead=6)
rapidfore


