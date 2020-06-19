install.packages("tseries")
install.packages("fracdiff")
install.packages("TTR")
install.packages("forecast")
install.packages("FitARMA")
install.packages("DMwR")
install.packages("quantmod")
install.packages("dlm")
install.packages("PerformanceAnalytics")
library(tseries)
library(fracdiff)
library(TTR)
library(forecast)
library(FitARMA)
library(DMwR)
library(quantmod)
library(dlm)
library(PerformanceAnalytics)

rapid <- read.csv('moc_transports.csv')
rapid$qyyyy <- paste(rapid$year,rapid$Quarter, sep='-')
rapidmean <- tapply(rapid$Overturning_Strength,rapid$qyyyy,mean)

# Convert it to a time series object
rapidmean.ts <- ts(as.vector(rapidmean),start=c(2004,2),frequency = 4)
rapidmean.ts
print(rapidmean.ts)
summary(rapidmean.ts)
ts.plot(rapidmean.ts) # This will plot the time series
# ACF help us determine what type of series we have, whether it is a White noise, Random walk, Auto regressive or Moving average
# to assess whether a time series is dependent on its past
acf(rapidmean.ts)
plot(ACF,type="h",ylim=c(-.8,1))
pacf(rapidmean.ts)
adf.test(rapidmean.ts)
auto.arima(rapidmean.ts)
#Fitting the AR Model to the time series
#AR model is an ARIMA(1, 0, 0) model
AR <- arima(rapidmean.ts, order = c(1,0,0))
print(AR)

#plotting the series along with the fitted values
ts.plot(rapidmean.ts)
AR_fit <- rapidmean.ts - residuals(AR)
points(AR_fit, type = "l", col = 2, lty = 2)

# Forecasting using AR model
# predict() function can be used to make forecasts from an estimated AR model
#the $pred value is the forceast, and the $se value is the standard error for the forceast
#Using predict() to make a 1-step forecast
predict_AR <- predict(AR)

#Obtaining the 1-step forecast using $pred[1]
predict_AR$pred[1]
#ALternatively Using predict to make 1-step through 6-step forecasts
predict(AR, n.ahead = 6)

#plotting the AirPassenger series plus the forecast and 95% prediction intervals
ts.plot(rapidmean.ts, xlim = c(2004, 2015))
AR_forecast <- predict(AR, n.ahead = 6)$pred
AR_forecast_se <- predict(AR, n.ahead = 6)$se
points(AR_forecast, type = "l", col = 2)
points(AR_forecast - 2*AR_forecast_se, type = "l", col = 2, lty = 2)
points(AR_forecast + 2*AR_forecast_se, type = "l", col = 2, lty = 2)


