require(dlm)

co2 <- read.table('co2_mm_mlo.txt',header=T)
#
# linear trend + monthly seasonal term

co2$average[co2$average< -90]=NA
plot(co2$date,co2$average)

#  make co2av a variable of class ts

co2av <- as.ts(co2$average)

dlmco2 <- dlmModPoly(order=3) + dlmModSeas(12)
buildFun <- function(x) {
  diag(W(dlmco2))[2:12] <- exp(x[1:11])
  V(dlmco2) <- exp(x[12])
  return(dlmco2)
}

# fit by MLE

(fit <- dlmMLE(co2av, parm = rep(0, 12), build = buildFun))$conv
dlmco2 <- buildFun(fit$par)
drop(V(dlmco2))
diag(W(dlmco2))[2:11]


# smoothing estimates of states
co2Smooth <- dlmSmooth(co2av, mod = dlmco2)
xs <- cbind(co2av, dropFirst(co2Smooth$s[,c(1,12)]))
colnames(xs) <- c("CO2", "Trend", "Seasonal")
plot(xs, type = 'o', main = "CO2 Concentration")

co2filt <- dlmFilter(co2av,mod=dlmco2)
co2fore <- dlmForecast(co2filt,nAhead=10)

sqrtR <- sapply(co2fore$R, function(x) sqrt(x[1,1]))
pl <- co2fore$a[,1] + qnorm(0.05, sd = sqrtR)
pu <- co2fore$a[,1] + qnorm(0.95, sd = sqrtR)
x <- ts.union(window(co2av, start = c(500, 1)), window(co2Smooth$s[,1], start = c(500, 1)),co2fore$a[,1], pl, pu) 
plot(x, plot.type = "single", type = 'o', pch = c(1, 0, 20, 3, 3),col = c("darkgrey", "darkgrey", "brown", "yellow", "yellow"), ylab = "CO2 concentration")
legend("bottomright", legend = c("Observed", "Smoothed (deseasonalized)", "Forecasted level", "90% probability limit"), bty = 'n', pch = c(1, 0, 20, 3, 3), lty = 1,col = c("darkgrey", "darkgrey", "brown", "yellow", "yellow"))

