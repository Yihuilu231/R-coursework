install.packages("geoR")
install.packages("mnormt")

require(geoR)
require(ggplot2)

set.seed(362)
#  Read data
data <- read.csv('GrandBanks_Dec_1997.csv')
gdata <- as.geodata(data,coords.col=2:3,data.col=6)
dup <-dup.coords(gdata) # locate duplicated coordinates
gdata2 <- jitterDupCoords(gdata,max=0.1,min=0.05)

par(mfrow=c(1,1))
plot(gdata2)
# The above belongs to my dear teacher.
summary(gdata2)

# Remove sst>28
which(data$sst>28)
data <- data[-366,]
which(data$sst>28)
data <- data[-663]
summary(data)
# Now do it again
gdata <- as.geodata(data,coords.col=2:3,data.col=6)
dup <-dup.coords(gdata) # locate duplicated coordinates
gdata2 <- jitterDupCoords(gdata,max=0.1,min=0.05)

par(mfrow=c(1,1))
plot(gdata2)
# Now it's better.

# Question 2
summary(gdata2)

# variog4 Computes Directional Variograms
var4 <- variog4(gdata2,max.dist = 1)
plot(var4)

# Question 3
# Maximum likelihood estimates are obtained by fitting the model to the data. 
# Models can be compared by measures of fit, such as the log likelihood.
ml <- likfit(gdata2,ini = c(1,0.15))
#summary(ml)
reml <- likfit(gdata2,trend="1st",ini = c(1,0.15))
#summary(reml)
ml2 <- likfit(gdata2,trend="2nd",ini = c(1,0.15))
logLik(ml)
logLik(reml)
logLik(ml2)

#Model prediction
#Finally, we start the spatial prediction defining a grid of points. 
#The kriging function by default performs ordinary kriging. 
#It minimaly requires the data, prediction locations and estimated 
#model parameters.
gdata2.gr <- expand.grid((0:100)/100,(0:100)/100)
gdata2.kc <- krige.conv(gdata2,locations=gdata2.gr,krige=krige.control(obj.model=ml2))
names(gdata2.kc)

args(krige.control)
args(output.control)
OC <- output.control(simulations = TRUE, n.pred = 1000, quantile = c(0.1, 0.25, 0.5, 0.75,0.9), threshold = 350)
#xv.ml <- xvalid(gdata2, model = ml2)
#par(mfcol = c(5, 2), mar = c(3, 3, 1, 0.5), mgp = c(1.5, 0.7, 0)) 
#plot(xv.ml)
#xvR.ml <- xvalid(gdata2, model = ml2, reest = TRUE)
#par(mfcol = c(5, 2), mar = c(3, 3, 1, 0.5), mgp = c(1.5, 0.7, 0))
#plot(xvR.ml)
#par(mfcol=c(1,1))
#plot(gdata2$coords, xlim = c(0, 1.), ylim = c(0,1.), xlab = "Coord X", ylab = "Coord Y") 
#loci <- matrix(c(0.2, 0.6, 0.2, 0.9, 0.2, 0.3, 0.8, 0.9), ncol = 2) 
#text(loci, as.character(1:4), col = "red") 
#polygon(x = c(0, 1, 1, 0), y = c(0, 0, 1, 1), lty = 2)

#  ordinary kriging
kc4 <- krige.conv(gdata2, locations = loci, krige = krige.control(obj.m = ml2))
kc4$predict
kc4$krige.var
pred.grid <- expand.grid(seq(0, 1, l = 51), seq(0, 1, l = 51)) 
kc <- krige.conv(gdata2, loc = pred.grid, krige = krige.control(obj.m = ml2))
image(kc, loc = pred.grid, col = gray(seq(1, 0.1, l = 30)), xlab = "Coord X", ylab = "Coord Y")
#par(mfrow = c(2, 2), mar = c(3, 3, 0.5, 0.5))
#image(gdata2.kc, col = terrain.colors(21), x.leg = c(500, 750), y.leg = c(0, 50))
#image(gdata2.kc, val = gdata2.kc$quantile[, 3], col = terrain.colors(21), x.leg = c(500, 750),y.leg = c(0, 50)) 
#image(gdata2.kc, val = gdata2.kc$simulation[, 1], col = terrain.colors(21), x.leg = c(500,750), y.leg = c(0, 50))
#image(gdata2.kc, val = 1 - gdata2.kc$prob, col = terrain.colors(21), x.leg = c(500, 750), y.leg = c(0,50)) 
                                                                                                    
# Bayesian Inference
args(krige.bayes)
gdata2.bayes <- krige.bayes(gdata2, loc = gdata2.gr, model = model.control(trend.d = "1st", 
                           trend.l = "1st"), prior = prior.control(phi.prior = "rec", phi.disc = seq(0, 150, by = 15)), 
                            output = OC)
names(gdata2.bayes)
names(gdata2.bayes$posterior)
par(mfrow = c(1, 1))
plot(gdata2.bayes)
names(gdata2.bayes$predictive)
par(mfrow = c(1, 3))
image(gdata2.bayes, col = terrain.colors(21))
image(gdata2.bayes, val = apply(gdata2.bayes$pred$simul, 1, quantile, prob = 0.9), col = terrain.colors(21))
hist(apply(gdata2.bayes$pred$simul, 2, median), main = "")








