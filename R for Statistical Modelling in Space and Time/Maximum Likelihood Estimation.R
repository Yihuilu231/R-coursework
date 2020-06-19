require(geoR)
set.seed(766)
# generating a simulated data-set
ex.data <- grf(70, cov.pars=c(10, .15), cov.model="matern", kappa=2)
#
# defining the grid of prediction locations:
ex.grid <- as.matrix(expand.grid(seq(0,1,l=21), seq(0,1,l=21)))

ml <- likfit(ex.data, ini=c(0.5, 0.5), fix.nug = TRUE, cov.model='matern',kappa=2)
ml
summary(ml)
reml <- likfit(ex.data, ini=c(0.5, 0.5), fix.nug = TRUE, cov.model='matern', kappa=2,
               lik.met = "REML")
summary(reml)
plot(variog(ex.data))
lines(ml)
lines(reml, lty = 2)

xv.ml <- xvalid(ex.data, model = ml)

par(mfcol = c(5, 2), mar = c(3, 3, 1, 0.5), mgp = c(1.5, 0.7, 0)) 
plot(xv.ml)

xvR.ml <- xvalid(ex.data, model = ml, reest = TRUE)
par(mfcol = c(5, 2), mar = c(3, 3, 1, 0.5), mgp = c(1.5, 0.7, 0))
plot(xvR.ml)

par(mfcol=c(1,1))

plot(ex.data$coords, xlim = c(0, 1.), ylim = c(0,1.), xlab = "Coord X", ylab = "Coord Y") 
loci <- matrix(c(0.2, 0.6, 0.2, 0.9, 0.2, 0.3, 0.8, 0.9), ncol = 2) 
text(loci, as.character(1:4), col = "red") 
polygon(x = c(0, 1, 1, 0), y = c(0, 0, 1, 1), lty = 2)

#  ordinary kriging
kc4 <- krige.conv(ex.data, locations = loci, krige = krige.control(obj.m = ml))

kc4$predict
kc4$krige.var

pred.grid <- expand.grid(seq(0, 1, l = 51), seq(0, 1, l = 51)) 
kc <- krige.conv(ex.data, loc = pred.grid, krige = krige.control(obj.m = ml))

image(kc, loc = pred.grid, col = gray(seq(1, 0.1, l = 30)), xlab = "Coord X", ylab = "Coord Y")