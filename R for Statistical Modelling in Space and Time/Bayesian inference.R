require(geoR)
# generating a simulated data-set
ex.data <- grf(70, cov.pars=c(10, .15), cov.model="matern", kappa=2)
#
# defining the grid of prediction locations:
ex.grid <- as.matrix(expand.grid(seq(0,1,l=21), seq(0,1,l=21)))
#
# computing posterior and predictive distributions
# (warning: the next command can be time demanding)
ex.bayes <- krige.bayes(ex.data, loc=ex.grid,
                        model = model.control(cov.m="matern", kappa=2),
                        prior = prior.control(phi.discrete=seq(0, 0.7, l=51),
                                              phi.prior="reciprocal"))
#
# Prior and posterior for the parameter phi
plot(ex.bayes, type="h", tausq.rel = FALSE, col=c("red", "blue"))
#   
# Plot histograms with samples from the posterior
par(mfrow=c(3,1))
hist(ex.bayes)
par(mfrow=c(1,1))

# Plotting empirical variograms and some Bayesian estimates:
# Empirical variogram
plot(variog(ex.data, max.dist = 1), ylim=c(0, 15))
# Since ex.data is a simulated data we can plot the line with the "true" model 
lines.variomodel(ex.data, lwd=2)
# adding lines with summaries of the posterior of the binned variogram
lines(ex.bayes, summ = mean, lwd=1, lty=2)
lines(ex.bayes, summ = median, lwd=2, lty=2)
# adding line with summary of the posterior of the parameters
lines(ex.bayes, summary = "mode", post = "parameters")

# Plotting again the empirical variogram
plot(variog(ex.data, max.dist=1), ylim=c(0, 15))
# and adding lines with median and quantiles estimates
my.summary <- function(x){quantile(x, prob = c(0.05, 0.5, 0.95))}
lines(ex.bayes, summ = my.summary, ty="l", lty=c(2,1,2), col=1)

# Plotting some prediction results
op <- par(no.readonly = TRUE)
par(mfrow=c(2,2), mar=c(4,4,2.5,0.5), mgp = c(2,1,0))
image(ex.bayes, main="predicted values")
image(ex.bayes, val="variance", main="prediction variance")
image(ex.bayes, val= "simulation", number.col=1,
      main="a simulation from the \npredictive distribution")
image(ex.bayes, val= "simulation", number.col=2,
      main="another simulation from \nthe predictive distribution")
#
par(op)
ex.bayes <- krige.bayes(ex.data, loc=ex.grid,
                        model = model.control(cov.m="matern", kappa=2),
                        prior = prior.control(phi.discrete=seq(0, 0.7, l=51),
                                              phi.prior="squared.reciprocal"))