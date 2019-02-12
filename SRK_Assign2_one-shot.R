# Install the package for getting kernel density estimation.
# documentation: https://cran.r-project.org/web/packages/kdensity/kdensity.pdf
install.packages("kdensity")
# Install the package for getting maximum penalized likelihood estimation.
# documentation: https://cran.r-project.org/web/packages/gss/gss.pdf
install.packages("gss")
# Load both of these packages
require(kdensity)  
require(gss)


## Simulate data from normal(2,3), exponential(10) and beta(3,1)
set.seed(1)
n <- 100
y1 <- rnorm(n,2,3)
y2 <- rexp(n,2)
y3 <- rbeta(n,3,1)

## Calculate the kernel density estimates
# Choose bandwidth as nrd0 and kernel as gaussian, since the data are symmetric.
f1Kde <- kdensity(y1,bw='nrd0',kernel='gaussian', support=c(range(y1)[1],range(y1)[2]))
# Choose bandwidth as ucv and kernel as gamma, since the data are asymmetric and
# the support of gamma kernel covers the range of simulated data.
f2Kde <- kdensity(y2,bw='ucv',kernel='gamma', support=c(range(y2)[1],range(y2)[2]))
# Choose bandwidth as ucv and kernel as beta since the data are asymmetric and
# the support of beta kernel covers the range of simulated data.
f3Kde <- kdensity(y3,bw='ucv',kernel='beta',support=c(range(y3)[1],range(y3)[2]))

## Calculate the maximum penalized likelihood density estimates; theoretically, 
## the true domain for the density is unknown, so set the domain based on the data simulated.
f1Mple <- ssden(~y1,domain=data.frame(y1=c(range(y1)[1],range(y1)[2])))
f2Mple <- ssden(~y2,domain=data.frame(y2=c(range(y2)[1],range(y2)[2])))
f3Mple <- ssden(~y3,domain=data.frame(y3=c(range(y3)[1],range(y3)[2])))

## Plot the true density function and estimated density function, the domain is 
## set based on the data we simulated as well.
# Plot for Normal(2,3) density and estimated densities.
yy1 <- seq(range(y1)[1],range(y1)[2],len=100)
plot(yy1, dnorm(yy1, 2, 3),type = "l", lwd = 2, ylim = c(0,0.16), xlim=c(range(y1)[1],range(y1)[2]),
     ylab = 'p.d.f',xlab = 'y1',bty = 'n')  # Plot the true density.
lines(yy1, f1Kde(yy1), lwd = 1, col = 'green')  # Plot kernel density estimate.
lines(yy1, dssden(f1Mple, yy1), lwd = 1, col = 'red')  # Plot maximum penalized likelihood estimate
legend(3.5,0.17,c("True density function","Kernel densit","Maximum penalized likelihood"), 
       col = c("black","green","red"), lty=1, cex = 0.3, box.lty=0)

# Plot for Exponential(2) density and estimated densities.
yy2 <- seq(range(y2)[1],range(y2)[2],len=100)
plot(yy2, dexp(yy2, 2),type="l", lwd=2, ylim=c(0,2.2), xlim=c(range(y2)[1],range(y2)[2]), ylab = 'p.d.f',
     xlab = 'y2', bty='n')  # Plot the true density.
lines(yy2, f2Kde(yy2), lwd = 1, col = 'green')  # Plot kernel density estimate.
lines(yy2, dssden(f2Mple,yy2), lwd = 1, col = 'red')  # Plot maximum penalized likelihood estimate.
legend(0.8,2,c("True density function","Kernel density","Maximum penalized likelihood"), 
       col = c("black","green","red"), lty=1, cex = 0.3, box.lty=0)

# Plot for Beta(3,1) density and estimated densities.
yy3 <- seq(range(y3)[1],range(y3)[2],len=100)
plot(yy3, dbeta(yy3,3,1),type="l", lwd=2, ylim=c(0,4), xlim=c(range(y3)[1],range(y3)[2]), ylab = 'p.d.f',
     xlab = 'y3', bty='n')  # Plot the true density.
lines(yy3, f3Kde(yy3), lwd = 1, col = 'green')  # Plot kernel density estimate.
lines(yy3, dssden(f3Mple, yy3), lwd = 1, col = 'red')  # Plot maximum penalized likelihood estimate.
legend(0.2,3.5,c("True density function","Kernel density","Maximum penalized likelihood"), 
       col = c("black","green","red"), lty=1, cex = 0.3, box.lty=0)


