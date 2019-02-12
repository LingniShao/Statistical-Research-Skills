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
x1 <- rnorm(n,2,3)
x2 <- rexp(n,2)
x3 <- rbeta(n,3,1)

## Calculate the kernel density estimates
# Choose bandwidth as nrd0 and kernel as gaussian, since the data are symmetric.
f1Kde <- kdensity(x1,bw='nrd0',kernel='gaussian', support=c(range(x1)[1],range(x1)[2]))
# Choose bandwidth as ucv and kernel as gamma, since the data are asymmetric and
# the support of gamma kernel covers the range of simulated data.
f2Kde <- kdensity(x2,bw='ucv',kernel='gamma', support=c(range(x2)[1],range(x2)[2]))
# Choose bandwidth as ucv and kernel as beta since the data are asymmetric and
# the support of beta kernel covers the range of simulated data.
f3Kde <- kdensity(x3,bw='ucv',kernel='beta',support=c(range(x3)[1],range(x3)[2]))

## Calculate the maximum penalized likelihood density estimates; theoretically, 
## the true domain for the density is unknown, so set the domain based on the data simulated.
f1Mple <- ssden(~x1,domain=data.frame(x1=c(range(x1)[1],range(x1)[2])))
f2Mple <- ssden(~x2,domain=data.frame(x2=c(range(x2)[1],range(x2)[2])))
f3Mple <- ssden(~x3,domain=data.frame(x3=c(range(x3)[1],range(x3)[2])))

## Plot the true density function and estimated density function, the domain is 
## set based on the data we simulated as well.
# Plot for Normal(2,3) density and estimated densities.
xx1 <- seq(range(x1)[1],range(x1)[2],len=100)
plot(xx1, dnorm(xx1, 2, 3),type = "l", lwd = 2, ylim = c(0,0.16), xlim=c(range(x1)[1],range(x1)[2]),
     ylab = 'p.d.f',xlab = 'x1 from Normal(2,3)',bty = 'n')  # Plot the true density.
lines(xx1, f1Kde(xx1), lwd = 1, col = 'green')  # Plot kernel density estimate.
lines(xx1, dssden(f1Mple, xx1), lwd = 1, col = 'red')  # Plot maximum penalized likelihood estimate
legend(3.5,0.17,c("True density function","Kernel density","Maximum penalized likelihood"), 
       col = c("black","green","red"), lty=1, cex = 0.3, box.lty=0)

# Plot for Exponential(2) density and estimated densities.
xx2 <- seq(range(x2)[1],range(x2)[2],len=100)
plot(xx2, dexp(xx2, 2),type="l", lwd=2, ylim=c(0,2.2), xlim=c(range(x2)[1],range(x2)[2]), ylab = 'p.d.f',
     xlab = 'x2 from Expnential(2)', bty='n')  # Plot the true density.
lines(xx2, f2Kde(xx2), lwd = 1, col = 'green')  # Plot kernel density estimate.
lines(xx2, dssden(f2Mple,xx2), lwd = 1, col = 'red')  # Plot maximum penalized likelihood estimate.
legend(0.8,2,c("True density function","Kernel density","Maximum penalized likelihood"), 
       col = c("black","green","red"), lty=1, cex = 0.3, box.lty=0)

# Plot for Beta(3,1) density and estimated densities.
xx3 <- seq(range(x3)[1],range(x3)[2],len=100)
plot(xx3, dbeta(xx3,3,1),type="l", lwd=2, ylim=c(0,4), xlim=c(range(x3)[1],range(x3)[2]), ylab = 'p.d.f',
     xlab = 'x3 from Beta(3,1)', bty='n')  # Plot the true density.
lines(xx3, f3Kde(xx3), lwd = 1, col = 'green')  # Plot kernel density estimate.
lines(xx3, dssden(f3Mple, xx3), lwd = 1, col = 'red')  # Plot maximum penalized likelihood estimate.
legend(0.2,3.5,c("True density function","Kernel density","Maximum penalized likelihood"), 
       col = c("black","green","red"), lty=1, cex = 0.3, box.lty=0)


