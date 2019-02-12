# Install the package for getting kernel density estimation.
# documentation: https://cran.r-project.org/web/packages/kdensity/kdensity.pdf
install.packages("kdensity")
# Install the package for getting maximum penalized likelihood estimation.
# documentation: https://cran.r-project.org/web/packages/gss/gss.pdf
install.packages("gss")
# Load both of these packages
require(kdensity)  
require(gss)


simu_norm <- function(n){
  # Given a data from Normal(2,3), plot true density and estimated densities
  # using kernel density and maximum penalized likelihood and print out integrated
  # squared error for both methods.
  
  x <- rnorm(n,2,3)  # simulated data from Normal(2,3).
  low <- min(x)
  upp <- max(x)
  
  # Get Kernel density estimate, with bandwidth as nrd0 and kernel as gaussian, 
  # since the data are symmetric.
  fKde <- kdensity(x,bw='nrd0',kernel='gaussian',support=c(low,upp))
  # Get maximum penalized likelihood estimate. Theoretically, the true domain for 
  # the density is unknown, so set the domain based on the data simulated.
  fMple <- ssden(~x,domain=data.frame(x=c(low,upp)))
  
  # Plot the true density and estimated densities, the domain is set based on the
  # data we simulated as well.
  x1 <- seq(low,upp,len=500)
  plot(x1, dnorm(x1, 2, 3),type = "l", lwd = 2, ylim = c(0,0.16), xlim=c(low,upp),
       ylab = 'p.d.f',xlab = 'x from Normal(2,3)',bty = 'n')  # Plot the true density.
  lines(x1, fKde(x1), lwd = 1, col = 'green')  # Plot kernel density estimate.
  lines(x1, dssden(fMple, x1), lwd = 1, col = 'red')  # Plot maximum penalized likelihood estimate
  legend(2,0.165,c("True density function","Kernel density","Maximum penalized likelihood"), 
         col = c("black","green","red"), lty=1, cex = 0.3, box.lty=0)
  
  # Get squared error functions.
  seKde <- function(x) (fKde(x)-dnorm(x,2,3))^2
  seMple <- function(x) (dssden(fMple, x)-dnorm(x,2,3))^2
  # Get the integrated squared errors over the domain based on the simulated data.
  iseKde <- integrate(seKde,lower=low,upper=upp)$value
  iseMple <- integrate(seMple,lower=low,upper=upp)$value
  sprintf('ise for kernel density estimate = %.7f, ise for maximum penalized likelihood estimate = %.7f', iseKde, iseMple)
}

simu_exp <- function(n){
  # Given a data from Exponential(2), plot true density and estimated densities
  # using kernel density and maximum penalized likelihood and print out integrated
  # squared error for both methods.
  
  x <- rexp(n,2)  # Simulated data from exponential(2)
  low <- min(x)
  upp <- max(x)
  
  # Get Kernel density estimate, with bandwidth as ucv and kernel as gamma, since 
  # the data are asymmetric, and the support of gamma kernel covers the range of simulated data.
  fKde <- kdensity(x,bw='ucv',kernel='gamma',support=c(low,upp))
  # Get maximum penalized likelihood estimate. Theoretically, the true domain for 
  # the density is unknown, so set the domain based on the data simulated.
  fMple <- ssden(~x,domain=data.frame(x=c(low,upp)))
  
  # Plot the true density and estimated densities, the domain is set based on the
  # data we simulated as well.
  x1 <- seq(low,upp,len=500)
  plot(x1, dexp(x1, 2),type="l", lwd=2, ylim=c(0,2.5), xlim=c(low,upp), ylab = 'p.d.f',
       xlab = 'x from Exponential(2)', bty='n')  # Plot the true density.
  lines(x1, fKde(x1), lwd = 1, col = 'green')  # Plot kernel density estimate.
  lines(x1, dssden(fMple,x1), lwd = 1, col = 'red')  # Plot maximum penalized likelihood estimate.
  legend(1.5,2,c("True density function","Kernel density","Maximum penalized likelihood"), 
         col = c("black","green","red"), lty=1, cex = 0.3, box.lty=0)
  
  # Get squared error functions.
  seKde <- function(x) (fKde(x)-dexp(x,2))^2
  seMple <- function(x) (dssden(fMple, x)-dexp(x,2))^2
  # Get the integrated squared errors over the domain based on the simulated data.
  iseKde <- integrate(seKde,lower=low,upper=upp)$value
  iseMple <- integrate(seMple,lower=low,upper=upp)$value
  sprintf('ise for kernel density estimate = %.7f, ise for maximum penalized likelihood estimate = %.7f', iseKde, iseMple)
}

simu_beta <- function(n){
  # Given a data from Beta(3,1), plot true density and estimated densities
  # using kernel density and maximum penalized likelihood and print out integrated
  # squared error for both methods.
  
  x <- rbeta(n,3,1)  # Simulated data from beta(3,1).
  low <- min(x)
  upp <- max(x)
  
  # Get Kernel density estimate, with bandwidth as ucv and kernel as beta, since 
  # the data are asymmetric, and the support of beta kernel covers the range of simulated data.
  fKde <- kdensity(x,bw='ucv',kernel='beta',support=c(low,upp))
  # Get maximum penalized likelihood estimate. Theoretically, the true domain for
  # the density is unknown, so set the domain based on the data simulated.
  fMple <- ssden(~x,domain=data.frame(x=c(low,upp)))
  
  # Plot the true density and estimated densities, the domain is set based on the
  # data we simulated as well.
  x1 <- seq(low,upp,len=500)
  plot(x1, dbeta(x1,3,1),type="l", lwd=2, ylim=c(0,4), xlim=c(low,upp), ylab = 'p.d.f',
       xlab = 'x from Beta(3,1)', bty='n')  # Plot the true density.
  lines(x1, fKde(x1), lwd = 1, col = 'green')  # Plot kernel density estimate.
  lines(x1, dssden(fMple, x1), lwd = 1, col = 'red')   # Plot maximum penalized likelihood estimate.
  legend(0.2,3.5,c("True density function","Kernel density","Maximum penalized likelihood"), 
         col = c("black","green","red"), lty=1, cex = 0.3, box.lty=0)
  
  # Get squared error functions.
  seKde <- function(x) (fKde(x)-dbeta(x,3,1))^2
  seMple <- function(x) (dssden(fMple, x)-dbeta(x,3,1))^2
  # Get the integrated squared errors over the domain based on the simulated data.
  iseKde <- integrate(seKde,lower=low,upper=upp)$value  
  iseMple <- integrate(seMple,lower=low,upper=upp)$value
  sprintf('ise for kernel density estimate = %.7f, ise for maximum penalized likelihood estimate = %.7f', iseKde, iseMple)
}

set.seed(1)
# Plot and print simulated results with 250 sample size.
simu_norm(250)
simu_exp(250)
simu_beta(250)


# Plot and print simulated results with 500 csample size.
simu_norm(500)
simu_exp(500)
simu_beta(500)


# Plot and print simulated results with 1000 sample size.
simu_norm(1000)
simu_exp(1000)
simu_beta(1000)






