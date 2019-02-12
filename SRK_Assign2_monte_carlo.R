# Install the package for getting kernel density estimation.
# documentation: https://cran.r-project.org/web/packages/kdensity/kdensity.pdf
install.packages("kdensity")
# Install the package for getting maximum penalized likelihood estimation.
# documentation: https://cran.r-project.org/web/packages/gss/gss.pdf
install.packages("gss")
# Load the packages we need.
require(kdensity)  
require(gss)


simu_norm <- function(n){
  # Given a data from Normal(2,3), get the integrated squared error for kernel density
  # estimation and maximum penalized likelihood estimation.
  
  x <- rnorm(n,2,3)  # simulated data from Normal(2,3).
  low <- min(x)
  upp <- max(x)
  
  # Get Kernel density estimate, with bandwidth as nrd0 and kernel as gaussian, 
  # since the data are symmetric.
  fKde <- kdensity(x,bw='nrd0',kernel='gaussian',support=c(low,upp))
  # Get maximum penalized likelihood estimate. Theoretically, the true domain for 
  # the density is unknown, so set the domain based on the data simulated.
  fMple <- ssden(~x,domain=data.frame(x=c(low,upp)))
  
  # Get squared error functions.
  seKde <- function(x) (fKde(x)-dnorm(x,2,3))^2
  seMple <- function(x) (dssden(fMple, x)-dnorm(x,2,3))^2
  # Get the integrated squared errors over the domain based on the simulated data.
  iseKde <- integrate(seKde,lower=low,upper=upp)$value
  iseMple <- integrate(seMple,lower=low,upper=upp)$value
  
  return(c(iseKde,iseMple))
}

simu_exp <- function(n){
  # Given a data from Exponential(2), get the integrated squared error for kernel density
  # estimation and maximum penalized likelihood estimation.
  
  x <- rexp(n,2)  # Simulated data from exponential(2)
  low <- min(x)
  upp <- max(x)
  
  # Get Kernel density estimate, with bandwidth as ucv and kernel as gamma, since 
  # the data are asymmetric, and the support of gamma kernel covers the range of simulated data.
  fKde <- kdensity(x,bw='ucv',kernel='gamma',support=c(low,upp))
  # Get maximum penalized likelihood estimate. Theoretically, the true domain for 
  # the density is unknown, so set the domain based on the data simulated.
  fMple <- ssden(~x,domain=data.frame(x=c(low,upp)))
  
  # Get squared error functions.
  seKde <- function(x) (fKde(x)-dexp(x,2))^2
  seMple <- function(x) (dssden(fMple, x)-dexp(x,2))^2
  # Get the integrated squared errors over the domain based on the simulated data.
  iseKde <- integrate(seKde,lower=low,upper=upp)$value
  iseMple <- integrate(seMple,lower=low,upper=upp)$value
  
  return(c(iseKde,iseMple))
}


getci <- function(y){
  # the function for getting confidence interval using bootstrap.
  sims <- 10000
  n <- length(y)
  means <- numeric(sims)
  for (s in 1:sims) {
    bs <- sample(y,n,replace=TRUE)
    means[s] <- mean(bs)
  }
  means <- sort(means)
  upp <- round(means[0.95*sims],7)
  return(upp)
}

set.seed(1)
R = 50
# Get integrated squared errors for 250 sample size.
isenorm250 <- matrix(0,nrow = R, ncol = 2)
for (i in 1:R) {
  isenorm250[i,1] <- simu_norm(250)[1]
  isenorm250[i,2] <- simu_norm(250)[2]
}
iseexp250 <- matrix(0,nrow = R, ncol = 2)
for (i in 1:R) {
  iseexp250[i,1] <- simu_exp(250)[1]
  iseexp250[i,2] <- simu_exp(250)[2]
}
MCisenorm250 <- apply(isenorm250,2,mean)
MCiseexp250 <- apply(iseexp250,2,mean)


# Get integrated squared errors for 500 csample size.
isenorm500 <- matrix(0,nrow = R, ncol = 2)
for (i in 1:R) {
  isenorm500[i,1] <- simu_norm(500)[1]
  isenorm500[i,2] <- simu_norm(500)[2]
}
iseexp500 <- matrix(0,nrow = R, ncol = 2)
for (i in 1:R) {
  iseexp500[i,1] <- simu_exp(500)[1]
  iseexp500[i,2] <- simu_exp(500)[2]
}
MCisenorm500 <- apply(isenorm500,2,mean)
MCiseexp500 <- apply(iseexp500,2,mean)

# Get integrated squared errors for 1000 sample size.
isenorm1000 <- matrix(0,nrow = R, ncol = 2)
for (i in 1:R) {
  isenorm1000[i,1] <- simu_norm(1000)[1]
  isenorm1000[i,2] <- simu_norm(1000)[2]
}
iseexp1000 <- matrix(0,nrow = R, ncol = 2)
for (i in 1:R) {
  iseexp1000[i,1] <- simu_exp(1000)[1]
  iseexp1000[i,2] <- simu_exp(1000)[2]
}
MCisenorm1000 <- apply(isenorm1000,2,mean)
MCiseexp1000 <- apply(iseexp1000,2,mean)


## Get the crirical value of integrated squared errors using bootstrap.
cinormkde250 <- getci(isenorm250[,1])
cinormkde500 <- getci(isenorm500[,1])
cinormkde1000 <- getci(isenorm1000[,1])
cinormmple250 <- getci(isenorm250[,2])
cinormmple500 <- getci(isenorm500[,2])
cinormmple1000 <- getci(isenorm1000[,2])

ciexpkde250 <- getci(iseexp250[,1])
ciexpkde500 <- getci(iseexp500[,1])
ciexpkde1000 <- getci(iseexp1000[,1])
ciexpmple250 <- getci(iseexp250[,2])
ciexpmple500 <- getci(iseexp500[,2])
ciexpmple1000 <- getci(iseexp1000[,2])


## Plot box plot for integrated squared errors.
boxplot(isenorm250~col(isenorm250),names=c('KDE,250','MPLE,250'), col=c('deepskyblue','gold'))
boxplot(isenorm500~col(isenorm500),names=c('KDE,500','MPLE,500'), col=c('deepskyblue','gold'))
boxplot(isenorm1000~col(isenorm1000),names=c('KDE,1000','MPLE,1000'), col=c('deepskyblue','gold'))

boxplot(iseexp250~col(iseexp250),names=c('KDE,250','MPLE,250'), col=c('deepskyblue','gold'))
boxplot(iseexp500~col(iseexp500),names=c('KDE,500','MPLE,500'), col=c('deepskyblue','gold'))
boxplot(iseexp1000~col(iseexp1000),names=c('KDE,1000','MPLE,1000'), col=c('deepskyblue','gold'))


