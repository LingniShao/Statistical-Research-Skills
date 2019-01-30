install.packages("COMPoissonReg")
library(COMPoissonReg)

da0 <- read.dta("http://www.stata-press.com/data/lf2/couart2.dta")
da <- data.frame(art = da0$art,phd = da0$phd)
View(da)

# Fit Poisson Regression Model
fit_poi <- glm(art ~ phd, data=da,
            family=poisson, na.action=na.exclude)
summary(fit_poi)

# Fit COM-Poisson Regression Model
fit_cmp <- glm.cmp(art ~ phd, data=da)
print(fit_cmp)


# Compute constant COM-Poisson and Poisson deviances
dev <- sum(deviance(fit_cmp))
dev1 <- sum(residuals.glm(fit_poi,type = 'deviance')^2)
dev;dev1

# Compute fitted values
y.hat <- predict(fit_cmp, newdata=da)

#Plot data and estimated model
plot(art~phd,da)
lines(da$phd,y.hat,col = 'blue')
lines(da$phd,predict(fit_poi,newdata = da),col = 'red')
legend('topright',inset = 0.05, c("COM-Poisson Regression", "Poisson Regression"),
       col=c('blue', 'red'),lty=1,cex = 0.5)



