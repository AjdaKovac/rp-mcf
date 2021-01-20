rm(list=ls())
library(vars)
library(tseries)
library(urca)
install.packages("ggplot2")
library(ggplot2)
library(gsl)

#setwd("~/Desktop/WU/Master/3_WS2020:21/Data Science and Machine Learning/Github/rp-mcf")
setwd("~/Documents/Github/rp-mcf/data")   # Ajda /Users/ajdakovac/Documents ...

library(readxl)
data1 <- read_excel("data.xlsx")
View(data1)  

##### Stationarity


y <- ts(data1$`Real GDP growth`, start = c(2004,1), frequency = 4)
plot.ts(y, col="red")

i <- ts(data1$Inflation, start=c(2004,1), frequency=4)
plot.ts(i, col="red")

sr <- ts(data1$`Shadow rate`, start = c(2004,1), frequency = 4)
plot.ts(sr, col="red")

lmargins <- ts(data1$`Lending margins`, start = c(2004,1), frequency = 4)
plot.ts(lmargins, col="red")

lstandards <- ts(data1$`Lending standards`, start = c(2004,1), frequency = 4)
plot.ts(lstandards, col="red")

ltl <- ts(data1$`Long term loans`, start = c(2004,1), frequency = 4)
plot.ts(ltl, col="red")

plot(cbind(y, lmargins, lstandards, i, sr, ltl), main = "Time Series Data")


################################ Dickey Fuller Tests for unit roots (different specifications)

library(stargazer)
### Y

y_norm <- adf.test(y) ### non-stationary
y_norm
pp.test(y)
y_trend <- ur.df(y, selectlags = "AIC", type = "trend")
summary(y_trend) ### has unit root => non-stationary
y_drift <- ur.df(y, selectlags = "AIC", type = "drift")
summary(y_drift) ### has unit root => non-stationary

### I
adf.test(i) ### non-stationary
pp.test(i) ### non-stationary
i_trend <- ur.df(i, selectlags = "AIC", type = "trend")
summary(i_trend) ### has unit root => non-stationary
i_drift <- ur.df(i, selectlags = "AIC", type = "drift")
summary(i_drift) ### non-stationary

### SR
adf.test(sr) ### non-stationary
pp.test(sr)  ### non-stationary
sr_trend <- ur.df(sr, selectlags = "AIC", type = "trend")
summary(sr_trend) ### has unit root => non-stationary
sr_drift <- ur.df(sr, selectlags = "AIC", type = "drift")
summary(sr_drift) ### has unit root => non-stationary
sr_none <- ur.df(sr, selectlags = "AIC", type = "none")
summary(sr_none) ### has unit root => non-stationary

### lmargins
adf.test(lmargins) ### non-stationary
pp.test(lmargins) ### non-stationary
lmargins_trend <- ur.df(lmargins, selectlags = "AIC", type = "trend")
summary(lmargins_trend) ### has unit root => non-stationary
lmargins_drift <- ur.df(lmargins, selectlags = "AIC", type = "drift")
summary(lmargins_drift) ### has unit root => non-stationary
lmargins_none <- ur.df(lmargins, selectlags = "AIC", type = "none")
summary(lmargins_none)  ### has unit root => non-stationary

### lstandards
adf.test(lstandards) ### non-stationary
pp.test(lstandards) ### non-stationary
lstandards_trend <- ur.df(lstandards, selectlags = "AIC", type = "trend")
summary(lstandards_trend) ### has unit root => non-stationary
lstandards_drift <- ur.df(lstandards, selectlags = "AIC", type = "drift")
summary(lstandards_drift) ### has unit root => non-stationary
lstandards_none <- ur.df(lstandards, selectlags = "AIC", type = "none")
summary(lstandards_none) ### has unit root => non-stationary

### long term loans
adf.test(ltl) ### non-stationary
pp.test(ltl) ### non-stationary
ltl_trend <- ur.df(ltl, selectlags = "AIC", type = "trend")
summary(ltl_trend) ### has unit root => non-stationary
ltl_drift <- ur.df(ltl, selectlags = "AIC", type = "drift")
summary(ltl_drift) ### has unit root => non-stationary
ltl_none <- ur.df(ltl, selectlags = "AIC", type = "none")
summary(ltl_none) ### has unit root => non-stationary

########################## First Differencing to obtain stationary data

### take first differences
deltay <- diff(y)
deltai <- diff(i)
deltasr <- diff(sr)
deltalmargins <- diff(lmargins)
deltalstandards <- diff(lstandards)
deltaltl <- diff(ltl)

### Y => first diff is stationary
adf.test(deltay) ## => stationary now => I(1)
pp.test(deltay) ## => stationary now => I(1)

### I
adf.test(deltai) ## => stationary now => I(1)

### SR
adf.test(deltasr) ## => stationary now at 5% => I(1) 

### lmargins

adf.test(deltalmargins) ## => stationary now at 5% => I(1) 

### lstandards

adf.test(deltalstandards) # => stationary now at 10% => I(1) 

### ltl

adf.test(deltaltl) # => stationary now at 10% => I(1) 


########################## Define vector x: specify order according to recursive identification scheme

#x<-cbind(deltay, deltai, deltalstandards, deltasr, deltalmargins) # is a covariance stationary process
x<-cbind(deltay, deltai, deltasr, deltalstandards, deltalmargins, deltaltl)
plot.ts(x, col="blue", main = "Covariance stationary vector x")
#x <- base::scale(x, center = TRUE, scale = TRUE)


### probably not needed
#C <- matrix(NA, nrow=6, ncol=6, dimnames=list(c("Y", "i", "p", "G", "U", "C"), c("Y", "i", "p", "G", "U", "C")))

#C[1,2] <- 0
#C[1,3] <- 0
#C[1,4] <- 0
#C


############################################### VAR-system

### select lags
x1 <- VARselect(x, lag.max = 16, type = "const")
VARselect(x) 
x1$selection # p=8 lags according to AIC as well as HQ   # upon adding ltl, p=7


VAR_1 <- VAR(x, p = 2, type = "trend", season = NULL, exog = NULL) #VAR
as.matrix(Bcoef(VAR_1))
VAR_1$varresult
summary(VAR_1) # this is for the whole timeperiod

plot.ts(resid(VAR_1))


### test for stability:
max(roots(VAR_1)) # max eigenvalue < 1 which implies that our system is stable and stationary!

normality.test(VAR_1)


x.serial <- serial.test(VAR_1, lags.pt = 12, type = "PT.asymptotic")
x.serial

### Cholesky decomposition of the Variance covariance matrix to get matrix B
#install.packages("svars")
library(svars)
B <- id.chol(VAR_1, order_k = c(1, 2, 3, 4, 5, 6)) 
B1 <- summary(B)
IRF <- irf(B, impulse = "deltasr", response= "deltag", n.ahead = 8, boot = TRUE, ortho = TRUE)
plot(IRF)
stargazer(B1)