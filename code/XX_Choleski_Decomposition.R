######################################################################
######################################################################
# VAR Models
######################################################################
######################################################################


#######################################################################
# Choleski decomposition in R
#######################################################################


m <- matrix(c(5,1,1,3),2,2)
m

tP <- chol(m)
tP

t(tP) %*% tP

t(solve(tP)) %*% m %*% solve(tP)

#######################################################################
# Set up a small VAR model
#######################################################################

pacman::p_load(vars)

Canada

var<- VAR(Canada[,c("e","rw")], p = 2, type = "const")
imp_irf <-irf(var, 
              boot=FALSE, 
              ortho=FALSE,
              n.ahead=8)
plot(imp_irf)

imp_oirf <- irf(var, 
                boot=FALSE, 
                ortho=TRUE,
                n.ahead=8)

######################################################################
# A second model:
######################################################################

pacman::p_load(vars)
pacman::p_load(texreg)
pacman::p_load(xtable)

# Download data
data <- read.table("http://www.jmulti.de/download/datasets/e1.dat", skip = 6, header = TRUE)

# Only use the first 76 observations so that there are 73 observations
# left for the estimated VAR(2) model after taking first differences.
data <- data[1:76, ]

# Convert to time series object
data <- ts(data, start = c(1960, 1), frequency = 4)

# Take logs and differences
data <- diff(log(data))

# Plot data
plot(data,  main = "Dataset E1 from Luetkepohl (2007)")

# Estimate model
model <- VAR(data, p = 2, type = "const")

# Look at summary statistics
model_summary <- summary(model)





model_summary$corres

feir <- irf(model, impulse = "income", response = "cons",
            n.ahead = 8, ortho = FALSE, runs = 1000)

plot(feir)

feir <- irf(model, 
            impulse = c("invest", "income", "cons"), 
            response = c("invest", "income", "cons"),
            n.ahead = 8, ortho = FALSE, runs = 1000)

plot(feir)

oir <- irf(model, 
           impulse = c("invest", "income", "cons"), 
           response = c("invest", "income", "cons"),
           n.ahead = 8, ortho = TRUE, runs = 1000)

t(chol(model_summary$corres))

plot(oir)

###################################################################
# Produce Latex Output:
###################################################################

Texreg_VAR <- list()

model.names <- colnames(data)

for (i0 in 1:length(model$varresult)){
 
  Texreg_VAR[[i0]] <-
    createTexreg(coef.names = as.character(rownames(model_summary$varresult[[i0]]$coefficients)),
                 coef = model_summary$varresult[[i0]]$coefficients[,1],
                 se = model_summary$varresult[[i0]]$coefficients[,2],
                 pvalues = model_summary$varresult[[i0]]$coefficients[,4],
                 ci.low = numeric(0),
                 ci.up = numeric(0),
                 gof.names = c("Adjusted R-Squared",
                               "No. of Obs."),
                 gof = c(round(model_summary$varresult[[i0]]$adj.r.squared, 2),
                         model_summary$varresult[[i0]]$df[2]+
                           model_summary$varresult[[i0]]$df[1]),
                 gof.decimal = c(TRUE, FALSE))
  
}

Table_BW_over_TA <-  texreg(Texreg_VAR,
                            dcolumn = TRUE, 
                            booktabs = TRUE,
                            custom.model.names = model.names, 
                            digits = 4,
                            use.packages = FALSE,
                            label = "table:BWoverTotalAssets",
                            caption = "Book value over total assets",
                            float.pos = "H")

Table_BW_over_TA

