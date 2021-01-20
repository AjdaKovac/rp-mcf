##################################################################################################
##################################################################################################
# Lecture Example Vector Autoregression models.
##################################################################################################
##################################################################################################

shift<-function(x,shift_by){
  stopifnot(is.numeric(shift_by))
  stopifnot(is.numeric(x))
  
  if (length(shift_by)>1)
    return(sapply(shift_by,shift, x=x))
  
  out<-NULL
  abs_shift_by=abs(shift_by)
  if (shift_by > 0 )
    out<-c(tail(x,-abs_shift_by),rep(NA,abs_shift_by))
  else if (shift_by < 0 )
    out<-c(rep(NA,abs_shift_by), head(x,-abs_shift_by))
  else
    out<-x
  out
}

# Load some required packages.
is.installed <- function(mypkg){
  is.element(mypkg, installed.packages()[,1])
} 

# check if package "forecast" is installed
if (!is.installed("openxlsx")){
  install.packages("openxlsx")
}

if (!is.installed("vars")){
  install.packages("vars")
}

require(vars)
require(openxlsx)
require(texreg)
require(xtable)

# Read in Data from xlsx-file "currencies.xls"
input.data.path <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/")
input.data <- "currencies.xlsx"
currencies <- read.xlsx(paste(input.data.path, input.data, sep=""), detectDates = TRUE)

# Generate lags in R. Use Shift Function above, the lag function in R does not work  for data.frames.
currencies$EUR_lag1 <- shift(currencies$EUR, -1)
currencies$GBP_lag1 <- shift(currencies$GBP, -1)
currencies$JPY_lag1 <- shift(currencies$JPY, -1)

currencies$EUR_lag2 <- shift(currencies$EUR, -2)

# Construct log return time series 
currencies$REUR <- ( log(currencies$EUR) - log(currencies$EUR_lag1) ) * 100
currencies$RGBP <- ( log(currencies$GBP) - log(currencies$GBP_lag1) ) * 100  
currencies$RJPY <- ( log(currencies$JPY) - log(currencies$JPY_lag1) ) * 100 

currencies$REUR_lag1 <- shift(currencies$REUR, -1)
                     
# Now calculate the VAR model:
# restricted to the following dates: 2002-07-10 to 2007-07-07. Osveration 4 to 1827.
var_eur_gbp_jpy <-  VAR(currencies[4:1827,c("REUR", "RGBP", "RJPY")], 
                        p = 1, type = c("const"),
                        season = NULL, 
                        exogen = NULL, 
                        lag.max = NULL)

var_eur_gbp_jpy_level <-  VAR(currencies[4:1827,c("EUR", "GBP", "JPY")], 
                        p = 1, type = c("const"),
                        season = NULL, 
                        exogen = NULL, 
                        lag.max = NULL)

# VAR model for the full sample:
var_eur_gbp_jpy_full <-  VAR(na.exclude(currencies[,c("REUR", "RGBP", "RJPY")]), p = 2, type = c("const"),
                        season = NULL, exogen = NULL, lag.max = NULL)

##
## Build growth rates for Canada e, prod and rw:
data("Canada")
prod_growth <- Canada[,"prod"]

prod_growth <-  diff( Canada[,"prod"], lag=4)/ lag( Canada[,"prod"], k=-4)
e_growth <- diff( Canada[,"e"], lag=4)/ lag( Canada[,"e"], k=-4)
u_growth <- diff( Canada[,"U"], lag=4)/ lag( Canada[,"U"], k=-4)
rw_growth <- diff( Canada[,"rw"], lag=4)/ lag( Canada[,"rw"], k=-4)

Canada_growth <- cbind(prod_growth,e_growth,u_growth, rw_growth)

test2v <- summary(VAR(Canada_growth, p = 1, type = "const"))                     

# why is the rw_growth equation so stange
                     
Canada_adjusted <- na.exclude(cbind(prod_growth,e_growth,Canada[,"U"]/100, rw_growth))

test3v <- summary(VAR(Canada_adjusted, p = 1, type = "const")) 



