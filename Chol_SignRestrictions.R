rm(list=ls())
require(bvarsv)
require(vars)
require(tseries)


# Some auxiliary function to make life easier

#function for building X matrix with lags
#X is Txm matrix
mlag <- function(X,lag)
{
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag)  
}



################## without packages from the scratch #######################

# dataset with inflation (inf), unemployment (une) and interest rate (tbi)
data(usmacro.update)
plot(usmacro.update)

#Lags
p <- 2
#Intercept
cons <- TRUE



# Construct variables for use in VAR with mlag-function
Yraw <- as.matrix(usmacro.update)
Y <- Yraw
X <- mlag(Y,p)

dim(X)
dim(Y)

#first p rows are zero
head(X)
Y <- Y[(p+1):nrow(Y),]
X <- X[(p+1):nrow(X),]

#for constant, ie if cons=TRUE
if(cons)
{
  X <- cbind(X,1)
}

M <- ncol(Y)
bigT <- nrow(Y)
K <- ncol(X)

#OLS estimator for the coeff. matrix
B.ols <- solve(crossprod(X))%*%crossprod(X,Y)
#or B.ols1<- lm(Y ~ X - 1)$coef

#OLS Sigma: Sum over T[(Yt - BXt)^2/T-K]
SSE <- crossprod(Y-X%*%B.ols)
SIGMA.ols <- SSE/(bigT-K)

#Fit of the model
yfit <- X%*%B.ols
resids <- Y-yfit


###----------------------------------- IRF Simulation --------------------------------------------------------

# Total simulation runs
ntot <- 1000 
set.seed(0)

#horizon: impact and next NHOR-1 quarters
nhor <- 13

# OLS variance of coefficients 
V.ols <- solve(crossprod(X))

# OLS mean
B.mean <- V.ols%*%crossprod(X,Y)

## Store objects
B.store <- array(0,c(ntot,K,M))

#IRFs -- dimensions: Number of simulations x Number of responses x Number of structural shocks x horizon
IRF.chol.store <- array(NA,c(ntot,M,M,nhor))        #for Cholesky
IRF.sign.store <- array(NA,c(ntot,M,M,nhor))        #for Sign restrictions

for(irep in 1:ntot)
{
  #Step 1: Assign OLS sigma
  SIGMA <- SIGMA.ols
  #Sigma (Kronecker) OLS Variance
  bigV <- kronecker(SIGMA, V.ols)
  
  #Cholesky decomposition of bigV
  cholV.ols <- t(chol(bigV))
  
  #Step 2: Construct a draw from the OLS coefficients with uncertainty (B.mean|SIGMA, Y from a Gaussian)
  B.sim <- as.vector(B.mean) + cholV.ols %*% rnorm(K*M) #vectorized
  B.sim <- matrix(B.sim,K,M) # reshape in matrix form 
  
  ##Calculate impulse responses
  #Compute companion matrix 
         B.comp <- matrix(0,M*p,M*p) 
  B.comp[1:M,] <- t(B.sim[1:(M*p),,drop=FALSE])   
  if(p > 1) 
    B.comp[(M+1):(M*p),1:(M*(p-1))] <- diag(M*(p-1))
  
  
  #Before computing IRFs, check stability --> real part of the eigenvalues should below one
  #If Eigenvalues > 1, transform data to something stationary
  if(max(Re(eigen(B.comp)$values))<1)
  {
    #Assign conditions for sign restrictions later on
    cond.overall <- 0
    counter <- 0
    
    #Identification via Choleski
    shock.chol <- t(chol(SIGMA))
    
    ##Normalise shocks to X units
    #unit <- 1
    #shock.chol <- shock.chol%*%diag(unit/diag(shock.chol))
    
    #Identification via Sign restrictions "while-loop" for sign restrictions
    
    #draw rotation matrices Q till you find a fitting one (while-loop)
    while(cond.overall == 0)
    {
      counter <- counter + 1
      
      #Define a rotation matrix with positive values on the main diagonal (ie. Householder Transformation, KL:p427)
         Stilde <- matrix(rnorm(M*M,0,1),M,M)
      qr.object <- qr(Stilde)
              Q <- qr.Q(qr.object)
              Q <- Q%*%diag((diag(Q)>0)-(diag(Q)<0))

      # #Define a rotation matrix using Givens rotation matrices (ie. Givens Rotation matrices, KL:p426)
      # th1   <- pi*runif(1)
      # th2   <- pi*runif(1)
      # th3   <- pi*runif(1)
      # qmat12 <- matrix(c(cos(th1), -sin(th1), 0,
      #                  sin(th1),  cos(th1), 0,
      #                  0      ,   0     , 1), nrow=3, byrow=T)
      # qmat13 <- matrix(c(cos(th2),   0,       -sin(th2),
      #                    0,   1,     0,
      #                    sin(th2),   0,      cos(th2) ), nrow=3, byrow=T)
      # qmat23 <- matrix(c(1,   0,       0,
      #                    0, cos(th3), -sin(th3),
      #                    0, sin(th3),  cos(th3) ), nrow=3, byrow=T)
      # qmat <- qmat12%*%qmat13%*%qmat23
      
      #shock is a full matrix
      shock.sign <- shock.chol%*%(Q)  #horizon=0

      # Price shock (e.g. Oil)
      cond.AS <- (shock.sign[1,1]>0)*(shock.sign[2,1]>0)*(shock.sign[3,1]<0)
      # Demand shock
      cond.AD <- (shock.sign[1,2]<0)*(shock.sign[2,2]>0)*(shock.sign[3,2]<0)
      # Monetary Policy shock
      cond.MP <- (shock.sign[1,3]<0)*(shock.sign[2,3]>0)*(shock.sign[3,3]>0)

      ## If you want to apply the restrictions not only for impulses on impact, 
      ## but also for a higher horizon, construct J-matrix
      ## and add conditions accordingly.
      #
      #       J <- matrix(0, M*p, M)
      # diag(J) <- 1
      # shock.sign1 <- t(J)%*%B.comp%*%J%*%shock.chol%*%Q #horizon=1
      #
      ## Demand shock effects in period t+1
      #
      # cond.AD1 <- (shock.sign1[1,2]<0)*(shock.sign1[2,2]>0)*(shock.sign1[3,2]<0)
      
      #Shocks have to be mutually exclusive (orthogonal)
      cond.overall <- cond.MP*cond.AS*cond.AD
      
      print(counter)
    }
    
    # Temporary objects for state space representation
    irf.mat.chol<- irf.mat.sign  <- array(0,c(M*p,M*p,nhor))
    
    # Impulse --> shock at t = 0:
    irf.mat.chol[1:M,1:M,1] <- shock.chol 
    irf.mat.sign[1:M,1:M,1] <- shock.sign
    
    #start at t = 1, as t = 0 is the impulse shock
    for(ihorz in 2:nhor)
    {
      irf.mat.chol[,,ihorz] <- B.comp%*%irf.mat.chol[,,ihorz-1]
      irf.mat.sign[,,ihorz] <- B.comp%*%irf.mat.sign[,,ihorz-1]
      
    }
    
    #Now Kp x Kp x nhor -> drop everything greater than M
    #clean up, as lagged values are not necessary
    irf.mat.chol <- irf.mat.chol[1:M,1:M,]
    irf.mat.sign <- irf.mat.sign[1:M,1:M,]
    
    #Store the IRFs
    IRF.chol.store[irep,,,] <- irf.mat.chol
    IRF.sign.store[irep,,,] <- irf.mat.sign
  }
}

#Quantiles over the first dimension ('confidence bounds')
#Under normality, the 16-th and 84-th percentiles correspond to the bounds of a one-SD confidence interval
IRF.chol.low <- apply(IRF.chol.store, c(2,3,4), quantile, 0.16,na.rm=TRUE) #16th percentile
IRF.chol.high <- apply(IRF.chol.store, c(2,3,4), quantile, 0.84,na.rm=TRUE)#84th percentile
IRF.chol.median <- apply(IRF.chol.store, c(2,3,4), median, na.rm=TRUE)

IRF.sign.low <- apply(IRF.sign.store, c(2,3,4), quantile, 0.16,na.rm=TRUE)
IRF.sign.high <- apply(IRF.sign.store, c(2,3,4), quantile, 0.84,na.rm=TRUE)
IRF.sign.median <- apply(IRF.sign.store, c(2,3,4), median, na.rm=TRUE)

#Start plotting the IRFs w.r.t the different shocks
#Cholesky identification
par(mfrow=c(3,3),mar=c(4,3,2,2))
for(ii in 1:M)
{
  for(jj in 1:M)
  {
    plot(IRF.chol.median[ii,jj,], ylab=paste("Response of",colnames(Y)[[ii]]), main=paste(colnames(Y)[[jj]],"(chol)"), ylim = c(min(IRF.chol.low[ii,jj,]) - 0.1, max(IRF.chol.high[ii,jj,]) + 0.1), type="l", xlab="Response in quarters") 
    lines(IRF.chol.low[ii,jj,], lty = 2)
    lines(IRF.chol.high[ii,jj,], lty = 2)
    abline(h=0,col="red")
  }
}

#Sign restrictions
par(mfrow=c(3,3),mar=c(4,3,2,2))
for(ii in 1:M)
{
  for(jj in 1:M)
  {
    plot(IRF.sign.median[ii,jj,], ylab=paste("Response of",colnames(Y)[[ii]]), main=paste(colnames(Y)[[jj]],"(sign)"), ylim = c(min(IRF.sign.low[ii,jj,]) - 0.1, max(IRF.sign.high[ii,jj,]) + 0.1), type="l", xlab="Response in quarters") 
    lines(IRF.sign.low[ii,jj,], lty = 2)
    lines(IRF.sign.high[ii,jj,], lty = 2)
    abline(h=0,col="red")
  }
}

# Plot the first 250 impulses that satisfy the assumed Signs
par(mfrow=c(3,3),mar=c(4,3,2,2))
for(ii in 1:M)
{
  for(jj in 1:M)
  {
    plot(IRF.sign.store[1,ii,jj,],ylab=colnames(Y)[[ii]],main=colnames(Y)[[jj]], ylim = c(min(na.omit(IRF.sign.store[,ii,jj,])), max(na.omit(IRF.sign.store[,ii,jj,]))), type = "l", col = "grey40", xlab="Response in quarters")
    for(kk in 2:250)
    {
      lines(IRF.sign.store[kk,ii,jj,], col = "grey50")
    }
    abline(h=0,col="red")
    lines(IRF.sign.median[ii,jj,], col = "green", lwd=2)
    lines(IRF.sign.low[ii,jj,], col = "blue", lty=2, lwd=2)
    lines(IRF.sign.high[ii,jj,], col = "blue", lty=2, lwd=2)
  }
}

