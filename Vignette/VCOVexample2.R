#-----------------------------------------------------------------
##################
# Setup
################## 
n <- 100

e <- runif(n)
rho <- .6
u <- sapply( seq(e), function(i) rho*e[i] )
v <- runif(n) + c( 0, head(u, -1) )

x <- 1:n
y <- 2*x + v
#-----------------------------------------------------------------
##################
# Regression and Output
################## 
reg <- lm(y~x)
X <- cbind(1,x)
r <- reg$residuals
#-----------------------------------------------------------------
##################
# Compare Covariance Matrices
################## 
# lag (m=2)

# Newey-West Sandwich Style
V0 <- sandwich::NeweyWest(reg, sandwich=F, prewhite=F, lag=2, adjust=F)

# Sparse Rolling Style
library(Matrix)
w <- as( toeplitz(c(1,2/3,1/3, rep(0,n-3))), "dgCMatrix")
V1 <- 1/n * t(X) %*% t( r*t(w*r) ) %*% X

# Differences < epsilon
(V0 - V1) < 1E3*.Machine$double.eps
#            (Intercept)    x
#(Intercept)        TRUE TRUE
#x                  TRUE TRUE


#-----------------------------------------------------------------
##################
# Compare Covariance Matrices
################## 
# lag (m=1)

# Newey-West Sandwich Style
V2 <- sandwich::NeweyWest(reg, sandwich=F, prewhite=F, lag=1, adjust=F)

# Sparse Rolling Style
library(Matrix)
w <- as( toeplitz(c(1,.5, rep(0,n-2))), "dgCMatrix")
V3 <- 1/n * t(X) %*% t( r*t(w*r) ) %*% X

# Differences < epsilon
(V2 - V3) < 1E3*.Machine$double.eps
#            (Intercept)    x
#(Intercept)        TRUE TRUE
#x                  TRUE TRUE
