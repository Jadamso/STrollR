```r
    #--------------------------------------------------------------
    # Setup
    n <- 100
    e <- runif(n)
    u <- sapply( seq(e), function(i) .6*e[i] )
    v <- runif(100) + c( 0, head(u, -1) )
    x <- 1:n
    y <- 2*x + v

    #--------------------------------------------------------------
    # Regression and Output
    reg <- lm(y~x)
    X <- cbind(1,x)
    r <- reg$residuals
    
    #--------------------------------------------------------------
    # Covariance Matrix `Meat`
    # Newey-West Sandwich Style
    V0 <- sandwich::NeweyWest(reg, sandwich=F, prewhite=F, lag=2, adjust=F)
    # Sparse Rolling Style
    library(Matrix)
    w <- as( toeplitz(c(1,2/3,1/3, rep(0,n-3))), "dgCMatrix")
    V2 <- 1/n * t(X) %*% t( r*t(w*r) ) %*% X

    #--------------------------------------------------------------
    # Compare
    (V0 - V2) < 1E3*.Machine$double.eps
    #            (Intercept)    x
    #(Intercept)        TRUE TRUE
    #x                  TRUE TRUE

```
