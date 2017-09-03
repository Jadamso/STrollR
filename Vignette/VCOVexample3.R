
##################
# CHUNK D
################## 

# u = (I - pM)^-1  e  = W e
# E(uu') = 1/n Omega = 1/n W e e' W'

# W <- solve( diag( nrow(wmat) ) - rho * wmat )
# image.plot(W, col=rainbow(10))


# WMAT <- toeplitz(.2^(0:4)); e <- 1:5
## If spatial and temporal autocorrelation was constant over time, 
## V(e) = Sigma_{N\timesN} \kronecker Sigma_{T}, where
## Sigma_{N\timesN} is (NxN)^2 matrix of spatial autocorrelation
## Sigma_{T} is (T)^2 matrix of temporal autocorrelation


#x <- cbind(1:5, 1:5*2)
#e <- runif(5)
#K <- toeplitz( .4^(0:4))

#V <- K * e %*% t(e) 
#XVX <- t(x) %*% V %*% x
#XVX2 <- (t(x) %*% (K * diag(e) ) ) %*% ( diag(e)  %*% x )

