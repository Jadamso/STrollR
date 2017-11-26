
#---------------------------------------------------------------------------
## Cluster Standard Errors

message("see multiwaycov::cluster.vcov")
# from http://www.r-bloggers.com/easy-clustered-standard-errors-in-r/
# also at https://thetarzan.wordpress.com/2011/05/28/heteroskedasticity-robust-and-clustered-standard-errors-in-r/
vcovCL <-function(model, cluster){
  require(sandwich, quietly = TRUE)
  require(lmtest, quietly = TRUE)
  
  #calculate degree of freedom adjustment
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- model$rank
  dfc <- (M/(M-1))*((N-1)/(N-K))
  
  #calculate the uj's
  uj  <- apply(estfun(model),2, function(x) tapply(x, cluster, sum))
  
  #use sandwich to get the var-covar matrix
  vcovCL <- dfc*sandwich(model, meat=crossprod(uj)/N)
  #cft    <- coeftest(model, vcovCL) 
  return(vcovCL)
}
vcovCL <- cmpfun(vcovCL)



## http://www.ne.su.se/polopoly_fs/1.216115.1426234213!/menu/standard/file/clustering1.pdf
clx <- function(fm, dfcw=1, cluster){
    M <- length(unique(cluster))
    N <- length(cluster)
    K <- fm$rank
    dfc <- (M/(M-1))*((N-1)/(N-K))

    if(class(fm)%in%c("lm","glm") ){
        ex  <- estfun(fm)
        u   <- apply(ex,2, function(x) tapply(x, cluster, sum))
        vcovCL <- dfc*sandwich(fm, meat=crossprod(u)/N)*dfcw

    } else {
        ex <- c(fm$residuals) * fm$cX
        b  <- solve( t(fm$cX)%*%fm$cX )

        u  <- apply(ex,2, function(x) tapply(x, cluster, sum))
        m  <- crossprod(u)/N
        #meat1 <- t(X) %*% diag(diag(crossprod(t( fm$residuals )))) %*% X
        vcovCL <- dfc*(b%*%m%*%b)*dfcw

    }

    coeftest(fm, vcovCL) 
}
clx <- cmpfun(clx)


mclx <- function(fm, dfcw, cluster1, cluster2){
    library(sandwich)
    library(lmtest)
    cluster12 = paste(cluster1,cluster2, sep="")
    M1  <- length(unique(cluster1))
    M2  <- length(unique(cluster2))
    M12 <- length(unique(cluster12))
    N   <- length(cluster1)
    K   <- fm$rank
    dfc1  <- (M1/(M1-1))*((N-1)/(N-K))
    dfc2  <- (M2/(M2-1))*((N-1)/(N-K))
    dfc12 <- (M12/(M12-1))*((N-1)/(N-K))
    u1   <- apply(estfun(fm), 2,
                  function(x) tapply(x, cluster1,  sum))
    u2   <- apply(estfun(fm), 2,
                  function(x) tapply(x, cluster2,  sum))
    u12  <- apply(estfun(fm), 2,
                  function(x) tapply(x, cluster12, sum))
    vc1   <-  dfc1*sandwich(fm, meat=crossprod(u1)/N )
    vc2   <-  dfc2*sandwich(fm, meat=crossprod(u2)/N )
    vc12  <- dfc12*sandwich(fm, meat=crossprod(u12)/N)
    vcovMCL <- (vc1 + vc2 - vc12)*dfcw
    coeftest(fm, vcovMCL)
}
mclx <- cmpfun(mclx)



#---------------------------------------------------------------------------
## Heteroskedasticity-robust standard error calculation.

summaryw <- function(model) {
    s <- summary(model)
    X <- model.matrix(model)
    u2 <- residuals(model)^2
    XDX <- 0
     
    ## Here one needs to calculate X'DX. But due to the fact that
    ## D is huge (NxN), it is better to do it with a cycle.
    for(i in 1:nrow(X)) {
    XDX <- XDX + u2[i]*X[i,]%*%t(X[i,])
    }
     
    # inverse(X'X)
    XX1 <- solve(t(X)%*%X)
     
    # Variance calculation (Bread x meat x Bread)
    varcovar <- XX1 %*% XDX %*% XX1
     
    # degrees of freedom adjustment
    dfc <- sqrt(nrow(X))/sqrt(nrow(X)-ncol(X))
     
    # Standard errors of the coefficient estimates are the
    # square roots of the diagonal elements
    stdh <- dfc*sqrt(diag(varcovar))
     
    t <- model$coefficients/stdh
    p <- 2*pnorm(-abs(t))
    results <- cbind(model$coefficients, stdh, t, p)
    dimnames(results) <- dimnames(s$coefficients)
    results
}


##sphet:::Omega
