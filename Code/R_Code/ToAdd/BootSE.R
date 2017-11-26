#---------------------------------------------------------------------------
## Bootstrapping
#library(boot)
#library(simpleboot)
#library(clusterSEs)  ## cluster.bs.glm
#library(multiwayvcov) ## cluster.boot

#---------------------------------------------------------------------------
## M out of N

#n <- nobs(DF)
#m <- sqrt(n)

# lim m log(log(n)) / n --> 0
# lim m / n --> 0

# Chernick & LaBudde - An Introduction to Bootstrap Methods with applications
# https://normaldeviate.wordpress.com/2013/01/27/bootstrapping-and-subsampling-part-ii/
#http://www.colorado.edu/Economics/seminars/SeminarArchive/2010-11/wilson.pdf
#On the choice of m in the m out of n bootstrap and confidence bounds for extrema

#---------------------------------------------------------------------------
## LOESS Standard Errors and BOOTSTRAP


## Bootstrap and Fit to X; library(simpleboot)
# Lo   <- loess(y ~ x, span= smoother, degree= deg)
# Lo.b <- simpleboot::loess.boot(Lo, R = 500, ngrid=100)

## Manual Bootstrapping
#boot.coef <- function(data, indices, ...){ #maxit=20
    #data <- data[indices,] # select obs. in bootstrap sample
    #predict(loess(LTRI10 ~ Start, data, span= smoother, degree= deg), se=FALSE)}
#BOOT <- boot(R_Battle, boot.coef, R=1000, parallel="multicore")

## Scatterplot bootstrap
#Lo.b <- scatboot(y, x, breps, .68, span=smoother, degree= deg)



#---------------------------------------------------------------------------
## SCATBOOT; Fit a loess curve to a scatterplot, and calculate a symmetric
## nonparametric bootstrap confidence region surrounding that curve.

#tryCatch( source("http://content.csbs.utah.edu/~rogers/datanal/R/scatboot.r"), 
#    error   = function(c) "scatboot failed!") #,warning = function(c) "warning", message = function(c) "message")


## Alan R. Rogers, 26 Feb 2011
# wrapper
#' Trim Whitespace from a String
#' @seealso \url{http://content.csbs.utah.edu/~rogers/datanal/R/scatboot.r}
#' @param x, y A vector of data
#' @return A list containing a confidence region and the input data/parameters
# @examples scatboot( )
#' @export
scatboot <- compiler::cmpfun( function(x, y, nreps=100, confidence=0.9,
                     degree=2, span=2/3, family="gaussian"){#

    # Put input data into a data frame, sorted by x, with no missing
    # values.
    dat <- na.omit(data.frame(x=x,y=y))
    if(nrow(dat) == 0) {
        print("Error: No data left after dropping NAs")
        print(dat)
        return(NULL)
    }
    ndx <- order(dat$x)
    dat$x <- dat$x[ndx]
    dat$y <- dat$y[ndx]#

    # Establish a series of x values for use in fits
    r <- range(dat$x, na.rm=T)
    x.out <- seq(r[1], r[2], length.out=40)

    # Fit curve to data
    f <- loess(y~x, data=dat, degree=degree, span=span,
               family=family)
    y.fit <- approx(f$x, fitted(f), x.out,rule=2)$y#

    # Generate bootstrap replicates
    len <- length(dat$x)
    mat <- matrix(0,nreps,length(x.out))
    for(i in seq(nreps)){
        ndx <- sample(len,replace=T)
        x.repl <- x[ndx]
        y.repl <- y[ndx]
        f <- loess(y.repl~x.repl, degree=degree,
                   span=span, family=family)
        mat[i,] <- predict(f, newdata=x.out)
    }#

    # Calculate limits and standard deviation
    n.na <- apply(is.na(mat), 2, sum)  # num NAs in each column
    nx <- ncol(mat)
    up.lim <- rep(NA, nx)
    low.lim <- rep(NA, nx)
    stddev <- rep(NA, nx)
    for(i in 1:nx) {
        if(n.na[i] > nreps*(1.0-confidence)) {
            # Too few good values to get estimate
            next
        }
        # effective confidence excluding NAs
        conf <- confidence*nreps/(nreps-n.na[i])
        pr <- 0.5*(1.0 - conf)
        up.lim[i] <- quantile(mat[,i], 1.0-pr, na.rm=T)
        low.lim[i] <- quantile(mat[,i], pr, na.rm=T)
        stddev[i] <- sd(mat[,i], na.rm=T)
    }
    ndx <- !is.na(up.lim)   # indices of good values#

    fit <- data.frame(x=x.out[ndx], y.fit=y.fit[ndx], up.lim=up.lim[ndx],
                      low.lim=low.lim[ndx], stddev=stddev[ndx])
    return(list(nreps=nreps, confidence=confidence,
                degree=degree, span=span,family=family,
                data=dat, fit=fit))
})

# Plot the object returned by scatboot.  Additional arguments are
# passed to xyplot.  Function returns a plot object.
#
# Alan R. Rogers, 26 Feb 2011
scatboot.plot <- function(sb,...) {
    require(lattice)
    require(grid)
    p <- xyplot(y~x, data=sb$data,
        panel=function(x,y,...) {
            panel.xyplot(x,y,...)
            panel.xyplot(sb$fit$x, sb$fit$y.fit, type="l",...)
            panel.xyplot(sb$fit$x, sb$fit$up.lim, type="l",
                col="gray",...)
            panel.xyplot(sb$fit$x, sb$fit$low.lim, type="l",
                col="gray",...)
            pg.x <- c(sb$fit$x, rev(sb$fit$x))
            pg.y <- c(sb$fit$up.lim, rev(sb$fit$low.lim))
            grid.polygon(pg.x, pg.y,
            gp=gpar(fill="pink",col="transparent",
                alpha=0.5),
                default.units="native")#
        }, ...)
    return(p)
}

# Test scatboot and scatboot.plot
#
# Alan R. Rogers, 26 Feb 2011
scatboot.test <- function(nreps=100, confidence=0.9,
                          span=2/3, degree=2,family="gaussian") {#

    x <- seq(0, 4, length.out=25)
    y <- sin(2*x)+ 0.5*x + rnorm(25, sd=0.5)
    sb <- scatboot(x, y, nreps=nreps, confidence=confidence,
        span=span, degree=degree, family=family)
    scatboot.plot(sb)
}

