<img src="https://raw.githubusercontent.com/Jadamso/STrollR/master/STsim/STvarX.gif"  align="right" width="150" height="150">




# `S'patio `T'emporal `roll`ed standed errors in `R'

A computationally efficient way to calculate covariance matrices that are corrected for spatial and temporal correlation using a method I call *rolling*. The figure to the right is a visual representation of one such process \{ VonNeumann(1, rho=.3), AR(1, rho=.8) \}


Huge spatiotemporal covariance matrices can be calculated using sparse matrix approaches with spam and spam64. To calculate large sparse spatial weights matrices, use spam::rdist


Some functions have been extended to take advantage of parallel processing. The correction is applied in a post-estimation style. Some of the math is layed out in the appendix to the third chapter of my  [dissertation](https://sites.google.com/a/g.clemson.edu/ja-resources/research/Adamson2017_Thesis.pdf?attredirects=0). 



<!-- ![Alt Text](https://raw.githubusercontent.com/Jadamso/STrollR/master/STsim/STvarX.gif)
-->



## To install latest release from Github: 

    #install.packages("devtools")
    devtools::install_github("Jadamso/STrollR",subdir="STrollR")
    citation("STrollR")


<!--
knit(input="README.rmd", output="README.md")
-->

## Example of usage
```r
    #--------------------------------------------------------------
    # Setup
    DDFF <- fake_data_traditional(50)

    #--------------------------------------------------------------
    # Prepare Covariance Inputs
    ## Do this early,
    ## then recycle for all your regressions
    XEUC      <- 1000  #km
    d_cutoff  <- XEUC*(360/(6378.388*2*pi))
    coord_mat <- as.matrix(DDFF[ DDFF$Time==1, c("x","y")])
    WMAT <- weight_mat(coord_mat, cutoff=d_cutoff)
    ## N.B. distances are calculated from x,y coordinates
    ## can use lon,lat if change both coord_mat and weight_mat 

    #--------------------------------------------------------------
    # Regression
    ## Main Covariates
    ## Time, Group Fixed Effects
    ## 0 Instruments
    ## Temporally Clustered on ID
    form <- as.formula( Y~X+Z | Group + Time | 0 | ID)
    REG <- lfe::felm( form, data=DDFF, keepCX=TRUE)

    #--------------------------------------------------------------
    # Covariance Estimation
    VCOV <- vcovST.loop(
        reg=REG,
        unit="ID",
        time="Time",
        sp_coords=c("x","y"),
        wmat=WMAT,
        DF=DDFF[, c("ID", "Time", "x", "y")],
        cores=1
    )

    #--------------------------------------------------------------
    # Comparison Standard Errors
    sqrt( diag( VCOV / REG$clustervcv ) )
    ## 1.438368 1.413572
    ## spatial+clustered SE's > clustered SE's
```
## Example of the approach

Here is an example to display the logic. The example is only for temporal correlation.

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


## Caveat Emptor
This is a pre-release to the public and not gauranteed to work correctly with your code and data.


Recent Name Changes:

    vcovSTsep --> vcovST
    vcovSCL --> vcovST.loop
    

## Things To Do 

* small sample / asymptotic adjustment options

* support lm objects fully

felm is supported in most functions
currently lm is supported in just in a few functions

* vcovST.loop non-default behaviour

If not passing DF to vcovST.loop
vcovST.format Assumes the second Fixed Effect is temporal
must augment the original formula

* eliminate data.tables 

sub_dat <- as.data.frame(sub_dat) in vcovST_meat 
should be eliminated but I cannot compile data.tables correctly

* vcovSTkron

