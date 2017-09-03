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
