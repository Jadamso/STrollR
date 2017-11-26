
#---------------------------------------------------------------------------
##################
#    CHUNK F     #
################## 

vcovSTkron.sep <- function(
    DT,
    LISTJ=NA,
    cores=1,
    tcores=NA,
    verbose=FALSE) {


    ## Input Data
    dt    <- DT$dt
    Xvars <- DT$Xvars
    #TSUB  <- DT$TSUB
    #IDSUB <- DT$IDSUB

    timeUnique  <- unique(dt[, time])
    panelUnique <- unique(dt[, unit])

    ##############################################
    # Prepare Covariance Matrix Inputs

    # Parameters
    n <- nrow(dt)
    k <- length(Xvars)

    # Bread Matrix
    X     <- as.matrix(dt[, eval(Xvars), with=FALSE])
    invXX <- solve(t(X) %*% X) * n

    # Residuals
    E <- dt$e

    ## Distance Matrix Data
    wmat <- LISTJ$wmat
    wid  <- LISTJ$wid
    d_cutoff <- LISTJ$d_cutoff


    ## Distance Matrix Data
    tmat <- LISTJ$tmat
    tid  <- LISTJ$tid
    t_cutoff <- LISTJ$t_cutoff

    ## Must setkey(dt, unit, time) in vcovST.format
    mt      <- order(dt$unit)
    ms      <- order(dt$time)

    ## Temporal Covariance
    It      <- Diagonal(length( panelUnique ))
    tmat0   <- tmat - diag(nrow(tmat))
    Sigma_t <- kronecker(It, tmat0)

    ## Spatial Covariance
    Is      <- Diagonal(length( timeUnique ))
    Sigma_s <- kronecker(Is, wmat)
    Sigma_s <- Sigma_s[order(ms),order(ms)]

    ## Spatial and Temporal Meat Matrix
    Sigma_st <- Sigma_s + Sigma_t
    XeeX     <- XOmegaX0(X, Sigma_st, E)/n

    ## Spatial and Temporal Sandwich Matrix
    V_spatial_HAC <- invXX %*% XeeX %*% invXX /n
    V_spatial_HAC <- (V_spatial_HAC + t(V_spatial_HAC)) / 2

    return(V_spatial_HAC)
}

vcovSTkron.sep <- cmpfun( vcovSTkron.sep )


#---------------------------------------------------------------------------
##################
#    CHUNK F     #
################## 

vcovSTkron <- function(
    DT,
    LISTJ=NA,
    cores=1,
    tcores=NA,
    verbose=FALSE) {


    ## Input Data
    dt    <- DT$dt
    Xvars <- DT$Xvars
    #TSUB  <- DT$TSUB
    #IDSUB <- DT$IDSUB

    timeUnique <- unique(dt[, time])
    panelUnique <- unique(dt[, unit])

    ##############################################
    # Prepare Covariance Matrix Inputs

    # Parameters
    n <- nrow(dt)
    k <- length(Xvars)

    # Bread Matrix
    X     <- as.matrix(dt[, eval(Xvars), with=FALSE])
    invXX <- solve(t(X) %*% X) * n

    # Residuals
    E <- dt$e

    ## Distance Matrix Data
    wmat <- LISTJ$wmat
    wid  <- LISTJ$wid
    d_cutoff <- LISTJ$d_cutoff


    ## Distance Matrix Data
    tmat <- LISTJ$tmat
    tid  <- LISTJ$tid
    t_cutoff <- LISTJ$t_cutoff

    ## Must setkey(dt, unit, time) in vcovST.format
    mt      <- order(dt$unit)
    ms      <- order(dt$time)

    ## Must setkey(dt, unit, time) in vcovST.format
    all_mat  <- kronecker(wmat, tmat)
    #all_mat2 <- all_mat[order(ms), order(ms)]
    # all_mat[1:15,1:15]

    XeeX    <- XOmegaX0(X, all_mat, dt$e) / n


    V_spatial_HAC <- invXX %*% XeeX %*% invXX /n
    V_spatial_HAC <- (V_spatial_HAC + t(V_spatial_HAC)) / 2

    return(V_spatial_HAC)
}

vcovSTkron <- cmpfun( vcovSTkron )


