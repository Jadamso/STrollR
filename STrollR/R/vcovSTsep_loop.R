#------------------------------------------------------------------
##################
# vcovSTsep, Parallel for Space and Time
##################
#source( paste0(sedir, "vcovST_meat.R") )
message("vcovSTsep_loop")

#------------------------------------------------------------------
##################
#' vcovSTsep with Parallel for Space and Time
################## 
#'
#' @param DAT list of regression objects from vcov.format
#' @param LISTJ list of weighting objects from listj
#' @param wmat weights matrix
#' @param t_double_count double count time?
#' @param verbose show messages
#' @param return_each for debugging, return only spatial and only temporal covariances
#' @param cores number of cores in spatial loop
#' @param tcores number of cores in temporal loop
#' @return covariance matrix w/ SHAC correction
#'
#  @examples
#'
#' @export


vcovSTsep.loop <- compiler::cmpfun( function(
    DAT,
    LISTJ=NA,
    wmat=LISTJ$wmat,
    t_double_count=FALSE,
    verbose=FALSE,
    return_each=FALSE, 
    cores=4,
    tcores=NA) {

    ## Input Data
    dat   <- DAT$dat
    Xvars <- DAT$Xvars
    #TSUB  <- DAT$TSUB
    #IDSUB <- DAT$IDSUB

    timeUnique  <- unique(dat[, time])
    panelUnique <- unique(dat[, unit])

    ##############################################
    # Prepare Covariance Matrix Inputs

    # Parameters
    n <- nrow(dat)
    k <- length(Xvars)

    # Bread Matrix
    X     <- as.matrix(dat[, eval(Xvars), with=FALSE])
    invXX <- solve(t(X) %*% X) * n

    # Residuals
    E <- dat$e


    ## Distance Matrix Data
    wmat <- LISTJ$wmat
    wid  <- LISTJ$wid
    d_cutoff <- LISTJ$d_cutoff
    TSUB     <- DAT$TSUB


    ## Distance Matrix Data
    tmat <- LISTJ$tmat
    tid  <- LISTJ$tid
    t_cutoff <- LISTJ$t_cutoff
    IDSUB    <- DAT$IDSUB


    ##############################################
    # Correct for spatial correlation:

    message("Loop over time periods...")

    if( names(LISTJ)[1]!="wmat") { 
        message("Error in listj")
        break
    }

    if(verbose){message("Using Jordan's Methods...")}

    XeeXhs <- parallel::mclapply(TSUB, function(tsub){
        iterateObsJSpatial(
            sub_dat=tsub,
            wmat=wmat,
            #wid=wid,
            Xvars=Xvars,
            verbose=verbose )
            }, mc.cores=cores)
    XeeX_spatial <- Reduce("+",  XeeXhs) / n

    if(verbose) {message("Computed Spatial VCOV.")}

    ##############################################
    # Correct for temporal correlation:
    # setkey(dat, time)

    if(is.na(tcores)) tcores <- cores
    message("Loop over units...")

    ## all(diff(mo)==1)
    ## Don't Double Count Space and Time Diagonal
    tmat0        <- tmat - diag(nrow(tmat))


    if(verbose) {message("Computed Serial VCOV.")}

    if(!t_double_count){
        ## Don't Double Count Space and Time
        diag(tmat) <- 0
        #cleanup(TMAT)  
    }

    XeeXhs <- parallel::mclapply(IDSUB, function(idsub) {
        iterateObsJTemporal(
            sub_dat=idsub,
            wmat=tmat,
            wid=tid,
            Xvars=Xvars,
            verbose=verbose)
        }, mc.cores=tcores
    )
    XeeX_serial <- Reduce("+",  XeeXhs) / n


    ##############################################
    # Correct for spatial and temporal correlation:
    XeeX <- (XeeX_spatial + XeeX_serial) 

    V_spatial_HAC <- invXX %*% XeeX %*% invXX /n
    V_spatial_HAC <- (V_spatial_HAC + t(V_spatial_HAC)) / 2

    if(return_each){

        # Generate VCE for only spatial correlation:
        V_spatial <- invXX %*% XeeX_spatial %*% invXX /n
        V_spatial <- (V_spatial + t(V_spatial)) / 2


        # Generate VCE for only temporal correlation:
        V_temporal <- invXX %*% XeeX_serial %*% invXX / n
        V_temporal <- (V_temporal + t(V_temporal)) / 2

        VCE <- list(
            V_temporal=V_temporal,
            V_spatial=V_spatial,
            V_spatial_HAC=V_spatial_HAC)

    } else {

        VCE <- V_spatial_HAC
    }


    return(VCE)
} )



#------------------------------------------------------------------
##################
#' vcovSpace with Parallel Approach
################## 
#'
#' @param DAT list of regression objects from vcov.format
#' @param LISTJ list of weighting objects from listj
#' @param wmat weights matrix
#' @param verbose show messages
#' @param cores number of cores in spatial loop
#' @return covariance matrix w/ spatial correction
#'
#  @examples
#' @export

vcovSpace.loop <- compiler::cmpfun( function(
    DAT,
    LISTJ=NA,
    wmat=LISTJ$wmat,
    verbose=FALSE,
    cores=4) {

    ## Input Data
    dat   <- DAT$dat
    Xvars <- DAT$Xvars

    ##############################################
    # Prepare Covariance Matrix Inputs

    # Parameters
    n <- nrow(dat)
    k <- length(Xvars)

    # Bread Matrix
    X     <- as.matrix(dat[, eval(Xvars), with=FALSE])
    invXX <- solve(t(X) %*% X) * n

    # Residuals
    E <- dat$e

    ## Distance Matrix Data
    TSUB <- DAT$TSUB


    ##############################################
    # Correct for spatial correlation:

    if(verbose){message("Loop over time periods...")}

    if( cores > 1){
        XeeXhs <- parallel::mclapply(TSUB, function(tsub){
            iterateObsJSpatial(
                sub_dat=tsub,
                wmat=wmat,
                #wid=wid,
                Xvars=Xvars,
                verbose=verbose )
                }, mc.cores=cores)
    } else {
        XeeXhs <- lapply(TSUB, function(tsub){
        iterateObsJSpatial(
            sub_dat=tsub,
            wmat=wmat,
            #wid=wid,
            Xvars=Xvars,
            verbose=verbose )
            })
    }

    XeeX_spatial <- Reduce("+",  XeeXhs) / n

    if(verbose) {message("Computed Spatial VCOV.")}

    # Generate VCE for only spatial correlation:
    V_spatial <- invXX %*% XeeX_spatial %*% invXX /n
    V_spatial <- (V_spatial + t(V_spatial)) / 2

    return(V_spatial)
} )


#------------------------------------------------------------------
##################
#' vcovSpace, Single Year Only
################## 
#'
#' @param DAT list of regression objects from vcov.format
#' @param LISTJ list of weighting objects from listj
#' @param wmat weights matrix
#' @param verbose show messages
#'
#' @return covariance matrix w/ spatial correction
#  @examples
#'
#' @export

vcovSpace.single <- compiler::cmpfun( function(
    DAT,
    LISTJ=NA,
    wmat=LISTJ$wmat,
    verbose=FALSE) {

    ## Input Data
    dat   <- DAT$dat
    Xvars <- DAT$Xvars

    ##############################################
    # Prepare Covariance Matrix Inputs

    # Parameters
    n <- nrow(dat)
    k <- length(Xvars)

    # Bread Matrix
    X     <- as.matrix(dat[, eval(Xvars), with=FALSE])
    invXX <- solve(t(X) %*% X) * n

    # Residuals
    E <- dat$e

    ##############################################
    # Correct for spatial correlation:

    if(verbose){message("Loop over time periods...")}

    XeeX_spatial <- XOmegaX0(X, wmat, E)


    if(verbose) {message("Computed Spatial VCOV.")}

    # Generate VCE for only spatial correlation:
    V_spatial <- invXX %*% XeeX_spatial %*% invXX /n
    V_spatial <- (V_spatial + t(V_spatial)) / 2

    return(V_spatial)
} )

