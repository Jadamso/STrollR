#------------------------------------------------------------------
##################
# Load the meat
##################
#sedir <- path.expand("~/Desktop/Common/R_Code/")
#source( paste0(sedir, "vcovST_meat.R") )
message("vcovSTsep")
#------------------------------------------------------------------
##################
#' vcovSTsep
################## 
#' @param DAT list of regression objects from vcov.format
#' @param LISTJ list of weighting objects from listj
#' @param t_double_count double count time?
#' @param verbose print output
#' @param return_each for debugging, return only spatial and only temporal covariances
#' @return covariance matrix w/ SHAC correction
#  @examples
#' @export


vcovSTsep <- compiler::cmpfun( function(
    DAT,
    LISTJ=NA,
    t_double_count=FALSE,
    verbose=FALSE,
    return_each=FALSE) {

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


    ## Temporal Matrix Data
    tmat <- LISTJ$tmat
    tid  <- LISTJ$tid
    t_cutoff <- LISTJ$t_cutoff

    ##############################################
    # Correct for spatial correlation:

    if( names(LISTJ)[1]!="wmat") { 
        message("Error in listj")
        break
    }

    if(verbose){message("Using Jordan's Methods...")}

    mo           <- order(dat$time)
    Is           <- Diagonal(length(unique(dat$time)))
    Sigma_s      <- kronecker(Is, wmat)
    XeeX_spatial <- XOmegaX0(X[mo,], Sigma_s, E[mo])/n

    if(verbose) {message("Computed Spatial VCOV.")}

    ##############################################
    # Correct for temporal correlation:
    # setkey(dat, time)

    ## all(diff(mo)==1)
    It           <- Diagonal(length(unique(dat$unit)))
    ## Don't Double Count Space and Time Diagonal
    tmat0        <- tmat - diag(nrow(tmat))
    Sigma_t      <- kronecker(It, tmat0)
    XeeX_serial  <- XOmegaX0(X, Sigma_t, E)/n

    if(verbose) {message("Computed Serial VCOV.")}

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

