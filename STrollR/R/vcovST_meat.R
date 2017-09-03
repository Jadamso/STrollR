#########################
# Calculating Weights Matrices
#########################
message("vcovST_meat")

#------------------------------------------------------------------
##################
#' Main `Meat` Matrix Calculation
################## 
#'
#' @param X design matrix 
#' @param WMAT weighting matrix (preferably sparse sparse)
#' @param e vector of residuals
#' @return object to be used in vcov* functions
#'
#  @examples
#' @export

XOmegaX0 <- compiler::cmpfun( function(X, WMAT, e){

    #requireNamespace("Matrix")
    
    WEE    <- Matrix::t(e * Matrix::t(WMAT * e))

    XWeeX  <- t(X) %*% WEE %*% X

    return( as.matrix(XWeeX) )
} )


## Potential Speedup if WEE is symmetric
## WEE=(e * t(WMAT * e))


#------------------------------------------------------------------
##################
#' Wrapper for Matrix Calculation for Space
################## 
#' @param sub_dat a dataframe object for 1 time period
#' @param Xvars RHS design matrix
#' @param wmat Weights Matrix
#' @param verbose print output
#' @param checkSym check if wmat is symmetric
#' @param XOmegaX meat function
#' @return object to be used in vcov* functions
#  @examples
#' @export


iterateObsJSpatial <- compiler::cmpfun( function(
    sub_dat,
    Xvars,
    wmat, #wid,
    verbose=TRUE, 
    checkSym=FALSE,
    XOmegaX=XOmegaX0,
    ...){




    # Make sure Distance Matrix matches with residuals
    # s0     <- order(sub_dat[,unit])
    # s1     <- order( match(sub_dat[,unit], wid) )
    # sub_dat <- sub_dat[s1,]

    if(verbose){
        message(paste("Starting on sub index:",
            sub_dat[,"time"][1]))
        if(checkSym) {
            message("WMAT symmetric: ", isSymmetric(WMAT))
        }
    }

    sub_dat <- as.data.frame(sub_dat)
    X   <- as.matrix(sub_dat[,Xvars])
    e   <- sub_dat[,"e"]
    
    XeeXhs <- XOmegaX(X, wmat, e)
    return(XeeXhs)
})

#------------------------------------------------------------------
##################
#' Wrapper for Matrix Calculation for Time
################## 
#' @param sub_dat a data frame object for one cell unit
#' @param Xvars RHS design matrix
#' @param wmat Weights Matrix
#' @param verbose print output
#' @param checkSym check if wmat is symmetric
#' @param XOmegaX meat function
#' @return object to be used in vcov* functions
#  @examples
#' @export


iterateObsJTemporal <- compiler::cmpfun( function(
    sub_dat,
    Xvars,
    wmat, #wid,
    verbose=TRUE, 
    checkSym=FALSE,
    XOmegaX=XOmegaX0,
    ...){

    # Make sure Distance Matrix matches with residuals
    # s0     <- order(sub_dat[,time])
    # s1     <- order( match(sub_dat[,time], wid) )
    # sub_dat <- sub_dat[s1,]

    if(verbose){
        message(paste("Starting on sub index:",
            sub_dat[,"unit"][1]))
        if(checkSym) {
            message("WMAT symmetric: ", isSymmetric(WMAT))
        }
    }

    sub_dat <- as.data.frame(sub_dat)
    X   <- as.matrix(sub_dat[,Xvars])
    e   <- sub_dat[, e]
    
    XeeXhs <- XOmegaX(X, wmat, e)
    return(XeeXhs)
} )

