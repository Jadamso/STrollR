
#------------------------------------------------------------------
##################
#' XOmegaX 
################## 
#' Meat Matrix Calculations
#' @param X design matrix 
#' @param WMAT weighting matrix (preferably sparse sparse)
#' @param loc_lat, loc_lon spatial-coordinate vectors
#' @param loc_y, loc_x spatial-coordinate vectors
#' @param loc_t time-coordinate vector
#' @param E vector of residuals
#' @param cutoff_s how far to extend bartlett kernel in spatial dimention (in km)
#' @param cutoff_t how far to extend bartlett kernel in time dimension (in time units)
#' @return object to be used in vcov* functions
#  @examples
#' @details 
#'  XOmegaX_rolled Fully is for Sparse-Matrix Meat-Matrix Calculation (Whole Matrix At Once)
#'  XOmegaX_semirolled is for Sparse-Matrix Meat-Matrix Calculation (One Row At A Time)
#'  XOmegaX_bruteforce_ll is for latlon data.
#'  XOmegaX_bruteforce_llt is for latlon and time data.
#'  XOmegaX_bruteforce_xy is for projected data.
#'  XOmegaX_bruteforce_xyt is for projected and time data.
#'  XOmegaX_bruteforce_*t weights the time-dimension according to bartlett kernel 
#' @name XOmegaX
NULL



#------------------------------------------------------------------
##################
# XOmegaX Rolled Weights Matrices
##################
#' @rdname XOmegaX 
# Potential Speedup if WEE is symmetric? WEE=(e * t(WMAT * e))
XOmegaX_rolled <- compiler::cmpfun( function(X, WMAT, E){
    require("Matrix")
    XWeeX  <- t(X) %*% t(E * t(WMAT * E)) %*% X
    return( as.matrix(XWeeX) )
} )

#' @rdname XOmegaX
XOmegaX_semirolled <- compiler::cmpfun( function(X, WMAT, E){
    require("Matrix")
    WEE <- sapply(1:nrow(WMAT), function(i){
        WMAT[i,]*E[i]*E
    })
    XWeeX  <- as.matrix(t(X) %*% WEE %*% X)
    return( XWeeX )
} )


#------------------------------------------------------------------
##################
# XOmegaX Brute Force Kernel Weighting for Space
##################
#' @rdname XOmegaX 
## Taken from Fiona Burlig
XOmegaX_bruteforce_ll <- compiler::cmpfun( function(X, E, loc_lat, loc_lon, cutoff_s, manual_dist=T){
    k <- ncol(X)
    n <- nrow(X)
    meatWeight <- lapply(1:length(E), function(i) {

        # Calculate Distances in km's
        if(manual_dist){
            ## Manual Calculation Turn longitude & latitude into KM
            ## Turn longitude & latitude into KMs
            ## 1 deg lat = 111 km, 1 deg lon = 111 km* cos(lat)
            lonscale <- cos(loc_lat[i]*pi / 180) * 111
            latscale <- 111
            disti <- as.numeric(sqrt((latscale*(loc_lat[i] - loc_lat))^2
            + (lonscale*(loc_lon[i] - loc_lon))^2))
        } else {
            ## Spam Calculation
            ll_mat <- cbind(loc_lat,loc_lon)
            disti <- c(spam::spam_rdist.earth( t(ll_mat[i,]), ll_mat, delta=cutoff_s, miles=FALSE))
        }

        ## Bartlett Kernel of distances
        wi <- (1 - disti/cutoff_s)*as.numeric(disti <= cutoff_s)

        # sum_j(X_iX_j'e_ie_j K(d_{ij})),
        #XWeeX <- ((t(t(X[i, ])) %*% matrix(1, 1, n) * E[i]) * (matrix(1, k, 1) %*% (t(E) * t(wi)))) %*% X
        XWeeX <- ((t(t(X[i, ])) %*% matrix(1, 1, n) * E[i]) * (matrix(1, k, 1) %*% t(E * wi))) %*% X
        return(XWeeX)
    })
    XWeeX <- as.matrix( Reduce("+", meatWeight) )

    return( XWeeX )
} )

#' @rdname XOmegaX
XOmegaX_bruteforce_xy <- compiler::cmpfun( function(X, E, loc_y, loc_x, cutoff_s, manual_dist=T){
    k <- ncol(X)
    n <- nrow(X)
    meatWeight <- lapply(1:length(E), function(i) {
        # Calculate Distances in km's
        if(manual_dist){
            ## Manual Calculation
            disti <- as.numeric(sqrt( (loc_x[i] - loc_x)^2 + (loc_y[i] - loc_y)^2))
        } else {
            ## Spam Calculation
            xy_mat <- cbind(loc_x,loc_y)
            disti  <- c(spam::spam_rdist( t(xy_mat[i,]), xy_mat, delta=cutoff_s))
        }
        wi <- (1 - disti/cutoff_s)*as.numeric(disti <= cutoff_s)
        #XWeeX <- ((t(t(X[i, ])) %*% matrix(1, 1, n) * E[i]) * (matrix(1, k, 1) %*% (t(E) * t(wi)))) %*% X
        XWeeX <- ((t(t(X[i, ])) %*% matrix(1, 1, n) * E[i]) * (matrix(1, k, 1) %*% t(E * wi))) %*% X
        return(XWeeX)
    })
    XWeeX <- as.matrix( Reduce("+", meatWeight) )
    return( XWeeX )
} )



#------------------------------------------------------------------
##################
# XOmegaX Brute Force Kernel Weighting for Space and Time
##################
#' @rdname XOmegaX
XOmegaX_bruteforce_llt <- compiler::cmpfun( function(
    X, E, loc_lat, loc_lon, loc_t, cutoff_s, cutoff_t, manual_dist=T){
    
    k <- ncol(X)
    n <- nrow(X)
    meatWeight <- lapply(1:length(E), function(i) {

        # Calculate Spatial Distances in km's
        if(manual_dist){
            ## Manual Calculation Turn longitude & latitude into KM
            ## Turn longitude & latitude into KMs
            ## 1 deg lat = 111 km, 1 deg lon = 111 km* cos(lat)
            lonscale <- cos(loc_lat[i]*pi / 180) * 111
            latscale <- 111
            disti_s <- as.numeric(sqrt((latscale*(loc_lat[i] - loc_lat))^2
            + (lonscale*(loc_lon[i] - loc_lon))^2))
        } else {
            ## Spam Calculation
            ll_mat <- cbind(loc_lat,loc_lon)
            disti_s <- c(spam::spam_rdist.earth( t(ll_mat[i,]), ll_mat, delta=cutoff_s, miles=FALSE))
        }

        # Calculate Time Distances
        disti_t <- as.numeric( abs(loc_t[i] - loc_t) )

        ## Bartlett Kernel Weights
        wi_s <- (1 - disti_s/cutoff_s)*as.numeric(disti_s <= cutoff_s)
        wi_t <- (1 - disti_t/cutoff_t)*as.numeric(disti_t <= cutoff_t)
        wi   <- wi_t*wi_s
        
        # sum_j(X_iX_j'e_ie_j K(d_{ij}))
        #XWeeX <- ((t(t(X[i, ])) %*% matrix(1, 1, n) * E[i]) * (matrix(1, k, 1) %*% (t(E) * t(wi)))) %*% X
        XWeeX <- ((t(t(X[i, ])) %*% matrix(1, 1, n) * E[i]) * (matrix(1, k, 1) %*% t(E * wi))) %*% X
        return(XWeeX)
    })
    XWeeX <- as.matrix( Reduce("+", meatWeight) )

    return( XWeeX )
} )

#' @rdname XOmegaX
XOmegaX_bruteforce_xyt <- compiler::cmpfun( function(
    X, E, loc_y, loc_x, loc_t, cutoff_s, cutoff_t, manual_dist=T){

    k <- ncol(X)
    n <- nrow(X)
    meatWeight <- lapply(1:length(E), function(i) {

        # Calculate Spatial Distances in km's
        if(manual_dist){
            ## Manual Calculation
            disti_s <- as.numeric(sqrt( (loc_x[i] - loc_x)^2 + (loc_y[i] - loc_y)^2))
        } else {
            ## Spam Calculation
            xy_mat <- cbind(loc_x,loc_y)
            disti  <- c(spam::spam_rdist( t(xy_mat[i,]), xy_mat, delta=cutoff_s))
        }

        # Calculate Time Distances
        disti_t <- as.numeric( abs(loc_t[i] - loc_t) )

        ## Bartlett Kernel Weights
        wi_s <- (1 - disti_s/cutoff_s)*as.numeric(disti_s <= cutoff_s)
        wi_t <- (1 - disti_t/cutoff_t)*as.numeric(disti_t <= cutoff_t)
        wi   <- wi_t*wi_s
        
        #XWeeX <- ((t(t(X[i, ])) %*% matrix(1, 1, n) * E[i]) * (matrix(1, k, 1) %*% (t(E) * t(wi)))) %*% X
        XWeeX <- ((t(t(X[i, ])) %*% matrix(1, 1, n) * E[i]) * (matrix(1, k, 1) %*% t(E * wi))) %*% X
        return(XWeeX)
    })
    
    XWeeX <- as.matrix( Reduce("+", meatWeight) )
    return( XWeeX )
} )



#------------------------------------------------------------------
##################
#' Model Frame with Regression Intercept
################## 
model.frame.i <- function(reg, check_int=T){

    Xvars <- names(reg$coefficients)
    XY    <- model.frame(reg)
    if( ('(Intercept)' %in% Xvars) & !( '(Intercept)' %in% names(XY) ) ){
        XY  <- cbind(XY, '(Intercept)'=1)
    }
    return(XY)
}

#------------------------------------------------------------------
##################
# vcovSHAC
################## 
#' Calculate a SHAC (Spatial Heteroskedastic and Autocorrelation Consistent) Variance Covariance Matrix
#' @param reg an `lm' object
#' @param wmat weights matrix
#' @param method c(`rolled', `semirolled', 'bruteforce_ll', 'bruteforce_xy')
#' @param cutoff_s include weights up to cutoff_s. Required to be in km if method='bruteforce_ll' or in map-units if 'bruteforce_xy'.
#' @param loc_lat,loc_lon required if method='bruteforce_ll'
#' @param loc_y,loc_x required if method='bruteforce_xy'
#' @param verbose show messages
#' @param options passed to sandwich
#' @param add_hc,add_cluster,add_hac logical for adding HC correction (default=T), Cluster correction (default=F), HAC correction (default=F)
#' @return covariance matrix
vcovSHAC <- compiler::cmpfun( function(
    reg,
    wmat,
    method='bruteforce_ll',
    cutoff_s,
    loc_lat, 
    loc_lon,
    loc_y,
    loc_x,
    add_hc=T,
    add_cluster=F,
    add_hac=F,
    verbose=FALSE,
    manual_dist=T,
    ...) {

    ##############################################
    # Prepare Covariance Matrix Inputs

    # Model
    Xvars <- names(coef(reg))
    XY <- model.frame.i(reg)

    X <- as.matrix(XY[,Xvars])

    # Residuals
    E <- resid(reg)
    
    # Parameters
    n <- length(E)
    #k <- length(Xvars)

    ##############################################
    # Correct for spatial correlation:
    
    ## SAC Meat Matrix
    if(verbose) {message("Computing XEEX.")}
    if(method=='rolled'){
        if(class(wmat)!='dgCMatrix'){ break }
        XeeX_space <- XOmegaX_rolled(X, wmat, E)/n
    } else if(method=='semirolled'){
        XeeX_space <- XOmegaX_semirolled(X, wmat, E)/n
    } else if(method=='bruteforce_ll'){
        if(verbose) {message("Requires lat, lon coordinates and cutoff (in km)")}
        XeeX_space <- XOmegaX_bruteforce_ll(X, E, loc_lat, loc_lon, cutoff_s, manual_dist)/n
    } else if(method=='bruteforce_xy'){
        if(verbose) {message("Requires projected x,y coordinates and cutoff (in map-projection units)")}
        XeeX_space <- XOmegaX_bruteforce_xy(X, E, loc_y, loc_x, cutoff_s, manual_dist)/n
    }
    
    ##############################################
    ## Calculate Consistent Sandwich
    if(verbose) {message("Computing Sandwhich.")}
    ## Correct For Heteroskedasticity
    if(add_hc) {
        XeeX0 <- sandwich::meatHC(reg,...)
    } else {
        XeeX0 <- 0
    }
    ## Correct For Clustering
    if(add_cluster){
        XeeXCL <- sandwich::meatCL(reg,...)
    } else {
        XeeXCL <- 0
    }
    ## Correct For Heteroskedasticity and Autocorrelation
    if(add_hac){
        XeeXHAC <- sandwich::meatHAC(reg,...)
    } else {
        XeeXHAC <- 0
    }
    
    invXX <- solve(t(X) %*% X) # Bread Matrix
    Omega <- ((XeeX0 + XeeX_space + XeeXCL + XeeXHAC)*n) # Meat Matrix
    V_SHAC <- t(invXX) %*% Omega %*% invXX # Sandwhich
    V_SHAC <- (V_SHAC + t(V_SHAC)) / 2

    return(V_SHAC)
} )


#------------------------------------------------------------------
##################
# vcovSHACsep
################## 
#
#' @describeIn vcovSHAC similar to vcovSTHAC, but treats spatial autocorrelation seperately within each time period
vcovSHACsep <- compiler::cmpfun( function(
    reg,
    wmat,
    method='bruteforce_llt',
    cutoff_s,
    cutoff_t,    
    loc_lat, 
    loc_lon,
    loc_y,
    loc_x,
    loc_t,
    add_hc=T,
    add_cluster=F,
    add_hac=F,
    verbose=FALSE,
    manual_dist=T,
    ...) {

    ##############################################
    # Prepare Covariance Matrix Inputs

    # Model
    Xvars <- names(coef(reg))
    XY <- model.frame.i(reg)

    X <- as.matrix(XY[,Xvars])

    # Residuals
    E <- resid(reg)
    
    # Parameters
    n <- length(E)
    #k <- length(Xvars)


    ##############################################
    # Correct for spatial correlation:
   
    ## STAC Meat Matrix
    if(verbose) {message("Computing XEEX.")}
    if(method=='rolled'){
        if(class(wmat)!='dgCMatrix'){ break }
        Et_list <- split(E, loc_t)
        Xt_list <- split(as.data.frame(X), loc_t)
        XeeX_space <- lapply( seq(Et_list), function(t_id){
            Et <- Et_list[[t_id]]
            Xt <- as.matrix(Xt_list[[t_id]])
            XOmegaX_rolled(Xt, wmat, Et)/n
        })
        XeeX_space <- Reduce("+",  XeeX_space)
        
    } else if(method=='semirolled'){
        Et_list <- split(E, loc_t)
        Xt_list <- split(as.data.frame(X), loc_t)
        XeeX_space <- lapply( seq(Et_list), function(t_id){
            Et <- Et_list[[t_id]]
            Xt <- as.matrix(Xt_list[[t_id]])
            XOmegaX_semirolled(Xt, wmat, Et)/n
        })
        XeeX_space <- Reduce("+",  XeeX_space)
    } else if(method=='bruteforce_ll'){
        if(verbose) {message("Requires lat, lon coordinates and cutoff (in km)")}
        Et_list <- split(E, loc_t)
        Xt_list <- split(as.data.frame(X), loc_t)
        loc_latt_list <- split(loc_lat, loc_t)
        loc_lont_list <- split(loc_lon, loc_t)
        XeeX_space <- lapply( seq(Et_list), function(t_id){
            Et <- Et_list[[t_id]]
            Xt <- as.matrix(Xt_list[[t_id]])
            loc_latt <- loc_latt_list[[t_id]]
            loc_lont <- loc_lont_list[[t_id]]            
            XOmegaX_bruteforce_ll(Xt, Et, loc_latt, loc_lont, cutoff_s, manual_dist)/n
        })
        XeeX_space <- Reduce("+",  XeeX_space)
    } else if(method=='bruteforce_xy'){
        if(verbose) {message("Requires projected x,y coordinates and cutoff (in map-projection units)")}
        Et_list <- split(E, loc_t)
        Xt_list <- split(as.data.frame(X), loc_t)
        loc_yt_list <- split(loc_y, loc_t)
        loc_xt_list <- split(loc_x, loc_t)
        XeeX_space <- lapply( seq(Et_list), function(t_id){
            Et <- Et_list[[t_id]]
            Xt <- as.matrix(Xt_list[[t_id]])
            loc_yt <- loc_yt_list[[t_id]]
            loc_xt <- loc_xt_list[[t_id]]            
            XOmegaX_bruteforce_xy(Xt, Et, loc_yt, loc_xt, cutoff_s, manual_dist)/n
        })
        XeeX_space <- Reduce("+",  XeeX_space)
    }


    ##############################################
    ## Calculate Consistent Sandwich
    if(verbose) {message("Computing Sandwhich.")}
    ## Correct For Heteroskedasticity
    if(add_hc) {
        XeeX0 <- sandwich::meatHC(reg,...)
    } else {
        XeeX0 <- 0
    }
    ## Correct For Clustering
    if(add_cluster){
        XeeXCL <- sandwich::meatCL(reg,...)
    } else {
        XeeXCL <- 0
    }
    ## Correct For Heteroskedasticity and Autocorrelation
    if(add_hac){
        XeeXHAC <- sandwich::meatHAC(reg,...)
    } else {
        XeeXHAC <- 0
    }
    
    invXX <- solve(t(X) %*% X) # Bread Matrix
    Omega <- ((XeeX0 + XeeX_space + XeeXCL + XeeXHAC)*n) # Meat Matrix
    V_SHAC <- t(invXX) %*% Omega %*% invXX # Sandwhich
    V_SHAC <- (V_SHAC + t(V_SHAC)) / 2

    return(V_SHAC)

} )


#------------------------------------------------------------------
##################
# vcovSTHAC
################## 
#' @describeIn vcovSHAC  Space and Time HAC. Also weights the time-dimension according to bartlett kernel (i.e., weight= K(space)*K(time)) for bartlett-kernel K.
vcovSTHAC <- compiler::cmpfun( function(
    reg,
    wmat,
    method='bruteforce_llt',
    cutoff_s,
    cutoff_t,    
    loc_lat, 
    loc_lon,
    loc_y,
    loc_x,
    loc_t,
    add_hc=T,
    add_cluster=F,
    add_hac=F,
    verbose=FALSE,
    manual_dist=T,
    ...) {

    ##############################################
    # Prepare Covariance Matrix Inputs

    # Model
    Xvars <- names(coef(reg))
    XY <- model.frame.i(reg)

    X <- as.matrix(XY[,Xvars])

    # Residuals
    E <- resid(reg)
    
    # Parameters
    n <- length(E)
    #k <- length(Xvars)

    ##############################################
    # Correct for spatial correlation:

    ## STAC Meat Matrix
    if(verbose) {message("Computing XEEX.")}
    if(method=='rolled'){
        if(class(wmat)!='dgCMatrix'){ break }
        XeeX_spacetime <- XOmegaX_rolled(X, wmat, E)/n
    } else if(method=='semirolled'){
        XeeX_spacetime <- XOmegaX_semirolled(X, wmat, E)/n
    } else if(method=='bruteforce_llt'){
        if(verbose) {message("Requires lat, lon coordinates and cutoff (in km)")}
        XeeX_spacetime <- XOmegaX_bruteforce_llt(X, E, loc_lat, loc_lon, loc_t, cutoff_s, cutoff_t, manual_dist)/n
    } else if(method=='bruteforce_xyt'){
        if(verbose) {message("Requires projected x,y coordinates and cutoff (in map-projection units)")}
        XeeX_spacetime <- XOmegaX_bruteforce_xyt(X, E, loc_y, loc_x, loc_t, cutoff_s, cutoff_t, manual_dist)/n
    }
    
    ##############################################
    ## Calculate Consistent Sandwich
    if(verbose) {message("Computing Sandwhich.")}
    ## Correct For Heteroskedasticity
    if(add_hc) {
        XeeX0 <- sandwich::meatHC(reg,...)
    } else {
        XeeX0 <- 0
    }
    ## Correct For Clustering
    if(add_cluster){
        XeeXCL <- sandwich::meatCL(reg,...)
    } else {
        XeeXCL <- 0
    }
    ## Correct For Heteroskedasticity and Autocorrelation
    if(add_hac){
        XeeXHAC <- sandwich::meatHAC(reg,...)
    } else {
        XeeXHAC <- 0
    }
    
    invXX <- solve(t(X) %*% X) # Bread Matrix
    Omega <- ((XeeX0 + XeeX_spacetime + XeeXCL + XeeXHAC)*n) # Meat Matrix
    V_SHAC <- t(invXX) %*% Omega %*% invXX # Sandwhich
    V_SHAC <- (V_SHAC + t(V_SHAC)) / 2

    return(V_SHAC)
} )

