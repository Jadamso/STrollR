message("SpaceTimeSE_format")
#------------------------------------------------------------------
##################
# Helper Functions
################## 
#' Converts Factor to Number
#' @param x numeric factor
#' @export
Fac2Num <- compiler::cmpfun( function(x) {
    as.numeric(as.character(x)) } )


#' checks if data table has lattice structure ?
#' @param DAT a data.table
#' @param TIME name of temporal column
#' @param ID name of cellular ID variable
#' @export
is.lattice <- compiler::cmpfun( function( DAT, TIME, ID) {

    requireNamespace("data.table")
    
    utime <- unique(DAT[,..TIME])
    uid   <- unique(DAT[,..ID])

    uit_i <- sapply( unlist(uid), function(uid_i){
        uid_i_row <- which( DAT[,..ID]==uid_i[[1]] )
        uid_i_t   <- DAT[ uid_i_row, ..TIME]
        all( utime %in% uid_i_t ) 
    })
    names(uit_i) <- unlist(uid)

    uit_i
})





#------------------------------------------------------------------
##################
#' Data Matrix Preparations
################## 
#' @param ddff formatted data.table from vcovST.format
#' @param wmat spatial weights matrix
#' @param tmat temporal weights matrix
#' @param t_cutoff temporal cutoff
#' @param rho_t unsupported vonneumann structure
#' @param d_cutoff distance cutoff
#' @param rho_sp unsupported vonneumann structure
#' @param latlon coordinates in lon,lat or x,y
#' @param convert_to_angles lon,lat to x,y?
#' @return list object to be passed to vcov* functions
#  @examples
#' @export

## Spatial Weights
#source( paste0(sedir, "SpatialWeights.R") )



listj <- compiler::cmpfun( function(
    ddff,
    wmat=NA,
    tmat=NA,
    t_cutoff=4,
    rho_t=NA,
    d_cutoff=1,
    rho_sp=NA,
    latlon=NA,
    convert_to_angles=TRUE){

    ##############################################
    ## Formatting

    which_numeric <- c( sapply( ddff, class) == "numeric")
    all_numeric <- all( which_numeric[
        c("unit", "time", "sp_coord1", "sp_coord2")])
    if( !all_numeric)  {
        message("Error: sp_coords, t_coords, id  must be numeric")
        break
    }

    TIME <- ddff[,time]
    UniqueTime <- as.matrix(unique(TIME))

    TIMEsub  <- ddff[TIME==UniqueTime[1],]
    ID       <- as.numeric( TIMEsub[, unit] )

    
    if( class(wmat)=="logical" ){

        ##############################################
        ## Compute Spatial Weights
        M <- as.matrix( cbind( TIMEsub[, sp_coord1],
            TIMEsub[, sp_coord2] ) )

        if( is.na(rho_sp) ){
            ## Using Bartlett Weights
            WMAT <- weight_mat(M, d_cutoff, latlon)
        } else {
        
            message("currently unsupported")
            WMAT <- weight_mat(M, d_cutoff, latlon)
        
            ## Using AR Weights
            #Mmat <- VonNeumann( SpatialPoints(M))
            #Wmat <- Diagonal(nrow(M)) -rho_sp*Mmat
            #Smat <- spdep::invIrW(Mmat, rho=rho_sp)
            #Smat <- spdep::powerWeights(Mmat,
            #    rho=rho_sp, order=2, X=Diagonal( nrow(Mmat)))
        }

        message("WMAT symmetric: ", isSymmetric(WMAT))

    } else {

        WMAT <- wmat
    }


    if( class(tmat)=="logical" ){

        ##############################################
        ## Compute Temporal Weights
        if(is.na(rho_t)){

            ## Using Bartlett Weights
            TMAT <- weight_mat( UniqueTime, t_cutoff)
            #TMAT <- TimeDistMat(UniqueTime, t_cutoff)
            #TMAT <- Matrix(TMAT)
        } else {

            ## Using AR Weights
            TMAT <- toeplitz( rho**seq(0, length(UniqueTime) -1) )
        }

        message("TMAT symmetric: ", isSymmetric(TMAT))
    } else {
        TMAT <- tmat
    }

    ##############################################
    return( list(
        wmat=WMAT, wid=ID,
        tmat=TMAT, tid=as.numeric(UniqueTime),
        t_cutoff=t_cutoff, d_cutoff=d_cutoff
        )
    )
})


#------------------------------------------------------------------
##################
#' Format Regression Output
################## 
#' @param reg felm object to be shaped
#' @param unit string name for cellular_id variable
#' @param time string name for time variable
#' @param sp_coords string name for coordinate variables
#' @param DF data.table to format, <NULL> formatts 
#' @param SPLIT_sp return list dataframe of DF for each cell
#' @param SPLIT_t return list dataframe of DF for each time
#' @param SPLIT SPLIT_sp=SPLIT_sp=TRUE
#' @return object to be used in vcov* functions
#  @examples
#' @export


vcovST.format <- function(
    reg,
    unit,
    time,
    sp_coords,
    DF=NA,
    SPLIT_sp=FALSE,
    SPLIT_t=FALSE,
    SPLIT=FALSE) {
  
    ##############################################
    # Format Data

    if( class(reg) == "felm") {
        #dat <- data.frame(reg$cY, reg$cX)
        dat <- as.data.frame( cbind(reg$cY, reg$cX) )
    } else if(class(reg) == "lm") {
        dat <- as.data.frame(reg$model)
    } else {
        message("Model class not recognized.")
        break
    }

    # Space and Time Clustering
    if( ! is.null( names(DF) ) ){
    
        dat[, "unit"] <- Fac2Num(DF[,unit]) 
        dat[, "time"] <- Fac2Num(DF[,time])
        dat[, "sp_coord1"] <- DF[,sp_coords[1]]
        dat[, "sp_coord2"] <- DF[,sp_coords[2]]
        
    } else {
    
        if( class(reg) == "felm") {
            message( paste0("Assumes FE1 is", unit) )
            message("sp_coords not yet supported")
            break
            
            ## data.table formatting
            dat <- dat[, unit:= Fac2Num(reg$fe[[1]])]
            dat <- dat[, time:= Fac2Num(reg$fe[[2]])]
         
        } else if (class(reg) == "lm") {
            message("must pass DF if reg is an lm object")
            break
        }
    }

    ## Residuals
    dat[, "e"] <- as.numeric(reg$residuals)


    ## Orders Data
    dat <- data.table::as.data.table(dat)
    data.table::setkey(dat, unit, time)

    ## Xvars
    Xvars <- rownames(reg$coefficients)


    if(SPLIT){ SPLIT_t <- SPLIT_sp <- TRUE }

    if(SPLIT_t) {
        ## List of Data Frames by TimePeriod
        TSUB <- split( dat , f = dat$time )
        #sub_dat <- dat[time == sub_index]
    } else {
        TSUB <- NA
    }

    if(SPLIT_sp) {
        ## List of Data Frames by ID
        IDSUB <- split( dat , f = dat$unit )
        #sub_dat <- dat[unit == sub_index]
    } else {
        IDSUB <- NA
    }

    # Return
    RETlist <- list( 
        Xvars=Xvars,
        dat=dat,
        TSUB=TSUB, 
        IDSUB=IDSUB
    )
    return(RETlist)
}

