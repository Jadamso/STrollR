#------------------------------------------------------------------
##################
#' make_space_data
################## 
#'
#' @param dimS spatial dimension 
#' @param K number of covariates
#' @param t_id Time ID
#' @param space_groups number of contiguous groups (e.g. countries)
#' @param xpars generate each Xvariable from uniform(xpar[1],xpar[2]) 
#' @param theta parameter vector for Y=X%*%theta
#' @param error_type c('focal', 'distance', 'spherical')
#' @param error_scale variance of errors (normally distributed) 
#' @param sar_factor parameter for error_type='focal'
#' @param wf_mat matrix specifing weights to smooth when error_type='focal'
#' @param xy_mat matrix specifing xy locations
#' @param ws_mat matrix specifing correlation structer of error_type='distance' 
#' 
#' @return dataframe
#' 
#' @description  Create Multivariate Spatial Data
#'
#' @details
#' 'focal' creates errors based on \eqn{e = u + \rho v}
#' 'distance' creates errors base on drawing e from a multivariate normal with covariance matrix ws_mat
#'  
#' @examples
#' set.seed(3)
#' DFs  <- make_space_data.raster(10,2,error_type='distance')
#' set.seed(3)
#' DFs2 <- make_space_data(10,2,error_type='distance')
#' @export
make_space_data <- compiler::cmpfun( function(
    dimS,
    K,
    t_id=NA,
    space_groups=round(sqrt(dimS)),
    xpars=c(0,1),
    theta=1:K,
    error_type='focal',
    error_scale=1,
    sar_factor=.5,
    wf_mat=rbind( rep(sar_factor,3), c(sar_factor,0,sar_factor),rep(sar_factor,3)),
    xy_mat=expand.grid(Loc_X=seq(0,1,length.out=dimS),Loc_Y=seq(1,0,length.out=dimS)),
    ws_mat=exp(-as.matrix(dist( xy_mat ) ))
){

    ## Independant Variables
    Xrast <- sapply(1:K, function(k){
        runif(dimS*dimS, min=xpars[1], max=xpars[2])
    })
    #Xrast <- mvtnorm::rmvnorm(K,sigma=diag(dimS*dimS) )
    colnames(Xrast) <- paste0('X', 1:K)

    ## Space ID
    Srast <- as.matrix(1:nrow(xy_mat))
    colnames(Srast) <- "Space_ID"
    
    ## Time ID
    Trast <- matrix(t_id, ncol=1, nrow=dimS*dimS)
    colnames(Trast) <- 'Time_ID'
        
    ## Space Groups
    g1 <- cut(xy_mat[,1], space_groups)
    g2 <- cut(xy_mat[,2], space_groups)
    g3 <- as.factor(paste0(g1, g2))
    for(i in seq(levels(g3))){ levels(g3)[i] <- i }
    SGrast <- as.matrix(as.numeric(g3))
    colnames(SGrast) <- "Space_Group"

    ## Errors
    if(error_type=='focal'){
        require(raster)
        Erast0 <- raster::raster( matrix( rnorm(dimS*dimS, mean=0, 1), ncol=dimS, nrow=dimS))
        Erast1 <- raster::focal(Erast0, wf_mat, sum, pad=TRUE, padValue=NA, na.rm=T)
        Erast  <- raster::raster( matrix( rnorm(dimS*dimS, mean=0, 1), ncol=dimS, nrow=dimS)) + Erast1
        Erast  <- as.matrix(c(as.matrix(Erast)))*error_scale
    } else if (error_type=='distance') {
        Erast0 <- mvtnorm::rmvnorm(1, mean=rep(0, nrow(ws_mat)), sigma=ws_mat)
        Erast <- t(Erast0)*error_scale
    } else if (error_type=='spherical') {
        Erast0 <- mvtnorm::rmvnorm(1, mean=rep(0, nrow(Xrast)), sigma=diag( nrow(Xrast)))
        Erast <- t(Erast0)*error_scale
    }
    colnames(Erast) <- 'E'

    ## Dependant Variables
    Yrast <- Xrast %*% as.matrix(theta) + Erast
    colnames(Yrast) <- 'Y'

    ## Raster Stack
    STACK <- cbind(
        xy_mat,
        Srast,
        Trast,
        Yrast,
        Xrast,
        SGrast,
        Erast)

    ## To Data Frame
    DFs <- as.data.frame(STACK)
    return(DFs)
})


#' @describeIn make_space_data make_space_data.df uses raster functions instead of matrices (primarily for transparent debugging)
make_space_data.raster <- compiler::cmpfun( function(
    dimS,
    K,
    t_id=NA,
    space_groups=round(sqrt(dimS)),
    xpars=c(0,1),
    theta=1:K,
    error_type='focal',
    error_scale=1,
    sar_factor=.5,
    wf_mat=rbind( rep(sar_factor,3), c(sar_factor,0,sar_factor),rep(sar_factor,3)),
    xy_mat=expand.grid(Loc_X=seq(0,1,length.out=dimS),Loc_Y=seq(1,0,length.out=dimS)),
    ws_mat=exp(-as.matrix(dist( xy_mat ) ))
){

    require(raster)
    base_rast <- raster::rasterFromXYZ(xy_mat)
    
    ## Independant Variables
    Xrast_list <- lapply(1:K, function(k){
        Xrast <- base_rast
        values(Xrast) <- runif(dimS*dimS, min=xpars[1], max=xpars[2])
        return(Xrast)
    })
    Xrast <- raster::stack(Xrast_list)
    names(Xrast) <- paste0('X', 1:K)
    
    ## Space ID
    Srast <- base_rast
    values(Srast) <- 1:nrow(xy_mat)
    names(Srast) <- "Space_ID"
    
    ## Time ID
    Trast <- base_rast
    values(Trast) <- t_id
    names(Trast) <- 'Time_ID'

    ## Space Groups
    SGrast <- base_rast
    g1 <- cut(xy_mat[,1], space_groups)
    g2 <- cut(xy_mat[,2], space_groups)
    g3 <- as.factor(paste0(g1, g2))
    for(i in seq(levels(g3))){ levels(g3)[i] <- i }
    raster::values(SGrast) <- as.numeric(g3)
    names(SGrast) <- "Space_Group"


    ## Errors
    if(error_type=='focal'){
        Erast0 <- Erast1 <- base_rast
        values(Erast0) <- rnorm(dimS*dimS, mean=0, 1)*error_scale
        values(Erast1) <- rnorm(dimS*dimS, mean=0, 1)*error_scale
        Erast <- focal(Erast0, wf_mat, sum, pad=TRUE, padValue=NA, na.rm=T)+ Erast1
    } else if (error_type=='distance') {
        Erast0 <- mvtnorm::rmvnorm(1, mean=rep(0, nrow(ws_mat)), sigma=ws_mat)
        Erast  <- raster::rasterFromXYZ(cbind(xy_mat, t(Erast0)))*error_scale
    } else if (error_type=='spherical') {
        Erast0 <- mvtnorm::rmvnorm(1, mean=rep(0,nrow(xy_mat)), sigma=diag(nrow(xy_mat)))
        Erast  <- raster::rasterFromXYZ(cbind(xy_mat, t(Erast0)))*error_scale
    }    
    names(Erast) <- 'E'


    ## Dependant Variables
    Yrast <- sum(raster::stack(Xrast_list)*theta) + Erast
    names(Yrast) <- 'Y'

    
    ## Raster Stack
    STACK <- raster::stack(
        Srast,
        Trast,
        Yrast,
        Xrast,
        SGrast,
        Erast)

    ## To Data Frame
    DFs <- raster::as.data.frame(STACK)
    DFs <- cbind(xy_mat,DFs)
    return(DFs)
})



#------------------------------------------------------------------
##################
#' Create Multivariate Temporal Data
################## 
#'
#' @param dimT temporal dimension 
#' @param K number of covariates
#' @param s_id Spatial ID
#' @param time_groups number of contiguous groups (e.g. time-regimes)
#' @param xpars generate each Xvariable from uniform(xpar[1],xpar[2]) 
#' @param theta parameter vector for Y=X%*%theta
#' @param error_type c('focal', 'distance', 'spherical')
#' @param error_scale variance of errors (normally distributed) 
#' @param ar_factor parameter for error_type='focal'
#' @param wt_mat matrix specifing correlation structer of error_type='distance'
#' @param xy_mat matrix specifing xy locations (defaults to NA)

#' 
#' @return dataframe
#'
#' @details 
#' 'focal' creates errors based on \eqn{e = u + \rho v}
#' 'distance' creates errors base on drawing e from a multivariate normal with covariance matrix wt_mat
#'
#' @examples DFt  <- make_time_data(6,2)
#'
#' @export


make_time_data <- compiler::cmpfun( function(
    dimT,
    K,
    s_id=NA,
    time_groups=round(sqrt(dimT)),
    xpars=c(0,1),
    theta=1:K,
    error_type='focal',
    error_scale=1,
    ar_factor=2/3,
    wt_mat=toeplitz(c(1,ar_factor,ar_factor^2, rep(0,dimT-3))),
    xy_mat=data.frame(Loc_X=NA,Loc_Y=NA)
){

    ## Independant Variables
    Xrast <- sapply(1:K, function(k){
        Xrast <- runif(dimT, min=xpars[1], max=xpars[2])
    })
    colnames(Xrast) <- paste0('X', 1:K)

    ## Space ID
    Srast <- matrix( s_id, ncol=1, nrow=dimT)
    colnames(Srast) <- 'Space_ID'

    ## Time ID
    Trast <- as.matrix(1:dimT)
    colnames(Trast) <- 'Time_ID'
          
    ## Time Groups
    TGrast <- cut( 1:dimT, time_groups)
    for(i in seq(levels(TGrast))){ levels(TGrast)[i] <- i }
    TGrast <- as.matrix(as.numeric(TGrast))
    colnames(TGrast) <- "Time_Group"

    ## Errors
    if(error_type=='focal'){
        Erast0 <- rnorm(dimT)*ar_factor
        Erast  <- as.matrix(rnorm(dimT) + c( 0, head(Erast0, -1) ))*error_scale
    } else if (error_type=='distance') {
        Erast0 <- mvtnorm::rmvnorm(1, mean=rep(0, nrow(wt_mat)), sigma=wt_mat)
        Erast <- t(Erast0)*error_scale
    } else if (error_type=='spherical') {
        Erast0 <- mvtnorm::rmvnorm(1, mean=rep(0, nrow(Xrast)), sigma=diag( nrow(Xrast)))
        Erast <- t(Erast0)*error_scale
    }
    colnames(Erast) <- 'E'

    ## Dependant Variables
    Yrast <- Xrast %*% as.matrix(theta) + Erast
    colnames(Yrast) <- 'Y'

    ## Raster Stack
    STACK <- cbind(
        xy_mat,
        Srast,
        Trast,
        Yrast,
        Xrast,
        TGrast,
        Erast)

    ## To Data Frame
    DFt <- as.data.frame(STACK)
    return(DFt)
})


##------------------------------------------------------------------
###################
##' Create Multivariate SpaceTime Data
################### 
#'
#' @param dimS, dimT spatial dimension and temporal dimension 
#' @param K number of covariates
#' @param space_groups, time_groups number of contiguous groups (e.g. countries,  time-regimes)
#' @param xpars generate each Xvariable from uniform(xpar[1],xpar[2]) 
#' @param theta parameter vector for Y=X%*%theta
#' @param error_type c('focal', 'distance')
#' @param error_scale variance of errors (normally distributed) 
#' @param xy_mat matrix specifing xy locations
#' @param ws_mat matrix specifing correlation structer of error_type='distance' 
#' @param wt_mat matrix specifing correlation structer of error_type='distance'
#' @param wst_mat=kronecker(wt_mat, ws_mat)
#' 
#' @return dataframe
#' 
#' @details 
#' 'focal' creates errors based on \eqn{e = u + \rho v}
#' 'distance' creates errors base on drawing e from a multivariate normal with covariance matrix wst_mat
#'
#' @examples 
#' DFst <- make_spacetime_data(11,6,2)
#' DFst[DFst$Space_ID==1,]
#'  
#' @export

make_spacetime_data <- compiler::cmpfun( function(
    dimS,
    dimT,
    K,
    xpars=c(0,1),
    theta=1:K,
    space_groups=round(sqrt(dimS)),
    time_groups=round(sqrt(dimT)),
    error_type='distance',
    error_scale=1,
    xy_mat=expand.grid(Loc_X=seq(0,1,length.out=dimS),Loc_Y=seq(1,0,length.out=dimS)),
    ws_mat=exp(-1 * as.matrix( dist( xy_mat )) ),
    wt_mat=toeplitz((2/3)^(0:(dimT-1))),
    wst_mat=kronecker(wt_mat, ws_mat),
    verbose=F)
{

    if(verbose){
        message("space-time lattice: ", dimS, " x ", dimS, " x ", dimT)
    }

    ## Independant Variables
    #Xrast <- mvtnorm::rmvnorm(K,sigma=diag(dimT*dimS*dimS) )
    Xrast <- sapply(1:K, function(k){ runif(dimT*dimS*dimS, min=xpars[1], max=xpars[2]) })
    colnames(Xrast) <- paste0('X', 1:K)

    ## Space ID
    Srast <- kronecker(rep(1,dimT), as.matrix(1:nrow(xy_mat)))
    Srast <- as.matrix(c(Srast))
    colnames(Srast) <- 'Space_ID'

    ## Time ID
    Trast <- kronecker( 1:dimT, c(matrix(1,dimS, dimS)))
    Trast <- as.matrix(c(Trast))
    colnames(Trast) <- 'Time_ID'

    ## Space Groups
    g1 <- cut(xy_mat[,1], space_groups)
    g2 <- cut(xy_mat[,2], space_groups)
    g3 <- as.factor(paste0(g1, g2))
    for(i in seq(levels(g3))){ levels(g3)[i] <- i }
    SGrast0 <- as.matrix(as.numeric(g3))
    SGrast <- kronecker(t(rep(1,dimT)), as.matrix(SGrast0))
    SGrast <- as.matrix(c(SGrast))
    colnames(SGrast) <- paste0("Space_Group")

    ## Time Groups
    TGrast0 <- cut( Trast, time_groups)
    for(i in seq(levels(TGrast0))){ levels(TGrast0)[i] <- i }
    TGrast <- as.matrix(as.numeric(TGrast0))
    #TGrast <- kronecker(TGrast0,  matrix(1,dimS, dimS))
    #TGrast <- as.matrix(c(TGrast))
    colnames(TGrast) <- "Time_Group"

    ## Errors
    if(error_type=='focal'){
        break
    } else if (error_type=='distance') {
        Erast0 <- mvtnorm::rmvnorm(1, mean=rep(0, nrow(wst_mat)), sigma=wst_mat)
        Erast <- t(Erast0)*error_scale
    } else if (error_type=='spherical') {
        Erast0 <- mvtnorm::rmvnorm(1, mean=rep(0, nrow(Xrast)), sigma=diag( nrow(Xrast)))
        Erast <- t(Erast0)*error_scale
    }
    colnames(Erast) <- 'E'

    ## Dependant Variables
    Yrast <- Xrast %*% as.matrix(theta) + Erast
    colnames(Yrast) <- 'Y'

    ## Raster Stack
    STACK <- cbind(
        xy_mat,
        Srast,
        Trast,
        Yrast,
        Xrast,
        SGrast,
        TGrast,
        Erast)

    ## To Data Frame
    DFst <- as.data.frame(STACK)
    return(DFst)

} )



# @param proj.1, proj.2 initial and projected coordinate systems
#    proj.1="+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
#    proj.2="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
#    xy_scale=1E5,
#    if( !is.na(xy_scale)){ DF[,c("x","y")] <-  DF[,c("x","y")]*xy_scale }
#    coords.xy <- sp::SpatialPoints( DF[,c("x","y")],
#        proj4string=sp::CRS(proj.1) )
#    coords.ll <- sp::spTransform(coords.xy, proj.2)
#    DF[,c("x2", "y2")] <- sp::coordinates(coords.ll)


# Create Spatial Data in a Loop
#    if(parallel){
#        RasterList <- parallel::mclapply(1:tf, make_space_data)
#    } else {
#        RasterList <- lapply(1:tf, make_space_data)
#    }
#    DF <- do.call("rbind", RasterList)

