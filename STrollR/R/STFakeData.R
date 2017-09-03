#------------------------------------------------------------------
##################
#' Create Space Time Lattice Data
################## 
#'
#' @param n spatial dimension 
#' @param tf temporal dimension
#' @param theta parameter vector for RHS
#' 
#' @return dataframe with (n1,n2,t) coordinates and variables Y,X,X, Country Time
#' 
# @details
# @examples
#'  
#' @export

fake_data_traditional <- compiler::cmpfun( function(
    n=10,
    tf=5,
    theta=c(5,1,1,1))
{

    message("space-time lattice: ",
        n, " x ", n, " x ", tf)

    RasterList <- parallel::mclapply(1:tf, function(t_id) {
    
        ## RHS vars
        Xrast <- round( raster::raster(
            matrix( runif(n*n), ncol=n, nrow=n ) ) *n )
        Zrast <- round( raster::raster(
            matrix( rnorm(n*n), ncol=n, nrow=n ) ) *n )
        erast <- raster::raster(
            matrix( rnbinom(n*n, mu=4, 1), ncol=n, nrow=n ) )
        Grast <- raster::raster(
            matrix( 1:5, ncol=n, nrow=n ) )
        Trast <- raster::raster(
            matrix( t_id, ncol=n, nrow=n ) )
            
        ## LHS vars
        Yrast <-theta[1]*Xrast +
                theta[2]*Grast +
                theta[3]*Zrast +
                theta[4]*Trast +
                erast

        ## Output Raster Stack
        STACK <- raster::stack(
            Yrast,
            Xrast,
            Zrast, 
            Grast,
            Trast)
        names(STACK) <- c("Y", "X", "Z", "Group", "Time")
        
        return(STACK)
    })

    ## Convert RasterStack to DF
    RasterDF <- parallel::mclapply(RasterList, function(si) {
        DFi   <- raster::as.data.frame( si, xy=TRUE)
        DFi$ID <- as.numeric(rownames(DFi))
        DFi
    } )
    DF <- do.call("rbind", RasterDF)

    ## Coordinates Information
    DF[,c("x","y")] <-  DF[,c("x","y")]*1E5
    proj.m <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    proj.w <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    coords.xy <- sp::SpatialPoints( DF[,c("x","y")],
        proj4string=sp::CRS(proj.m) )
    coords.ll <- sp::spTransform(coords.xy, proj.w)
    DF[,c("lon", "lat")] <- sp::coordinates(coords.ll)

    ## Return
    return(DF)
    
} )


