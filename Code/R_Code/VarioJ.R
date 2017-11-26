#-------------------------------------------------------------------
##################
#' Variogram Calculation
################## 
#' @param coords coordinate matrix
#' @param cutt cutoff from which to calculate variogram
#' @param residu vector of values (i.e. OLS residuals) associated coords
#' @param latlon coordinates are lon,lat or x,y
#' @param indices return indices?
#' @param clean unused currently
#' @return data.frame of dij and (ei-ej)^2
#  @examples
#' @export

varioJ <-  compiler::cmpfun( function(
    coords,
    cutt,
    residu,
    latlon=FALSE,
    indices=FALSE,
    clean=FALSE) {

    requireNamespace("spam")
    requireNamespace("spam64")


	#### Distance Matrix Calculation
    message("Distance Matrix...")
	if(latlon){
        cutt_ll <- ((cutt/1000)*360)/(6378.388*2*pi)
		DIST <- spam::spam_rdist.earth(
            coords,
            coords,
            delta=cutt_ll,
            miles=FALSE)
        ##system.time() 20 mins
	 } else {
		DIST <- spam::spam_rdist(coords, coords, delta=cutt)
	}
	#diag(DIST) <- 0

	## Calculate Squared Difference Between Residuals
    message("Residual Differences Squared...")
	DM   <- (DIST<=cutt)*1
	E    <- diag.spam(residu)
	EE   <- ( (DM %*% E)-(E %*% DM) )**2

	#### Sparse Matrix to Data Frame
    message("Sparse Matrix DF...")
	VarioList <- triplet(EE)
	DistList  <- triplet(DIST)$values
	VDList    <- as.data.frame( cbind(VarioList$values, DistList) )
	names(VDList) <- c("Vario", "Dist")

	if(indices) {
        VDList <- as.data.frame( cbind(VDList, VarioList$indices) )
    }

	return(VDList)
} )
 

