#-------------------------------------------------------------------
##################
#' Variogram Calculation
################## 
#' @param coords coordinate matrix
#' @param cuttoff cutoff from which to calculate variogram
#' @param E vector of values (i.e. OLS residuals) associated coords
#' @param latlon coordinates are lon,lat or x,y
#' @param indices return indices?
#' @param clean unused currently
#' @param verbose 
#' @return data.frame of dij and (ei-ej)^2
#  @examples
#' @export

varioJ <-  compiler::cmpfun( function(
    coords,
    cuttoff,
    E,
    latlon=FALSE,
    indices=FALSE,
    clean=FALSE, 
    verbose=FALSE) {

    requireNamespace("spam")
    requireNamespace("spam64")


	#### Distance Matrix Calculation
    if(verbose){ message("Distance Matrix...")}
	if(latlon){
        cutt_ll <- ((cuttoff/1000)*360)/(6378.388*2*pi)
		DIST <- spam::spam_rdist.earth(
            coords,
            coords,
            delta=cutt_ll,
            miles=FALSE)
        ##system.time() 20 mins
	 } else {
		DIST <- spam::spam_rdist(coords, coords, delta=cuttoff)
	}
	#diag(DIST) <- 0

	## Calculate Squared Difference Between Residuals
   if(verbose){message("Residual Differences Squared...")}
	DME  <-  1*(DIST<=cuttoff) %*% spam::diag.spam(E)
	EE   <- ( DME- t(DME) )^2
    
	#### Sparse Matrix to Data Frame
    if(verbose){message("Sparse Matrix DF...")}
	if(indices) {
    	VDList <- as.data.frame( cbind(
	        Vario=spam::triplet(EE)$values,
	        Dist=spam::triplet(DIST)$values,
	        VarioList$indices) )
    } else {
    	VDList <- as.data.frame( cbind(
	        Vario=spam::triplet(EE)$values,
	        Dist=spam::triplet(DIST)$values ) )
    }

	return(VDList)
} )
 

