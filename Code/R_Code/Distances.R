# Compare Speed and Accuracy of distance calculation
#------------------------------------------------------------------
##################
# Load Packages
################## 

for(pack in c("fields", "sp", "geosphere", "spdep", "spam")) {
    library(pack, character.only=TRUE) }


#------------------------------------------------------------------
##################
# Distance Calculation Functions
################## 

dfun0   <-  dist

## fields::rdist.earth
dfun1 <- function(fordist){ message("rdist.earth")
    distances <- fields::rdist.earth(fordist, miles=FALSE)
    diag(distances) <- 0
    distances
}
dfun1 <- cmpfun(dfun1)

## sp::spDists
dfun2 <- function(fordist){ message("spDists")
    distances <- sp::spDists(fordist, longlat=TRUE)
    diag(distances) <- 0
    distances
}
dfun2 <- cmpfun(dfun2)

## geosphere::distHaversine
dfun3 <- function(fordist){ message("distHaversine")
    distances <- lapply(1:nrow(fordist), function(k) geosphere::distHaversine(fordist[k,],as.matrix(fordist)))
    distances <- do.call(cbind, distances)/1000
    diag(distances) <- 0
    distances
}
dfun3 <- cmpfun(dfun3)

## spam::spam_rdist.earth & spam::nearest.dist, Fortran, Fast
dfun4 <- function(fordist, delta1=180 ){ message("spam_rdist.earth")
    # to restrict to distances smaller than delta.km
    # one has to specify delta=delta.km*(360/(6378.388*2*pi))
    distances <- spam::spam_rdist.earth(fordist, fordist, miles=FALSE, delta=180)
    diag(distances) <- 0
    as.matrix(distances) 
}
dfun4 <- cmpfun(dfun4)

## spdep::nbdist
dfun5 <- function(fordist, cutoff=40075){ message("nbdists") # 40,075 km is earth circumference
    coord_nb  <- spdep::dnearneigh(fordist, d1=0, d2=cutoff, longlat=TRUE) 
    distances <- spdep::nbdists(coord_nb, fordist, longlat=TRUE)
    nld       <- length(distances)
    distances <- lapply(1:nld, function(r) {
        dlist <- distances[[r]]
        dlist <- c(dlist[1:r-1], 0, dlist[r:nld] )
        as.vector(na.omit(dlist))
        } )
    distances <- do.call(cbind, distances)
    diag(distances) <- 0
    distances
}
dfun5 <- cmpfun(dfun4)


## distances for knearest neighbours
dfun6 <- function(fordist){ library(spdep); message("B")
    #coord_nb  <- dnearneigh(fordist, d1=0, d2=cutoff, longlat=TRUE)
    coord_nb  <- knn2nb(knearneigh(fordist, k=(nrow(fordist)-1), longlat=TRUE))
    ## Spatial Weights Matrix
    listW     <- nb2listw(coord_nb, style="B")
    W         <- as_dgRMatrix_listw(listW) 
    WW        <- as(W, "CsparseMatrix")
    #trMat   <- trW(W, type="mult", p=32)
}
dfun6 <- cmpfun(dfun6)




#------------------------------------------------------------------
##################
# Unused Distances
################## 

## SpatialTools::dist1(x)  ;; slow
## fossil::earth.dist
## spdep::autocov_dist
## raster::distanceFromPoints, only raster data can use


#------------------------------------------------------------------
##################
# Unused Functions
################## 


# Function to be called by main loop
abt.fun <- cmpfun( function(e,r) {  r*acos( ifelse( abs(e)>1, 1*sign(e), e ) ) } )


rdist.earth.big <- function (x1, miles = TRUE, R = NULL) {
    # Miles or KiloMeters
    if (is.null(R)) {
        if (miles) 
             R <- 3963.34
        else R <- 6378.388
    }
    # Calculations
    coslat1 <- ( cos((x1[, 2] * pi)/180) )
    sinlat1 <- ( sin((x1[, 2] * pi)/180) )
    coslon1 <- ( cos((x1[, 1] * pi)/180) )
    sinlon1 <- ( sin((x1[, 1] * pi)/180) )
    # Memory Eaters
    pp  <- cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1) %*%
           t(cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1))
    gc()
    for(i in 1:nrow(pp) ){ 
        pp[i,] <- abt.fun(pp[i,], R) #apply( as.matrix(PP[i,]), 1, abt.fun ) 
        #   if(i%in%SQ==TRUE) gc()
    }
    gc()
    return( pp )
}
#rdist.earth.big <- cmpfun(rdist.earth.big)
