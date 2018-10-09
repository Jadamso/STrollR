#########################
# Calculating Weights Matrices
#########################
message("SpatialWeights")
library(spdep)
library(sphet)

## Highly Suggest library(spam); library(spam64)

#------------------------------------------------------------------
##################
#' K nearest neighbours
################## 
## Nearest Neighbours & Spatial Weights Matrix
#' Calculate the number of neighbours within a neighbourhood.
#' 
#' @param w number of neighbours wide (east-west).
#' @param h number of neighbours long (north-south).
#' @param type type of neighbourhood; "Moore" or "VonNeumann"
#' @return the number of nearest neighbours
#' @examples KNN( 4 )
#' @export


KNN <- function(w,h=w, type="Moore") {
    if(type=="Moore"){
        neighs <- (2*w+1)*(2*h+1)-1
    } else if(type=="VonNeumann") {
        # neighs <- (2*w+1)*(2*h+1)-1
    }
}
KNN <- compiler::cmpfun(KNN)

#------------------------------------------------------------------
##################
#' Write a listw object as a GWT file
##################
#' 
#' @param listw returned from spdep::listw
#' @param dgwt_outfile name of file 
#' @return dgwt_outfile
#  @examples 
#    coord_mat <- expand.grid( list(x=1:10, y=1:10))
#    listw_dist <- mat2listw( dist(coord_mat) )
#    listw2gwt( listw_dist ) )
#' @export

write.listw2gwt <- function(
    listw,
    dgwt_outfile=paste0( tempdir(),"dgwt.GWT") )
{
    dgwt  <- spdep::listw2sn(listw) 
    #distW <- as(dgwt, "distance")
    spdep::write.sn2gwt(dgwt, dgwt_outfile)
    return(dgwt_outfile)
}
write.listw2gwt <- compiler::cmpfun(write.listw2gwt)

#------------------------------------------------------------------
##################
#' Calculate the weights objects used in {spdep} {sphet} 
##################
#' @param coord_sp matrix of coordinates or a SpatialPoints object
#' @param neigh number of neighbours to use in calculation
#' @param knn calculate weights using knn approach
#' @param adj calculate vonneumann weights (see VonNeumann)
#' @param dnn dnn approach unsupported
#' @param rast raster approach unsupported
#' @param vario is coord_sp a weights matrix? 
#' @param sphet create objects used in sphet?
#' @param tracer create trace matrix objects?
#' @param tr_type type of trace matrix
#' @param tr_m trace matrix m
#' @param tr_p trace matrix p
#' @param symm make weights symmetric
#' @param symm_check check for symmetric weights matrix
#' @param SAVE filename to save to, NA <default> returns as object
#'
#' @return filename of saved objects, or returns objects if SAVE=NA
#  @examples 
#    coord_mat <- sp::SpatialPoints(expand.grid( list(x=1:10, y=1:10)))
#    gridded(coord_mat) <- TRUE
#    NEIGH(coord_sp, neigh=2, sphet=FALSE, trace=FALSE)
#'
#' @export

NEIGH <- compiler::cmpfun( function(
    coord_sp,
    neigh=1,
    knn=TRUE,
    adj=FALSE, 
    dnn=FALSE,
    rast=FALSE,
    vario=FALSE,
    sphet=FALSE,
    tracer=TRUE,
    tr_type="mult", #"moments""
    tr_m=20,
    tr_p=16,
    symm=TRUE,
    symm_check=TRUE,
    SAVE=NA){
    
    requireNamespace("spdep")
    requireNamespace("sphet")

    ##############################################
    # Create listw
    
    ## Nearest Neighbour Weights
    if(knn){
        neigh_i <- KNN(neigh)
        message("Nearest Neighbours: k=", neigh_i)

        requireNamespace("RANN")
        system.time({
            coord_nn <- knearneigh(
                coord_sp, k=neigh_i, RANN=TRUE)
        })
        ## 18mins

        coord_nb <- knn2nb( coord_nn )
        listw <- nb2listw( coord_nb )
        ## 1min

    } else if(adj) {
        neigh_i <- 4
        message("Rook Neighbours: directions=", neigh_i)

        wmatspars <- VonNeumann(coord_sp)
        system.time({
            listw <- mat2listw( as.matrix(wmatspars)) 
        })
        #system.time( coord_nb <- tri2nb(nmat) )

    } else if(dnn){
        message("DNN unsupported")
        break
        ## coord_spc <- coordinates(coord_sp)
        ## coord_dn <- dnearneigh(coord_spc, 0, dst)
    } else if(rast){
        message("rast unsupported")
        break
        ## cell2nb(7, 7, type="rook", torus=TRUE)

    } else if(vario){
        message("coord_sp must be dgCMatrix of weights")
        neigh_i <- "vario"
        coord_nb <- NA
        system.time({
            listw <- mat2listw( coord_sp )
            # Default: unknown style="M"
        })
    }



    ##############################################
    # Create dgwt

    if( sphet | (vario & symm) ){
        message("Writing GWT")
        system.time({ 
            dgwt_outfile <- paste0(tempdir(),"dgwt",neigh,".GWT")
            write.listw2gwt(listw, dgwt_outfile)
        })
    }


    ##############################################
    # Check and create symmetric weights

    if( vario & symm ){
        message("Reading GWT for Symmetry")
        symm_check <- FALSE
        coord_nb <- read.gwt2nb(dgwt_outfile)
    }

    if( symm_check ) {
        issymm <- is.symmetric.nb(coord_nb, force=TRUE)
        message(paste0("Is Symmetric: ", issymm ) )
    }

    ## Make Symmetric
    if( symm ) {
        message("Make Symmetric")
        coord_nb <- make.sym.nb(coord_nb)
        #1min
        listw <- nb2listw(coord_nb)
        #10seconds
    }



    ##############################################
    # Create Trace Matrices

    ## Trace Matrices
    if( tracer ){
        message("Make Trace")
        W <- as(as_dgRMatrix_listw(listw),
                "CsparseMatrix")
        system.time({
            trMatb <- trW(W,
                listw=listw,
                type=tr_type,
                m=tr_m,
                p=tr_p)
        })
        # 15mins (type="mult",m=20,p=16)
    } else {
        trMatb <- W <- NA
    }

    ## Sparse Distance Matrix  For SPHET HAC estimation
    if( sphet ){ 
        message("Make Sparse")
        distW <- read.gwt2dist(file=dgwt_outfile)
        # 20 seconds
    } else {
        distW <- NA
    }

    ##############################################
    # Saving

    if( is.na(SAVE) ){
        return( list(coord_nb, listw, W, trMatb, distW) )
    } else {
        message("Saving")
        save( coord_nb, listw, W, trMatb, distW,
            file=paste0(SAVE,"NB",neigh_i,".Rdata") )
    }
})



#------------------------------------------------------------------
##################
#' Weighting Kernel
##################
#' 
#' @param d returned from spdep::listw
#' @param dmax name of file 
#' @return bartlett weight
#  @examples coord_mat <- expand.grid( list(x=1:10, y=1:10))
#' @export

bartlettSparse <- compiler::cmpfun( function(d,dmax){1 - d/dmax })
    # sandwich::kweights(.,kernel="Bartlett")
    # stats::kernel( )

#------------------------------------------------------------------
##################
#' Compute Sparse Spatial Weights Matrix
##################
#' 
#' @param M matrix of coordinates
#' @param cutoff use distances up cutoff
#' @param latlon are the rows (lat,lon) coordinates?
#' @param convert_to_angles convert cutoff from km to angles?
#' @return the number of nearest neighbours
#' @examples weight_mat(expand.grid( list(x=1:10, y=1:10)), cutoff=.5)
#' @export


weight_mat <- compiler::cmpfun( function(
    M,
    cutoff,
    latlon=NA,
    convert_to_angles=TRUE)
{

    requireNamespace("spam")
    requireNamespace("spam64")
    
    if( is.na(latlon) | (latlon==FALSE) ) {
        message("Euclidean distances in km")

        WMAT <- spam::spam_rdist(M,M, delta=cutoff)

        bart_cut <- cutoff

    } else {

        if(latlon=="SH") {

            message("Not Dooable for Sparse Matrix")
            message("try DistMat(.,'SH')")
            # (distance_1,2)^2 = 
            # 111*(lat1 - lat2)**2) + 
            # (cos(lat1*3.141593/180)*111*(long1-long2)**2

        } else {

            a_scale  <- (360/(6378.388*2*pi)) 

            if(convert_to_angles){
                message("converts cutoff from km to angles")
                a_cutoff <- cutoff*a_scale
                bart_cut <- cutoff

            } else {
                message("Assumes cutoff is in degrees, i.e.")
                message("     delta=delta.km*360/(6378.388*2*pi)")
                a_cutoff <- cutoff
                bart_cut <- a_cutoff/a_scale
            }

            WMAT <- spam::spam_rdist.earth(M,M,
                delta=a_cutoff, miles=FALSE)

            ## Cleanup
            #diag(WMAT) <- 0
            #WMAT <- spam::cleanup(WMAT)
        }
    }



    ## Convert to Matric class used by vcovST
    WMAT     <- spam::as.dgCMatrix.spam( WMAT )
    ## https://github.com/rstats-gsoc/gsoc2017/wiki/Sparse-matrix-automatic-conversion-in-RcppArmadillo

    ## If class(WMAT) =="spam", use WMAT@entries 

    bartlettSparse <- cmpfun( function(d,dmax){1 - d/dmax } )
    WMAT@x  <- bartlettSparse( WMAT@x, bart_cut)

    ## Clean
    diag(WMAT) <- 1
    WMAT <- WMAT * (WMAT > .Machine$double.eps)
    WMAT <- Matrix::drop0(WMAT)

    return(WMAT)
})



#------------------------------------------------------------------
##################
#' Compute VonNeumann Neighbours
##################
#' 
#' @param coord_sp SpatialPoints object or coordinate-matrix
#' @param directions see adjacent
#' @return sparse weights matrix
# @examples VonNeumann( sp::SpatialPoints(expand.grid( list(x=1:10, y=1:10))) )
#' @export


VonNeumann <- function(coord_sp, directions=4){
        
        if( class(coord_sp) %in% c('data.frame','matrix')){
            message('converting to SpatialPoints')
            coord_sp <- sp::SpatialPoints( coord_sp )
        }
        sp::gridded(coord_sp) <- TRUE
        rast <- raster(coord_sp)

        vv   <- coord_sp@grid.index

        nmat <- adjacent(rast,
            cells=vv,
            sorted=T, id=TRUE,
            directions=directions)

        nmat_1 <- tapply(nmat[,"from"], nmat[,"id"], unique)

        nmat_2 <- match(nmat[,"to"], nmat_1)


        nmat_id1 <- c(nmat[,"id"])

        nmat_id2 <- as.numeric(rownames( nmat_1[nmat_2] ) )

        nmat_id <- na.omit( data.frame(
            from=nmat_id1,
            to=nmat_id2
        ) )

        wmatsparse <- Matrix::sparseMatrix(
            i=nmat_id[,"from"],
            j=nmat_id[,"to"],
            x=rep(1,nrow(nmat_id)))

        #wmatsparse2 <- spa,::as.spam(
        #    list(indices=nmat[,c("from", "to")],
        #    matrix(rep(1), nrow(nmat))) )
        
        return(wmatsparse)
}
VonNeumann <- compiler::cmpfun(VonNeumann)

## Note that this creates VonNeumann Neighbours (SpatW)
## To spatially weight the residuals according to correlatio rho, need to use
## solve(I - rho*SpatW) %*% solve(I - rho*t(SpatW) ) 

## See  http://www.sciencedirect.com.libproxy.clemson.edu/science/article/pii/S0304407698000840
# https://www.jstor.org/stable/pdf/2648817.pdf
