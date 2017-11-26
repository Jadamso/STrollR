#########################
# Calculating Weights Matrices
#########################

#------------------------------------------------------------------
##################
# Covariance Matrices
##################
#pkgs <- c("data.table", "parallel", "lfe") 
#invisible(sapply(pkgs, require, character.only=TRUE))

#sedir <- path.expand("~/Desktop/Common/R_Code/")

#source( paste0(sedir, "vcovSTsep.R") )
#source( paste0(sedir, "vcovSTsep_loop.R") )
#source( paste0(sedir, "/Experimental/vcovSTkron.R") )
#message("vcovSTkron")

## Wrapper
#source( paste0(sedir, "SpaceTimeSE_format.R") )

#------------------------------------------------------------------
##################
# Future Considerations
##################

#sourceCpp( paste0(sedir, "SparseMatMult.cpp") )
#lltLS <- cppFunction( "")

#http://casualinference.org/2016/07/07/big_fast_lm/
#RcppEigen::fastLm
#RcppArmadillo::fastLm
#RcppGSL::fastLm

#------------------------------------------------------------------
##################
# Sources
##################
# http://freigeist.devmag.net/economics/936-spatial-temporal-autocorrelation-correction-with-r-conley.html
# http://www.trfetzer.com/using-r-to-estimate-spatial-hac-errors-per-conley/
# http://darinchristensen.github.io/resources/
# https://dl.dropboxusercontent.com/u/16041354/ConleySEs.zip

# http://www.sciencedirect.com.libproxy.clemson.edu/science/article/pii/S0304407698000840
# https://www.jstor.org/stable/pdf/2648817.pdf
# Solomon Hsiang (2010)


# https://gist.github.com/fabian-s/9665846#file-2014-03-20-dynamic-dispatch-for-sparse-matrices-rmd
# http://www.stata.com/meeting/sandiego12/materials/sd12_jeanty.pdf


#------------------------------------------------------------------
##################
#' calculate vcovSTsep for felm object
##################
#'
#' @param vcovfun function for which type of SHAC correction
#' @param reg an lm or felm object
#' @param DF data.frame with unit, time, sp_coords
#' @param unit cellular id
#' @param time time id
#' @param sp_coords coordinate id
#' @param t_cutoff cutoff for considering time correlation
#' @param d_cutoff cutoff for considering space correlation
#' @param wmat spatial weights matrix
#' @param tmat temporal weights matrix
#' @param latlon are coordinates in lon,lat or x,y form
#' @param convert_to_angles convert lon,lat to x,y
#' @param verbose show messages
#' @param SPLIT_sp,SPLIT_t,SPLIT see vcovST.format
#'
#' @return covariance matrix w/ SHAC correction
#
#  @examples
#
#' @export

vcovST <- compiler::cmpfun( function(
    vcovfun=vcovSTsep,
    reg,
    DF,
    unit,
    time,
    sp_coords,
    t_cutoff=NA,
    d_cutoff=NA,
    wmat=NA,
    tmat=NA,
    latlon=TRUE,
    convert_to_angles=TRUE,
    verbose=FALSE,
    SPLIT_sp=FALSE,
    SPLIT_t=FALSE,
    SPLIT=FALSE
    ){

    if(verbose) {message("Formatting Data")}

    DAT <- vcovST.format(
        reg=reg,
        DF=DF, 
        unit=unit,
        time=time, 
        sp_coords=sp_coords,
        SPLIT_sp=SPLIT_sp,
        SPLIT_t=SPLIT_t,
        SPLIT=SPLIT
    )

    if(verbose) {message("Creating LISTJ")}


    LISTJ <- listj(
        ddff=DAT$dat,
        wmat=wmat,
        tmat=tmat,
        t_cutoff=t_cutoff,
        d_cutoff=d_cutoff,
        latlon=latlon,
        convert_to_angles=convert_to_angles
    )

    if(verbose) {message("Calculating SEs")}
    SEs <-vcovfun(
        DAT=DAT,
        LISTJ=LISTJ,
        verbose=FALSE
    )

    return(SEs)

} )


#------------------------------------------------------------------
##################
#' vcovST for seperable covariances using split approach
##################
#'
#' @describeIn vcovST
#'
#' @param vcvtime add temporal clustering
#' @param cores number of cores to use
#'
#' @details if vcvtime==FALSE, wmat should have diagonal elements
#'
#' @export

# combine rolling method with parallel loop method

vcovST.loop <- compiler::cmpfun( function(
    vcovfun=vcovSpace.loop,
    reg,
    DF,
    unit,
    time,
    sp_coords,
    t_cutoff=NA,
    d_cutoff=NA,
    wmat=NA,
    tmat=NA,
    latlon=TRUE,
    convert_to_angles=TRUE,
    verbose=FALSE,
    cores=4,
    vcvtime=TRUE
    ){

    if(vcvtime){
        VCOV_Time <- reg$clustervcv
        #cluster.vcov(reg, ALL$ID)
        #vcovCL(reg, cluster=DF$ID)
    } else {
        VCOV_Time <- reg$vcv
        #vcovHC(reg, type="HC3")
    }

    if(verbose) {message("Formatting Data")}

    if( !is.null(reg$na.action) ){
        DF <- DF[-reg$na.action, ]
    } 

    DAT <- vcovST.format(
        reg=reg,
        DF=DF[, c(unit, time, sp_coords)], 
        unit=unit,
        time=time,
        sp_coords=sp_coords,
        SPLIT_t=TRUE
    )

    ## vcovSpace.loop doesn't require LISTJ
    #if(verbose) {message("Creating LISTJ")}
    #LISTJ <- listj(
    #    ddff=DAT$dat,
    #    wmat=wmat,
    #    tmat=tmat,
    #    t_cutoff=t_cutoff,
    #    d_cutoff=d_cutoff,
    #    latlon=latlon,
    #    convert_to_angles=convert_to_angles)
    
    if(verbose) {message("Calculating SEs")}

    VCOV_Space <- vcovfun(
        DAT=DAT,
        wmat=wmat,
        verbose=verbose,
        cores=cores
    )


    VCV <- (VCOV_Space + VCOV_Time)
    return(VCV)

} )

#------------------------------------------------------------------
##################
# vcovST.loop for space only
##################
#' @rdname vcovST.loop 

vcovS.loop <- compiler::cmpfun( function(
    vcovfun=vcovSTsep.loop,
    reg,
    DF,
    unit,
    time,
    sp_coords,
    wmat=NA,
    verbose=FALSE,
    cores=4)
{
    
    vcovST.loop(
        vcovfun,
        reg,
        DF,
        unit,
        time,
        sp_coords,
        wmat,
        verbose,
        cores,
        vcvtime=FALSE)
})



