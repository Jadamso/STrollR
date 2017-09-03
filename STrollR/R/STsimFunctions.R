#------------------------------------------------------------------
##################
#' cleanup spam::mvtnorm
################## 
#' 
#' @param m matrix of simulations (each row a realization)
#' @param nsim number of simulations
#' @param xyt lattice coordinates
#' 
#' @return dataframe
#' 
# @examples
#' @export

spt_cleanup <-  compiler::cmpfun( function(m, nsim, xyt){
    m <- t(m)
    colnames(m) <- paste0("Sim", 1:nsim)
    DF <- as.data.frame( cbind(xyt, m ) )
    return(DF)
})

#------------------------------------------------------------------
##################
#' Convert Dataframe with 1 variable to raster for one realization
################## 
#' 
#' @param df_i dataframe
#' @param sim_i which simulation
#' 
#' @return raster
#' 
# @examples
#' @export

var2stack <- compiler::cmpfun( function(df_i, sim_i) {
    coordinates(df_i) <- ~ x + y
    gridded(df_i) <- TRUE
    rast <- raster::raster( df_i[sim_i] )
    return(rast)
})

#------------------------------------------------------------------
##################
#' Convert list of dataframes with to rasterstack
################## 
#' 
#' @param DF list of dataframe
#' @param sim_i which simulation
#' 
#' @return rasterstack
#' 
# @examples
#' @export

df2stack <- compiler::cmpfun( function(sim_i, DF){

    DFi <- DF[ c("Time","x","y",sim_i)]
    DFit <- split(DFi, DFi$Time)
    
    STACKi <- lapply(DFit, var2stack, sim_i=sim_i)
    names(STACKi) <- paste0("Time", 1:ti)
    STACKi <- raster::stack(STACKi)
    
    return(STACKi)

})

#------------------------------------------------------------------
##################
#' Convert simulation to rasterstack
################## 
#' 
#' @param e_spt matrix of draws from spam.mvtnorm
#' @param number of simulations
#' @param xyt lattice structure
#' 
#' @return list of rasterstacks
#' 
# @examples
#' @export

sim2stack <- compiler::cmpfun( function(e_spt, nsim, xyt){

    DF <- spt_cleanup(e_spt, nsim, xyt)

    DFlist <- parallel::mclapply(
        paste0("Sim", 1:nsim), 
        df2stack,
        DF)
    
    return(DFlist)
})

#------------------------------------------------------------------
##################
# Old Functions
################## 

## Convert to STSDF
#DF      <- DF[order(DF$Year),] 
#stDATA  <- DF[order(DF$ID),] 
## Convert to Dates
#DFi     <- stDATA[stDATA$ID==stDATA$ID[1], ]
#stfTIME <- as.POSIXct( paste0(DFi$Year, "-01", "-01") )
## Convert to Space
#DFy     <- stDATA[ stDATA$Time==1, ]
#stfSP   <- SpatialPoints(DFy[,c("x","y")])
## Convert to STFDF
#stfDF  <- STFDF(sp=stfSP, time=stfTIME, data=stDATA)




#------------------------------------------------------------------
##################
# Create Gif Plots
################## 
#' 
#' @rdname mkGif
#' @export

DynPlot <- compiler::cmpfun( function(
    DFlist,
    ti,
    pname="STvarX",
    ind=1){

    for (i in 1:ti) {

        pdfname <- paste0(fdir, pname, i, ".pdf")    
        
        pdf(pdfname, width=12, height=12)
        par( mar=c(0,0,1,0), bg="white")
        
        image( DFlist[[ind]][[i]],
            zlim=val_range,
            col=rev(rainbow(100)))

        title( paste0("Time ", i))   
        dev.off()
    }
})



#' @rdname mkGif
#' @export
DynGif <- compiler::cmpfun( function( pname, vw=FALSE){

    mkGif <- paste0('convert -delay 100 ',
        '-size 1000x1000 ',
        '-loop 0 ',
        '~/Desktop/Packages/STrollR/STsim/', pname,'*.pdf ',
        '~/Desktop/Packages/STrollR/STsim/', pname,'.gif')
        
    clGif <- paste0('rm ~/Desktop/Packages/STrollR/STsim/',
        pname, '*.pdf')
        
    vwGif <- paste0('eog ~/Desktop/Packages/STrollR/STsim/',
        pname,'.gif')
        
    system(mkGif)  
    system(clGif)
    if(vw) { system(vwGif)  } 
} )


#------------------------------------------------------------------
##################
# Create Gif Plots
################## 
#' 
#' @param DFlist 
#' @param ti number of time periods
#' @param pname name of file
#' @param ind which simulation
#' @param vw view output
#' 
#' @return list of rasterstacks
#' 
# @examples
#' @export

mkGif <-  compiler::cmpfun( function(
    DFlist,
    ti,
    pname="STvarX",
    ind=1,
    vw=FALSE){
    
    DynPlot(DFlist, ti, pname, ind)
    DynGif(pname, vw)
})

