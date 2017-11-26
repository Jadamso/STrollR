
#------------------------------------------------------------------
##################
# Create Gif Plots
################## 
#
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
#' Create Gif Plots
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

