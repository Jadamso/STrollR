#------------------------------------------------------------------
##################
#' Run multiple regressions on the same dataset
################## 
#' 
#' @param FORMS list of regression formula
#' @param ... args passed to mfx_fun
#' @param mfx_fun what type of regression
#' @param mc.cores number of cores if parallel=TRUE
#' @param parallel use parallel processing?
#'
#' @return List of Regressions Summaries
#' 
# @details
# @examples
#'  
#' @export

mfxall <- compiler::cmpfun(function(
    FORMS,
    ...,
    mfx_fun=mfxi,
    mc.cores=as.numeric(system("nproc", intern=TRUE)),
    parallel=TRUE){
    
    if(parallel){
        MFX <- parallel::mclapply(FORMS, mc.cores=mc.cores,
            mfxi, ...)
    } else {
        MFX <- lapply(FORMS, mfxi, ...)
    }
    
    return(MFX)
})


