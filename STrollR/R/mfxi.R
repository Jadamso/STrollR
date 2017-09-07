#------------------------------------------------------------------
##################
#' Run a Regression
################## 
#' 
#' @param formi regression formula
#' @param datai data for regression
#' @param scl vcov correct for spatial+temporal covariance?
#' @param wmat0 spatial weights matrix passed to vcovSCL
#' @param unit_id,time_id,coord_id passed to vcovSCL
#'
#' @return summary table
#' 
# @details
# @examples
#'  
#' @export

#library(lmtest) #library(lfe) #library(sandwich)#library(stargazer)
#source("~/Desktop/Common/R_Code/SpaceTimeSE.R") # vcovSCL
#source("~/Desktop/Common/R_Code/TableMaker_trim.R")


mfxi <- compiler::cmpfun(function(
    formi,
    datai,
    scl=TRUE,
    wmat0=get("WMAT0"),
    unit_id="ID",
    time_id="Year",
    coord_id=c("x","y")
){

    requireNamespace("PrettyR")
    requireNamespace("STrollR")
    
        ## Regression
        reg <- lfe::felm(
            formula(formi),
            data=datai,
            keepCX=TRUE)
            
        ## Standard Errors
        if(scl){
            VCOV <- vcovST.loop(
                reg=reg,
                unit=unit_id,
                time=time_id,
                sp_coords=coord_id,
                wmat=wmat0,
                DF=datai[,c(unit_id, time_id, coord_id)],
                cores=1)
        } else {
            VCOV <- vcov(reg)
        }
        
        ## Coefficient and Standard Errors
        b.s  <- lmtest::coeftest(reg, vcov=VCOV )
        return(b.s)
})
    
