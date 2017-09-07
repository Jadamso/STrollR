#------------------------------------------------------------------
##################
#' Run a Regression 
################## 
#'
#' @param formi regression formula
#' @param datai data for regression
#' @param vcv type of covariance correction
#' 
#' @return summary table
#' 
# @details
# @examples
#'  
#' @export

mfxi.lm <- compiler::cmpfun(function(
    formi,
    datai,
    vcv="standard"
){
        ## Regression
        REG <- lm( formula(formi), data=datai)
            
        ## Standard Errors
        if(vcv=="standard"){
            VCOV <- vcov(reg)
        } else {
            message("unsupported")
        }
        
        ## Coefficient and Standard Errors
        SUMMARY <- lmtest::coeftest(REG, vcov=VCOV )
        return(SUMMARY)
})
