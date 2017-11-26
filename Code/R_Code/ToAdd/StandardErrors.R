#### Jordan  Adamson ####
#########################
#       R Commands      #
#########################
#source("~/Desktop/Common/R_Code/StandardErrors.R")

library(sandwich) ## vcovHAC( lmobject..,order.by )
library(lmtest)
library(compiler)

message("heteroskedastic: sandwich::vcovHC")
message("heteroskedastic autocorrelation: sandwich::vcovHAC")
message("panel heteroskedastic autocorrelation: plm::vcovHC")
#message("heteroskedastic clustering: vcovCL")
#message("heteroskedastic clustering temporal autocorrelation: multiwaycov::cluster.vcov")
message("spatial and temporal autocorrelation: vcovST")

message("'~/Desktop/Common/R_Code/StandardErrors.R'")
#?plm:::vcovHC.plm
#?plm::vcovNW

#---------------------------------------------------------------------------
## Plot Standard Errors


## Polygons
#reg <- lm(y~x)
#D   <- cbind(x, reg$fit-1*reg$s, Loregfit+1*reg$s)
#polygon( c( D[,1],rev(D[,1]) ), c( D[,2],rev(D[,3]) ), border=NA)

## Loess with Defaul Standard Error Lines and Polygons
#Lo   <- predict(loess(y ~ x, data, span=smoother, degree=deg), se=TRUE)
    # to predict for points along the x-axis set
    # data <- data.frame(x = seq(a,b,length) )
#lines(Lo$fit$x, Lo$fit$up.lim); lines(Lo$fit$x, Lo$fit$low.lim)

    
#---------------------------------------------------------------------------
##################
# Bootstrap SE 
################## 

source("~/Desktop/Common/R_Code/BootSE.R")
source("~/Desktop/Common/R_Code/SpaghettiPlot.R")

#-------------------------------------------------------------------
##################
# Clustered SE 
################## 

source("~/Desktop/Common/R_Code/ClusterSE.R")

# plm::vcovHC(..., cluster = "group")
#plm::vcovHC(reg, type = "HC1", cluster = "ID")


# clubSandwich::vcovCR.lm


## library(clusterSEs)
## library(clubSandwich)
#-------------------------------------------------------------------
##################
# Spatial HAC 
################## 

source("~/Desktop/Common/R_Code/SpaceTimeSE.R")

## Other Conoley post-estimation options:
## splm::vcov.splm

## Driscoal and Kray correction
# plm::vcovSCC


## Other
# nlme::ls( formula, correlation = corGaus(form = ~x+y, nugget = TRUE) + corAR1(0.8, form = ~t), data =dat )
# lme4::glmer.nb
## spdep::ME 



