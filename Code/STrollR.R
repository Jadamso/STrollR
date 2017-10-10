#-------------------------------------------------------------------
##################
# Setup 
##################
library(devtools)
library(roxygen2)
library(MiscUtils)

hmdir <- path.expand("~/Desktop/Common/R_Code/")
pmdir <- path.expand("~/Desktop/Packages/")

#-------------------------------------------------------------------
##################
# Make 
##################
pack  <- "STrollR"
pdir  <- paste0(pmdir, pack,"/")
packg <- paste0(pdir, pack)


Version <- numeric_version("0.1.3")

# Setup R Package
source(paste0(pdir,"Code/PackageSetup.R") )

# Create R Package Contents
source(paste0(pdir,"Code/CodeSetup.R") )

pack_up(pdir)

message(" Name Changes from Paper2:
    vcovSTsep --> vcovST
    vcovSCL --> vcovST.loop")
    
#-------------------------------------------------------------------
##################
# Install 
##################
devtools::install(packg) ## Locally Works

devtools::install_github( paste0("Jadamso/",pack), subdir=pack)
## Public Package From Github Fails Often

citation(pack)

print("Done")

## importFrom(data.table, parallel)

## source("~/Desktop/Packages/STrollR/Code/STrollR.R")

## R CMD BATCH Code/STrollR.R && rm STrollR.Rout .Rdata
