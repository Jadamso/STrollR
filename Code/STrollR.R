
pack  <- "STrollR"
pdir  <- paste0(path.expand("~/Desktop/Packages/"),pack,"/")
packg <- paste0(pdir, pack)


Version <- numeric_version("0.1.2")

#-------------------------------------------------------------------
##################
# Setup R Package
################## 
source(paste0(pdir,"Code/PackageSetup.R") )

#-------------------------------------------------------------------
##################
# Create R Package Contents
################## 
source(paste0(pdir,"Code/CodeSetup.R") )


message(" Name Changes from Paper2:
    vcovSTsep --> vcovST
    vcovSCL --> vcovST.loop")


#-------------------------------------------------------------------
##################
# Make and Upload R Package
################## 
source(paste0(pdir,"Code/PackageMake.R") )


## importFrom(data.table, parallel)

## source("~/Desktop/Packages/STrollR/Code/STrollR.R")

## R CMD BATCH STrollR.R
