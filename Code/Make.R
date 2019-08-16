#-------------------------------------------------------------------
##################
# Setup 
##################
library(devtools)
library(roxygen2)
library(MiscUtils)

 
pmdir <- path.expand("~/Desktop/Packages/")

#-------------------------------------------------------------------
##################
# Make 
##################
pack  <- "STrollR"
pdir  <- paste0(pmdir, pack,"/")
packg <- paste0(pdir, pack)


Version <- numeric_version("0.1.6")
#usethis::use_version()

# Setup R Package
source(paste0(pdir,"Code/PackageSetup.R") )

# Create R Package Contents
source(paste0(pdir,"Code/CodeSetup.R") )

message(" Name Changes from Paper2:
    vcovSTsep --> vcovST
    vcovSCL --> vcovST.loop")

devtools::build(pkg=packg, manual=TRUE)

system( paste0('R CMD Rd2pdf ', packg) )
#devtools::check_built(manual=TRUE, path=packg)

#-------------------------------------------------------------------
##################
# Upload and Install 
##################

MiscUtils::pack_up(pdir)
## https://www.rstudio.com/wp-content/uploads/2015/03/devtools-cheatsheet.pdf

devtools::install(packg) ## Locally Works

devtools::install_github( paste0("Jadamso/",pack), subdir=pack)

citation(pack)

print("Done")

## R CMD BATCH --no-save Code/Make.R Code/Make.Rout
