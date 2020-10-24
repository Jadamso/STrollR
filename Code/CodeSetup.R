#-------------------------------------------------------------------
##################
# Libraries
##################


## Imported Packages: (.packages())
rfiles <- c("Matrix", "spam", "spam64", "raster", "sp", "spdep", 'mvtnorm')
#"parallel", "data.table", "lfe",
#"sphet", "sf", 'spacetime'

for( i in rfiles) { usethis::use_package( i ) }

# devtools::use_package( i, "Suggests", pkg=pdir)}

#install.packages()
#devtools::install_github("thk686/rgdal2")


#-------------------------------------------------------------------
##################
# Which Codes
################## 

rfile <- c(
    "MakeData.R",
    "SpatialWeights.R",
    "STsimFunctions.R",
    "vcovSHAC.R",
    "VarioJ.R")

#"TwoWayFE.R"
#"mfxall.R"
#"mfxi.R"
#"vcovST_meat.R"
#"vcovSTsep.R"
#"vcovSTsep_loop.R"
#"SpaceTimeSE_format.R"
#"SpaceTimeSE.R"

rfiles <- paste0(pdir,"Code/R_Code/",rfile)

# Move Code
file.copy(rfiles, rdir, overwrite=T )
devtools::load_all( rdir )

# Create Code Documentation
devtools::document( pkg=packg)

