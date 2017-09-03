#-------------------------------------------------------------------
##################
# Libraries
##################

## Imported Packages: (.packages())
rfiles <- c("Matrix","spam","spam64",
    "parallel","data.table","lfe", "spdep", "sphet")
for( i in rfiles) {
    devtools::use_package( i, pkg=packg)
}

# devtools::use_package( i, "Suggests", pkg=pdir)}


#install.packages("sf")
#devtools::install_github("thk686/rgdal2")
# library(spacetime)

#-------------------------------------------------------------------
##################
# CHUNK A
################## 
# Add Codes


# Which Codes
rfile0 <- c(
    "VarioJ.R",
    #"TwoWayFE.R",
    "SpatialWeights.R",
    "vcovST_meat.R",
    "vcovSTsep.R",
    "vcovSTsep_loop.R",
    "SpaceTimeSE_format.R",  
    "SpaceTimeSE.R",
    "STFakeData.R",
    "STsimFunctions.R"
)

codedire <- path.expand("~/Desktop/Common/R_Code/")
rfiles <- paste0(codedire, rfile0)

# Move Code
file.copy(rfiles, rdir, overwrite=T )
devtools::load_all( rdir )

# Create Code Documentation
devtools::document( pkg=packg)

