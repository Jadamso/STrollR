#-------------------------------------------------------------------
##################
# Create Data
##################
library(STrollR)

DDFF <- fake_data_traditional(50)

#-------------------------------------------------------------------
##################
# Prepare Covariance Inputs
##################

# Do this early,
# then recycle for all your regressions

XEUC <- 1000  #km
XHAV <- XEUC*(360/(6378.388*2*pi))
d_cutoff <- XHAV
t_cutoff <- 10

coord_mat <- as.matrix(DDFF[,c("x","y")])
coord_mat_u <- coord_mat[DDFF$Time==1,]
WMAT <- weight_mat(coord_mat_u, cutoff=d_cutoff)

# N.B. distances are calculated from x,y coordinates
# can use lon,lat if change both coord_mat and weight_mat 

#-------------------------------------------------------------------
##################
# Regression
##################

## Main Covariates
## Time, Group Fixed Effects
## 0 Instruments
## Temporally Clustered on ID
form <- as.formula( Y~X+Z | Group + Time | 0 | ID)

REG <- lfe::felm( form, data=DDFF, keepCX=TRUE)

#-------------------------------------------------------------------
##################
# Standard Errors
##################

VCOV <- vcovST.loop(
    reg=REG,
    unit="ID",
    time="Time",
    sp_coords=c("x","y"),
    wmat=WMAT,
    DF=DDFF[, c("ID", "Time", "x", "y")],
    cores=1
)

## N.B. if not passing DF to vcovST.loop
## vcovST.format Assumes the second Fixed Effect is temporal
## must augment the original formula

#-------------------------------------------------------------------
##################
# Comparison
##################

## How much larger are spatial SE's than clustered?

diag( VCOV / REG$clustervcv )




