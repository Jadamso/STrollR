

# STrollR

Calculate covariance matrices that are corrected for spatial and temporal correlation. Brute-force options are available by default, as are **S**patio **T**emporal **roll**ed standed errors (in **R**).

<p  align="center">
<figure>
  <img src="https://raw.githubusercontent.com/Jadamso/STrollR/master/STsim/STvarX.gif"  align="center" width="150" height="150">
  <figcaption align="center">*VonNeumann(1, rho=.3), AR(1, rho=.8)*  </figcaption>
</figure>
<p>



To install latest release from Github:
```r
#install.packages("devtools")
devtools::install_github("Jadamso/STrollR",subdir="STrollR")
citation("STrollR")
```

<!--
knit(input="README.rmd", output="README.md")
-->


Huge spatiotemporal covariance matrices can be efficiently calculated by using sparse matrix approaches with spam and spam64, specifically spam::rdist. The correction is applied in a post-estimation style and some functions have been extended to take advantage of parallel processing. See [Theory.pdf](https://github.com/Jadamso/STrollR/blob/master/Theory.pdf) for a little theory behind the *rolled* calculations. Rolling is not only fast, but also transparent---eventually helping to correct bugs in the original codes ([stata update](http://www.fight-entropy.com/2010/06/standard-error-adjustment-ols-for.html), [stata code](http://www.globalpolicy.science/code/)) ([R update](https://darinchristensen.com/post/conley-correction/), [R code](http://www.trfetzer.com/using-r-to-estimate-spatial-hac-errors-per-conley/)).



## Usage Examples

#### Toy DataSets and Diagnostics

Create Datasets with Space, Time, or Space-Time attributes
```r
library(raster)
set.seed(33)
## Make Data
DFst <- make_spacetime_data(11,6,2)
DFs <- make_space_data(21,3)
```

Plot a Time Series of Maps of the Error Terms (correlated by default)
```r
Elist <- lapply(split(DFst, DFst$Time_ID), function(dft){
    raster::rasterFromXYZ(dft)$E
})
E  <- raster::stack( Elist)
rasterVis::levelplot(E, names.attr=paste0('Time=', names(Elist)))
```

Spatial AutoCorrelation Diagnotistics (Semi-Variogram)
```r
### Smooth Estimates
vario1 <- varioJ(DFs[,c('Loc_X', 'Loc_Y')], cutt=2, E=DFs$E)
vario1 <- vario1[order(vario1$Dist),]
vgram1 <- predict(loess(Vario~Dist, data=vario1, deg=1, span=.3),
    control=loess.control(statistics='none', trace.hat='approximate'))
### Binned Estimates
vgram2 <- geoR:::variog(coords=DFs[,c('Loc_X', 'Loc_Y')], data=DFs$E, bin.cloud=TRUE)

par(mfrow=c(1,2))
ylm <- c(0,max(vgram1, vgram2$v))
plot.new()
plot.window(xlim=range(vario1$Dist), ylim=ylm)
points(Vario~Dist, data=vario1, pch=16, col=grey(0.5,.05))
lines(vario1$Dist, vgram1)
axis(1)
axis(2)
mtext('distance', 1, line=2)
mtext("semivariance",2, line=2)
title('Smooth Estimates')
box()

geoR:::plot.variogram(vgram2, ylim=ylm)
title('Binned Estimates')
title("Empirical Semi-Variogram (Semi-Variance=|e_xy - e_x'y'|^2)", outer=T)
```

<!--
<p  align="center">
<figure>
  <img src="https://raw.githubusercontent.com/Jadamso/STrollR/master/STsim/STvarX.gif"  align="center" width="150" height="150">
  <figcaption align="center">
  Spatiotemporal Process Visualized

  *VonNeumann(1, rho=.3), AR(1, rho=.8)* </figcaption>
</figure>
<p>
--->



#### SHAC Covariance-Matrices for lm's

Run A Regression and Correct for Spatial Errors
```r
DFs <- make_space_data.raster(32,1, xpars=c(0,.1), theta=rep(4,1),
     sar_factor=.5, error_type='distance', error_scale=2)
reg <- lm(Y~X1, data=DFs)

vcv_list <- list(
    IID=vcov(reg),
    HC=sandwich::vcovHC(reg),
    SHAC=vcovSHAC(reg, method='bruteforce_xy',
        cutoff_s=.5, loc_x=DFs$Loc_X, loc_y=DFs$Loc_Y))
se_list <- lapply(vcv_list, function(vcv){ sqrt(diag(vcv)) })
reg_list <- lapply(seq(se_list), function(i) reg)
## Compare se's: Normal vs. SHAC
stargazer::stargazer( reg_list, se=se_list, type='text', omit.stat='all')
```



Run 400 Regressions and Summarize SE's: SHAC vs. IID
```r
## Simulate 400 Regressions
se_sims <- lapply(1:400, function(b){
    DFs <- make_space_data(21,1, xpars=c(0,.1), theta=rep(4,1),
         sar_factor=.5, error_type='focal', error_scale=2)
    reg <- lm(Y~X1, data=DFs)
    vcv_list <- list(
        IID=vcov(reg),
        SHAC=vcovSHAC(reg, method='bruteforce_xy',
            cutoff_s=.5, loc_x=DFs$Loc_X, loc_y=DFs$Loc_Y))
    se_list <- sapply(vcv_list, function(vcv){ sqrt(diag(vcv)) })
    se_list
})

## Summarize SE's
## Notice the standard errors were not *always* larger
ses.x1 <- t(sapply(se_sims, function(j){ j['X1', ] }))
bks <- seq(min(ses.x1), max(ses.x1), length.out=50)
h_shac <- hist(ses.x1[,'SHAC'], breaks=bks, plot=F)
h_iid <- hist(ses.x1[,'IID'], breaks=bks)
plot(h_shac, col=rgb(1,0,0,.25),
    xlim=range(ses.x1), ylim=range(h_shac$counts,h_iid$counts),
    main='Standard Error Estimates for X1', xlab='se')
plot(h_iid, col=rgb(0,0,1,.25), add=T)
legend('topright', legend=c('SHAC', 'IID'), cex=1,
    pch=15, col=c(rgb(1,0,0,.25),rgb(0,0,1,.25)))
```


Note that the standard errors do not even need to be systematically larger. See e.g.,
* http://www.mostlyharmlesseconometrics.com/2010/12/heteroskedasticity-and-standard-errors-big-and-small/
* https://stats.stackexchange.com/questions/412460/my-standard-errors-are-smaller-when-i-use-newey-west-help

In the above example, neighbouring cells have errors that do not contribute as much new information compared to if all errors were independant draws. In the extreme, if they were perfectly positively correctly we should really consider the entire grid a single draw (i.e., the coefficient estimates are based on a single data point so the se's should be larger).


Note that more any type of weights-matrix can be passed to vcovSHAC with method in c('rolled', 'semirolled'). This means it could be used for other corrections; e.g., a weights matrix that indicates clusters. *The general principal is that you should weight the residuals according to how strong their interdependance is.*

#### Standard Error Corrections for Spatio-Temporal Data

SE correction for `lm` objects
```r
DFst <- make_spacetime_data(21,6,2)
reg <- lm(Y~X1+X2, data=DFst)

## Compute se's
vcv_list <- list(
    STHAC=vcovSTHAC(reg, method='bruteforce_xyt',
        cutoff_s=.6, loc_x=DFst$Loc_X, loc_y=DFst$Loc_Y,
        cutoff_t=3, loc_t=DFst$Time_ID),
    SHACsepHC=vcovSHACsep(reg, method='bruteforce_xy',
        cutoff_s=.6, loc_x=DFst$Loc_X, loc_y=DFst$Loc_X, loc_t=DFst$Time_ID),
    SHACsepHAC=vcovSHACsep(reg, method='bruteforce_xy',
        cutoff_s=.6, loc_x=DFst$Loc_X, loc_y=DFst$Loc_X, loc_t=DFst$Time_ID,
        add_hc=F, add_hac=T),
    SHACcl=vcovSHAC(reg, method='bruteforce_xy',
        cutoff_s=.6, loc_x=DFst$Loc_X, loc_y=DFst$Loc_Y,
        add_hc=F, add_cl=T, cluster=DFst[,c('Time_ID')]),
    PL=sandwich::vcovPL(reg, cluster=DFst[,c('Space_ID','Time_ID')]),
    PC=sandwich::vcovPC(reg, cluster=DFst[,c('Space_ID','Time_ID')]),
    2CCE=sandwich::vcovCL(reg, cluster=DFst[,c('Space_ID','Time_ID')]),
    HC=sandwich::vcovHC(reg),
    IID=vcov(reg))
se_list <- lapply(vcv_list, function(vcv){ sqrt(diag(vcv)) })
reg_list <- lapply(seq(se_list), function(i) reg)

## Compare se's
stargazer::stargazer( reg_list, se=se_list, type='text', omit.stat='all')
```
STHAC (a.k.a. TS-SAR) accounts for the ``off-diagonal'' correlations, which often produces the largest se's in the above example (although SHACsepHC and SHACsepHAC are both close to STHAC). STHAC explicitly accounts for *cor(E(s,t),E(s',t'))!=0*, while SHACsep* treats spatial correlation seperately. Off diagonal correlations occur because *E(s,t)* spatially correlates with *E(s',t)* which temporally correlates with *E(s',t')*. Using `SHACsepHAC(method='rolled')` can be very efficient for long time series.


SE Correction for `lfe::felm` objects
(Not all heterkedasticity *robust* corrections are working)
```r
DFst <- make_spacetime_data(21, 6, 2)

## Add Clusters for Space and Time ID
DFst2 <-dummies::dummy.data.frame(DFst,
    names=c("Space_Group", "Time_Group"), sep='.')
Space_Group_FE <- grep('Space_Group', names(DFst2), value=T)
Space_Group_Coef <- as.numeric(gsub('Space_Group\\.','',Spatial_Group_FE))
Time_Group_FE <- grep('Time_Group', names(DFst2), value=T)
Time_Group_Coef <- as.numeric(gsub('Time_Group\\.','',Time_Group_FE))
Space_FE <- as.matrix(DFst2[,Space_Group_FE])%*%as.matrix(Space_Group_Coef)
Time_FE <- as.matrix(DFst2[,Time_Group_FE])%*%as.matrix(Time_Group_Coef)
## Add Fixed Effects to Outcome
DFst$Y <- DFst$Y + Space_FE + Time_FE

## Regression (Main Covariates | Fixed Effects | Instruments | Clusters)
reg <- lfe::felm( formula( Y~X1+X2 | Space_Group + Time_Group ), data=DFst)

## Compute SE's
vcv_list <- list(
    STHAC=vcovSTHAC(reg, method='bruteforce_xyt',
        cutoff_s=.6, loc_x=DFst$Loc_X, loc_y=DFst$Loc_Y,
        cutoff_t=3, loc_t=DFst$Time_ID,
        add_hc=T, type="HC0"),
    SHACsepHAC=vcovSHACsep(reg, method='bruteforce_xy',
        cutoff_s=.6, loc_x=DFst$Loc_X, loc_y=DFst$Loc_X, loc_t=DFst$Time_ID,
        add_hc=F, add_hac=T),
    SHACcl=vcovSHAC(reg, method='bruteforce_xy',
        cutoff_s=.6, loc_x=DFst$Loc_X, loc_y=DFst$Loc_Y,
        add_hc=F, add_cl=T, cluster=DFst[,c('Time_ID')]),
    PL=sandwich::vcovPL(reg, cluster=DFst[,c('Space_ID','Time_ID')]),
    PC=sandwich::vcovPC(reg, cluster=DFst[,c('Space_ID','Time_ID')]),
    2CCE=sandwich::vcovCL(reg, cluster=DFst[,c('Space_ID','Time_ID')]),
    HC2=vcov(reg))

se_list <- lapply(vcv_list, function(vcv){ sqrt(diag(vcv)) })
reg_list <- lapply(seq(se_list), function(i) reg)

## Compute Default se's (reghdfe is Stata's Clustering)
vcv_list2<- lapply( c("cgm", "cgm2", "reghdfe"), function(c_meth){
    reg2 <- lfe::felm( formula(Y~X1+X2|Space_Group+Time_Group|0|Space_ID+Time_ID),
        data=DFst, cmethod=c_meth)
    vcov(reg2)
})
se_list2 <- lapply(vcv_list2, function(vcv){ sqrt(diag(vcv)) })
reg_list2 <- lapply(seq(se_list2), function(i) reg)

## Compare se's
stargazer::stargazer( c(reg_list, reg_list2), se=c(se_list, se_list2), type='text', omit.stat='all')
```

        
<!--
#### Space-Time Standard Error Corrections with a Real Dataset
```r
```
--->


## Some of the Nitty-Gritty

#### Intuition of the approach for Time-Series Data

Here is an example to display the logic, applied to temporal correlations for easy comparison.
```r
# Setup
n <- 100
e <- runif(n)
v <- runif(100) + c( 0, head(e, -1) )*.6
X <- cbind(1,1:n)
Y <- X%*%as.matrix(c(1,2)) + v

# Regression and Output
reg <- lm(Y~X)
r <- reg$residuals

# Covariance Matrix `Meat`
# Newey-West Sandwich Style
V0 <- sandwich::NeweyWest(reg, sandwich=F, prewhite=F, lag=2, adjust=F)
# Sparse Rolling Style
library(Matrix)
w <- as( toeplitz(c(1,2/3,1/3, rep(0,n-3))), "dgCMatrix")
V2 <- 1/n * t(X) %*% t( r*t(w*r) ) %*% X

# Compare
(V0 - V2) < 1E3*.Machine$double.eps
#            (Intercept)    x
#(Intercept)        TRUE TRUE
#x                  TRUE TRUE
```


#### Comparison of Internal Methods

Compare Meat Matrix Calculations by Computational Approach (for data with X,Y coordinates)
```r
DFs <- make_space_data(21,2)
reg <- lm(Y~X1+X2, data=DFs)
E <- resid(reg)
Xvars <- names(coef(reg))
X <- as.matrix(model.frame.i(reg)[,Xvars])

## Compute Meat Matrices
cutoff_s <- .5  #generic map units
system.time(XWeeX1 <- XOmegaX_rolled(X, weight_mat.df(DFs, cutoff_s=cutoff_s), E)) ##.066
system.time(XWeeX2 <- XOmegaX_semirolled(X, weight_mat.df(DFs, cutoff_s=cutoff_s), E)) ##.186
system.time(XWeeX3 <- XOmegaX_bruteforce_xy(X, E, DFs$Loc_Y, DFs$Loc_X, cutoff_s)) ##.027

## See that any differences are negligable
sum((XWeeX3-XWeeX2)^2)
sum((XWeeX3-XWeeX1)^2)
sum((XWeeX2-XWeeX1)^2)
```

Compare Brute-Force Meat Matrix Calculations for X,Y and Lat,Lon Data
```r
DFs <- make_space_data(21,2)
reg <- lm(Y~X1+X2, data=DFs)
E <- resid(reg)
Xvars <- names(coef(reg))
X <- as.matrix(model.frame.i(reg)[,Xvars])
XY_Coords <- DFs[,c('Loc_X', 'Loc_Y')]
cutoff_s <- .5  #generic map units
WMAT  <- weight_mat.df(XY_Coords, cutoff_s=cutoff_s)

## Convert from Generic Coordinates to LatLon
## (assume generic coordinates were centered on Rome with Lambert Azimuthal Equal Area projection)
library(sp)
LL_Coords   <- sp::SpatialPoints( XY_Coords,
    proj4string=sp::CRS('+proj=laea +lat_0=41.891775 +lon_0=12.486137'))
LL_Coords <- as.data.frame(sp::spTransform(LL_Coords, sp::CRS("+init=epsg:4326") ))

## Compute Meat Matrix via lat,lon coordinates
## LAEA map units are in meters, so transform to km before converting to angles
cutoff_km <- cutoff_s/1000
XWeeX3 <- XOmegaX_bruteforce_xy(X, E, DFs$Loc_Y, DFs$Loc_X, cutoff_s)
XWeeX4 <- XOmegaX_bruteforce_ll(X, E, LL_Coords$Loc_Y, LL_Coords$Loc_X, cutoff_km)
XWeeX5 <- XOmegaX_bruteforce_ll(X, E, LL_Coords$Loc_Y, LL_Coords$Loc_X, cutoff_km, manual_dist=F)
## See that the differences are very small, but not zero
(XWeeX4-XWeeX3)/(XWeeX3)*100
(XWeeX5-XWeeX3)/(XWeeX3)*100

## Plot Weights Matrices
WMAT2 <- weight_mat.df(LL_Coords, cutoff_s=cutoff_km, latlon=T, cutoff_km2angles=T)
gridExtra::grid.arrange(image(WMAT),image(WMAT2), ncol=2)
## notice they are very similar, but not identical
```





## Some Resources I found useful along the way
* Solomon Hsiang (2010), Kelejian and Prucha (2006), Kelejian and Prucha (1999), Conley (1999)
* https://msu.edu/~tjv/ferobustse-working-paper.pdf
* https://editorialexpress.com/cgi-bin/conference/download.cgi?db_name=IAAE2017&paper_id=612
* https://spatial.uchicago.edu/sites/spatial.uchicago.edu/files/11_spatial2sls_slides.pdf
* http://freigeist.devmag.net/economics/936-spatial-temporal-autocorrelation-correction-with-r-conley.html, http://www.trfetzer.com/using-r-to-estimate-spatial-hac-errors-per-conley/, https://darinchristensen.com/post/conley-correction/
* http://www.stata.com/meeting/sandiego12/materials/sd12_jeanty.pdf
<!---
* https://gist.github.com/fabian-s/9665846#file-2014-03-20-dynamic-dispatch-for-sparse-matrices-rmd

--->


## Caveat Emptor
This is a pre-release to the public and not gauranteed to work correctly with your code and data.

