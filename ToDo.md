# <p align=center> Things to Do </p>


## Better Make_*_Data Options

add `focal` error_type for `make_spacetime_data`
* include option `@param sar_factor, ar_factor parameter for error_type='focal'`

add better weights_matrix functionality for `make_time_data`
* also, use stats::arima.sim, stats::filter

`make_space_data` has two types of "focal" errors: SAR, SMA

## Better vcov* Options
* different kernels (currently only Bartlett)
* small sample adjustments
* add cluster correction via `sandwhich::cluster`
* option to calculate space weights-matrix
* option to calculate time weights-matrix
* option to calculate space-time kronecker weights-matrix 



## Better Support for non-lm objects

`lfe::felm`
* SHACsep with heteroskedasticityy **FAILS** for `felm` objects with type='HC3'
* `SHACsepHC=vcovSHACsep(reg, method='bruteforce_xy',cutoff_s=.6, loc_x=DFst$Loc_X, loc_y=DFst$Loc_X, loc_t=DFst$Time_ID)`

`plm`


## Support Bootstrapping
This is where pre-computing and recycling a weights matrix is most useful


## More Documentation (README or Vignette)
provide real example (e.g., stata example)
* http://freigeist.devmag.net/economics/936-spatial-temporal-autocorrelation-correction-with-r-conley.html
* https://github.com/darinchristensen/conley-se

example of bootstrapping efficiency with method='rolled'

## Recent Changes

Eliminated data.tables dependancy

Eliminated pre-formatting data-matrices

Name Changes
* vcovSTsep -> vcovST -> vcovSHACsep(add_hc=F, add_hac=T)
* vcovSCL -> vcovST.loop -> vcovSHACsep(add_hc=F, add_cl=T)




<!-- 
## RBLOCK_SHACtrials1 VCOVspatial works with DistMat
    ## Write Time series autocorrelation in RBLOCK1 as massive sparse matrix via kronecker.

    ## note: weight_mat like DistMat only for Eclidean Distances
    ## XOmegaX faster with MatMult.cpp?
    ## Write TimeDist as Sparse Mat
    
<details>
  <summary>click to expand</summary>
  ...
</details>


<!-- 
Name Changes from Paper2:
    vcovSTsep -> vcovST
    vcovSCL -> vcovST.loop

<!-- 

#------------------------------------------------------------------
##################
# Future Considerations
##################

#sourceCpp( paste0(sedir, "SparseMatMult.cpp") )
#lltLS <- cppFunction( "")

#http://casualinference.org/2016/07/07/big_fast_lm/
#RcppEigen::fastLm
#RcppArmadillo::fastLm
#RcppGSL::fastLm


<!-- 
#------------------------------------------------------------------
##################
# Another Example
##################
# u = (I - pM)^-1  e  = W e
# E(uu') = 1/n Omega = 1/n W e e' W'

# W <- solve( diag( nrow(wmat) ) - rho * wmat )
# image.plot(W, col=rainbow(10))

# WMAT <- toeplitz(.2^(0:4)); e <- 1:5
## If spatial and temporal autocorrelation was constant over time, 
## V(e) = Sigma_{N\timesN} \kronecker Sigma_{T}, where
## Sigma_{N\timesN} is (NxN)^2 matrix of spatial autocorrelation
## Sigma_{T} is (T)^2 matrix of temporal autocorrelation


#x <- cbind(1:5, 1:5*2)
#e <- runif(5)
#K <- toeplitz( .4^(0:4))

#V <- K * e %*% t(e) 
#XVX <- t(x) %*% V %*% x
#XVX2 <- (t(x) %*% (K * diag(e) ) ) %*% ( diag(e)  %*% x )

-->

