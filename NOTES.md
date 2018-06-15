

Recent Name Changes:
* vcovSTsep -> vcovST
* vcovSCL -> vcovST.loop
 

## Things To Do 

* `vcovSTkron`

* small sample / asymptotic adjustment options

* support `lm` objects fully

`felm` is supported in most functions
currently `lm` is supported in just in a few functions

* `vcovST.loop` non-default behaviour

If not passing DF to vcovST.loop
`vcovST.format` Assumes the second Fixed Effect is temporal
must augment the original formula

* eliminate data.tables dependancy

`sub_dat <- as.data.frame(sub_dat)` in `vcovST_meat` 
should be eliminated but I cannot compile data.tables correctly



<!-- vcovSTkron

## RBLOCK_SHACtrials1 VCOVspatial works with DistMat
    ## Write Time series autocorrelation in RBLOCK1 as massive sparse matrix via kronecker.

    ## note: weight_mat like DistMat only for Eclidean Distances
    ## XOmegaX faster with MatMult.cpp?
    ## Write TimeDist as Sparse Mat
-->
