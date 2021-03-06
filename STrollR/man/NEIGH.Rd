% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialWeights.R
\name{NEIGH}
\alias{NEIGH}
\title{Calculate the weights objects}
\usage{
NEIGH(
  coord_sp,
  neigh = 1,
  knn = TRUE,
  adj = FALSE,
  dnn = FALSE,
  rast = FALSE,
  vario = FALSE,
  tracer = TRUE,
  tr_type = "mult",
  tr_m = 20,
  tr_p = 16,
  symm = TRUE,
  symm_check = TRUE,
  SAVE = NA,
  write_gwt = F
)
}
\arguments{
\item{coord_sp}{matrix of coordinates or a SpatialPoints object}

\item{neigh}{number of neighbours to use in calculation}

\item{knn}{calculate weights using knn approach}

\item{adj}{calculate vonneumann weights (see VonNeumann)}

\item{dnn}{dnn approach unsupported}

\item{rast}{raster approach unsupported}

\item{vario}{is coord_sp a weights matrix?}

\item{tracer}{create trace matrix objects?}

\item{tr_type}{type of trace matrix}

\item{tr_m}{trace matrix m}

\item{tr_p}{trace matrix p}

\item{symm}{make weights symmetric}

\item{symm_check}{check for symmetric weights matrix}

\item{SAVE}{filename to save to, NA (default) returns as object}

\item{write_gwt}{create GWT objects used in spdep or sphet}
}
\value{
filename of saved objects, or returns objects if SAVE=NA
}
\description{
Calculate the weights objects
}
