% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L3Analyses.R
\name{GetRetainedCatchStats}
\alias{GetRetainedCatchStats}
\title{Get retained catch summary statistics}
\usage{
GetRetainedCatchStats(
  ObsRetCatchFreqAtLen,
  ObsTotCatchAtLen,
  ObsLenClRetCatchMidPt_Fem,
  ObsLenClRetCatchMidPt_Mal,
  ObsAgeClRetCatch_Fem,
  ObsAgeClRetCatch_Mal
)
}
\arguments{
\item{ObsRetCatchFreqAtLen}{randomly-generated length frequency data for retained fish}

\item{ObsTotCatchAtLen}{randomly-generated length frequency data for retained and discarded fish}

\item{ObsLenClRetCatchMidPt_Fem}{randomly-generated length data for females}

\item{ObsLenClRetCatchMidPt_Mal}{randomly-generated length data for males}

\item{ObsAgeClRetCatch_Fem}{randomly-generated age class data for females}

\item{ObsAgeClRetCatch_Mal}{randomly-generated length class data for females}
}
\value{
ObsPropRetFish, ObsPropRetFish_sd, ObsMeanLenRetFish_Fem, ObsMeanLenRetFish_Mal, ObsMeanLenRetFish,
ObsLenRetFish_sd_Fem, ObsLenRetFish_sd_Mal, ObsLenRetFish_sd, ObsMeanLenRetFish_sd_Fem, ObsMeanLenRetFish_sd_Mal,
ObsMeanLenRetFish_sd, ObsMeanAgeRetFish_Fem, ObsMeanAgeRetFish_Mal, ObsMeanAgeRetFish, ObsAgeRetFish_sd_Fem,
ObsAgeRetFish_sd_Mal, ObsAgeRetFish_sd, ObsMeanAgeRetFish_sd_Fem, ObsMeanAgeRetFish_sd_Mal, ObsMeanAgeRetFish_sd
}
\description{
Calculate prop retained and mean length of retained fish
}
\keyword{internal}
