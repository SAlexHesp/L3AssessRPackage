% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L3Analyses.R
\name{GetAgeLenVectors_AgeLengthCatchCurve_Growth}
\alias{GetAgeLenVectors_AgeLengthCatchCurve_Growth}
\title{Get vectors of observed ages and lengths for plotting}
\usage{
GetAgeLenVectors_AgeLengthCatchCurve_Growth(
  nTimeSteps,
  TimeStep,
  nLenCl,
  midpt,
  ObsRetCatchFreqAtLen
)
}
\arguments{
\item{nTimeSteps}{c(Linf,vbK) von Bertalanffy growth parameters}

\item{TimeStep}{observed frequencies in length classes}

\item{nLenCl}{minimum frequency to include}

\item{midpt}{lower bounds of length classes}

\item{ObsRetCatchFreqAtLen}{upper bounds of length classes}
}
\value{
ObsAge, ObsAge_F, ObsAge_M, ObsLenClRetCatchMidPt, ObsLenClRetCatchMidPt_F, ObsLenClRetCatchMidPt_M
}
\description{
Get vectors of observed ages and lengths for plotting of growth for age and length based catch curve
}
\keyword{internal}
