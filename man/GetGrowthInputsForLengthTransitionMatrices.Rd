% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L3Analyses.R
\name{GetGrowthInputsForLengthTransitionMatrices}
\alias{GetGrowthInputsForLengthTransitionMatrices}
\title{Get inputs for length transition matrices}
\usage{
GetGrowthInputsForLengthTransitionMatrices(
  MaxAge,
  TimeStep,
  nLenCl,
  midpt,
  GrowthCurveType,
  GrowthModelType,
  GrowthParams,
  RefnceAges,
  SelectivityType
)
}
\arguments{
\item{MaxAge}{maximum age considered in model}

\item{TimeStep}{model timestep (e.g. 1 = annual, 1/12 = monthly)}

\item{nLenCl}{number of length classes}

\item{midpt}{mid points of length classes}

\item{GrowthCurveType}{1=von Bertalanffy, 2=Schnute}

\item{GrowthModelType}{1=von Bertalanffy - 1 sex, 2=von Bertalanffy - 2 sexes,
Schnute - 1 sex, 2=Schnute - 2 sexes}

\item{GrowthParams}{growth parameters (used for GrowthModelType 2 or 4)}

\item{RefnceAges}{Schnute reference ages (used for GrowthModelType 4)}

\item{SelectivityType}{1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve}
}
\value{
mean size at age from growth curve (MeanSizeAtAge), estimated length after one year for each length class
mid point, given growth curve (MeanEndingLength) and associated change in length after one year of growth (TimestepGrowthSizeInc)
}
\description{
Get key inputs required to calculate length transition matrices, for models with either one or two sexes,
and using von Bertalanffy or Schnute growth model
}
\keyword{internal}
