% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L3Analyses.R
\name{GetParamRes_AgeAndLengthBasedCatchCurve}
\alias{GetParamRes_AgeAndLengthBasedCatchCurve}
\title{Get estimates of model parameters and associated 95 percent confidence intervals}
\usage{
GetParamRes_AgeAndLengthBasedCatchCurve(
  GrowthCurveType,
  SelectivityType,
  params,
  nlmb,
  ses
)
}
\arguments{
\item{GrowthCurveType}{1=von Bertalaffy, 2=Schnute}

\item{SelectivityType}{1=selectivity inputted as vector, 2=logistic}

\item{params}{estimated parameters, varying depending on growth model type, selectivity type and number of sexes}

\item{nlmb}{output from nlminb optimisation}

\item{ses}{estimates of standard errors for estimated parameters}
}
\value{
ParamEst
}
\description{
Get estimates of model parameters and associated 95 percent confidence intervals. Outputs vary
depending on growth model type, selectivity type and number of sexes
}
\keyword{internal}
