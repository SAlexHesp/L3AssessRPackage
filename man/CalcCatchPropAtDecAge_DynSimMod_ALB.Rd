% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L3Analyses.R
\name{CalcCatchPropAtDecAge_DynSimMod_ALB}
\alias{CalcCatchPropAtDecAge_DynSimMod_ALB}
\title{Calculate expected annual retained and discarded catches and catch proportions at age, for dynamic simulation model}
\usage{
CalcCatchPropAtDecAge_DynSimMod_ALB(
  nYears,
  nTimeSteps,
  Ages,
  RetCatch_Fem_yal,
  RetCatch_Mal_yal,
  DiscCatch_Fem_yal,
  DiscCatch_Mal_yal
)
}
\arguments{
\item{nYears}{number of length classes}

\item{nTimeSteps}{number of model time steps per year}

\item{Ages}{ages considered by model}

\item{RetCatch_Fem_yal}{expected retained catch numbers of females in each year, at age and and length class}

\item{RetCatch_Mal_yal}{expected retained catch numbers of males in each year, at age and and length class}

\item{DiscCatch_Fem_yal}{expected discarded catch numbers of females in each year, at age and and length class}

\item{DiscCatch_Mal_yal}{expected discarded catch numbers of males in each year, at age and and length class}
}
\value{
RetCatchAtDecAge_Fem, RetCatchAtDecAge_Mal, DiscCatchAtDecAge_Fem, DiscCatchAtDecAge_Mal,
ExpRetCatchPropAtDecAge_Fem, ExpRetCatchPropAtDecAge_Mal, ExpDiscCatchPropAtDecAge_Fem,
ExpDiscCatchPropAtDecAge_Mal
}
\description{
Calculates expected annual retained and discarded catches and catch proportions at at age. Outputs may be used for generating random
catches at age from a from a dynamic simulation model.
}
\keyword{internal}
