% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L3Analyses.R
\name{CalcLogisticSelOrReten}
\alias{CalcLogisticSelOrReten}
\title{Calculate logistic length-based selectivity (or retention)}
\usage{
CalcLogisticSelOrReten(L50, L95, midpt)
}
\arguments{
\item{L50}{length at which 50 percent of fish are selected into the fishery}

\item{L95}{length at which 95 percent of fish are selected into the fishery}

\item{midpt}{mid points of length classes}

\item{nLenCl}{number of length classes}
}
\value{
Selectivity at length (SelAtLength)
}
\description{
Calculate logistic length-based selectivity or retention (asymptotic curve)
}
\keyword{internal}
