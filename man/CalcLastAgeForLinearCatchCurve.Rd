% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L3Analyses.R
\name{CalcLastAgeForLinearCatchCurve}
\alias{CalcLastAgeForLinearCatchCurve}
\title{Determine age at full recruitment into the fishery (assuming knife-edge selection)}
\usage{
CalcLastAgeForLinearCatchCurve(MinFreq, RecAge, Ages, ObsAgeFreq)
}
\arguments{
\item{MinFreq}{minimum frequency of fish for including data for old fish}

\item{RecAge}{recruitment age}

\item{Ages}{ages for analysis}

\item{ObsAgeFreq}{observed age frequency}
}
\value{
oldest age to be included for linear catch curve analysis (LastAgeForLinearCC)
}
\description{
This function calculates the oldest age to be included in data for a linear catch curve analysis,
according to a specified minimum frequency
}
\examples{
# simulate data
set.seed(123)
MinAge = 1
MaxAge = 40
Ages = MinAge:MaxAge
NatMort <- exp(1.46 - (1.01 * (log(MaxAge)))) # i.e. Hoenig's (1983) eqn for fish
FMort = 0.1
ZMort = FMort + NatMort
SelA50 = 6
SelA95 = 8
SampleSize = 1000 # required sample size. For 2 sex model, same sample size generated for each sex.
Res=SimAgeFreqData_EqMod(SampleSize, MinAge, MaxAge, SelA50, SelA95, NatMort, FMort)
ObsAgeFreq = Res$CatchSample
RecAge = CalcRecruitmentAge(RecAssump=0, Ages, ObsAgeFreq)
CalcLastAgeForLinearCatchCurve(MinFreq=1, RecAge, Ages, ObsAgeFreq)
}
\keyword{internal}
