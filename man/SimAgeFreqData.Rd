% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L3Analyses.R
\name{SimAgeFreqData}
\alias{SimAgeFreqData}
\title{Simulate age frequency data}
\usage{
SimAgeFreqData(SampleSize, MinAge, MaxAge, SelA50, SelA95, NatMort, FMort)
}
\arguments{
\item{SampleSize}{required sample size}

\item{MinAge}{minimum age}

\item{MaxAge}{maximum age}

\item{SelA50}{age at 50 percent selectivity}

\item{SelA95}{age at 95 percent selectivity}

\item{NatMort}{natural mortality}

\item{FMort}{fully-selected fishing mortality}
}
\value{
Age classes (Ages), selectivity at age (SelAtAge), expected proportions at age (PropAtAge), observed
catch frequency at age (CatchSample)
}
\description{
Simulate age frequency data with specified selectivity and mortality parameters
}
\examples{
set.seed(123)
MinAge = 1
MaxAge = 40
Ages = MinAge:MaxAge
NatMort <- exp(1.46 - (1.01 * (log(MaxAge)))) # i.e. Hoenig's (1983) eqn for fish
FMort = 0.1
ZMort = FMort + NatMort
SelA50 = 6
SelA95 = 8
SampleSize = 1000 # required number of fish for age sample
Res=SimAgeFreqData(SampleSize, MinAge, MaxAge, SelA50, SelA95, NatMort, FMort)
ObsAgeFreq = unlist(as.vector(Res$CatchSample))
}
