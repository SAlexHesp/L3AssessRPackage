% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L3Analyses.R
\name{GetChapmanRobsonMortalityResults}
\alias{GetChapmanRobsonMortalityResults}
\title{Estimate total mortality applying the Chapman & Robson (1960) mortality estimator}
\usage{
GetChapmanRobsonMortalityResults(
  RecAssump,
  SpecRecAge,
  MinAge,
  MaxAge,
  ObsAgeFreq
)
}
\arguments{
\item{RecAssump}{0=age at peak frequency, 1=age at peak frequency + 1, 2=specified age at full recruitment}

\item{SpecRecAge}{specified at at full recruitment, set to NA when RecAssump is set to 0 or 1}

\item{MinAge}{minimum age}

\item{MaxAge}{maximum age}

\item{ObsAgeFreq}{observed age frequency}
}
\value{
total mortality and approximate lower and upper 95 percent confidence
limits, (EstZMort), approximate standard error for total mortality (EstZMort_se),
median estimate and lower and upper 95 percent confidence limits for Z from resampling (EstZMort_resamp),
variables for plotting, including expected frequencies at age with associated confidence limits, calculated from resampling
(EstFreq, EstFreq_Zup, EstFreq_Zlow), recruitment age (RecAge),  maximum age in sample (MaxAgeInSample),
relative ages used in analysis (CRAges), sample size for relative ages (n) and additional variable
calculated in analysis (CR_T)
}
\description{
This function provides an estimate of total mortality, applying the Chapman & Robson (1960) mortality
estimator, with associated 95 percent confidence limits. Additional variables are returned for plotting
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
ObsAgeFreq = as.vector(unlist(Res$CatchSample))
res=GetChapmanRobsonMortalityResults(RecAssump=1, SpecRecAge=NA, MinAge, MaxAge, ObsAgeFreq)
}
