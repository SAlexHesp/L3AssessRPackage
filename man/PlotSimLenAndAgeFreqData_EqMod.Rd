% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L3Analyses.R
\name{PlotSimLenAndAgeFreqData_EqMod}
\alias{PlotSimLenAndAgeFreqData_EqMod}
\title{Produce plots of simulated age and length data from equilibrium population simulation model}
\usage{
PlotSimLenAndAgeFreqData_EqMod(MaxAge, MaxLen, SimRes, PlotOpt)
}
\arguments{
\item{MaxAge}{Maximum age}

\item{MaxLen}{Maximum length}

\item{SimRes}{Outputs of SimLenAndAgeFreqData_EqMod function}

\item{PlotOpt}{# 0=all plots, 1=retained lengths at age, 2=retained plus discarded lengths at age, 3=length frequency, 4=age frequency}
}
\value{
plots of simulated age and length data
}
\description{
Produce plots of simulated age and length data from equilibrium population simulation model
}
\examples{
# Simulate data
set.seed(123)
SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
MaxAge = 30
TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = 4.22/MaxAge
FishMort = 0.2
MaxLen = 1200
LenInc = 20
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityAtLen = NA # selectivity vector
SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
RetenParams = c(400, 50) # L50, L95-L50 for retention
DiscMort = 0 # proportion of fish that die due to natural mortality
# single sex, von Bertalanffy
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 800
vbK = 0.2
CVSizeAtAge = 0.08
GrowthParams = c(Linf, vbK)
RefnceAges = NA
SimRes=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                            SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
# plot data
PlotOpt=0 # 0=all plots, 1=retained lengths at age and growth curves, 2=retained plus discard lengths at age
# and growth curves, 3=retained plus discard length frequency, 4=retained plus discard age frequency,
# 5=sex specific retained lengths at age and growth curves, 6=retained plus discard lengths at age and growth curves,
# 7=female length frequency, 8=male length frequency, 9=female age frequency, 10=male age frequency,
# 11=selectivity/retention, 12=F-at-age reten + disc, 13=F-at-age reten, 14=F-at-age disc
PlotSimLenAndAgeFreqData_EqMod(MaxAge, MaxLen, SimRes, PlotOpt)
}
