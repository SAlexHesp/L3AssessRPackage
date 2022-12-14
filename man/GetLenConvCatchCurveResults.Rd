% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L3Analyses.R
\name{GetLenConvCatchCurveResults}
\alias{GetLenConvCatchCurveResults}
\title{Fit length-converted catch model and return results}
\usage{
GetLenConvCatchCurveResults(
  ModelType,
  GrowthParams,
  RefnceAges,
  ObsCatchFreqAtLen,
  MinFreq,
  lbnd,
  midpt,
  ubnd
)
}
\arguments{
\item{ModelType}{1=von Bertalanffy growth, 2=Schnute growth}

\item{GrowthParams}{c(Linf, vbK) von Bertalanffy growth or c(y1, y2, a, b) Schnute growth}

\item{RefnceAges}{Schnute growth curve reference ages}

\item{ObsCatchFreqAtLen}{observed frequencies in length classes}

\item{MinFreq}{minimum frequency to include}

\item{lbnd}{lower bounds of length classes}

\item{midpt}{mid points of length classes}

\item{ubnd}{upper bounds of length classes}
}
\value{
sample size for observed data (SampleSize), Estimates of parameters and associated 95 percent
confidence limits (ParamEst), length class with greatest frequency of fish (PeakLencl), last length
class considered in analysis (LastLenCl), observed frequency data to which catch curve is fitted (ObsCatchFreqAtLen2),
relative ages at mid points of length classes (Age_midptlencl), natural logarithms of observed and expected numbers in length
classes divided by time interval taken to grow through length class (Obs_ln_n_dt, Est_ln_n_dt), lower and upper
95 percent confidence limits for (Est_ln_n_dtlow, Est_ln_n_dtup)
}
\description{
Fit specified length-converted catch curve model with growth described either by a von
Bertalanffy growth model or Schunute growth model (with a and b not equal to zero)
}
\examples{
# Simulate data
SampleSize=1000
set.seed(123)
MaxAge = 30
TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = 4.22/MaxAge
FishMort = 0.2
MaxLen = 1500
LenInc = 50
MLL = NA
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
L50 = 300 # selectivity
L95 = 400 # selectivity
SelectivityVec = NA
DiscMort = 0
# # von Bertalanffy
# GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
# Linf = 800
# vbK = 0.2
# CVSizeAtAge = 0.08
# GrowthParams = c(Linf, vbK)
# RefnceAges = NA
# Schnute
GrowthCurveType = 2 # 1 = von Bertalanffy, 2 = Schnute
t1 = 1 # growth - Schnute
t2 = 10 # growth - Schnute
y1 = 400 # growth - Schnute
y2 = 1000 # growth - Schnute
a = 0.1 # growth - Schnute
b = 2.0 # growth - Schnute
CVSizeAtAge = 0.08
GrowthParams = c(y1, y2, a, b)
RefnceAges = c(t1,t2)
Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                         L50, L95, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)

ObsCatchFreqAtLen = as.vector(Res$ObsCatchFreqAtLen)
MinFreq = 20 # set minimum frequency for larger lengths for analysis
# note, this needs to be high enough so that data for ln(n/dt) vs relative age for essentially straight
# line - if not, Z will be biased!!!
midpt=Res$midpt
lbnd=Res$lbnd
ubnd=Res$ubnd
ModelType = 2 # 1 = von Bertalanffy growth curve (Pauly), 2 = length-converted catch curve - Schnute growth curve)
res=GetLenConvCatchCurveResults(ModelType, GrowthParams, RefnceAges, ObsCatchFreqAtLen,
                                MinFreq, lbnd, midpt, ubnd)
}
