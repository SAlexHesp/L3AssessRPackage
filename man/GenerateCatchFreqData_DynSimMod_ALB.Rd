% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L3Analyses.R
\name{GenerateCatchFreqData_DynSimMod_ALB}
\alias{GenerateCatchFreqData_DynSimMod_ALB}
\title{Get age, length, and condition age-at-length frequency distributions and individuals ages and lengths}
\usage{
GenerateCatchFreqData_DynSimMod_ALB(
  nYears,
  TimeStep,
  nTimeSteps,
  midpt,
  nLenCl,
  ret_L50,
  SampleSize_Fem,
  SampleSize_Mal,
  SampleSize,
  DiscSampleSize_Fem,
  DiscSampleSize_Mal,
  DiscSampleSize,
  ObsRetCatchFreqAtLen_Fem,
  ExpRetCatchPropDecAgeGivenLength_Fem,
  ObsRetCatchFreqAtLen_Mal,
  ExpRetCatchPropDecAgeGivenLength_Mal,
  ObsDiscCatchFreqAtLen_Fem,
  ExpDiscCatchPropDecAgeGivenLength_Fem,
  ObsDiscCatchFreqAtLen_Mal,
  ExpDiscCatchPropDecAgeGivenLength_Mal
)
}
\arguments{
\item{nYears}{number of length classes}

\item{TimeStep}{model timestep}

\item{nTimeSteps}{number of mode time steps per year}

\item{midpt}{mid points of length classes}

\item{nLenCl}{number of length classes}

\item{ret_L50}{user input retention parameters}

\item{SampleSize_Fem}{annual sample sizes for retained females}

\item{SampleSize_Mal}{annual sample sizes for retained males}

\item{SampleSize}{annual sample sizes for retained females and males}

\item{DiscSampleSize_Fem}{annual sample sizes for discarded females}

\item{DiscSampleSize_Mal}{annual sample sizes for discarded males}

\item{DiscSampleSize}{annual sample sizes for discarded females and males}

\item{ObsRetCatchFreqAtLen_Fem}{randomly-generated female retained catch length frequency distributions}

\item{ExpRetCatchPropDecAgeGivenLength_Fem}{expected retained female catch age proportions at length, from dynamic simulation model}

\item{ObsRetCatchFreqAtLen_Mal}{randomly-generated male retained catch length frequency distributions}

\item{ExpRetCatchPropDecAgeGivenLength_Mal}{expected retained male catch age proportions at length, from dynamic simulation model}

\item{ObsDiscCatchFreqAtLen_Fem}{randomly-generated female discarded catch length frequency distributions}

\item{ExpDiscCatchPropDecAgeGivenLength_Fem}{expected discarded female catch age proportions at length, from dynamic simulation model}

\item{ObsDiscCatchFreqAtLen_Mal}{randomly-generated male discarded catch length frequency distributions}

\item{ExpDiscCatchPropDecAgeGivenLength_Mal}{expected discarded male catch age proportions at length, from dynamic simulation model}
}
\value{
ObsRetCatchFreqAtLengthAndDecAge_Fem, ObsRetCatchFreqAtLengthAndDecAge_Mal, ObsRetCatchFreqAtLengthAndDecAge,
ObsDiscCatchFreqAtLengthAndDecAge_Fem, ObsDiscCatchFreqAtLengthAndDecAge_Mal, ObsDiscCatchFreqAtLengthAndDecAge,
ObsDecAgeRetCatch_Fem, ObsDecAgeRetCatch_Mal, ObsDecAgeRetCatch, ObsDecAgeDiscCatch_Fem, ObsDecAgeDiscCatch_Mal,
ObsDecAgeDiscCatch, ObsAgeClRetCatch_Fem, ObsAgeClRetCatch_Mal, ObsAgeClRetCatch, ObsAgeClDiscCatch_Fem,
ObsAgeClDiscCatch_Mal, ObsAgeClDiscCatch, ObsLenClRetCatchMidPt_Fem, ObsLenClRetCatchMidPt_Mal, ObsLenClRetCatchMidPt,
ObsLenClDiscCatchMidPt_Fem, ObsLenClDiscCatchMidPt_Mal, ObsLenClDiscCatchMidPt, ObsRetCatchFreqAtDecAge_Fem,
ObsRetCatchFreqAtDecAge_Mal, ObsRetCatchFreqAtDecAge, ObsDiscCatchFreqAtDecAge_Fem, ObsDiscCatchFreqAtDecAge_Mal,
ObsDiscCatchFreqAtDecAge
}
\description{
Gets age, length, and condition age-at-length frequency distributions and individuals ages and lengths,
generated from dynamic simulation model
}
\keyword{internal}
