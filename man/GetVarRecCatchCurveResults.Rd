% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L3Analyses.R
\name{GetVarRecCatchCurveResults}
\alias{GetVarRecCatchCurveResults}
\title{Fit multi-year catch curve allowing for variable recruitment and get results}
\usage{
GetVarRecCatchCurveResults(dat, params)
}
\arguments{
\item{dat}{input data list}

\item{params}{model parameter list (fishing mortality and recruitment deviations)}
}
\value{
SummaryRes (convergence, NLL, sdr, param_est, params_sd, EstFMort, EstRecDevs),
RTMBReport (ExpCatchPropAtAge, EstRecDevYrs, RecruitYrs, FirstSampYr, t1_yrcl, RecDevs)
}
\description{
This function fits a  multi-year catch curve allowing for variable recruitment.
The method estimated fishing mortality and recruitment deviations. As currently
implemented, knife-edge recruitment is assumed.
}
\examples{
library(RTMB)
library(L3Assess)
# Specify exploitation history
nPeriods <- 3 # number of periods of mortality trends
InitYr <- 1950 # first year of mortality
PeriodEndYr <- c(1995,2000,2020) # periods when mortality changes
# Specifying stable F
NatMort = 4.22/40 # natural mortality  (year-1)
InitYr_FMort <- 0.67*NatMort # initial mortality, at beginning of first period
PeriodEndYr_FMort <- c(0.67*NatMort,0.67*NatMort,0.67*NatMort)
# # Specifying initial in increase F, reduction in F then stable F period
#InitYr_FMort <- 1.5*NatMort # initial mortality, at beginning of first period
#PeriodEndYr_FMort <- c(1.5*NatMort,1.5*NatMort,1.5*NatMort)
autocorr <- 0
randerr_sd <- 0
set.seed(1975)
res=SimExploitationHistory_DynMod(nPeriods, InitYr, PeriodEndYr, InitYr_FMort,
                                  PeriodEndYr_FMort, autocorr, randerr_sd)
Years = res$Years
nYears <- res$nYears
FMortByYear <- res$rand_FMort
# par(mfrow=c(1,1))
# Years = res$Years
# plot(Years, res$rand_FMort, ylim=c(0,0.5))
# Generate random age data
# Simulate age data from a dynamic age-based population model
InitRecruit <- 1 # Initial recruitment
MaxModelAge <- 50 # maximum age considered by model, years
TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be approprriate to use a smaller time step)
Linf <- c(1000, 1000) # mm - von Bertalanffy growth model parameters - Females, males
vbK <- c(0.1, 0.1) # year-1 - von Bertalanffy growth model parameters - Females, males
tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
EstLenAtAge <- data.frame(EstFemLenAtAge=NA, EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
ln_lenwt_a <- NA # for log-log relationship
lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
WLrel_Type <- 1 # 1=power, 2=log-log relationship
EstWtAtAge <- data.frame(EstFemWtAtAge=NA, EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
InitRatioFem <- 0.5 # Ratio of females to males at age zero
FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
FinalSex_A95 <- NA # Logistic sex change relationship parameters (95\% of max probability)
EstSexRatioAtAge <- NA  # sex ratio at age (from age 0) inputted as vector
EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, ~0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
Gear_sel_A50 <- NA # females, males - Logistic age selectivity relationship parameters
Gear_sel_A95 <- NA # females, males - Logistic age selectivity relationship parameters
EstGearSelAtAge <- data.frame(FemGearSelAtAge=NA, MalGearSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
Land_sel_A50 <- NA # females, males - Logistic age selectivity relationship parameters
Land_sel_A95 <- NA # females, males - Logistic age selectivity relationship parameters
FemSelLandAtAge <- rep(1,MaxModelAge+1) # assuming knife edge selectivity (to simplify)
FemSelLandAtAge[1:6] = 0
MalSelLandAtAge <- FemSelLandAtAge
EstLandSelAtAge <- data.frame(FemSelLandAtAge=FemSelLandAtAge, MalSelLandAtAge=MalSelLandAtAge) # gear selectivity at age (from age 0), inputted as values in data frame
ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
EstRetenAtAge <- data.frame(FemRetProbAtAge=NA, MalRetProbAtAge=NA) # fish retention at age (from age 0), inputted as values in data frame
DiscMort <- 0.0 # discard mortality (e.g. 50\% released fish die = 0.5)
Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
SimAnnSampSize <- 1000
lnSigmaR <- 0.6
autocorr <- 0
InitRec <- 1000
res=SimAgeFreqData_DynMod_AB(SimAnnSampSize, nYears, lnSigmaR, autocorr, InitRec, MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                             ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge,
                             EggFertParam, mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95,
                             EstGearSelAtAge, Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax,
                             ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, FMortByYear)
# Specify inputs to catch curve model
nSampYrs = 2
ObsFreqAtAge <- t(res$RandObsCatchAgeFreq_CombSex[seq(nYears-nSampYrs+1,nYears,1),])
FirstSampYr = Years[nYears] - nSampYrs + 1
MaxModelAge = 50
Ages = 0:MaxModelAge
nAges = length(Ages)
AgeAtFullRecruit = 6 # knife edge selectivity (for now)
NatMort = 4.22/40
FirstRecDevYr = 1993
LastSampYr = FirstSampYr + nSampYrs - 1
LastRecDevYr = LastSampYr - AgeAtFullRecruit - (nSampYrs - 2) # assuming 2 years to reliably estimate a rec dev
nRecDevs = LastRecDevYr - FirstRecDevYr + 1
SigmaR = 0.6
Ages = 0:MaxModelAge
nAges = length(Ages)
LastSampYr = FirstSampYr + nSampYrs - 1
LastRecDevYr = LastSampYr - AgeAtFullRecruit - (nSampYrs - 2) # assuming 2 years to reliably estimate a rec dev
nRecDevs = LastRecDevYr - FirstRecDevYr + 1
# RTMB data input list
dat <- list(nSampYrs=nSampYrs,
            FirstSampYr=FirstSampYr,
            ObsFreqAtAge=ObsFreqAtAge,
            AgeAtFullRecruit=AgeAtFullRecruit,
            MaxModelAge=MaxModelAge,
            Ages=Ages,
            nAges=nAges,
            NatMort=NatMort,
            FirstRecDevYr=FirstRecDevYr,
            LastRecDevYr=LastRecDevYr,
            nRecDevs=nRecDevs,
            SigmaR=SigmaR)
# RTMB param input list
params <- list(lnFMort=log(0.1),
               lnRecDevs=rep(0,nRecDevs))
# Fit catch curve model with RTMB
FittedRes=GetVarRecCatchCurveResults(dat, params)
FMortByYear[nYears] # True F
FittedRes$SummaryRes$EstFMort # Est F
# plot fit of catch curve to observed data
par(mfcol=c(nSampYrs,1), mar=c(4,3,2,2))
for (t in 1:nSampYrs) {
  plot(Ages, ObsFreqAtAge[,t]/sum(ObsFreqAtAge[,t]),"o", main=FirstSampYr+t-1, ylab="Prop")
  lines(Ages, FittedRes$RTMBReport$ExpCatchPropAtAge[,t], col="blue")
}
# plot simulated vs estimted rec devs
par(mfrow=c(1,1))
plot(c(Years,2021), res$random_dev, ylim=c(0,6), "l", col="blue")
lines(FittedRes$RTMBReport$EstRecDevYrs, FittedRes$SummaryRes$EstRecDevs[,1], "l", col="red")
lines(FittedRes$RTMBReport$EstRecDevYrs, FittedRes$SummaryRes$EstRecDevs[,2], "l", col="pink")
lines(FittedRes$RTMBReport$EstRecDevYrs, FittedRes$SummaryRes$EstRecDevs[,3], "l", col="pink")
}
