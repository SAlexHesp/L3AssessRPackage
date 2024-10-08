% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L3Analyses.R
\name{SimLenFreqData_DynMod_LB}
\alias{SimLenFreqData_DynMod_LB}
\title{Simulate length composition data from a length-based dynamic model with recruitment variation}
\usage{
SimLenFreqData_DynMod_LB(
  SimAnnSampSize,
  nYears,
  lnSigmaR,
  autocorr,
  InitRec,
  MaxModelAge,
  TimeStep,
  lbnd,
  ubnd,
  midpt,
  nLenCl,
  GrowthCurveType,
  GrowthParams,
  RefnceAges,
  CVSizeAtAge,
  lenwt_a,
  ln_lenwt_a,
  lenwt_b,
  WLrel_Type,
  EstWtAtLen,
  ReprodScale,
  ReprodPattern,
  InitRatioFem,
  FinalSex_L50,
  FinalSex_L95,
  EstSexRatioAtLen,
  EggFertParam,
  mat_L50,
  mat_L95,
  EstMatAtLen,
  sel_L50,
  sel_L95,
  EstGearSelAtLen,
  ret_Pmax,
  ret_L50,
  ret_L95,
  EstRetenAtLen,
  DiscMort,
  Steepness,
  SRrel_Type,
  NatMort,
  FMortByYear
)
}
\arguments{
\item{SimAnnSampSize}{specified sample size (for both sexes)}

\item{nYears}{number of years of deviations required}

\item{lnSigmaR}{natural logarithm of standard deviation of annual recruitment deviations}

\item{autocorr}{level of autocorrelation of annual recruitment deviations}

\item{InitRec}{initial recruitment (R0)}

\item{MaxModelAge}{maximum age considered by model}

\item{TimeStep}{model time step (in y)}

\item{lbnd}{lower bounds of length classes}

\item{ubnd}{upper bounds of length classes}

\item{midpt}{mid points bounds of length classes}

\item{nLenCl}{number of length classes}

\item{GrowthCurveType}{1=von Bertalanffy, 2=Schnute}

\item{GrowthParams}{growth parameters of either von Bertalanffy or Schnute model}

\item{RefnceAges}{reference ages for Schnute model, set to NA if using von Bertalanffy model}

\item{CVSizeAtAge}{coefficient of variation for size at age}

\item{lenwt_a}{weight-length parameter}

\item{ln_lenwt_a}{weight-length parameter}

\item{lenwt_b}{weight-length parameter}

\item{WLrel_Type}{1=power, 2=log-log}

\item{EstWtAtLen}{user-specified weights at lengths}

\item{ReprodScale}{<- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)}

\item{ReprodPattern}{reproductive pattern, 1=gonochoristic, 2=protogynous (female to male sex change) hermaphroditism,
3=protandrous hermaphroditism (male to female sex change)}

\item{InitRatioFem}{proportion of fish that are females at hatching}

\item{FinalSex_L50}{logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)}

\item{FinalSex_L95}{logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)}

\item{EstSexRatioAtLen}{NA  # sex ratio at length, inputted as vector}

\item{EggFertParam}{<- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects}

\item{mat_L50}{logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)}

\item{mat_L95}{logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)}

\item{EstMatAtLen}{vector of proportion mature at length (set to NA if using age at maturity parameters)}

\item{sel_L50}{logistic parameter for gear selectivity curve}

\item{sel_L95}{logistic parameter for gear selectivity curve}

\item{EstGearSelAtLen}{gear selectivity curve inputted as vector}

\item{ret_Pmax}{logistic parameter for fish retention curve}

\item{ret_L50}{logistic parameter for fish retention curve}

\item{ret_L95}{logistic parameter for fish retention curve}

\item{EstRetenAtLen}{retention curve inputted as vector}

\item{DiscMort}{proportion of fish that die following capture and release}

\item{Steepness}{steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship}

\item{SRrel_Type}{1 = Beverton-Holt, 2=Ricker stock-recruitment relationship}

\item{NatMort}{natural mortality}

\item{FMortByYear}{annual fully-selected fishing mortality}
}
\value{
catches in numbers for females, males and combined sexes (Catch_Fem, Catch_Mal, Catch_CombSex),
catch biomass for combined sexes (Catch_Biom), female, male and combined sex spawning biomass (FemSpBiom,
MalSpBiom, CombSexSpBiom), annual recruitment (AnnRecruit), recruitment deviations (random_dev),
simulated sample sizes for females and males, given combined sex imput and model assumptions (AnnSimSampSize_Fem,
AnnSimSampSize_Mal), randomly-generated length frequencies, by year, for females, males and combined sexes
(RandObsCatchLenFreq_Fem, RandObsCatchLenFreq_Mal, RandObsCatchLenFreq_CombSex), mean and 95 percent
confidence limits for each year, for simulated length data (MeanLengthStats)
}
\description{
This function provides simulated length composition data from a length-based dynamic model
with specified annual fishing and (natural mortality), annual recruitment with specified process
error and autocorrelation.
}
\examples{
# Simulate length data from dynamic model, from specified biology and exploitation history
MaxModelAge <- 60 # maximum age considered by model, years
TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
MaxLen <- 1400
LenInc <- 20
lbnd <- seq(0,MaxLen - LenInc, LenInc)
ubnd <- lbnd + LenInc
midpt <- lbnd + (LenInc/2)
nLenCl <- length(midpt)
GrowthCurveType <- 1 # 1 = von Bertalanffy, 2 = Schnute
Linf <- c(1000, 1000) # mm - von Bertalanffy growth model parameters - Females, males
vbK <- c(0.12, 0.12) # year-1 - von Bertalanffy growth model parameters - Females, males
tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
GrowthParams <- data.frame(Linf=Linf, vbK=vbK, tzero=tzero)
RefnceAges <- NA
CVSizeAtAge = c(0.06,0.06)
lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
ln_lenwt_a <- NA # for log-log relationship
lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
WLrel_Type <- 1 # 1=power, 2=log-log relationship
EstWtAtLen <- data.frame(EstFemWtAtLen=NA,
                         EstMalWtAtLen=NA) # weight at length, inputted as values in data frame
ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
InitRatioFem <- 0.5 # Ratio of females to males at recruitment age/length
FinalSex_L50 <- NA # Logistic sex change relationship parameters (inflection point)
FinalSex_L95 <- NA # Logistic sex change relationship parameters (95\% of max probability)
EstSexRatioAtLen <- NA  # sex ratio at length inputted as vector
EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, ~0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
mat_L50 <- c(350, 350) # females, males - Logistic length (mm) at maturity relationship parameters
mat_L95 <- c(400, 400) # females, males - Logistic length (mm) at maturity relationship parameters
EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
                          EstMalMatAtLen=NA) # maturity at length, inputted as values in data frame
sel_L50 <- c(350, 350) # females, males - Logistic length selectivity relationship parameters
sel_L95 <- c(400, 400) # females, males - Logistic length selectivity relationship parameters
EstGearSelAtLen <- data.frame(EstFemGearSelAtLen=NA,
                              EstMalGearSelAtLen=NA)
ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
ret_L50 <- c(490, 490) # females, males - Logistic fish retention at length parameters
ret_L95 <- c(500, 500) # females, males - Logistic fish retention at length parameters
EstRetenAtLen <- data.frame(EstFemRetenAtLen=NA,
                            EstMalRetenAtLen=NA)
DiscMort <- 0.5 # discard mortality (e.g. 50\% released fish die = 0.5)
Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
NatMort <- 4.22 / 41 # natural mortality  (year-1)
# example exploitation history - initial in increase F, reduction in F then stable F period
nPeriods <- 3 # number of periods of mortality trends
InitYr <- 1970 # first year of mortality
PeriodEndYr <- c(1995,2000,2020) # periods when mortality changes
InitYr_FMort <- 0.33*NatMort # initial mortality, at beginning of first period
PeriodEndYr_FMort <- c(3*NatMort,0.67*NatMort,0.67*NatMort)
autocorr <- 0.2
randerr_sd <- 0.1
set.seed(123)
res=SimExploitationHistory_DynMod(nPeriods, InitYr, PeriodEndYr, InitYr_FMort,
                                  PeriodEndYr_FMort, autocorr, randerr_sd)
# # example exploitation history -  stable F with no error
# nPeriods = 1 # number of periods of mortality trends
# InitYr = 1970 # first year of mortality
# InitYr_FMort = NatMort # initial mortality, at beginning of first period
# PeriodEndYr_FMort = NatMort
# FMort_SD = 0
# autocorr = 0
# randerr_sd = 0
# set.seed(123)
res=SimExploitationHistory_DynMod(nPeriods, InitYr, PeriodEndYr, InitYr_FMort,
                                  PeriodEndYr_FMort, autocorr, randerr_sd)
# Get random length data
nYears <- res$nYears
FMortByYear <- res$rand_FMort
SimAnnSampSize <- 1000
lnSigmaR <- 0.6
autocorr <- 0.3
InitRec <- 1000
Res=SimLenFreqData_DynMod_LB(SimAnnSampSize, nYears, lnSigmaR, autocorr, InitRec,
                                MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                                ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                                EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                                ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, FMortByYear)
}
