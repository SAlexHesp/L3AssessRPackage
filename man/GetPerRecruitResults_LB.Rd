% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L3Analyses.R
\name{GetPerRecruitResults_LB}
\alias{GetPerRecruitResults_LB}
\title{Get outputs from length-based per recruit analysis across a range of fishing mortality values}
\usage{
GetPerRecruitResults_LB(
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
  Current_F,
  Output_Opt
)
}
\arguments{
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

\item{ReprodScale}{<- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)'}

\item{ReprodPattern}{reproductive pattern, 1=gonochoristic, 2=protogynous (female to male sex change) hermaphroditism,
3=protandrous hermaphroditism (male to female sex change)}

\item{InitRatioFem}{proportion of fish that are females at hatching}

\item{FinalSex_L50}{logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)}

\item{FinalSex_L95}{logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)}

\item{EstSexRatioAtLen}{NA  # sex ratio at length, inputted as vector}

\item{EggFertParam}{<- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects}

\item{mat_L50}{logistic length at maturity parameter (set to NA if directly inputting proportion mature at length)}

\item{mat_L95}{logistic length at maturity parameter (set to NA if directly inputting proportion mature at length)}

\item{EstMatAtLen}{vector of proportion mature at length (set to NA if using length at maturity parameters)}

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

\item{Current_F}{estimated current fishing mortality}

\item{Output_Opt}{1 = standard, 2=include expected length and weight distributions at F (slower)}
}
\value{
yield per recruit for combined sexes (YPR), female, male and combined sex spawning potential ratio
(Fem_SPR, Mal_SPR, CombSex_SPR), female and male length at age (FemLenAtAge, MalLenAtAge),
female and male weight at length (FemWtAtLen, MalWtAtLen), female and male proportion mature at length
(FemPropMatAtLen, MalPropMatAtLen), unfished female and male numbers per recruit at age and length
(Unfish_FemNPerRec, Unfish_MalNPerRec), fished female and male numbers per recruit at age and length (Fish_FemNPerRec,
Fish_MalNPerRec), unfished female and male biomass per recruit at age and length (Unfish_FemSpBiomPerRecAtAge,
Unfish_MalSpBiomPerRecAtAge), fished female and male biomass per recruit at age and length (Fish_FemSpBiomPerRecAtAge,
Fish_MalSpBiomPerRecAtAge), unfished female and male biomass at age (Unfish_FemBiomAtAge, Unfish_MalBiomAtAge),
fished female and male biomass at age (Fish_FemBiomAtAge, Fish_MalBiomAtAge), female and male catch biomass
(FemCatchBiom, MalCatchBiom), equilibrium recruitment for either Beverton-Holt or Ricker relationship (Eq_Rec),
equilibrium catch (Eq_Catch), female and male and spawning biomass (Eq_FemSpBiom, Eq_MalSpBiom), equilibrium relative
female, male and combined sex spawning biomass (Eq_FemRelSpBiom, Eq_MalRelSpBiom, Eq_CombSexRelSpBiom),
female and male gear selectivity at length (FemGearSelAtLen, MalGearSelAtLen), female and male retention at length
(FemRetProbAtLen, MalRetProbAtLen), female and male selectivity of landings at length (FemSelLandAtLen, MalSelLandAtLen),
female and male selectivity of discards at length (FemSelDiscAtLen, MalSelDiscAtLen),
female and male discard mortality at length (FemDiscFAtLen, MalDiscFAtLen), female and male mortality associated with landings
at length (FemLandFAtLen, MalLandFAtLen), female and male fishing mortality at length (FemFAtLen, MalFAtLen),
female and male total mortality at length (FemZAtLen, MalZAtLen), fishing mortality values evaluated (FishMort = seq(0,2,0.01)),
maximum yield per recruit (YPRmax), fishing mortality associated with maxeqCatch (F_MSY or FmaxeqCatch),
biomass target at 1.2maxeqCatch, (BMSY_Targ), biomass threshold at maxeqCatch (BMSY_Thresh), biomass threshold at 0.5maxeqCatch (BMSY_Lim),
YPR vs FishMort (YPRResults), Eq_Catch vs FishMort (Eq_CatchResults), Eq_FemRelSpBiom vs FishMort (Fem_SPRResults),
Eq_MalRelSpBiom vs FishMort (Mal_SPRResults), Eq_FemRelSpBiom vs FishMort (Eq_FemRelSpBiomResults),
Eq_MalRelSpBiom vs FishMort (Eq_MalRelSpBiomResults), Eq_CombSexRelSpBiom vs FishMort (Eq_CombSexRelSpBiomResults),
Eq_Rec vs FishMort (Eq_RecResults)
}
\description{
This function provides outputs associated with per recruit analysis, and an
extended form of this analysis with a Beverton-Holt stock recruitment relationship to account
for potential impacts of fishing on recruitment. Outputs are provided for a range of
fishing mortality values, including the current, estimated value.
}
\examples{
# Non-hermaphroditic species
MaxModelAge <- 20 # maximum age considered by model, years
TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
MaxLen = 800
LenInc = 20
lbnd = seq(0,MaxLen - LenInc, LenInc)
ubnd = lbnd + LenInc
midpt = lbnd + (LenInc/2)
nLenCl = length(midpt)
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
vbK <- c(0.2, 0.2) # year-1 - von Bertalanffy growth model parameters - Females, males
tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
GrowthParams = data.frame(Linf=Linf, vbK=vbK, tzero=tzero)
RefnceAges = NA
# GrowthParams c(Linf, vbK, tzero) single sex von Bertalanffy, or data.frame(Linf=Linf, vbK=vbK, tzero=tzero),
#' both sexes von Bertalanffy, or c(y1, y2, a, b) single sex Schnute, or data.frame(y1=y1, y2=y2, a=a, b=b), both sexes Schnute
CVSizeAtAge = c(0.05,0.05)
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
mat_L50 <- c(250, 250) # females, males - Logistic length (mm) at maturity relationship parameters
mat_L95 <- c(300, 300) # females, males - Logistic length (mm) at maturity relationship parameters
EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
                          EstMalMatAtLen=NA) # maturity at length, inputted as values in data frame
sel_L50 <- c(250, 250) # females, males - Logistic length selectivity relationship parameters
sel_L95 <- c(300, 300) # females, males - Logistic length selectivity relationship parameters
EstGearSelAtLen <- data.frame(EstFemGearSelAtLen=NA,
                              EstMalGearSelAtLen=NA)
ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
ret_L50 <- c(250, 250) # females, males - Logistic fish retention at length parameters
ret_L95 <- c(300, 300) # females, males - Logistic fish retention at length parameters
EstRetenAtLen <- data.frame(EstFemRetenAtLen=NA,
                            EstMalRetenAtLen=NA)
DiscMort <- 0.25 # discard mortality (e.g. 50\% released fish die = 0.5)
Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
NatMort = 4.22 / MaxModelAge # natural mortality  (year-1)
Current_F = 0.2
Output_Opt = 1 # 1=standard output, 2=with added length and weight outputs (slower)
Res=GetPerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                            RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                            ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                            EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                            ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, Current_F, Output_Opt)
# Hermaphroditic species
MaxModelAge <- 100 # maximum age considered by model, years
TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
MaxLen = 1400
LenInc = 20
lbnd = seq(0,MaxLen - LenInc, LenInc)
ubnd = lbnd + LenInc
midpt = lbnd + (LenInc/2)
nLenCl = length(midpt)
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf <- c(682, 982) # mm - von Bertalanffy growth model parameters - Females, males
vbK <- c(0.14, 0.08) # year-1 - von Bertalanffy growth model parameters - Females, males
tzero <- c(0.06, -0.48) # years - von Bertalanffy growth model parameters - Females, males
GrowthParams = data.frame(Linf=Linf, vbK=vbK, tzero=tzero)
RefnceAges = NA
# GrowthParams c(Linf, vbK, tzero) single sex von Bertalanffy, or data.frame(Linf=Linf, vbK=vbK, tzero=tzero),
#' both sexes von Bertalanffy, or c(y1, y2, a, b) single sex Schnute, or data.frame(y1=y1, y2=y2, a=a, b=b), both sexes Schnute
CVSizeAtAge = c(0.05,0.05)
lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters
ln_lenwt_a <- -11.017 # for log-log relationship
lenwt_b <- 3.041 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
WLrel_Type <- 2 # 1=power, 2=log-log relationship
EstWtAtLen <- data.frame(EstFemWtAtLen=NA,
                         EstMalWtAtLen=NA) # weight at length, inputted as values in data frame
ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
ReprodPattern <- 2 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
InitRatioFem <- 1.0 # Ratio of females to males at recruitment age
FinalSex_L50 <- 821 # Logistic sex change relationship parameters (inflection point)
FinalSex_L95 <- 930 # Logistic sex change relationship parameters (95\% of max probability)
EstSexRatioAtLen <- NA  # sex ratio at length inputted as vector
EggFertParam <- 0.5 # (NA or from ~0.2-1) NA = no effect, ~0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
mat_L50 <- c(653, 0) # females, males - Logistic length (mm) at maturity relationship parameters
mat_L95 <- c(950, 1) # females, males - Logistic length (mm) at maturity relationship parameters
EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
                          EstMalMatAtLen=NA) # maturity at length, inputted as values in data frame
sel_L50 <- c(500, 500) # females, males - Logistic length selectivity relationship parameters
sel_L95 <- c(600, 600) # females, males - Logistic length selectivity relationship parameters
EstGearSelAtLen <- data.frame(EstFemGearSelAtLen=NA,
                              EstMalGearSelAtLen=NA)
ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
ret_L50 <- c(499, 499) # females, males - Logistic fish retention at length parameters
ret_L95 <- c(500, 500) # females, males - Logistic fish retention at length parameters
EstRetenAtLen <- data.frame(EstFemRetenAtLen=NA,
                            EstMalRetenAtLen=NA)
DiscMort <- 0.25 # discard mortality (e.g. 50\% released fish die = 0.5)
Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
NatMort <- 4.22 / 70 # natural mortality  (year-1)
RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
PlotOpt <- 0 # 0=all plots, 1=len at-age, 2=wt at length, 3=fem mat/sel/ret at length, 4=mal mat/sel/ret at length,
# 5=fem F at length, 6=mal F at length, 7=fem rel surv, 8=mal rel surv, 9=fem biom at age, 10=fem biom at age,
# 11=ypr/eq catch, 12=fem SPR/Brel, 13=mal SPR/Brel, 14=eq recruit
Current_F = 0.05
Output_Opt = 1 # 1=standard output, 2=with added length and weight outputs (slower)
Res=GetPerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                            RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                            ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                            EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                            ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, Current_F, Output_Opt)
}
