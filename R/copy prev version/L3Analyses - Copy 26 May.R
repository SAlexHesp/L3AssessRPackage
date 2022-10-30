#' @useDynLib L3Assess, .registration = TRUE
#' @importFrom Rcpp sourceCpp
 NULL
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#
# This does something
# @export
# hello <- function() {
#   print("Hello, world!")
# }

# seasonal growth curve
# @export
# CalcGrowth <- function(Linf, vbK, tzero, C, tc, Ages, nAges) {
#
#   EstLengths = CalcSeasonalGrowth_cpp(Linf, vbK, tzero, C, tc, Ages, nAges)
#
#   return(EstLengths)
# }

# Alex Hesp, July 2020
# Catch curve and per recruit analysis package

# **************************
# General plotting functions
# **************************


#******************************************
# Catch curve analyses - size-based methods
#******************************************

# Alex Hesp 16/9/2020. Length-based catch curve method allowing for two alternative selectivity functions.
# This method applies a length-transition matrix to calculate growth and associated uncertainty. Length
# composition data are generated from length-based per recruit calculations (i.e. applying survival function and
# Baranov catch equation). Model is fitted assuming the length composition conforms to a multinomial distribution.

#' Calculate gillnet selectivity at length, by mesh and overall (assuming equal fishing intensity
#' among meshes), using the method of Kirkwood and Walker (1986) and associated model parameters
#' (theta1 and theta2).
#'
#' @param params theta1 (Number), theta2 (Number), MeshSize_mm (vector), nMeshes (Number),
#' nLenCl (Number), midpt (Vector)
#' @return Selectivity at length for each mesh, SelAtLengthForMesh (data frame), overall length-based
#' selectivity SelAtLength (vector)
#' @examples
#' MaxLen = 1000
#' LenInc = 20
#' lbnd = seq(0,MaxLen - LenInc, LenInc)
#' midpt = lbnd + (LenInc/2)
#' nLenCl = length(midpt)
#' MeshSize_mm = c(115, 127)
#' theta1 = 6
#' theta2 = 3000
#' nMeshes = length(MeshSize_mm)
#' CalcGillnetSelectivity(theta1, theta2, MeshSize_mm, nMeshes, nLenCl, midpt)
#' @export
CalcGillnetSelectivity <- function(theta1, theta2, MeshSize_mm, nMeshes, nLenCl, midpt) {

  # calculate alpha and beta parameters for selectivity schedule
  alpha_beta = rep(NA,nMeshes)
  alpha_beta = MeshSize_mm * theta1
  beta = -0.5 * (alpha_beta - sqrt(alpha_beta * alpha_beta + 4 * theta2))
  alpha = alpha_beta / beta

  SelAtLengthForMesh = data.frame(matrix(nrow=nMeshes, ncol=nLenCl))
  colnames(SelAtLengthForMesh) = midpt

  SumSelAtLengthForMesh = rep(0,nLenCl)
  for (m in 1:nMeshes) {
    SelAtLengthForMesh[m,] = ((midpt / alpha_beta[m]) ^ alpha[m]) *
      exp(alpha[m] - midpt / beta[m])

    SumSelAtLengthForMesh = SumSelAtLengthForMesh + SelAtLengthForMesh[m,]
  }
  SelAtLength = SumSelAtLengthForMesh / max(SumSelAtLengthForMesh)

  SelectResults = list(SelAtLengthForMesh = SelAtLengthForMesh,
                       SelAtLength = SelAtLength)

  return(SelectResults)

}

#' Calculate selectivity at length using an asymptotic, logistic curve.
#'
#' @keywords internal
#' @param params L50 (Number), L95 (Number), nLenCl (Number), midpt (Vector)
#' @return Selectivity at length SelAtLength (vector)
CalcLogisticSelectivity <- function(L50, L95, nLenCl, midpt) {

  SelAtLength = rep(0,nLenCl)
  SelAtLength = 1 / (1 + exp(-log(19) * (midpt - L95) / (L95-L50)))

  SelectResults = list(SelAtLength = SelAtLength)

  return(SelectResults)

}

#' Calculate size distribution of 1+ year old recruits, given mean length
#' at age 1 and specified cv for lengths at age, assuming a normal distribution.
#'
#' @keywords internal
#' @param params MeanSizeAtAge (Vector), CVSizeAtAge (Number)
#' @return Selectivity at length SelAtLength (vector)
CalcSizeDistOfRecruits <- function(MeanSizeAtAge, CVSizeAtAge, lbnd, ubnd, nLenCl) {

  # Mean length and SD of 0+ recruits, at 1 year of age
  MeanLenRec <- MeanSizeAtAge[1]
  SDAgeOneRecruits = MeanLenRec * CVSizeAtAge

  RecLenDist = rep(0,nLenCl)
  RecLenDist = pnorm(ubnd, mean=MeanLenRec, sd=SDAgeOneRecruits, lower.tail = T) -
    pnorm(lbnd, mean=MeanLenRec, sd=SDAgeOneRecruits,lower.tail = T)

  RecLenDist = RecLenDist / sum(RecLenDist)

  return(RecLenDist)

}

#' Visualise length at age data, when applying a length transition matrix and specified
#' von Bertalanffy growth parameters and associated variation (CV for length at age),
#' i.e. this function simulates growth trajectories of individual fish
#'
#' @param params Number of fish for which to simulate growth, nFish (Number),
#' MeanSizeAtAge (Vector), AnnGrowthSizeInc (Vector), CVSizeAtAge (Number), midpt (Vector), nLenCl (Number)
#' @return plot of mean growth curve and simulated, individual fish growth trajectories
#' @examples
#' set.seed(123)
#' MaxAge = 20
#' Ages = 1:MaxAge
#' MaxLen = 1000
#' LenInc = 20
#' lbnd = seq(0,MaxLen - LenInc, LenInc)
#' ubnd = lbnd + LenInc
#' midpt = lbnd + (LenInc/2)
#' nLenCl = length(midpt)
#' Linf = 800
#' vbK = 0.2
#' tzero = 0
#' CVSizeAtAge = 0.1
#' nFish = 100
#' MeanSizeAtAge = Linf * (1 - exp (-vbK * (Ages - tzero))) # e.g. as calculated from a von Bertalanffy growth curve
#' MeanEndingLength = midpt + (Linf - midpt) * (1 - exp(-vbK)) # e.g. as calculated from a von Bertalanffy growth curve
#' AnnGrowthSizeInc = MeanEndingLength-midpt # amount of annual growth with respect to initial length
#' VisualiseGrowthApplyingLTM(nFish, MeanSizeAtAge, AnnGrowthSizeInc, CVSizeAtAge, midpt, nLenCl)
#' @export
VisualiseGrowthApplyingLTM <- function (nFish, MeanSizeAtAge, AnnGrowthSizeInc, CVSizeAtAge, midpt, nLenCl) {

  # a bit of code to check growth transition matrix is likely to
  # give sensible results - not part of catch curve model

  # note, using the normal distribution allows for the possibility
  # of negative growth of fish, though on average, this won't occur

  # LTM = CalcLTM(AnnGrowthSizeInc, CVSizeAtAge, lbnd, midpt, ubnd, nLenCl)
  LTM = CalcLTM_cpp(AnnGrowthSizeInc, CVSizeAtAge, lbnd, midpt, ubnd, nLenCl)
  plot(Ages, MeanSizeAtAge,  "l", cex.main=1, frame.plot=F, xlim=c(0,20), ylim=c(0,1200),
       xlab = "Age", ylab = "Length")

  for (j in 1:nFish) { # fish
    FishLenAtAge = rep(0,MaxAge)
    for (k in 1:MaxAge) { # age

      if (k ==1) {
        rndlen = rnorm(1,MeanSizeAtAge[1],CVSizeAtAge[1]*CVSizeAtAge)
        # rndlen = rnorm(1,MeanSizeAtAge[1],CVSizeAtAge)
        FishLenAtAge[k] = rndlen
        rndlc = which(lbnd<rndlen & ubnd>rndlen)
      } else {

        # get cum. distn for final length
        cdistn = rep(0,nLenCl)
        for (i in 1:nLenCl) {
          if (i == 1) {
            cdistn[i] = LTM[i,rndlc]
          } else {
            cdistn[i] = cdistn[i-1] + LTM[i,rndlc]
          }
        }

        # get random number from uniform distn, between 0 and 1
        rndnum = runif(1,0,1)

        rndlc = min(which(cdistn>=rndnum))

        FishLenAtAge[k] = midpt[rndlc]
      } # else
    } # age
    points(Ages, FishLenAtAge, pch=16, cex=0.6, col=j)
  } # fish
}


#' Get statistical outputs from a length-based catch curve
#'
#' This function fits a length-based catch curve with length-based selectivity, based on
#' either 1) specified gillnet selectity parameters from the Kirkwood and Walker (1986) model
#' or 2) estimated selectivity parameters for an asymptotic logistic selectivity curve. The model is
#' to a sample of fish length frequency data, by minimising the negative log-likelihood associated
#' with the parameters and data, using nlminb. It provides various statistical outputs in include
#' convergence statistics, parameter estimates and associated 95% confidence limits and associated
#' variance-covariance matrix, calculated using the MASS package.
#'
#' @param params Vector of model parameters in log space (params) to be estimated, including several or all of
#' the following depending on selectivity option: initial fishing mortality, InitFishMort, logistic selectivity
#' parameters InitL50, InitL95, von Bertalanffy growth parameters, InitLinf, InitvbK, and coefficient of
#' variation, i.e. CV, for size at age (InitCVSizeAtAge), (params), Selectivity model type, 1=gillnet selectivity
#' model with fixed parameters (Kirkwood and Walker, 1986), 2=asymptotic logistic selectivity curve with estimated
#' parameters (SelectivityType), vector of observed catches in each length class (ObsCatchFreqAtLen), vectors
#' of lower bounds, upper bounds and midpoints of each length class (lbnd, ubnd, midpt), values for gillnet selectivity parameters
#' mesh sizes and number of meshes  (set to NA if for SelectivityType=2), values of fixed parameters when applying
#' gillnet selectivity including L50, L95, Linf, vbK and CVSizeAtAge, maximum age in analysis (must be the same or
#' greater than observed maximum age), natural mortality, NatMort.
#' @return negative log-likelihood (nll), nlminb convergence diagnostic (convergence)
#' sample size (SampleSize), growth parameter estimates with lower and upper 95%
#' confidence limits (ParamEst), point estimates for estimated parameters (params)
#' and variance-covariance matrix (vcov.Params), selectivity at length (SelAtLength),
#' growth curve (MeanSizeAtAge), midpoint of each length class (midpt), mean length after 1 year from
#' growth curve, given initial length (MeanEndingLength), mean change in length after 1 year,
#' from initial length - note, assuming normal a distribution allows for possibility of negative growth
#' if above asyptotic length (AnnGrowthSizeInc), length distribution of 1+ year old recruits (RecLenDist),
#' expected catches, at length (ExpCatchAtLen), proportion of catch at each length (ExpCatchPropInLenClass),
#' observed catch at length data (ObsCatchFreqAtLen)
#' @examples
#' SampleSize=10000
#' MaxAge = 20
#' NatMort = 0.2
#' FishMort = 0.2
#' MaxLen = 1000
#' LenInc = 20
#' lbnd = seq(0,MaxLen - LenInc, LenInc)
#' ubnd = lbnd + LenInc
#' midpt = lbnd + (LenInc/2)
#' SelectivityType=2 # 1=gillnet selectivity model (Kirkwood and Walker, 1986), 2=asymptotic logistic selectivity curve
#' theta1=NA # selectivity
#' theta2=NA # selectivity
#' MeshSize_mm=NA # selectivity
#' nMeshes=NA # selectivity
#' L50 = 300 # selectivity
#' L95 = 400 # selectivity
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.08
#' MeanSizeAtAge = Linf * (1 - exp (-vbK * (1:MaxAge))) # e.g. as calculated from a von Bertalanffy growth curve
#' MeanEndingLength = midpt + (Linf - midpt) * (1 - exp(-vbK)) # e.g. as calculated from a von Bertalanffy growth curve
#' AnnGrowthSizeInc = MeanEndingLength-midpt # amount of annual growth with respect to initial length
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, NatMort, FishMort, MaxLen, LenInc, SelectivityType, theta1, theta2,
#'                          MeshSize_mm, nMeshes, L50, L95, MeanSizeAtAge, AnnGrowthSizeInc, CVSizeAtAge)
#' lbnd=Res$lbnd
#' midpt=Res$midpt
#' ubnd=Res$ubnd
#' ObsCatchFreqAtLen=Res$ObsCatchFreqAtLen
#' InitFishMort = 0.5 # specify starting parameters
#' InitL50 = 300
#' InitL95 = 500
#' params = log(c(InitFishMort, InitL50, InitL95))
#' res=GetLengthBasedCatchCurveResults(params, SelectivityType, ObsCatchFreqAtLen,
#'                                     lbnd, ubnd, midpt, theta1, theta2, MeshSize_mm, nMeshes,
#'                                     Linf, vbK, CVSizeAtAge, MaxAge, NatMort)
#' @export
GetLengthBasedCatchCurveResults <- function (params, SelectivityType, ObsCatchFreqAtLen,
                                             lbnd, ubnd, midpt, theta1, theta2, MeshSize_mm, nMeshes,
                                             Linf, vbK, CVSizeAtAge, MaxAge, NatMort)
{
  nlmb <- nlminb(params, CalcObjFunc_LengthBasedCatchCurve,
                 gradient = NULL, hessian = TRUE)
  (hess.out = optimHess(nlmb$par, CalcObjFunc_LengthBasedCatchCurve))
  (vcov.Params = solve(hess.out))
  (ses = sqrt(diag(vcov.Params)))
  EstFMort = exp(c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) *
                     ses[1]))
  if (SelectivityType == 1) {
    ParamEst = t(data.frame(FMort = round(EstFMort, 2)))
    colnames(ParamEst) = c("Estimate", "lw_95%CL",
                           "up_95%CL")
  }
  if (SelectivityType == 2) {
    EstL50 = exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) *
                     ses[2]))
    EstL95 = exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) *
                     ses[3]))
    ParamEst = t(data.frame(FMort = round(EstFMort, 3), SelL50 = round(EstL50, 3), SelL95 = round(EstL95, 3)))
    colnames(ParamEst) = c("Estimate", "lw_95%CL",
                           "up_95%CL")
  }

  # store some diagnostic outputs from model
  params = nlmb$par
  L50=NA
  L95=NA
  CatchCurveType=1 #1=length-based, 2=age and length based
  res = LengthAndAgeBasedCatchCurvesCalcs(params, CatchCurveType, SelectivityType, lbnd, ubnd, midpt,
                                          theta1, theta2, MeshSize_mm, nMeshes, L50, L95,
                                          Linf, vbK, CVSizeAtAge, MaxAge, NatMort)

  SampleSize = sum(ObsCatchFreqAtLen)
  nll = nlmb$objective
  convergence = nlmb$convergence

  Results = list(nll = nll,
                 convergence = convergence,
                 SampleSize = SampleSize,
                 ParamEst = ParamEst,
                 params = nlmb$par,
                 vcov.Params = vcov.Params,
                 SelAtLength=res$SelAtLength,
                 MeanSizeAtAge=res$MeanSizeAtAge,
                 midpt=res$midpt,
                 MeanEndingLength=res$MeanEndingLength,
                 AnnGrowthSizeInc=res$AnnGrowthSizeInc,
                 RecLenDist=res$RecLenDist,
                 ExpCatchAtLen=res$ExpCatchAtLen,
                 ExpCatchPropInLenClass=res$ExpCatchPropInLenClass,
                 ObsCatchFreqAtLen=res$ObsCatchFreqAtLen)

  return(Results)
}

# CalcCatches_AgeAndLengthBasedCatchCurves <- function(params, NatMort, RecLenDist, MaxAge, nLenCl, midpt,
#                                                      SelAtLength, LTM) {
#
#   FishMort = exp(params[1])
#
#   # per recruit numbers surviving after natural mortality
#   Fish_FemNPerRec <- rep(0,nLenCl)
#   Fish_FemNPerRecAtAge <- data.frame(matrix(nrow = MaxAge, ncol = nLenCl))
#   colnames(Fish_FemNPerRecAtAge) <- midpt
#   Fish_FemNPerRecAtAge[1:length(1:MaxAge),1:nLenCl] = 0
#   Fish_FemSurvPerRecAtAge = Fish_FemNPerRecAtAge
#   Catch = Fish_FemNPerRecAtAge
#   CatchAtLen = rep(0,nLenCl)
#
#   Fish_FemNPerRecAtAge[1,] = RecLenDist
#   Fish_FemNPerRec = Fish_FemNPerRec + Fish_FemNPerRecAtAge[1,]
#
#   # fishing mortality, calcuated based on selectivity associated with research fishing
#   FAtLen = SelAtLength * FishMort
#
#   # total mortality experienced by the stock, due to commercial fishing
#   ZAtLen = NatMort + FAtLen
#
#   # calculate catch at length for timestep (in numbers)
#   Catch[1,] = Fish_FemNPerRecAtAge[1,] * (FAtLen / ZAtLen) * (1 - exp(-ZAtLen))
#   CatchAtLen = CatchAtLen + Catch[1,]
#
#   # apply mortality to calculate survival
#   for (t in 2:MaxAge) {
#     if (t < MaxAge) {
#       Fish_FemSurvPerRecAtAge[t,] = Fish_FemNPerRecAtAge[t-1,] * exp(-ZAtLen)
#     } else {
#       Fish_FemSurvPerRecAtAge[t,] = Fish_FemNPerRecAtAge[t-1,] * exp(-ZAtLen) /
#         (1 - exp(-ZAtLen))
#     }
#
#     # apply growth - looping - slower
#     # for (ii in 1:nLenCl) {  # starting length class
#     #   for (i in 1:nLenCl) { # ending length class
#     #     Fish_FemNPerRecAtAge[t,i] = Fish_FemNPerRecAtAge[t,i] + (Fish_FemSurvPerRecAtAge[t,ii] * LTM[i,ii])
#     #   }
#     Fish_FemSurvPerRecAtAge_temp = as.vector(Fish_FemSurvPerRecAtAge[t,])
#     temp = GrowFishForTimeStep_cpp(LTM, Fish_FemSurvPerRecAtAge_temp, nLenCl)
#     Fish_FemNPerRecAtAge[t,] = temp
#
#     # apply growth - matrix multiplication - faster
#     # M1=t(t(Fish_FemSurvPerRecAtAge[t,]))
#     # M2=t(LTM[,])
#     # Fish_FemNPerRecAtAge[t,] = as.numeric(M1 %*% M2)
#     # Fish_FemNPerRec = Fish_FemNPerRec + Fish_FemNPerRecAtAge[t,]
#
#     # calculate catch at length for timestep (in numbers)
#     Catch[t,] = Fish_FemNPerRecAtAge[t,] * (FAtLen / ZAtLen) * (1 - exp(-ZAtLen))
#
#     CatchAtLen = CatchAtLen + Catch[t,]
#
#   } # t
#
#   if (!is.numeric(sum(Catch))) {
#     cat("Problem - OperatingModel. sum(Catch) not numeric")
#   }
#
#   Res=list(Catch=Catch,
#            CatchAtLen=CatchAtLen,
#            Fish_FemNPerRec = Fish_FemNPerRec)
#
#   return(Res)
#
# }

#' Get NLL for length based catch curve, for optimisation
#'
#' @keywords internal
#' @param params (Vector)
#' @return negative log-likelihood (NLL)
CalcObjFunc_LengthBasedCatchCurve <- function(params) {
  # get NLL for length based catch curve, for optimisation

  CatchCurveType=1 #1=length-based, 2=age and length based
  L50=NA
  L95=NA
  Res = LengthAndAgeBasedCatchCurvesCalcs(params, CatchCurveType, SelectivityType, lbnd, ubnd, midpt,
                                          theta1, theta2, MeshSize_mm, nMeshes, L50, L95,
                                          Linf, vbK, CVSizeAtAge, MaxAge, NatMort)

  ExpCatchPropInLenClass = Res$ExpCatchPropInLenClass
  NLL = CalcNLLMargLengthComposition(ObsCatchFreqAtLen, ExpCatchPropInLenClass)
  cat("NLL", NLL, "exp(params)", exp(params), "\n")

  return(NLL)

}

#' Get NLL for length and age-based catch curve, for optimisation
#'
#' @keywords internal
#' @param params (Vector)
#' @return negative log-likelihood (NLL)
CalcObjFunc_LengthAndAgeBasedCatchCurve <- function(params) {
  # get NLL for length based catch curve, for optimisation

  L50=NA
  L95=NA
  Linf=NA
  vbK=NA
  CVSizeAtAge=NA
  CatchCurveType=2
  Res = LengthAndAgeBasedCatchCurvesCalcs(params, CatchCurveType, SelectivityType, lbnd, ubnd, midpt,
                                          theta1, theta2, MeshSize_mm, nMeshes, L50, L95,
                                          Linf, vbK, CVSizeAtAge, MaxAge, NatMort)

  L50_Pen = Res$L50_Pen
  L95_Pen = Res$L95_Pen

  # get NLL for marginal length composition
  ExpCatchPropInLenClass = Res$ExpCatchPropInLenClass
  Length_NLL = CalcNLLMargLengthComposition(ObsCatchFreqAtLen, ExpCatchPropInLenClass)

  # get NLL for age at length observations
  ExpCatchPropAtAge = as.matrix(Res$ExpCatchPropAtAge)
  ExpCatchPropLengthGivenAge = as.matrix(Res$ExpCatchPropLengthGivenAge)

  nLenCl = length(midpt)
  CondAgeAtLengthNLL = CalcNLLCondAgeAtLength_cpp(nLenCl, MaxAge, ObsCatchFreqAtLengthAndAge, ExpCatchPropLengthGivenAge, ExpCatchPropAtAge)
  NLL = Length_NLL + CondAgeAtLengthNLL + L50_Pen + L95_Pen

  cat("NLL", NLL, "Length_NLL", Length_NLL, "CondAgeAtLengthNLL", CondAgeAtLengthNLL, "L50_Pen", L50_Pen,
      "L95_Pen",L95_Pen, "\n")
  cat("exp(params)", exp(params), "\n")
  cat("", "\n")

  return(NLL)

}

#' Get NLL for marginal length composition data, for optimisation
#'
#' @keywords internal
#' @param ObsCatchFreqAtLen (Vector), ExpCatchPropInLenClass (Vector)
#' @return negative log-likelihood (Length_NLL)
CalcNLLMargLengthComposition <- function(ObsCatchFreqAtLen, ExpCatchPropInLenClass) {

  Length_NLL = -sum(ObsCatchFreqAtLen * log(ExpCatchPropInLenClass + 1E-4))

  return(Length_NLL)
}

#' Do length and length and age-based catch curve calcs
#'
#' @keywords internal
#' @param params params (Vector), SelectivityType (Number), ObsCatchFreqAtLen (Vector),
#' lbnd (Vector), ubnd (Vector), midpt (Vector), theta1 (Number), theta2 (Number),
#' MeshSize_mm (Number), nMeshes (Number), L50 (Number), L95 (Number),
#' Linf (Number), vbK (Number), CVSizeAtAge (Number), MaxAge (Number), NatMort (Number)
#' @return NLL (Number), selectivity at length (SelAtLength),
#' von Bertalanffy growth curve (MeanSizeAtAge), midpt, mean length after 1 year from
#' growth curve, given initial length (MeanEndingLength), mean change in length after 1 year,
#' from initial length - note, assuming normal a distribution allows for possibility of negative growth
#' if above asyptotic length (AnnGrowthSizeInc), length distribution of 1+ year old recruits (RecLenDist),
#' expected catches, at length (ExpCatchAtLen), proportion of catch at each length (ExpCatchPropInLenClass),
#' expected catches, at age (ExpCatchAtAge), proportion of catch at each age (ExpCatchPropAtAge),
#' relative catches at length and age (Catch), expected relative catches at each length and age
#' (ExpCatchPropLengthGivenAge), penalty values for L50 and L95 selectivity parameters (L50_Pen, L95_Pen)
#' @export
LengthAndAgeBasedCatchCurvesCalcs <- function (params, CatchCurveType, SelectivityType, lbnd, ubnd, midpt,
                                               theta1, theta2, MeshSize_mm, nMeshes, L50, L95,
                                               Linf, vbK, CVSizeAtAge, MaxAge, NatMort)
{

  nLenCl = length(midpt)

  if (SelectivityType == 1) { # gillnet selectivity
    if (CatchCurveType == 2) { # age and length based catch curve
      Linf = exp(params[2])
      vbK = exp(params[3])
      CVSizeAtAge = exp(params[4])
    }
  }

  if (SelectivityType == 2) { # logistic selectivity
    L50 = exp(params[2])
    L95 = exp(params[3])
    if (CatchCurveType == 2) { # age and length based catch curve
      Linf = exp(params[4])
      vbK = exp(params[5])
      CVSizeAtAge = exp(params[6])
    }
  }

    # calculate L95 penalty
    #cat("L95",L95,"max(ubnd)",max(ubnd),'\n')
    L95_Pen = 0
    L50_Pen = 0
    if (SelectivityType == 2) { # logistic selectivity
      if (L95 > Linf) {
        L95_Pen = 100.0 * (L95 - Linf)^2
        L95 = Linf
      }
      # calculate L50 penalty
      if (L50 > (L95 - 2)) {
        L50_Pen = 100.0 * ((L95 - 2) - L50)^2
        L50 = L95 - 2
      }
    }

  if (SelectivityType == 1) { # Kirkwood & Walker method
    SelAtLengthResults = CalcGillnetSelectivity(theta1, theta2,
                                                MeshSize_mm, nMeshes, nLenCl, midpt)
    SelAtLength = unlist(SelAtLengthResults$SelAtLength)
  }
  if (SelectivityType == 2) { # logistic selectivity (asymptotic)
    SelAtLengthResults = CalcLogisticSelectivity(L50, L95,
                                                 nLenCl, midpt)
    SelAtLength = SelAtLengthResults$SelAtLength
  }

  MeanSizeAtAge = Linf * (1 - exp (-vbK * (1:MaxAge))) # e.g. as calculated from a von Bertalanffy growth curve

  MeanEndingLength = midpt + (Linf - midpt) * (1 - exp(-vbK)) # e.g. as calculated from a von Bertalanffy growth curve

  AnnGrowthSizeInc = MeanEndingLength-midpt # amount of annual growth with respect to initial length

  RecLenDist = CalcSizeDistOfRecruits(MeanSizeAtAge, CVSizeAtAge,
                                      lbnd, ubnd, nLenCl)

  # LTM = CalcLTM(AnnGrowthSizeInc, CVSizeAtAge, lbnd, midpt, ubnd, nLenCl)
  LTM = CalcLTM_cpp(AnnGrowthSizeInc, CVSizeAtAge, lbnd, midpt, ubnd, nLenCl)

  #CatchCurveResults = CalcCatches_AgeAndLengthBasedCatchCurves(params, NatMort, RecLenDist, MaxAge, nLenCl,
  #                                                             midpt, SelAtLength, LTM)
  CatchCurveResults = CalcCatches_AgeAndLengthBasedCatchCurves_cpp(params, NatMort, RecLenDist, MaxAge, nLenCl,
                                                                   midpt, SelAtLength, LTM)

  # expected prop at length
  CatchAtLen = CatchCurveResults$CatchAtLen
  ExpCatchAtLen = CatchAtLen/sum(CatchAtLen)
  ExpCatchPropInLenClass = ExpCatchAtLen / sum(ExpCatchAtLen)

  # expected prop at age
  CatchAtLengthAndAge = CatchCurveResults$Catch
  ExpCatchAtAge = rowSums(CatchAtLengthAndAge)
  ExpCatchPropAtAge = ExpCatchAtAge / sum(ExpCatchAtAge)

  # expected prop at age, given length
  Catch = CatchCurveResults$Catch
  ExpCatchPropLengthGivenAge <- data.frame(matrix(nrow = MaxAge, ncol = nLenCl))
  colnames(ExpCatchPropLengthGivenAge) <- midpt
  for (i in 1:MaxAge) {
    ExpCatchPropLengthGivenAge[i,] = Catch[i,] / sum(Catch[i,])
  }

  Results = list(SelAtLength=SelAtLength,
                 MeanSizeAtAge=MeanSizeAtAge,
                 midpt=midpt,
                 MeanEndingLength=MeanEndingLength,
                 AnnGrowthSizeInc=AnnGrowthSizeInc,
                 RecLenDist=RecLenDist,
                 ExpCatchAtLen=ExpCatchAtLen,
                 ExpCatchPropInLenClass=ExpCatchPropInLenClass,
                 ExpCatchAtAge=ExpCatchAtAge,
                 ExpCatchPropAtAge=ExpCatchPropAtAge,
                 Catch=Catch,
                 ExpCatchPropLengthGivenAge=ExpCatchPropLengthGivenAge,
                 L50_Pen=L50_Pen,
                 L95_Pen=L95_Pen)

  return(Results)
}

#' Get statistical outputs from a length and age-based catch curve
#'
#' This function fits a length-based catch curve with length-based selectivity, based on
#' either 1) specified gillnet selectity parameters from the Kirkwood and Walker (1986) model
#' or 2) estimated selectivity parameters for an asymptotic logistic selectivity curve. The model is
#' fitted to a sample of fish length and age data, by minimising the overall negative log-likelihood,
#' including the NLL associated with the marginal length composition and a conditional age at length NLL,
#' given the parameters (selectivity, growth and mortality) and data, using nlminb.
#' It provides various statistical outputs in include convergence statistics, parameter estimates
#' and associated 95% confidence limits and associated variance-covariance matrix, calculated using
#' the MASS package.
#'
#' @param Vector of model parameters in log space (params) to be estimated, including several or all of
#' the following depending on selectivity option: initial fishing mortality, InitFishMort, logistic selectivity
#' parameters InitL50, InitL95, von Bertalanffy growth parameters, InitLinf, InitvbK, and coefficient of
#' variation, i.e. CV, for size at age (InitCVSizeAtAge), (params), Selectivity model type, 1=gillnet selectivity
#' model with fixed parameters (Kirkwood and Walker, 1986), 2=asymptotic logistic selectivity curve with estimated
#' parameters (SelectivityType), vector of observed catches in each length class (ObsCatchFreqAtLen),
#' matrix of observed catches in each age class (from 1:MaxAge) and length class (ObsCatchFreqAtLengthAndAge), vectors
#' of lower bounds, upper bounds and midpoints of each length class (lbnd, ubnd, midpt), values for gillnet selectivity parameters
#' mesh sizes and number of meshes  (set to NA if for SelectivityType=2), values of fixed parameters when applying
#' gillnet selectivity including L50, L95, Linf, vbK and CVSizeAtAge, maximum age in analysis (must be the same or
#' greater than observed maximum age), natural mortality, NatMort.
#' @return negative log-likelihood (nll), nlminb convergence diagnostic (convergence)
#' sample size (SampleSize), growth parameter estimates with lower and upper 95%
#' confidence limits (ParamEst), point estimates for estimated parameters (params)
#' and variance-covariance matrix (vcov.Params), selectivity at length (SelAtLength),
#' growth curve (MeanSizeAtAge), midpoint of each length class (midpt), mean length after 1 year from
#' growth curve, given initial length (MeanEndingLength), mean change in length after 1 year,
#' from initial length - note, assuming normal a distribution allows for possibility of negative growth
#' if above asyptotic length (AnnGrowthSizeInc), length distribution of 1+ year old recruits (RecLenDist),
#' expected catches, at length (ExpCatchAtLen), proportion of catch at each length (ExpCatchPropInLenClass),
#' expected catches, at age (ExpCatchAtAge), expected relative catch at length and age (Catch), proportion of
#' catch at each age and length (ExpCatchPropLengthGivenAge), observed catch at length data (ObsCatchFreqAtLen),
#' observed catch at length and age data (ObsCatchFreqAtLen)
#' @examples
#' set.seed(123)
#' SampleSize=10000
#' MaxAge = 20
#' NatMort = 0.2
#' FishMort = 0.2
#' MaxLen = 1000
#' LenInc = 20
#' midpt = seq(0,MaxLen - LenInc, LenInc) + (LenInc/2)
#' SelectivityType=2 # 1=gillnet selectivity model (Kirkwood and Walker, 1986), 2=asymptotic logistic selectivity curve
#' theta1=NA # selectivity
#' theta2=NA # selectivity
#' MeshSize_mm=NA # selectivity
#' nMeshes=NA # selectivity
#' L50 = 300 # selectivity
#' L95 = 500 # selectivity
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.08
#' MeanSizeAtAge = Linf * (1 - exp (-vbK * (1:MaxAge))) # e.g. as calculated from a von Bertalanffy growth curve
#' MeanEndingLength = midpt + (Linf - midpt) * (1 - exp(-vbK)) # e.g. as calculated from a von Bertalanffy growth curve
#' AnnGrowthSizeInc = MeanEndingLength-midpt # amount of annual growth with respect to initial length
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, NatMort, FishMort, MaxLen, LenInc,
#'                          SelectivityType, theta1, theta2, MeshSize_mm, nMeshes, L50, L95,
#'                          MeanSizeAtAge, AnnGrowthSizeInc, CVSizeAtAge)
#' lbnd=Res$lbnd
#' midpt=Res$midpt
#' ubnd=Res$ubnd
#' ObsCatchFreqAtLen = Res$ObsCatchFreqAtLen
#' ObsCatchFreqAtLengthAndAge = as.matrix(Res$ObsCatchFreqAtLengthAndAge)
#' CatchCurveType=2 # age and length based
#' InitFishMort = 0.5 # specify starting parameters
#' InitL50 = 300
#' InitL95 = 500
#' InitLinf = 800
#' InitvbK = 0.25
#' InitCVSizeAtAge = 0.1
#' params = log(c(InitFishMort, InitL50, InitL95, InitLinf, InitvbK, InitCVSizeAtAge))
#' GetAgeAndLengthBasedCatchCurveResults(params, SelectivityType, ObsCatchFreqAtLen, ObsCatchFreqAtLengthAndAge,
#'                                       lbnd, ubnd, midpt, theta1, theta2, MeshSize_mm, nMeshes, L50, L95,
#'                                       Linf, vbK, CVSizeAtAge, MaxAge, NatMort)
#' @export
GetAgeAndLengthBasedCatchCurveResults <- function (params, SelectivityType, ObsCatchFreqAtLen, ObsCatchFreqAtLengthAndAge,
                                                   lbnd, ubnd, midpt, theta1, theta2, MeshSize_mm, nMeshes, L50, L95,
                                                   Linf, vbK, CVSizeAtAge, MaxAge, NatMort)
{


  nlmb <- nlminb(params, CalcObjFunc_LengthAndAgeBasedCatchCurve,
                 gradient = NULL, hessian = TRUE)

  (hess.out = optimHess(nlmb$par, CalcObjFunc_LengthAndAgeBasedCatchCurve))
  (vcov.Params = solve(hess.out))
  (ses = sqrt(diag(vcov.Params)))
  EstFMort = exp(c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) *
                     ses[1]))

  if (SelectivityType == 1) { # gillnet selectivity (specified parameters)
    EstLinf = exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) *
                      ses[2]))
    Estk = exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) *
                   ses[3]))
    EstCV = exp(c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) *
                    ses[4]))
    ParamEst = t(data.frame(FMort = round(EstFMort, 3), Linf = round(EstLinf, 3),
                            vbK = round(Estk, 3), CV = round(EstCV, 3)))
    colnames(ParamEst) = c("Estimate", "lw_95%CL", "up_95%CL")
    }
  if (SelectivityType == 2) { # logistic selectivity (estimated parameters)
    EstL50 = exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) *
                     ses[2]))
    EstL95 = exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) *
                     ses[3]))
    EstLinf = exp(c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) *
                      ses[4]))
    Estk = exp(c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) *
                   ses[5]))
    EstCV = exp(c(nlmb$par[6], nlmb$par[6] + c(-1.96, 1.96) *
                    ses[6]))
    ParamEst = t(data.frame(FMort = round(EstFMort, 3), SelL50 = round(EstL50, 3), SelL95 = round(EstL95, 3),
                            Linf = round(EstLinf, 3), vbK = round(Estk, 3), CV = round(EstCV, 3)))
    colnames(ParamEst) = c("Estimate", "lw_95%CL", "up_95%CL")
  }

  # store some diagnostic outputs from model
  L50=NA
  L95=NA
  Linf=NA
  vbK=NA
  CVSizeAtAge=NA
  CatchCurveType=2
  params = nlmb$par
  res = LengthAndAgeBasedCatchCurvesCalcs(params, CatchCurveType, SelectivityType, lbnd, ubnd, midpt,
                                          theta1, theta2, MeshSize_mm, nMeshes, L50, L95,
                                          Linf, vbK, CVSizeAtAge, MaxAge, NatMort)
  SampleSize = sum(ObsCatchFreqAtLen)
  nll = nlmb$objective
  convergence = nlmb$convergence
  names(res)
  Results = list(nll = nll,
                 convergence = convergence,
                 SampleSize = SampleSize,
                 ParamEst = ParamEst,
                 params = nlmb$par,
                 vcov.Params = vcov.Params,
                 SelAtLength=res$SelAtLength,
                 MeanSizeAtAge=res$MeanSizeAtAge,
                 midpt=res$midpt,
                 MeanEndingLength=res$MeanEndingLength,
                 AnnGrowthSizeInc=res$AnnGrowthSizeInc,
                 RecLenDist=res$RecLenDist,
                 ExpCatchAtLen=res$ExpCatchAtLen,
                 ExpCatchPropInLenClass=res$ExpCatchPropInLenClass,
                 ExpCatchAtAge=res$ExpCatchAtAge,
                 ExpCatchPropAtAge=res$ExpCatchPropAtAge,
                 Catch=res$Catch,
                 ExpCatchPropLengthGivenAge=res$ExpCatchPropLengthGivenAge,
                 ObsCatchFreqAtLen=ObsCatchFreqAtLen,
                 ObsCatchFreqAtLengthAndAge=ObsCatchFreqAtLengthAndAge)

  return(Results)
}

#' Simulate some age and length frequency data with specified selectivity, mortality and growth parameters.
#'
#' @param params SampleSize (Number), MaxAge (Number), FishMort (Number), MaxLen (Number),
#' LenInc (Number), SelectivityType (Number) 1=gillnet selectivity, 2=asymptotic logistic
#' selectivity curve, theta1 (Number), theta2 (Number), MeshSize_mm (Number), nMeshes
#' (Number), L50 (Number), L95 (Number), MeanSizeAtAge (Vector), AnnGrowthSizeInc (Vector),
#' CVSizeAtAge (Number)
#' @return Observed catch frequency at age, ObsCatchFreqAtAge (Vector), Observed catch proportions at
#' age, ObsRelCatchAtAge (Vectors), Observed catch frequency at length, ObsCatchFreqAtLen (Vector),
#' observed catch proportions at length, ObsRelCatchAtLen (Vectors), Observed catch frequency at length
#' and age, ObsCatchFreqAtLengthAndAge (Matrix), age classes, ObsAgeCl (Vector), length class
#' mid points, ObslenClMidPt (Vector)
#' @examples
#' set.seed(123)
#' SampleSize=10000
#' MaxAge = 20
#' NatMort = 0.2
#' FishMort = 0.2
#' MaxLen = 1000
#' LenInc = 20
#' midpt = seq(0,MaxLen - LenInc, LenInc) + (LenInc/2)
#' SelectivityType=2 # 1=gillnet selectivity model (Kirkwood and Walker, 1986), 2=asymptotic logistic selectivity curve
#' theta1=NA # selectivity
#' theta2=NA # selectivity
#' MeshSize_mm=NA # selectivity
#' nMeshes=NA # selectivity
#' L50 = 300 # selectivity
#' L95 = 500 # selectivity
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.08
#' MeanSizeAtAge = Linf * (1 - exp (-vbK * (1:MaxAge))) # e.g. as calculated from a von Bertalanffy growth curve
#' MeanEndingLength = midpt + (Linf - midpt) * (1 - exp(-vbK)) # e.g. as calculated from a von Bertalanffy growth curve
#' AnnGrowthSizeInc = MeanEndingLength-midpt # amount of annual growth with respect to initial length
#' SimLenAndAgeFreqData(SampleSize, MaxAge, NatMort, FishMort, MaxLen, LenInc,
#'                      SelectivityType, theta1, theta2, MeshSize_mm, nMeshes, L50, L95,
#'                      MeanSizeAtAge, AnnGrowthSizeInc, CVSizeAtAge)
#' @export
SimLenAndAgeFreqData <- function(SampleSize, MaxAge, NatMort, FishMort, MaxLen, LenInc,
                                 SelectivityType, theta1, theta2, MeshSize_mm, nMeshes, L50, L95,
                                 MeanSizeAtAge, AnnGrowthSizeInc, CVSizeAtAge) {

  lbnd = seq(0,MaxLen - LenInc, LenInc)
  ubnd = lbnd + LenInc
  midpt = lbnd + (LenInc/2)
  nLenCl = length(midpt)

  # Simulate some observed catch at length data from catch curve model
  # (assuming data should conform to a multinomial distribution)

  # selectivity
  if (SelectivityType == 1) {
    SelAtLengthResults = CalcGillnetSelectivity(theta1, theta2, MeshSize_mm, nMeshes, nLenCl, midpt)
    SelAtLength = unlist(SelAtLengthResults$SelAtLength)
    SelAtLengthForMesh = unlist(SelAtLengthResults$SelAtLengthForMesh)
  }
  if (SelectivityType == 2) {
    SelAtLengthResults = CalcLogisticSelectivity(L50, L95, nLenCl, midpt)
    SelAtLength = SelAtLengthResults$SelAtLength
  }

  # size distribution of 1+ recruits
  RecLenDist = CalcSizeDistOfRecruits(MeanSizeAtAge, CVSizeAtAge, lbnd, ubnd, nLenCl)

  # LTM exploration results
  # LTM = CalcLTM(AnnGrowthSizeInc, CVSizeAtAge, lbnd, midpt, ubnd, nLenCl)
  LTM = CalcLTM_cpp(AnnGrowthSizeInc, CVSizeAtAge, lbnd, midpt, ubnd, nLenCl)

  # survival at length
  params=log(FishMort)
  # CatchCurveResults = CalcCatches_AgeAndLengthBasedCatchCurves(params, NatMort, RecLenDist, MaxAge,
  #                                                              nLenCl, midpt, SelAtLength, LTM)
  CatchCurveResults = CalcCatches_AgeAndLengthBasedCatchCurves_cpp(params, NatMort, RecLenDist,
                                                                   MaxAge, nLenCl, midpt, SelAtLength, LTM)

  # generate catches at age
  CatchAtLengthAndAge = CatchCurveResults$Catch
  CatchAtAge = rowSums(CatchAtLengthAndAge)
  ExpCatchPropAtAge = CatchAtAge / sum(CatchAtAge)
  ObsCatchFreqAtAge = as.vector(rmultinom(1, SampleSize, ExpCatchPropAtAge))
  ObsRelCatchAtAge = ObsCatchFreqAtAge/sum(ObsCatchFreqAtAge)

  # generate size at age data
  ObsCatchFreqAtLengthAndAge <- data.frame(matrix(nrow = MaxAge, ncol = nLenCl))
  colnames(ObsCatchFreqAtLengthAndAge) <- midpt
  for (i in 1:MaxAge) {
    ObsCatchFreqAtLengthAndAge[i,] = rmultinom(1, ObsCatchFreqAtAge[i], CatchAtLengthAndAge[i,])
  }

  # generate data for individual fish (age classes and mid points of length classes)
  ObsAgeCl = rep(NA, SampleSize)
  ObslenClMidPt = rep(NA, SampleSize)
  strt=1; fnsh=0
  for (i in 1:MaxAge) {
    for (j in 1:nLenCl) {

      x=ObsCatchFreqAtLengthAndAge[i,j] # number of fish in current length and age class
      if(x>0) {
        fnsh=strt+x-1
        # cat("i",i,"j",j,"strt",strt,"fnsh",fnsh,"#",ObsCatchFreqAtLengthAndAge[i,j],'\n')
        ObsAgeCl[strt:fnsh]=i
        ObslenClMidPt[strt:fnsh]=midpt[j]
        strt=strt+x
      }
    }
  }

  ObsCatchFreqAtLen = colSums(ObsCatchFreqAtLengthAndAge)
  ObsRelCatchAtLen = ObsCatchFreqAtLen/sum(ObsCatchFreqAtLen)

  Results = list(lbnd=lbnd,
                 midpt=midpt,
                 ubnd=ubnd,
                 ObsCatchFreqAtAge = ObsCatchFreqAtAge,
                 ObsRelCatchAtAge = ObsRelCatchAtAge,
                 ObsCatchFreqAtLen = ObsCatchFreqAtLen,
                 ObsRelCatchAtLen = ObsRelCatchAtLen,
                 ObsCatchFreqAtLengthAndAge = ObsCatchFreqAtLengthAndAge,
                 ObsAgeCl=ObsAgeCl,
                 ObslenClMidPt=ObslenClMidPt)

  return(Results)

}


#' Produce plot of fitted length-based catch curve to length composition data
#'
#' @param params params (Vector), SelectivityType (Number), ObsCatchFreqAtLen (Vector),
#' lbnd (Vector), ubnd (Vector), midpt (Vector), theta1 (Number), theta2 (Number),
#' MeshSize_mm (Number), nMeshes (Number), L50 (Number), L95 (Number),
#' Linf (Number), vbK (Number), CVSizeAtAge (Number), MaxAge (Number), NatMort (Number),
#' Plotting options including main plot label, set to NA for defaults (MainLabel), x and y
#' axis labels (xaxis_lab, yaxis_lab), max and intervals for x axis (xmax,xint) and y axis
#' (ymax, yint), option to plot 95 percent confidence limits for estiamted line (PlotCLs)
#' @return plot of fitted length-based catch curve to length composition data
#' @examples
#' set.seed(123)
#' SampleSize=10000
#' MaxAge = 20
#' NatMort = 0.2
#' FishMort = 0.2
#' MaxLen = 1000
#' LenInc = 20
#' midpt = seq(0,MaxLen - LenInc, LenInc) + (LenInc/2)
#' SelectivityType=2 # 1=gillnet selectivity model (Kirkwood and Walker, 1986), 2=asymptotic logistic selectivity curve
#' theta1=NA # selectivity
#' theta2=NA # selectivity
#' MeshSize_mm=NA # selectivity
#' nMeshes=NA # selectivity
#' L50 = 300 # selectivity
#' L95 = 500 # selectivity
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.08
#' MeanSizeAtAge = Linf * (1 - exp (-vbK * (1:MaxAge))) # e.g. as calculated from a von Bertalanffy growth curve
#' MeanEndingLength = midpt + (Linf - midpt) * (1 - exp(-vbK)) # e.g. as calculated from a von Bertalanffy growth curve
#' AnnGrowthSizeInc = MeanEndingLength-midpt # amount of annual growth with respect to initial length
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, NatMort, FishMort, MaxLen, LenInc,
#'                          SelectivityType, theta1, theta2, MeshSize_mm, nMeshes, L50, L95,
#'                          MeanSizeAtAge, AnnGrowthSizeInc, CVSizeAtAge)
#' lbnd=Res$lbnd
#' midpt=Res$midpt
#' ubnd=Res$ubnd
#' ObsCatchFreqAtLen=Res$ObsCatchFreqAtLen
#' InitFishMort = 0.5 # specify starting parameters
#' InitL50 = 300
#' InitL95 = 500
#' params = log(c(InitFishMort, InitL50, InitL95))
#' PlotLengthBasedCatchCurveResults(params, SelectivityType, ObsCatchFreqAtLen, lbnd, ubnd, midpt,
#'                                  theta1, theta2, MeshSize_mm, nMeshes, L50, L95,
#'                                  Linf, vbK, CVSizeAtAge, MaxAge, NatMort, MainLabel=NA,
#'                                  xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA,
#'                                  ymax=NA, yint=NA, PlotCLs=TRUE)
#' @export
PlotLengthBasedCatchCurveResults <- function(params, SelectivityType, ObsCatchFreqAtLen, lbnd, ubnd, midpt,
                                             theta1, theta2, MeshSize_mm, nMeshes, L50, L95,
                                             Linf, vbK, CVSizeAtAge, MaxAge, NatMort, MainLabel,
                                             xaxis_lab, yaxis_lab, xmax, xint,
                                             ymax, yint, PlotCLs) {


  res=GetLengthBasedCatchCurveResults(params, SelectivityType, ObsCatchFreqAtLen,
                                      lbnd, ubnd, midpt, theta1, theta2, MeshSize_mm, nMeshes,
                                      Linf, vbK, CVSizeAtAge, MaxAge, NatMort)
  res$ParamEst
  params = res$params
  vcov.params = res$vcov.Params
  ExpCatchAtLen = res$ExpCatchAtLen
  ObsRelCatchAtLen = ObsCatchFreqAtLen/sum(ObsCatchFreqAtLen)

  set.seed(123)
  sims = data.frame(MASS::mvrnorm(n = 200, params, vcov.params))
  EstPropAtLen.sim = data.frame(matrix(nrow = 200, ncol = length(midpt)))

  for (j in 1:200) {
    params = unlist(sims[j,])
    CatchCurveType=1 #1=length-based, 2=age and length based
    Res=LengthAndAgeBasedCatchCurvesCalcs(params, CatchCurveType, SelectivityType, lbnd, ubnd, midpt,
                                          theta1, theta2, MeshSize_mm, nMeshes, L50, L95,
                                          Linf, vbK, CVSizeAtAge, MaxAge, NatMort)

    EstPropAtLen.sim[j,] = unlist(Res$ExpCatchPropInLenClass)
    cat("j",j,'\n')
  }

  EstProp.sim = apply(EstPropAtLen.sim, 2, median)
  EstProp.sim_low = apply(EstPropAtLen.sim, 2, quantile, probs = 0.025)
  EstProp.sim_up = apply(EstPropAtLen.sim, 2, quantile, probs = 0.975)

  if (is.na(xaxis_lab)) xaxis_lab = "Length (mm)"
  if (is.na(yaxis_lab)) yaxis_lab = "Proportion"
  xlims = Get_xaxis_scale(ubnd)
  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint
  ylims = Get_yaxis_scale(ObsRelCatchAtLen)
  if (is.na(ymax)) ymax = ylims$ymax
  if (is.na(yint)) yint = ylims$yint

  plot(midpt, ObsRelCatchAtLen, "p", main=MainLabel, cex.main=1.0, pch=16, cex=0.8, xaxt = "n", yaxt = "n",
       xlab=xaxis_lab,ylab=yaxis_lab, frame=F, xlim=c(0,xmax), ylim=c(0,ymax))
  if (PlotCLs == TRUE) {
    x = c(Res$midpt,rev(Res$midpt)) # using shading for 95% CLs
    y = c(EstProp.sim_low, rev(EstProp.sim_up))
    polygon(x,y,col="pink",border=NA)
  }
  points(midpt, ObsRelCatchAtLen, col="black", pch=16, cex=0.8)
  points(midpt, ExpCatchAtLen, col="red", pch=1, cex=0.8)
  axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
  axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
  axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0,
       cex.axis = 0.8, las = 1)
  axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0,
       cex.axis = 0.8, las = 1)
  Fval = round(exp(params[1]),2)
  Fest = bquote("F =" ~ .(Fval) ~ y^-1)
  if (SelectivityType==1) {
    legend("topright", pch=-1, legend=as.expression(Fest),
           lty="solid",col="black",
           bty='n', cex=0.8,lwd=-1, y.intersp=1.2, adj=0)
  }
  if (SelectivityType==2) {
    L50est=paste("L50 =",round(exp(params[2]),0),"mm")
    L95est=paste("L95 =",round(exp(params[3]),0),"mm")
    legend("topright", pch=-1, legend=c(as.expression(Fest), L50est, L95est),
           lty="solid",col="black",
           bty='n', cex=0.8,lwd=-1, y.intersp=1.2)
  }
  legend("topleft", legend=c("Observed","Estimated"), y.intersp = 1.0, inset=c(0.13,0),
         lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16,1), col=c("black","red"))
}

#' Produce plot of fitted length-based catch curve to age and length composition data
#'
#' This function fits a length-based catch curve with length-based selectivity, based on
#' either 1) specified gillnet selectity parameters from the Kirkwood and Walker (1986) model
#' or 2) estimated selectivity parameters for an asymptotic logistic selectivity curve. The model is
#' fitted to a sample of fish length and age data, by minimising the overall negative log-likelihood,
#' including the NLL associated with the marginal length composition and a conditional age at length NLL,
#' given the parameters (selectivity, growth and mortality) and data, using nlminb.
#' It provides various statistical outputs in include convergence statistics, parameter estimates
#' and associated 95% confidence limits and associated variance-covariance matrix, calculated using
#' the MASS package.
#' @param params params (Vector), SelectivityType (Number), ObsCatchFreqAtLen (Vector),
#' ObsCatchFreqAtLengthAndAge (data frame), lbnd (Vector), ubnd (Vector), midpt (Vector),
#' theta1 (Number), theta2 (Number), MeshSize_mm (Number), nMeshes (Number), L50 (Number), L95 (Number),
#' Linf (Number), vbK (Number), CVSizeAtAge (Number), MaxAge (Number), NatMort (Number)
#' @return Plot of fitted length-based catch curve to age and length composition data
#' @examples
#' set.seed(123)
#' SampleSize=10000
#' MaxAge = 20
#' NatMort = 0.2
#' FishMort = 0.2
#' MaxLen = 1000
#' LenInc = 20
#' midpt = seq(0,MaxLen - LenInc, LenInc) + (LenInc/2)
#' SelectivityType=2 # 1=gillnet selectivity model (Kirkwood and Walker, 1986), 2=asymptotic logistic selectivity curve
#' theta1=NA # selectivity
#' theta2=NA # selectivity
#' MeshSize_mm=NA # selectivity
#' nMeshes=NA # selectivity
#' L50 = 300 # selectivity
#' L95 = 500 # selectivity
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.08
#' MeanSizeAtAge = Linf * (1 - exp (-vbK * (1:MaxAge))) # e.g. as calculated from a von Bertalanffy growth curve
#' MeanEndingLength = midpt + (Linf - midpt) * (1 - exp(-vbK)) # e.g. as calculated from a von Bertalanffy growth curve
#' AnnGrowthSizeInc = MeanEndingLength-midpt # amount of annual growth with respect to initial length
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, NatMort, FishMort, MaxLen, LenInc,
#'                          SelectivityType, theta1, theta2, MeshSize_mm, nMeshes, L50, L95,
#'                          MeanSizeAtAge, AnnGrowthSizeInc, CVSizeAtAge)
#' lbnd=Res$lbnd
#' midpt=Res$midpt
#' ubnd=Res$ubnd
#' ObsCatchFreqAtLen = Res$ObsCatchFreqAtLen
#' ObsCatchFreqAtLengthAndAge = as.matrix(Res$ObsCatchFreqAtLengthAndAge)
#' CatchCurveType=2 # age and length based
#' InitFishMort = 0.5 # specify starting parameters
#' InitL50 = 300
#' InitL95 = 500
#' InitLinf = 800
#' InitvbK = 0.25
#' InitCVSizeAtAge = 0.1
#' params = log(c(InitFishMort, InitL50, InitL95, InitLinf, InitvbK, InitCVSizeAtAge))
#'
#' PlotAgeAndLengthBasedCatchCurveResults(params, SelectivityType, ObsCatchFreqAtLen, ObsCatchFreqAtLengthAndAge,
#'                                        lbnd, ubnd, midpt, theta1, theta2, MeshSize_mm, nMeshes, L50=NA, L95=NA,
#'                                        Linf=NA, vbK=NA, CVSizeAtAge=NA, MaxAge, NatMort)
#' @export
PlotAgeAndLengthBasedCatchCurveResults <- function(params, SelectivityType, ObsCatchFreqAtLen, ObsCatchFreqAtLengthAndAge,
                                                   lbnd, ubnd, midpt, theta1, theta2, MeshSize_mm, nMeshes, L50, L95,
                                                   Linf, vbK, CVSizeAtAge, MaxAge, NatMort) {

  res=GetAgeAndLengthBasedCatchCurveResults(params, SelectivityType, ObsCatchFreqAtLen, ObsCatchFreqAtLengthAndAge,
                                            lbnd, ubnd, midpt, theta1, theta2, MeshSize_mm, nMeshes, L50, L95,
                                            Linf, vbK, CVSizeAtAge, MaxAge, NatMort)

  params = res$params
  ExpCatchAtLen = res$ExpCatchAtLen

  ObsRelCatchAtLen = ObsCatchFreqAtLen/sum(ObsCatchFreqAtLen)
  ymax = 1.2 * (max(ObsRelCatchAtLen))

  par(mfrow=c(1,1))

  # marginal length distribution
  plot(midpt, ObsRelCatchAtLen, xlab="Length", ylab="Prob - catch at length", ylim=c(0,ymax), cex.main=0.8)
  lines(midpt, ExpCatchAtLen, col="blue")
  Fval = round(exp(params[1]),2)
  Fest = bquote("F =" ~ .(Fval) ~ y^-1)

  if (SelectivityType==1) { # gillnet
    Linfest=paste("Linf =",round(exp(params[2]),0),"mm")
    Kval = round(exp(params[3]),3)
    Kest = bquote("k =" ~ .(Kval) ~ y^-1)
    CVest=paste("CV =",round(exp(params[4]),3))
    legend("topleft", pch=-1, legend=c(as.expression(Fest),Linfest,
                                       as.expression(Kest), CVest), lty="solid",col="black",
           bty='n', cex=0.8,lwd=-1, y.intersp=1.2, adj=0)
  }
  if (SelectivityType==2) { # logistic selectivity
    L50est=paste("L50 =",round(exp(params[2]),0),"mm")
    L95est=paste("L95 =",round(exp(params[3]),0),"mm")
    Linfest=paste("Linf =",round(exp(params[4]),0),"mm")
    Kval = round(exp(params[5]),3)
    Kest = bquote("k =" ~ .(Kval) ~ y^-1)
    CVest=paste("CV =",round(exp(params[6]),3))
    legend("topleft", pch=-1, legend=c(as.expression(Fest), L50est, L95est, Linfest,
                                       as.expression(Kest), CVest), lty="solid",col="black", bty='n', cex=0.8,lwd=-1, y.intersp=1.0)

    # plot selectivity
    plot(midpt, res$SelAtLength, xlab="Length", ylab="Selectivity", "l", ylim=c(0,1), cex.main=0.8)
    legend("topleft", pch=-1, legend=c(L50est, L95est), lty="solid",col="black", bty='n',
           cex=0.8,lwd=-1, y.intersp=1.0)
  }


  # predicted conditional proportion of age given length
  nLenCl = length(midpt)
  ExpCatchPropAtAge=res$ExpCatchPropAtAge
  ExpCatchPropLengthGivenAge=res$ExpCatchPropLengthGivenAge
  ObsCatchFreqAtLengthAndAge=res$ObsCatchFreqAtLengthAndAge
  FractDenom = rep(0,nLenCl)
  ExpCatchPropAgeGivenLength <- data.frame(matrix(nrow = MaxAge, ncol = nLenCl))
  colnames(ExpCatchPropAgeGivenLength) <- midpt
  ObsCatchPropAgeAtLength = ExpCatchPropAgeGivenLength
  for (i in 1:nLenCl) {
    for (j in 1:MaxAge) {
      FractDenom[i] = FractDenom[i] + (ExpCatchPropLengthGivenAge[j,i] * ExpCatchPropAtAge[j])
      # cat("i",i,"j",j,"FractDenom[i]",FractDenom[i],'\n')
    }

    for (j in 1:MaxAge) {
      ExpCatchPropAgeGivenLength[j,i] = (ExpCatchPropLengthGivenAge[j,i] * ExpCatchPropAtAge[j]) / FractDenom[i]
    }
  }

  # plot observed and expected age proportions for each length class
  par(mfcol=c(3,3), mar=c(4,4,2,2), oma=c(1,1,1,1))
  for (i in 1:nLenCl) {
    plot(1:MaxAge, ExpCatchPropAgeGivenLength[,i], "l", xlab="Age, y", ylab="Prop.", ylim=c(0,1),
         col="blue", cex.main=0.8, main=paste0("Length bin = ",lbnd[i],"-",ubnd[i]," mm"))
    ObsCatchPropAgeAtLength[,i] = ObsCatchFreqAtLengthAndAge[,i] / sum(ObsCatchFreqAtLengthAndAge[,i])
    lines(1:MaxAge, ObsCatchPropAgeAtLength[,i], col="black")
    legend("topleft", pch=-1, legend=c("Exp","Obs"), lty="solid", col=c("blue","black"), bty='n',
           cex=0.8,lwd=1, y.intersp=1.0)
  }

  # plot growth curve over age and length data

  # generate data for individual fish (age classes and mid points of length classes)
  SampleSize = sum(ObsCatchFreqAtLen)
  ObsAgeCl = rep(NA, SampleSize)
  ObslenClMidPt = rep(NA, SampleSize)
  strt=1; fnsh=0
  for (i in 1:MaxAge) {
    for (j in 1:nLenCl) {

      x=ObsCatchFreqAtLengthAndAge[i,j] # number of fish in current length and age class
      if(x>0) {
        fnsh=strt+x-1
        # cat("i",i,"j",j,"strt",strt,"fnsh",fnsh,"#",ObsCatchFreqAtLengthAndAge[i,j],'\n')
        ObsAgeCl[strt:fnsh]=i
        ObslenClMidPt[strt:fnsh]=midpt[j]
        strt=strt+x
      }
    }
  }

  par(mfrow=c(1,1))
  plot(ObsAgeCl, ObslenClMidPt, xlim=c(0,MaxAge), ylim=c(0,max(ubnd)),pch=16, cex=0.6)
  if (SelectivityType == 1) {
    EstCurve = res$ParamEst[2,1] * (1 - exp(-res$ParamEst[3,1]*(0:MaxAge)))
  }
  if (SelectivityType == 2) {
    EstCurve = res$ParamEst[4,1] * (1 - exp(-res$ParamEst[5,1]*(0:MaxAge)))
  }
  lines(0:MaxAge,EstCurve,col="blue")

}



# Pauly's length-converted catch curve 16/9/2020

#' Return outputs associated with Pauly's (1990) length converted catch curve method, required to fit the model
#' @keywords internal
#' @param ObsCatchFreqAtLen (Vector), MinFreq (Number), lbnd (Vector), midpt (Vector), ubnd (Vector),
#' vbK (Number), Linf (Number)
#' @return PeakLencl (Number), LastLenCl (Number), Observed data to which catch curve is fitted,
#' ObsCatchFreqAtLen2 (Vector), MeanAge_lencl (Vector), ln_n_dt (Vector)
Calcs_PaulyLenConvertCatchCurve <- function(ObsCatchFreqAtLen, MinFreq, lbnd, midpt, ubnd, vbK, Linf) {

  # determine length class at peak frequency. If there is more than one length class
  # with the same "peak" frequency, choose the larger one
  PeakFrequency <- max(ObsCatchFreqAtLen)
  PeakLencl <- which(ObsCatchFreqAtLen==PeakFrequency)
  if (length(PeakLencl)>1) PeakLencl = max(PeakLencl)

  # Select all length classes with > specified number of required observations
  # where length category is below Linf. (Need to remove categories for which
  # ln(n*dl/dt) becomes negative)
  nLenCl = length(lbnd)
  k = PeakLencl
  while (ObsCatchFreqAtLen[k] >= MinFreq & ubnd[k] < Linf)  {
    LastLenCl <- k
    k = k + 1
    if (k > nLenCl)
      stop("problem with while loop")
  }

  # number of length categories for analysis
  nLenCl_CC <- length(ObsCatchFreqAtLen[PeakLencl:LastLenCl])

  # data to which catch curve is fitted
  ObsCatchFreqAtLen2 <- rep(NA,nLenCl)
  ObsCatchFreqAtLen2[PeakLencl:LastLenCl] <- ObsCatchFreqAtLen[PeakLencl:LastLenCl]

  # calc MeanAge_lencl, delta_t and ln_n_dt
  Age_lbndlencl <- rep(0,nLenCl_CC)
  Age_ubndlencl <- rep(0,nLenCl_CC)
  MeanAge_lencl <- rep(0,nLenCl_CC)
  DeltaT_yrs <- rep(0,nLenCl_CC)
  ln_n_dt <- rep(0,nLenCl_CC)
  i=0
  for (j in PeakLencl:LastLenCl) {
    i=i+1
    Age_lbndlencl[i] = (-1 / vbK * log(1 - lbnd[j] / Linf))
    Age_ubndlencl[i] = (-1 / vbK * log(1 - ubnd[j] / Linf))
    MeanAge_lencl[i] = (-1 / vbK * log(1 - midpt[j] / Linf))
    DeltaT_yrs[i] = Age_ubndlencl[i] - Age_lbndlencl[i]
    ln_n_dt[i] = log(ObsCatchFreqAtLen[j] / DeltaT_yrs[i])
  }

  results = list(PeakLencl = PeakLencl,
                 LastLenCl = LastLenCl,
                 ObsCatchFreqAtLen2 = ObsCatchFreqAtLen2,
                 MeanAge_lencl = MeanAge_lencl,
                 ln_n_dt = ln_n_dt)

  return(results)
}

#' Return objective function for Pauly's (1990) length converted catch curve method
#' @keywords internal
#' @param params (Vector)
#' @return NLL (Number)
CalcNLL_PaulyLenConvCatchCurve <- function(params) {

  res=Calcs_PaulyLenConvertCatchCurve(ObsCatchFreqAtLen, MinFreq, lbnd, midpt, ubnd, vbK, Linf)

  nLenCl_CC <- length(res$PeakLencl:res$LastLenCl)
  Est_ln_n_dt <- rep(0,nLenCl_CC)
  SqRes <- rep(0,nLenCl_CC)

  NLL = 0
  i=0
  for (j in res$PeakLencl:res$LastLenCl) {
    i=i+1

    Est_ln_n_dt[i] = (-exp(params[1]) * res$MeanAge_lencl[i]) + params[2]

    SqRes[i] = (res$ln_n_dt[i] - Est_ln_n_dt[i]) *
      (res$ln_n_dt[i] - Est_ln_n_dt[i])

  } # j

  SumSq = sum(SqRes)
  Stdev = sqrt(SumSq/nLenCl_CC)
  NLL = (nLenCl_CC/2.) * log((2. * pi) + (2 * log(Stdev) + 1))

  return(NLL)

} # end function


#' Fit Pauly's (1990) length converted catch curve method and get results
#'
#' @param ObsCatchFreqAtLen (Vector), MinFreq (Number), lbnd (Vector), midpt (Vector), ubnd (Vector),
#' vbK (Number), Linf (Number)
#' @return nll (Number), convergence (Number), SampleSize (number), ParamEst (Data frame),
#' params (Vector), vcov.Params (Data frame), PeakLencl(Number), LastLenCl(Number),
#' data to which catch curve is fitted, ObsCatchFreqAtLen2 (Vector), MeanAge_lencl (Vector),
#' ln_n_dt (Vector), Est_ln_dt (Vector)
#' @examples
#' SampleSize=10000
#' MaxAge = 20
#' NatMort = 0.2
#' FishMort = 0.2
#' MaxLen = 1000
#' LenInc = 20
#' lbnd = seq(0,MaxLen - LenInc, LenInc)
#' ubnd = lbnd + LenInc
#' midpt = lbnd + (LenInc/2)
#' SelectivityType=2 # 1=gillnet selectivity model (Kirkwood and Walker, 1986), 2=asymptotic logistic selectivity curve
#' theta1=NA # selectivity
#' theta2=NA # selectivity
#' MeshSize_mm=NA # selectivity
#' nMeshes=NA # selectivity
#' L50 = 300 # selectivity
#' L95 = 400 # selectivity
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.08
#' MeanSizeAtAge = Linf * (1 - exp (-vbK * (1:MaxAge))) # e.g. as calculated from a von Bertalanffy growth curve
#' MeanEndingLength = midpt + (Linf - midpt) * (1 - exp(-vbK)) # e.g. as calculated from a von Bertalanffy growth curve
#' AnnGrowthSizeInc = MeanEndingLength-midpt # amount of annual growth with respect to initial length
#' res=SimLenAndAgeFreqData(SampleSize, MaxAge, NatMort, FishMort, MaxLen, LenInc, SelectivityType, theta1, theta2,
#'                          MeshSize_mm, nMeshes, L50, L95, MeanSizeAtAge, AnnGrowthSizeInc, CVSizeAtAge)
#' ObsCatchFreqAtLen = as.vector(res$ObsCatchFreqAtLen)
#' MinFreq = 1
#' lbnd = seq(0,MaxLen - LenInc, LenInc)
#' ubnd = lbnd + LenInc
#' nLenCl = length(midpt)
#' Pauly_Z <- 0.5
#' Pauly_yint <- 100
#' params = c(log(Pauly_Z),Pauly_yint)
#' GetPaulyLenConvCatchCurveResults(ObsCatchFreqAtLen, MinFreq, lbnd, midpt, ubnd, vbK, Linf)
#' @export
GetPaulyLenConvCatchCurveResults <- function(ObsCatchFreqAtLen, MinFreq, lbnd, midpt, ubnd, vbK, Linf) {

  # fit model with nlminb optimizer
  nlmb <- nlminb(params, CalcNLL_PaulyLenConvCatchCurve, gradient = NULL, hessian = TRUE)

  # calculate uncertianty for parameter estimates by getting variance-covariance matrix,
  # from fitted model, to get standard errors
  (hess.out = optimHess(nlmb$par, CalcNLL_PaulyLenConvCatchCurve))
  (vcov.Params = solve(hess.out))
  (ses = sqrt(diag(vcov.Params))) # asymptotic standard errors of parameter estimates

  # get approximate 95% confidence limits for the parameter estimates
  EstPauly_Z = exp(c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
  EstPauly_yint = c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2])

  # store results in data frame
  ParamEst = t(data.frame(ZMort=round(EstPauly_Z,2), YInt=round(EstPauly_yint,2)))
  colnames(ParamEst) = c("Estimate","lw_95%CL","up_95%CL")

  # store sample size
  SampleSize = sum(ObsCatchFreqAtLen)

  # store value of objective function
  nll = nlmb$objective

  # store convergence value
  convergence = nlmb$convergence

  # get additional info from fitted model
  res=Calcs_PaulyLenConvertCatchCurve(ObsCatchFreqAtLen, MinFreq, lbnd, midpt, ubnd, vbK, Linf)
  nLenCl_CC <- length(res$PeakLencl:res$LastLenCl)
  Est_ln_n_dt <- rep(0,nLenCl_CC)
  Est_ln_n_dt = (-exp(nlmb$par[1]) * res$MeanAge_lencl) + nlmb$par[2]

  results = list(nll = nll,
                 convergence = convergence,
                 SampleSize = SampleSize,
                 ParamEst = ParamEst,
                 params = nlmb$par,
                 vcov.Params = vcov.Params,
                 PeakLencl = res$PeakLencl,
                 LastLenCl = res$LastLenCl,
                 ObsCatchFreqAtLen2 = res$ObsCatchFreqAtLen2,
                 MeanAge_lencl = res$MeanAge_lencl,
                 ln_n_dt = res$ln_n_dt,
                 Est_ln_n_dt = Est_ln_n_dt)

  return(results)
}

#' Plot data and fit of Pauly's (1990) length converted catch curve method
#'
#' @param ObsCatchFreqAtLen (Vector), MinFreq (Number), lbnd (Vector), midpt (Vector), ubnd (Vector),
#' vbK (Number), Linf (Number)
#' @return plots relative to data, growth and catch curve fit
#' @examples
#' SampleSize=10000
#' MaxAge = 20
#' NatMort = 0.2
#' FishMort = 0.2
#' MaxLen = 1000
#' LenInc = 20
#' lbnd = seq(0,MaxLen - LenInc, LenInc)
#' ubnd = lbnd + LenInc
#' midpt = lbnd + (LenInc/2)
#' SelectivityType=2 # 1=gillnet selectivity model (Kirkwood and Walker, 1986), 2=asymptotic logistic selectivity curve
#' theta1=NA # selectivity
#' theta2=NA # selectivity
#' MeshSize_mm=NA # selectivity
#' nMeshes=NA # selectivity
#' L50 = 300 # selectivity
#' L95 = 400 # selectivity
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.08
#' MeanSizeAtAge = Linf * (1 - exp (-vbK * (1:MaxAge))) # e.g. as calculated from a von Bertalanffy growth curve
#' MeanEndingLength = midpt + (Linf - midpt) * (1 - exp(-vbK)) # e.g. as calculated from a von Bertalanffy growth curve
#' AnnGrowthSizeInc = MeanEndingLength-midpt # amount of annual growth with respect to initial length
#' res=SimLenAndAgeFreqData(SampleSize, MaxAge, NatMort, FishMort, MaxLen, LenInc, SelectivityType, theta1, theta2,
#'                          MeshSize_mm, nMeshes, L50, L95, MeanSizeAtAge, AnnGrowthSizeInc, CVSizeAtAge)
#' ObsCatchFreqAtLen = as.vector(res$ObsCatchFreqAtLen)
#' MinFreq = 1
#' lbnd = seq(0,MaxLen - LenInc, LenInc)
#' ubnd = lbnd + LenInc
#' nLenCl = length(midpt)
#' Pauly_Z <- 0.5
#' Pauly_yint <- 100
#' params = c(log(Pauly_Z),Pauly_yint)
#' PlotPaulyLenConvCatchCurveResults(Ages, MeanSizeAtAge, ObsCatchFreqAtLen, MinFreq, lbnd, midpt, ubnd, vbK, Linf)
#' @export
PlotPaulyLenConvCatchCurveResults <- function(Ages, MeanSizeAtAge, ObsCatchFreqAtLen, MinFreq, lbnd, midpt, ubnd, vbK, Linf) {


  res=GetPaulyLenConvCatchCurveResults(ObsCatchFreqAtLen, MinFreq, lbnd, midpt, ubnd, vbK, Linf)

  # plot data that can be used for the catch curve,
  # relative to full set of length data
  par(mfrow=c(2,2), tck=-0.03)
  plot(midpt, ObsCatchFreqAtLen, ylab="Frequency", xlab="Length", frame.plot=F,
       xlim=c(0,MaxLen))
  points(midpt, res$ObsCatchFreqAtLen2, pch=16, col="blue")
  legend('topleft', col=c("grey","blue"),legend=c("all data","c.crv. data"),
         bty='n', cex=0.8,lwd=1.75)

  # plot growth curve, and overlay range of lengths available for catch curve analysis
  MeanSizeAtAge = Linf * (1 - exp (-vbK * (Ages))) # von Bertalanffy growth curve,
  plot(Ages,MeanSizeAtAge,"l",frame.plot=F, xlim=c(min(Ages),max(Ages)),
       ylim=c(0,MaxLen),col="blue", ylab="Length",xlab="Age")
  abline(h=lbnd[res$PeakLencl],lty="solid")
  abline(h=ubnd[res$LastLenCl],lty="dotted")
  legend('bottomright', col=c("black","black"),lty=c("solid","dotted"),
         legend=c("start of catch curve","end of catch curve"),bty='n', cex=0.8,lwd=1.75)

  # plot catch curve
  plot(res$MeanAge_lencl, res$ln_n_dt,frame.plot=F, xlim=c(0,MaxAge), xlab="Mean Age", ylab="ln_n/dt")
  yvals <- (-exp(res$params[1]) * res$MeanAge_lencl) + res$params[2]
  lines(res$MeanAge_lencl,yvals)
  Zval = round(exp(res$params[1]),2)
  Zest = bquote("Z =" ~ .(Zval) ~ y^-1)
  legend("bottomleft", pch=-1, legend=as.expression(Zest),
         lty="solid",col="black",
         bty='n', cex=0.8,lwd=-1, y.intersp=1.2, adj=0)
}


#*****************************************
# Catch curve analyses - age-based methods
#*****************************************

#' Determine age at full recruitment into the fishery (assuming knife-edge selection).
#'
#' This function calculates the oldest age to be included in data for a linear catch curve analysis,
#' according to a specified minimum frequency.
#' @keywords internal
#' @param params MinFreq, RecAge (Numbers) ObsAgeFreq (Numeric Vector)
#' @return LastAgeForLinearCC (oldest age to be included for linear catch curve analysis)
#' @examples
#' set.seed(123)
#' MaxAge = 40
#' Ages = 1:MaxAge
#' NatMort <- exp(1.46 - (1.01 * (log(MaxAge)))) # i.e. Hoenig's (1983) eqn for fish
#' FMort = 0.1
#' ZMort = FMort + NatMort
#' SelA50 = 6
#' SelA95 = 8
#' SampleSize = 500 # required number of fish for age sample
#' SelAtAge = 1 / (1 + exp(-log(19) * (Ages - SelA50) / (SelA95 - SelA50)))
#' FAtAge = SelAtAge * FMort
#' ZAtAge = NatMort + FAtAge
#' N = numeric(MaxAge)
#' N[1] = 1
#' for (i in 2:MaxAge) {
#'   if (i < MaxAge) {
#'     N[i] = N[i-1] * exp(-ZAtAge[i-1])
#'   } else {
#'     N[i] = N[i-1] * exp(-ZAtAge[i-1] / (1 - exp(-ZAtAge[i])))
#'   }
#' }
#' CatchAtAge = N * (FAtAge / ZAtAge) * (1 - exp(-ZAtAge)) # catch at age
#' PropAtAge = CatchAtAge / sum(CatchAtAge)
#' CatchSample = rmultinom(n=1, size=SampleSize, prob=PropAtAge)
#' ObsAgeFreq = CatchSample[,1]
#' CalcLastAgeForLinearCatchCurve(MinFreq=1, RecAge=8, Ages, ObsAgeFreq)
CalcLastAgeForLinearCatchCurve <- function (MinFreq, RecAge, Ages, ObsAgeFreq)
{
  x = which(Ages == RecAge)
  xx = length(Ages)
  AgesEqualToAndAboveRecAge = Ages[x:xx]
  AgeFreqEqualToAndAboveRecAge = ObsAgeFreq[x:xx]
  if (min(AgeFreqEqualToAndAboveRecAge) >= MinFreq) {
    LastAgeForLinearCC = Ages[xx]
  } else {
    xxx = min(which(AgeFreqEqualToAndAboveRecAge < MinFreq)) - 1
    LastAgeForLinearCC = AgesEqualToAndAboveRecAge[xxx]
  }
  return(LastAgeForLinearCC)
}

#' Determine age at full recruitment into the fishery (assuming knife-edge selection).
#'
#' This function determines the age at full recruitment, based on the age at
#' peak frequency of fish. If RecAssump=0, RecAge is set to the age at
#' peak frequency, but if RecAssump=1, RecAge is set to one year above the age at
#' peak frequency.
#' @keywords internal
#' @param params RecAssump - must be 0 (peak age) or 1 (peak age + 1), Ages, ObsAgeFreq (Numeric Vectors)
#' @return RecAge (age at full recruitment into the fishery)
#' @examples
#' set.seed(123)
#' MaxAge = 40
#' Ages = 1:MaxAge
#' NatMort <- exp(1.46 - (1.01 * (log(MaxAge)))) # i.e. Hoenig's (1983) eqn for fish
#' FMort = 0.1
#' ZMort = FMort + NatMort
#' SelA50 = 6
#' SelA95 = 8
#' SampleSize = 500 # required number of fish for age sample
#' SelAtAge = 1 / (1 + exp(-log(19) * (Ages - SelA50) / (SelA95 - SelA50)))
#' FAtAge = SelAtAge * FMort
#' ZAtAge = NatMort + FAtAge
#' N = numeric(MaxAge)
#' N[1] = 1
#' for (i in 2:MaxAge) {
#'   if (i < MaxAge) {
#'     N[i] = N[i-1] * exp(-ZAtAge[i-1])
#'   } else {
#'     N[i] = N[i-1] * exp(-ZAtAge[i-1] / (1 - exp(-ZAtAge[i])))
#'   }
#' }
#' CatchAtAge = N * (FAtAge / ZAtAge) * (1 - exp(-ZAtAge)) # catch at age
#' PropAtAge = CatchAtAge / sum(CatchAtAge)
#' CatchSample = rmultinom(n=1, size=SampleSize, prob=PropAtAge)
#' ObsAgeFreq = CatchSample[,1]
#' CalcRecruitmentAge(RecAssump=0, Ages, ObsAgeFreq)
CalcRecruitmentAge <- function(RecAssump, Ages, ObsAgeFreq) {

  # determine peak frequency
  PeakFreq = max(ObsAgeFreq)

  # determine age at which peak frequency occurs
  x = which(ObsAgeFreq == PeakFreq)
  x = max(x) # in case there are more than one ages with same 'peak' frequency

  #RecAssump: 0=peak age, 1=peak age + 1
  if (RecAssump == 0) {
    RecAge = Ages[x]
  }
  if (RecAssump == 1) {
    RecAge = Ages[x+1]
  }
  Result = RecAge

  return(Result)
}

#' Fit linear catch curve and get results.
#'
#' This function fits a linear catch curve to fish age composition data.
#' The age at full recruitment is determined using on the age at
#' peak frequency of fish. If RecAssump=0, RecAge is set to the age at
#' peak frequency, but if RecAssump=1, RecAge is set to one year above the age at
#' peak frequency. The oldest age that should be considered for analysis with
#' the linear catch curve analysis, is calcualted given a specified minimum frequency
#' (must be >0). The catch curve is fitted to the natural logarithms of the frequency
#' data using the lm function.
#' @param params RecAssump - must be 0 (peak age) or 1 (peak age + 1), MinFreq (Numbers), Ages, ObsAgeFreq (Numeric Vectors)
#' @return Est_yint (y intercept with 95% CLs), Est_ZMort (total mortality with 95% CLs)
#' (Numeric Vectors), RecAge, LastAgeForLinearCC (age range for analysis) (Numbers),
#' EstlnFreq (estimated natural logarithms of frequencies at age) (Numeric Vector).
#' @examples
#' set.seed(123)
#' MaxAge = 40
#' Ages = 1:MaxAge
#' NatMort <- exp(1.46 - (1.01 * (log(MaxAge)))) # i.e. Hoenig's (1983) eqn for fish
#' FMort = 0.1
#' ZMort = FMort + NatMort
#' SelA50 = 6
#' SelA95 = 8
#' SampleSize = 500 # required number of fish for age sample
#' SelAtAge = 1 / (1 + exp(-log(19) * (Ages - SelA50) / (SelA95 - SelA50)))
#' FAtAge = SelAtAge * FMort
#' ZAtAge = NatMort + FAtAge
#' N = numeric(MaxAge)
#' N[1] = 1
#' for (i in 2:MaxAge) {
#'   if (i < MaxAge) {
#'     N[i] = N[i-1] * exp(-ZAtAge[i-1])
#'   } else {
#'     N[i] = N[i-1] * exp(-ZAtAge[i-1] / (1 - exp(-ZAtAge[i])))
#'   }
#' }
#' CatchAtAge = N * (FAtAge / ZAtAge) * (1 - exp(-ZAtAge)) # catch at age
#' PropAtAge = CatchAtAge / sum(CatchAtAge)
#' CatchSample = rmultinom(n=1, size=SampleSize, prob=PropAtAge)
#' ObsAgeFreq = CatchSample[,1]
#' GetLinearCatchCurveResults(RecAssump=0, MinFreq=1, Ages, ObsAgeFreq)
#' @export
GetLinearCatchCurveResults <- function(RecAssump, MinFreq, Ages, ObsAgeFreq) {

  # get recruitment age, given recruitment assumption
  RecAge=CalcRecruitmentAge(RecAssump, Ages, ObsAgeFreq)

  # calculate last age to be considered in analysis, given minimum frequency
  LastAgeForLinearCC=CalcLastAgeForLinearCatchCurve(MinFreq, RecAge, Ages, ObsAgeFreq)

  # get frequencies in the required age range
  x=which(Ages==RecAge) # RecAge position
  xx=which(Ages==LastAgeForLinearCC) # position of last age
  ObsAgeFreqForLinearCC = ObsAgeFreq[x:xx]

  # get corresponding ages
  AgeClassesForLinearCC = RecAge:LastAgeForLinearCC

  # take the natural logarithms of the observed frequencies
  lnObsAgeFreqForLinearCC <- log(ObsAgeFreqForLinearCC)

  # run the linear catch curve model.
  mod1 <- lm(lnObsAgeFreqForLinearCC~AgeClassesForLinearCC)

  # get linear model parameters
  yint <- as.numeric(mod1$coefficients[1])
  ZMort <- - as.numeric(mod1$coefficients[2])

  # get standard errors for the parameters
  yint_se <- as.numeric(summary(mod1)$coef[3])
  ZMort_se <- as.numeric(summary(mod1)$coef[4])

  # get approx 95% confidence limits for parameters
  yint_low <- yint - (1.96 * yint_se)
  yint_up <- yint + (1.96 * yint_se)
  ZMort_low <- ZMort - (1.96 * ZMort_se)
  ZMort_up <- ZMort + (1.96 * ZMort_se)

  # note, confidence limits reported here are based on normal distribution,
  # rather than t distribution. May change in future.
  # compare confint vs confint.default
  # confint(mod1) # based on t distribution
  # confint.default(mod1) # based on normal distribution

  Est_yint = c(yint, yint_low, yint_up)
  Est_ZMort = c(ZMort, ZMort_low, ZMort_up)

  # get predicted frequencies (in log space) at age.
  EstlnFreq <- predict(lm(lnObsAgeFreqForLinearCC~AgeClassesForLinearCC))

  # get 95% confidence limits
  newdat <- data.frame(AgeClassesForLinearCC=AgeClassesForLinearCC)
  conf_int = predict(mod1, newdata = newdat, interval = 'confidence')

  EstlnFreq = conf_int[,1]
  EstlnFreq_Zlow = conf_int[,2]
  EstlnFreq_Zup = exp(conf_int[,3])
  EstFreq = exp(conf_int[,1])
  EstFreq_Zlow = exp(conf_int[,2])
  EstFreq_Zup = exp(conf_int[,3])

  results = list(yint = round(Est_yint,3),
                 yint_se = round(yint_se,3),
                 ZMort = round(Est_ZMort,3),
                 ZMort_se = round(ZMort_se,3),
                 RecAge = RecAge,
                 LastAgeForLinearCC = LastAgeForLinearCC,
                 AgeClassesForLinearCC = AgeClassesForLinearCC,
                 lnObsAgeFreqForLinearCC = lnObsAgeFreqForLinearCC,
                 EstlnFreq = EstlnFreq,
                 EstlnFreq_Zlow = EstlnFreq_Zlow,
                 EstlnFreq_Zup = EstlnFreq_Zup,
                 Ages=Ages,
                 EstFreq = EstFreq,
                 EstFreq_Zlow=EstFreq_Zlow,
                 EstFreq_Zup = EstFreq_Zup)

  return(results)
}

#' Estimate total mortality applying the Chapman & Robson (1960) mortality estimator
#'
#' This function provides an estimate of total mortality, applying the Chapman & Robson (1960) mortality
#' estimator, with associated 95 percent confidence limits. Additional variables are returned for plotting.
#' @param params RecAssump - must be 0 (peak age) or 1 (peak age + 1), MinAge, MaxAge (integers), ObsAgeFreq for Ages=1:MaxAge (Numeric Vector)
#' @return total mortality, Z (ZMort), standard error and approximate lower and upper 95 percent confidence
#' limits for Z, (ZMort_se, ZMort_low, ZMort_up), median estimate and lower and upper 95 percent confidence
#' limits for Z from resampling (ZMort_resamp, ZMort_lowresamp, ZMort_upresamp), variable for plotting,
#' including expected frequencies at age with associated confidence limits, calculated from resampling
#'(EstFreq, EstFreq_Zup, EstFreq_Zlow), recruitment age (RecAge),  maximum age in sample (MaxAgeInSample),
#' relative ages used in analysis (CRAges), sample size for relative ages (n) and additional variable
#' calculated in analysis (CR_T)
#' @examples
#' set.seed(123)
#' MinAge = 1
#' MaxAge = 40
#' Ages = MinAge:MaxAge
#' NatMort <- exp(1.46 - (1.01 * (log(MaxAge)))) # i.e. Hoenig's (1983) eqn for fish
#' FMort = 0.1
#' ZMort = FMort + NatMort
#' SelA50 = 6
#' SelA95 = 8
#' SampleSize = 500 # required number of fish for age sample
#' SelAtAge = 1 / (1 + exp(-log(19) * (Ages - SelA50) / (SelA95 - SelA50)))
#' FAtAge = SelAtAge * FMort
#' ZAtAge = NatMort + FAtAge
#' N = numeric(MaxAge)
#' N[1] = 1
#' for (i in 2:MaxAge) {
#'   if (i < MaxAge) {
#'     N[i] = N[i-1] * exp(-ZAtAge[i-1])
#'   } else {
#'     N[i] = N[i-1] * exp(-ZAtAge[i-1] / (1 - exp(-ZAtAge[i])))
#'   }
#' }
#' CatchAtAge = N * (FAtAge / ZAtAge) * (1 - exp(-ZAtAge)) # catch at age
#' PropAtAge = CatchAtAge / sum(CatchAtAge)
#' CatchSample = rmultinom(n=1, size=SampleSize, prob=PropAtAge)
#' ObsAgeFreq = CatchSample[,1]
#' GetChapmanRobsonMortalityResults(RecAssump=0, MinAge, MaxAge, ObsAgeFreq)
#' @export
GetChapmanRobsonMortalityResults <- function(RecAssump, MinAge, MaxAge, ObsAgeFreq)
{
  Ages = MinAge:MaxAge
  MaxAgeInSample = max(Ages)

  FishAges <- rep(Ages, ObsAgeFreq)
  RecAge = CalcRecruitmentAge(RecAssump, Ages, ObsAgeFreq)
  x = which(Ages == RecAge)
  Ages_minus_RecAge <- FishAges - RecAge
  CRAges <- Ages_minus_RecAge[which(Ages_minus_RecAge >= 0)]
  n <- length(CRAges)
  AveAge <- mean(CRAges)
  CR_T <- n * AveAge
  ZMort <- log(1 + AveAge - (n^-1)) - log(AveAge) -
    (((n - 1) * (n - 2)) * (n * (CR_T + 1) * (n + CR_T - 1))^-1)
  # analytical solution
  Var_ZMort <- ((1 - exp(-ZMort))^2)/(n * exp(-ZMort))
  ZMort_se <- sqrt(Var_ZMort)
  ZMort_up <- ZMort + (1.96 * ZMort_se)
  ZMort_low <- ZMort - (1.96 * ZMort_se)
  # resampling
  set.seed(123)
  RandZMort = rep(0,500)
  ArrSize = (MaxAgeInSample-RecAge) + 1
  RandEstFreq = data.frame(matrix(nrow=200, ncol=ArrSize))
  colnames(RandEstFreq) <- seq(0,MaxAgeInSample-RecAge,1)
  for (i in 1:500) {
    # resampling only ages at or above recruitment age (i.e. fixed RecAge)
    RandCRAges <- sample(CRAges, replace = T)
    RandCRAgeFreq = as.vector(table(factor(RandCRAges, levels=seq(0,MaxAgeInSample-RecAge,1))))
    RandAveAge <- mean(RandCRAges)
    RandCR_T <- n * RandAveAge
    RandZMort[i] <- log(1 + RandAveAge - (n^-1)) - log(RandAveAge) -
      (((n - 1) * (n - 2)) * (n * (RandCR_T + 1) * (n + RandCR_T - 1))^-1)
    RandEstFreq[i,] = RandCRAgeFreq[1] * exp(-RandZMort[i] * seq(0, ArrSize - 1,1))
  }

  EstFreq = apply(RandEstFreq,2,median)
  EstFreq_Zup = apply(RandEstFreq,2,quantile,probs=0.025)
  EstFreq_Zlow = apply(RandEstFreq,2,quantile,probs=0.975)
  ZMort_resamp = median(RandZMort)
  ZMort_lowresamp = quantile(RandZMort, 0.025)
  ZMort_upresamp = quantile(RandZMort, 0.975)

  results = list(ZMort = round(ZMort, 3),
                 ZMort_se = round(ZMort_se, 3),
                 ZMort_low = round(ZMort_low, 3),
                 ZMort_up = round(ZMort_up, 3),
                 ZMort_resamp = round(ZMort_resamp, 3),
                 ZMort_lowresamp = round(ZMort_lowresamp, 3),
                 ZMort_upresamp = round(ZMort_upresamp, 3),
                 EstFreq = EstFreq,
                 Ages = Ages,
                 EstFreq = EstFreq,
                 EstFreq_Zup = EstFreq_Zup,
                 EstFreq_Zlow = EstFreq_Zlow,
                 RecAge = RecAge,
                 MaxAgeInSample = MaxAgeInSample,
                 CRAges = CRAges,
                 n = n,
                 CR_T = CR_T)

  return(results)
}



#' Plot age based catch curve results in normal space
#'
#' This function produces a plots of outputs of age-based catch curve analyses in normal space
#' @param params RecAssump - must be 0 (peak age) or 1 (peak age + 1) or NA, MinFreq (integer or NA),
#' MaxAge (integers), ObsAgeFreq for Ages=MinAge:MaxAge (Numeric Vector), CatchCurveModel
#' (1=C&R, 2=Linear, 3=logist. sel.) various plotting options (use NA for defaults), including
#' specifying plot label, MainLabel, x and y axis labels xaxis_lab, yaxis_lab, x and y axis
#' maximum limits and intervals, xmax, xint, ymax, yint, option to plot 95 percent confidence limits, PlotCLs
#' @return plot of expected vs observed proportions at age
#' @examples
#' set.seed(123)
#' MinAge = 1
#' MaxAge = 40
#' Ages = MinAge:MaxAge
#' NatMort <- exp(1.46 - (1.01 * (log(MaxAge)))) # i.e. Hoenig's (1983) eqn for fish
#' FMort = 0.1
#' ZMort = FMort + NatMort
#' SelA50 = 6
#' SelA95 = 8
#' SampleSize = 500 # required number of fish for age sample
#' SelAtAge = 1 / (1 + exp(-log(19) * (Ages - SelA50) / (SelA95 - SelA50)))
#' FAtAge = SelAtAge * FMort
#' ZAtAge = NatMort + FAtAge
#' N = numeric(MaxAge)
#' N[1] = 1
#' for (i in 2:MaxAge) {
#'   if (i < MaxAge) {
#'     N[i] = N[i-1] * exp(-ZAtAge[i-1])
#'   } else {
#'     N[i] = N[i-1] * exp(-ZAtAge[i-1] / (1 - exp(-ZAtAge[i])))
#'   }
#' }
#' CatchAtAge = N * (FAtAge / ZAtAge) * (1 - exp(-ZAtAge)) # catch at age
#' PropAtAge = CatchAtAge / sum(CatchAtAge)
#' CatchSample = rmultinom(n=1, size=SampleSize, prob=PropAtAge)
#' ObsAgeFreq = CatchSample[,1]
#' # CatchCurveModel = 1 # Chapman Robson
#' # NatMort = NA # Chapman Robson
#' # RecAssump = 1 # Chapman Robson
#' # MinFreq = NA # Chapman Robson
#' # CatchCurveModel = 2 # Linear
#' # NatMort = NA # Linear
#' # RecAssump = 1 # Linear
#' # MinFreq = 1 # Linear
#' CatchCurveModel = 3 # Logistic selectivity
#' RecAssump = NA # Logistic selectivity
#' NatMort = 0.104 # Logistic selectivity
#' MinFreq = NA # # Logistic selectivity
#' ln_params = log(c(FMort, SelA50, SelA95))
#' PlotAgeBasedCatchCurveResults_NormalSpace(RecAssump, MinFreq, MinAge, MaxAge, NatMort,
#'                                           ObsAgeFreq, CatchCurveModel, MainLabel=NA,
#'                                           xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA,
#'                                          ymax=NA, yint=NA, PlotCLs=T)
#' @export
PlotAgeBasedCatchCurveResults_NormalSpace <- function(RecAssump, MinFreq, MinAge, MaxAge, NatMort,
                                                      ObsAgeFreq, CatchCurveModel, MainLabel,
                                                      xaxis_lab, yaxis_lab, xmax, xint,
                                                      ymax, yint, PlotCLs) {

  if (is.na(xaxis_lab)) xaxis_lab = "Age class"
  if (is.na(yaxis_lab)) yaxis_lab = "Frequency"
  ylims = Get_yaxis_scale(ObsAgeFreq)
  if (is.na(ymax)) ymax = ylims$ymax
  if (is.na(yint)) yint = ylims$yint

  # Ages
  Ages = MinAge:MaxAge
  xlims = Get_xaxis_scale(Ages)
  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint

  # Chapman-Robson
  if (CatchCurveModel == 1) {
    Res = GetChapmanRobsonMortalityResults(RecAssump, MinAge, MaxAge, ObsAgeFreq)
    if (is.na(MainLabel)) MainLabel = "Chapman & Robson"
  }
  # Linear
  if (CatchCurveModel == 2) {
    Res = GetLinearCatchCurveResults(RecAssump, MinFreq, Ages, ObsAgeFreq)
    if (is.na(MainLabel)) MainLabel = "Linear"
  }
  # logistic age-based selectivity
  if (CatchCurveModel == 3) {
    Res = GetLogisticCatchCurveResults(ln_params, NatMort, Ages, ObsAgeFreq)
    if (is.na(MainLabel)) MainLabel = "Logistic selectivity"
  }


  # Chapman-Robson
  if (CatchCurveModel == 1) {
    Z_value = round(Res$ZMort,digits=3)
    x=which(Ages==Res$RecAge) # RecAge position
    xx=length(which(Res$EstFreq_Zup>0)) # position of last age
    xxx=length(which(Res$EstFreq_Zlow>0)) # position of last age
    xxxx=length(which(Res$EstFreq>0)) # position of last age
  }
  # Linear
  if (CatchCurveModel == 2) {
    Z_value = round(Res$ZMort[1],digits=3)
    x=which(Ages==Res$RecAge) # RecAge position
    xx=which(Ages==Res$LastAgeForLinearCC) # position of last age
    xxx=which(Ages==Res$LastAgeForLinearCC) #  # position of last age
    xxxx=which(Ages==Res$LastAgeForLinearCC) #  # position of last age
  }
  if (CatchCurveModel == 3) {
    Z_value = round(Res$EstZMort,digits=3)
    x=which(Ages==min(Ages)) # RecAge position
    xx=which(Ages==max(Ages)) # position of last age
    xxx=which(Ages==max(Ages)) #  # position of last age
    xxxx=which(Ages==max(Ages)) #  # position of last age
  }

  j = seq(x,xx,1) # up
  jj = seq(x,xxx,1) # low
  jjj = seq(x,xxxx,1) # est

  # plot catch curve over age frequency data, in normal space
  plot(Ages, ObsAgeFreq, "p", main=MainLabel, cex.main=1.0, pch=16, cex=0.8, xaxt = "n", yaxt = "n",
       xlab=xaxis_lab,ylab=yaxis_lab, frame=F, xlim=c(0,xmax), ylim=c(0,ymax)) # observed data (normal space)
  axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
  axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
  axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0,
       cex.axis = 0.8, las = 1)
  axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0,
       cex.axis = 0.8, las = 1)
  if (PlotCLs == TRUE) {
    x = c(Ages[j],rev(Ages[j])) # using shading for 95% CLs
    y = c(Res$EstFreq_Zlow[1:length(jj)], rev(Res$EstFreq_Zup[1:length(j)]))
    polygon(x,y,col="pink",border=NA)
    # lines(Ages[j], Res$EstFreq_Zup[1:length(j)] , col="black", lty="dotted")
    # lines(Ages[jj], Res$EstFreq_Zlow[1:length(jj)] , col="black", lty="dotted")
  }
  # lines(Ages[jjj], Res$EstFreq[1:length(jjj)], col="black")
  points(Ages, ObsAgeFreq, pch=16, cex=0.8)
  points(Ages[j], Res$EstFreq[1:length(j)], col="red", pch=1, cex=0.8)

  legend("topright", legend=bquote(paste("Z = ", .(Z_value), " ",y^-1)), y.intersp = 1.5, inset=c(0.13,0),
         lty=1, cex = 1, bty="n",seg.len = 0)
  legend("topleft", legend=c("Observed","Estimated"), y.intersp = 1.0, inset=c(0.13,0),
         lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16,1), col=c("black","red"))
}

#' Plot age based catch curve results in log space
#'
#' This function produces a plots of outputs of age-based catch curve analyses in log space
#' @param RecAssump - must be 0 (peak age) or 1 (peak age + 1) or NA, MinFreq (integer or NA),
#' MaxAge (integers), ObsAgeFreq for Ages=MinAge:MaxAge (Numeric Vector), CatchCurveModel
#' (1=C&R, 2=Linear, 3=logist. sel.) various plotting options (use NA for defaults), including
#' specifying plot label, MainLabel, x and y axis labels xaxis_lab, yaxis_lab, x and y axis
#' maximum limits and intervals, xmax, xint, ymax, yint, option to plot 95 percent confidence limits, PlotCLs
#' @return plot of expected vs observed proportions at age
#' @examples
#' set.seed(123)
#' MinAge = 1
#' MaxAge = 40
#' Ages = MinAge:MaxAge
#' NatMort <- exp(1.46 - (1.01 * (log(MaxAge)))) # i.e. Hoenig's (1983) eqn for fish
#' FMort = 0.1
#' ZMort = FMort + NatMort
#' SelA50 = 6
#' SelA95 = 8
#' SampleSize = 500 # required number of fish for age sample
#' SelAtAge = 1 / (1 + exp(-log(19) * (Ages - SelA50) / (SelA95 - SelA50)))
#' FAtAge = SelAtAge * FMort
#' ZAtAge = NatMort + FAtAge
#' N = numeric(MaxAge)
#' N[1] = 1
#' for (i in 2:MaxAge) {
#'   if (i < MaxAge) {
#'     N[i] = N[i-1] * exp(-ZAtAge[i-1])
#'   } else {
#'     N[i] = N[i-1] * exp(-ZAtAge[i-1] / (1 - exp(-ZAtAge[i])))
#'   }
#' }
#' CatchAtAge = N * (FAtAge / ZAtAge) * (1 - exp(-ZAtAge)) # catch at age
#' PropAtAge = CatchAtAge / sum(CatchAtAge)
#' CatchSample = rmultinom(n=1, size=SampleSize, prob=PropAtAge)
#' ObsAgeFreq = CatchSample[,1]
#' # CatchCurveModel = 1 # Chapman Robson
#' # NatMort = NA # Chapman Robson
#' # RecAssump = 1 # Chapman Robson
#' # MinFreq = NA # Chapman Robson
#' # CatchCurveModel = 2 # Linear
#' # NatMort = NA # Linear
#' # RecAssump = 1 # Linear
#' # MinFreq = 1 # Linear
#' CatchCurveModel = 3 # Logistic selectivity
#' RecAssump = NA # Logistic selectivity
#' NatMort = 0.104 # Logistic selectivity
#' MinFreq = NA # # Logistic selectivity
#' ln_params = log(c(FMort, SelA50, SelA95))
#' PlotAgeBasedCatchCurveResults_LogSpace(RecAssump, MinFreq, MinAge, MaxAge, NatMort,
#'                                        ObsAgeFreq, CatchCurveModel, MainLabel=NA,
#'                                        xaxis_lab=NA, yaxis_lab=NA, ymin=NA, xmax=NA, xint=NA,
#'                                        ymax=NA, yint=NA, PlotCLs=T)
#' @export
PlotAgeBasedCatchCurveResults_LogSpace <- function(RecAssump, MinFreq, MinAge, MaxAge, NatMort,
                                                   ObsAgeFreq, CatchCurveModel, MainLabel,
                                                   xaxis_lab, yaxis_lab, ymin, xmax, xint,
                                                   ymax, yint, PlotCLs) {

  if (is.na(xaxis_lab)) xaxis_lab = "Age class"
  if (is.na(yaxis_lab)) yaxis_lab = "ln Frequency"
  ylims = Get_yaxis_scale(log(ObsAgeFreq))
  if (is.na(ymax)) ymax = ylims$ymax
  if (is.na(yint)) yint = ylims$yint
  if (is.na(ymin)) ymin = -2

  # Ages
  Ages = MinAge:MaxAge
  xlims = Get_xaxis_scale(Ages)
  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint

  # Chapman-Robson
  if (CatchCurveModel == 1) {
    Res = GetChapmanRobsonMortalityResults(RecAssump, MinAge, MaxAge, ObsAgeFreq)
    if (is.na(MainLabel)) MainLabel = "Chapman & Robson"
  }
  # Linear
  if (CatchCurveModel == 2) {
    Res = GetLinearCatchCurveResults(RecAssump, MinFreq, Ages, ObsAgeFreq)
    if (is.na(MainLabel)) MainLabel = "Linear"
  }
  # logistic age-based selectivity
  if (CatchCurveModel == 3) {
    Res = GetLogisticCatchCurveResults(ln_params, NatMort, Ages, ObsAgeFreq)
    if (is.na(MainLabel)) MainLabel = "Logistic selectivity"
  }

  # Chapman-Robson
  if (CatchCurveModel == 1) {
    Z_value = round(Res$ZMort,digits=3)
    x=which(Ages==Res$RecAge) # RecAge position
    xx=max(which(log(Res$EstFreq_Zup)>0)) # position of last age
    xxx=max(which(log(Res$EstFreq_Zlow)>0)) # position of last age
    xxxx=max(which(log(Res$EstFreq)>0)) # position of last age
    # ages
    j = seq(x,xx+x-1,1) # up
    jj = seq(x,xxx+x-1,1) # low
    jjj = seq(x,xxxx+x-1,1) # est
    # estimates above rec age
    k = seq(1,length(j),1) # up
    kk = seq(1,length(jj),1) # low
    kkk = seq(1,length(jjj),1) # est
  }
  # Linear
  if (CatchCurveModel == 2) {
    Z_value = round(Res$ZMort[1],digits=3)
    x=which(Ages==Res$RecAge) # RecAge position
    xx=which(Ages==Res$LastAgeForLinearCC) # position of last age
    xxx=which(Ages==Res$LastAgeForLinearCC) #  # position of last age
    xxxx=which(Ages==Res$LastAgeForLinearCC) #  # position of last age
    j = seq(x,xx,1) # up
    jj = seq(x,xxx,1) # low
    jjj = seq(x,xxxx,1) # est
  }
  # logistic selectivity
  if (CatchCurveModel == 3) {
    Z_value = round(Res$EstZMort,digits=3)
    x=which(Ages==min(which(log(Res$EstFreq_Zlow)>0))) # RecAge position
    xx=length(which(log(Res$EstFreq_Zlow) > -1)) # position of last age
    j = seq(x,xx,1) # up
  }

  plot(Ages, log(ObsAgeFreq), "p", main=MainLabel, pch=16, cex=0.8, cex.main=1.0, xaxt = "n", yaxt = "n",
       xlab=xaxis_lab,ylab=yaxis_lab, frame=F, xlim=c(0,xmax), ylim=c(ymin,ymax)) # observed data (normal space)
  axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
  axis(2, at = seq(ymin, ymax, yint), line = 0.2, labels = F)
  axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0,
       cex.axis = 0.8, las = 1)
  axis(2, at = seq(ymin, ymax, yint), lwd = 0, labels = T, line = 0,
       cex.axis = 0.8, las = 1)
  # Chap-Rob
  if (CatchCurveModel == 1) {
    if (PlotCLs == TRUE) {
      x = c(Ages[j],rev(Ages[jj])) # using shading for 95% CLs
      y = c(log(Res$EstFreq_Zup[k]), rev(log(Res$EstFreq_Zlow[kk])))
      polygon(x,y,col="pink",border=NA)
      # lines(Ages[j], log(Res$EstFreq_Zup[k]), col="black", lty="dotted")
      # lines(Ages[jj], log(Res$EstFreq_Zlow[kk]), col="black", lty="dotted")
    }
    # lines(Ages[jjj], log(Res$EstFreq[kkk]), col="red")
    points(Ages[jjj], log(Res$EstFreq[kkk]), pch=1, col="red", cex=0.6)
  }
  if (CatchCurveModel == 2) {
    if (PlotCLs == TRUE) {
      x = c(Ages[j],rev(Ages[jj])) # using shading for 95% CLs
      y = c(log(Res$EstFreq_Zup[1:length(j)]), rev(log(Res$EstFreq_Zlow[1:length(jj)])))
      polygon(x,y,col= "pink",border=NA)
      # lines(Ages[j], log(Res$EstFreq_Zup[1:length(j)]), col="black", lty="dotted")
      # lines(Ages[jj], log(Res$EstFreq_Zlow[1:length(jj)]), col="black", lty="dotted")
    }
    # lines(Ages[jjj], log(Res$EstFreq[1:length(jjj)]), col="red")
    points(Ages[jjj], log(Res$EstFreq[1:length(jjj)]), pch=1, col="red", cex=0.6)
  }
  if (CatchCurveModel == 3) {
    if (PlotCLs == TRUE) {
      # arrows(Ages[j], log(Res$EstFreq_Zup[j]), Ages[j], log(Res$EstFreq_Zlow[j]), length=0.02, angle=90, code=3,
      #        col="dark grey")
      x = c(Ages[j],rev(Ages[j])) # using shading for 95% CLs
      y = c(log(Res$EstFreq_Zup[j]), rev(log(Res$EstFreq_Zlow[j])))
      polygon(x,y,col="pink",border=NA)
      # lines(Ages[j], log(Res$EstFreq_Zup[j]), col="pink", lty="dotted")
      # lines(Ages[j], log(Res$EstFreq_Zlow[j]), col="pink", lty="dotted")
    }
    #lines(Ages[j], log(Res$EstFreq[j]), col="red")
    points(Ages[j], log(Res$EstFreq[j]), pch=1, col="red", cex=0.6)
  }

  points(Ages, log(ObsAgeFreq), pch=16, cex=0.8)
  legend("topright", legend=bquote(paste("Z = ", .(Z_value), " ",y^-1)), y.intersp = 1.5, inset=c(0.13,0),
         lty=1, cex = 1, bty="n",seg.len = 0)
  legend("topleft", legend=c("Observed","Estimated"), y.intersp = 1.0, inset=c(0.13,0),
         lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16,1), col=c("black","red"))
}

#' Return negative log-likelihood for a catch curve with age-based, logistic selectivity
#'
#' This function returns the negative log-likelihood, given age frequency data, a value total mortality,
#' and age-based selectivity parameters, for a catch curve with age-based, logistic selectivity.
#' Function requires an estimate of natural mortality (NatMort), a value for maximum age (MaxAge)
#' and age frequency data (stored in memory in R).
#' @keywords internal
#' @param params FMort, SelA50, SelA95 (Numbers)
#' @return NLL (Number)
#' @examples
#' set.seed(123)
#' MaxAge = 40
#' Ages = 1:MaxAge
#' NatMort <- exp(1.46 - (1.01 * (log(MaxAge)))) # i.e. Hoenig's (1983) eqn for fish
#' FMort = 0.1
#' ZMort = FMort + NatMort
#' SelA50 = 6
#' SelA95 = 8
#' SampleSize = 500 # required number of fish for age sample
#' SelAtAge = 1 / (1 + exp(-log(19) * (Ages - SelA50) / (SelA95 - SelA50)))
#' FAtAge = SelAtAge * FMort
#' ZAtAge = NatMort + FAtAge
#' N = numeric(MaxAge)
#' N[1] = 1
#' for (i in 2:MaxAge) {
#'   if (i < MaxAge) {
#'     N[i] = N[i-1] * exp(-ZAtAge[i-1])
#'   } else {
#'     N[i] = N[i-1] * exp(-ZAtAge[i-1] / (1 - exp(-ZAtAge[i])))
#'   }
#' }
#' CatchAtAge = N * (FAtAge / ZAtAge) * (1 - exp(-ZAtAge)) # catch at age
#' PropAtAge = CatchAtAge / sum(CatchAtAge)
#' CatchSample = rmultinom(n=1, size=SampleSize, prob=PropAtAge)
#' ObsAgeFreq = CatchSample[,1]
#' Calculate_NLL_LogisticCatchCurve(ln_params)
Calculate_NLL_LogisticCatchCurve <- function(ln_params) {

  FMort = exp(ln_params[1])
  SelA50 = exp(ln_params[2])
  SelA95 = exp(ln_params[3])

  # selectivity
  SelAtAge = rep(0,length(Ages))
  SelAtAge = 1 / (1 + exp(-log(19) * (Ages - SelA50) / (SelA95 - SelA50)))

  # fishing mortality at age
  FAtAge = SelAtAge * FMort

  # total mortality at age
  ZAtAge = NatMort + FAtAge

  # popn. numbers at age
  N = numeric(length(Ages))
  N[1] = 1
  k=1
  MinAge = min(Ages)
  MaxAge = max(Ages)
  for (i in seq(MinAge+1,MaxAge,1)) {
    k=k+1
    if (i < MaxAge) {
      N[k] = N[k-1] * exp(-ZAtAge[k])
    } else {
      N[k] = N[k-1] * exp(-ZAtAge[k-1]) / (1 - exp(-ZAtAge[k]))
    }
  }

  # catch at age
  CatchAtAge = N * (FAtAge / ZAtAge) * (1 - exp(-ZAtAge))

  # calculate expected catch proportion at age
  ExpPropAtAge = CatchAtAge / sum(CatchAtAge)

  # calculate F penalty
  F_Pen = 0
  if (FMort > 2.0) {
    F_Pen = 1000 * (FMort - 2.0)^2
  }

  # calculate multinonmial negative log-likelihood
  NLL = -sum((ObsAgeFreq * log(ExpPropAtAge + 1E-4))) + F_Pen

  return(NLL)

}

#' Get statistical outputs from a fitted catch curve with age-based, logistic selectivity
#'
#' This function fits a catch curve with an asymptotic, age-based logistic selectivity curve,
#' to a sample of fish age frequency data, by minimising the negative log-likelihood associated
#' with the parameters and data, using nlminb. It provides various statistical outputs in include
#' convergence statistics, parameter estimates and associated 95 percent confidence limits and associated
#' variance-covariance matrix, calculated using the MASS package.
#' Function requires an estimate of natural mortality (NatMort), a value for maximum age (MaxAge)
#' and age frequency data (stored in memory in R).
#' (Numbers and a Numeric Vectors)
#' @param params Fishing mortality (FMort), ages at which 50 percent (SelA50) and
#'     95 percent of fish are selected by the gear (SelA95), Numbers.
#' @return negative log-likelihood (nll), nlminb convergence diagnostic (convergence)
#' sample size (SampleSize), growth parameter estimates with lower and upper 95 percent
#' confidence limits (ParamEst), point estimates for growth parameters (params)
#' and variance-covariance matrix (vcov.Params)
#' @examples
#' set.seed(123)
#' MaxAge = 40
#' Ages = 1:MaxAge
#' NatMort <- exp(1.46 - (1.01 * (log(MaxAge)))) # i.e. Hoenig's (1983) eqn for fish
#' FMort = 0.1
#' ZMort = FMort + NatMort
#' SelA50 = 6
#' SelA95 = 8
#' SampleSize = 500 # required number of fish for age sample
#' SelAtAge = 1 / (1 + exp(-log(19) * (Ages - SelA50) / (SelA95 - SelA50)))
#' FAtAge = SelAtAge * FMort
#' ZAtAge = NatMort + FAtAge
#' N = numeric(MaxAge)
#' N[1] = 1
#' for (i in 2:MaxAge) {
#'   if (i < MaxAge) {
#'     N[i] = N[i-1] * exp(-ZAtAge[i-1])
#'   } else {
#'     N[i] = N[i-1] * exp(-ZAtAge[i-1] / (1 - exp(-ZAtAge[i])))
#'   }
#' }
#' CatchAtAge = N * (FAtAge / ZAtAge) * (1 - exp(-ZAtAge)) # catch at age
#' PropAtAge = CatchAtAge / sum(CatchAtAge)
#' CatchSample = rmultinom(n=1, size=SampleSize, prob=PropAtAge)
#' ObsAgeFreq = CatchSample[,1]
#' NatMort = 0.104
#' ln_params = log(c(FMort, SelA50, SelA95))
#' GetLogisticCatchCurveResults(ln_params, NatMort, Ages, ObsAgeFreq)
#' @export
GetLogisticCatchCurveResults <- function (ln_params, NatMort, Ages, ObsAgeFreq)
{
  nlmb = nlminb(ln_params, Calculate_NLL_LogisticCatchCurve)
  nlmb$objective
  nlmb$convergence
  nlmb$par
  (hess.out = optimHess(nlmb$par, Calculate_NLL_LogisticCatchCurve))
  (vcov.Params = solve(hess.out))
  (ses = sqrt(diag(vcov.Params)))
  EstFMort = c(exp(nlmb$par[1]), exp(nlmb$par[1] + c(-1.96,
                                                     1.96) * ses[1]))
  EstSelA50 = c(exp(nlmb$par[2]), exp(nlmb$par[2] + c(-1.96,
                                                      1.96) * ses[2]))
  EstSelA95 = c(exp(nlmb$par[3]), exp(nlmb$par[3] + c(-1.96,
                                                      1.96) * ses[3]))
  SelAtAge = rep(0, length(Ages))
  SelAtAge = 1/(1 + exp(-log(19) * (Ages - EstSelA50[1])/(EstSelA95[1] -
                                                            EstSelA50[1])))
  FAtAge = SelAtAge * EstFMort[1]
  ZAtAge = NatMort + FAtAge
  ParamEst = t(data.frame(FMort = round(EstFMort, 3), SelA50 = round(EstSelA50,3),
                          SelA95 = round(EstSelA95, 3)))
  colnames(ParamEst) = c("Estimate", "lw_95%CL",
                         "up_95%CL")
  SampleSize = sum(ObsAgeFreq)
  nll = nlmb$objective
  convergence = nlmb$convergence
  set.seed(123)
  params = nlmb$par
  vcov.params = vcov.Params
  sims = data.frame(MASS::mvrnorm(n = 500, params, vcov.params))
  EstFreq.sim = data.frame(matrix(nrow = 500, ncol = length(Ages)))
  colnames(EstFreq.sim) <- Ages
  SelA50.sim = rep(0,500)
  SelA95.sim = rep(0,500)
  FMort.sim = rep(0,500)
  EstZMort.sim = rep(0,500)
  for (j in 1:500) {
    FMort.sim[j] = exp(sims[j, 1])
    SelA50.sim[j] = exp(sims[j, 2])
    SelA95.sim[j] = exp(sims[j, 3])
    SelAtAge.sim = rep(0, length(Ages))
    SelAtAge.sim = 1/(1 + exp(-log(19) * (Ages - SelA50.sim[j])/(SelA95.sim[j] -
                                                                SelA50.sim[j])))
    FAtAge.sim = SelAtAge.sim * FMort.sim[j]
    ZAtAge.sim = NatMort + FAtAge.sim
    N.sim = numeric(length(Ages))
    MinAge = min(Ages)
    MaxAge = max(Ages)
    N.sim[1] = 1
    k = 1
    for (i in seq(MinAge + 1, MaxAge, 1)) {
      k = k + 1
      if (i < MaxAge) {
        N.sim[k] = N.sim[k - 1] * exp(-ZAtAge.sim[k - 1])
      }
      else {
        N.sim[k] = N.sim[k - 1] * exp(-ZAtAge.sim[k - 1])/(1 - exp(-ZAtAge.sim[k]))
      }
    }
    CatchAtAge.sim = N.sim * (FAtAge.sim/ZAtAge.sim) * (1 - exp(-ZAtAge.sim))
    ExpPropAtAge.sim = CatchAtAge.sim/sum(CatchAtAge.sim)
    EstFreq.sim[j, ] = SampleSize * ExpPropAtAge.sim
    EstZMort.sim[j] = FMort.sim[j] + NatMort
  }
  EstFreq = apply(EstFreq.sim, 2, median)
  EstFreq_Zlow = apply(EstFreq.sim, 2, quantile, probs = 0.025)
  EstFreq_Zup = apply(EstFreq.sim, 2, quantile, probs = 0.975)

  results = list(nll = nll,
                 convergence = convergence,
                 SampleSize = SampleSize,
                 ParamEst = ParamEst,
                 Estlnparams = nlmb$par,
                 lnparams.sims = sims,
                 lnEstFMort_se = ses[1],
                 lnEstSelA50_se = ses[2],
                 lnEstFSelA95_se = ses[3],
                 vcov.Params = vcov.Params,
                 SelAtAge = SelAtAge,
                 FAtAge = FAtAge,
                 EstFreq.sim = EstFreq.sim,
                 EstFreq = EstFreq,
                 EstFreq_Zlow = EstFreq_Zlow,
                 EstFreq_Zup = EstFreq_Zup,
                 SelA50.sim = SelA50.sim,
                 SelA95.sim = SelA95.sim,
                 FMort.sim = FMort.sim,
                 EstZMort.sim = EstZMort.sim)
  return(results)
}

#*********************
# Per recruit analyses
#*********************

#' Get outputs from per recruit analysis for a specified value of fishing mortality
#'
#' This function provides outputs associated with per recruit analysis, and an
#' extended form of this analysis with a Beverton-Holt stock recruitment relationship to account
#' for potential impacts of fishing on recruitment
#'
#' @param params maximum age considered in model (MaxModelAge), model timestep, e.g. 1 = annual, 1/12 = monthly (TimeStep),
#' von Bertalanffy growth parameters (Linf, vbK, tzero), Female and male estimated lengths-at-age (set to NA if
#' using von Bertalanffy growth equation, or set the von Bertalanffy growth parameters to NA and input as data frame
#' (columns for sex) from any growth model (EstLenAtAge), weight-length parameters, using a power or log-log relationship
#' (lenwt_a, ln_lenwt_a, lenwt_b), weight-length relationship type (WLrel_Type), Female and male estimated
#' weights-at-ages  (set to NA if weight-length relationship parameters, or set the weight-length parameters
#' to NA and input these as data frame (columns for sex), (EstWtAtAge), reproductive pattern, 1=gonochoristic,
#' 2=protogynous, 3=protandrous (ReprodPattern), initial sex ratio of juvenile recruits (InitRatioFem),
#' Maximum proportion of final sex for hermaphroditic species (set to NA for gonochoristic species) (FinalSex_Pmax),
#' logistic sex change parameters for hermaphroditic species (set to NA for gonochoristic species) (FinalSex_A50, FinalSex_A95),
#' logistic length at maturity parameters (mat_L50, mat_L95), (set to NA if using age at maturity parameters, or set the
#' age at maturity parameters to NA and input these as data fra,e (columns for sex) (EstMatAtAge), logistic gear selectivity
#' parameters (sel_A50, sel_A95),logistic retention curve parameters (ret_Pmax, ret_A50, ret_A95), discard mortality rate
#' (DiscMort), steepness parameter of the stock-recruitment relationship (Steepness), stock-recruitment relationship type,
#' 1 = Beverton-Holt, 2=Ricker (SRrel_Type), natural mortality (NatMort), fishing mortality (FMort). Number.
#' @return Ages considered in analysis (Ages), yield per recruit for combined sexes (YPR),
#' female, male and combined sex spawning potential ratio (Fem_SPR, Mal_SPR, CombSex_SPR), unfished female
#' survival at age (UnfishFemSurvAtAge, UnfishMalSurvAtAge), unfished female and male biomass mature biomass at age
#' (UnfishFemBiomAtAge, UnfishMalBiomAtAge), fished female survival at age (FishFemSurvAtAge, FishMalSurvAtAge),
#' fished female and male biomass mature biomass at age (FishFemBiomAtAge, FishMalBiomAtAge), female and male
#' catch at age in numbers (FemCatchAtAgeNum, MalCatchAtAgeNum), female and male catch at age in biomass
#' (FemCatchAtAge, MalCatchAtAge), derived Beverton-Holt stock recruitment parameters and associated equilibrium
#' recruitment derived from those parameters (BH_SRRa, BH_SRRb, BH_Equil_Rec), equilibrium recruitment for
#' either Beverton-Holt or Ricker relationship (Equil_Rec), equilibrium catch (Equil_Catch), equilibrium
#' female and male and spawning biomass (Equil_FemSpBiom, Equil_MalSpBiom), equilibrium relative female,
#' male and combined sex spawning biomass (Equilmod_SPR, Equilmod_MalRelBiom, Equilmod_CombSexRelBiom),
#' female and male length at age (FemLenAtAge, MalLenAtAge), female and male weight at age
#' (FemWtAtAge, MalWtAtAge), female and male proportion mature at age (FemPropMatAtAge, MalPropMatAtAge), female and
#' male gear selectivity at age (FemSelAtAge, MalSelAtAge), female and male retention at age probabilities,
#' (FemRetProbAtAge, MalRetProbAtAge), selectivity of female and male fish landings (FemSelLandAtAge, MalSelLandAtAge),
#' female and male selectivity of discards (FemSelDiscAtAge, MalSelDiscAtAge), female and male fishing mortality
#' associated with discarding (FemDiscFAtAge, MalDiscFAtAge), female and male fishing mortality associated with
#' fish landings (FemLandFAtAge, MalLandFAtAge), female and male total fishing mortality at age (FemFAtAge, MalFAtAge),
#' female and male total mortality at age (FemZAtAge, MalZAtAge).
#' @examples
#' # Example 1. Non-hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 20 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be approprriate to use a smaller time step)
#' Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.5, 0.5) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstFemLenAtAge=NA,
#'                           EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- NA # for log-log relationship
#' lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA,
#'                          EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#'
#' mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA,
#'                           EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- c(2.5, 2.5) # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- c(3.5, 3.5) # females, males - Logistic age fish retention at age parameters
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 0.2 # natural mortality  (year-1)
#' FMort <- 0.4 # estimate of fishing mortality, e.g. from catch curve analysis
#' Res=CalcYPRAndSPRForFMort(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                           lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
#'                           ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
#'                           mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
#'                           DiscMort, Steepness, SRrel_Type, NatMort, FMort)
#' # Example 2: hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 100 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be approprriate to use a smaller time step)
#' Linf <- c(1000, 1000) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.1, 0.1) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstMalLenAtAge=NA,
#'                           EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- -11.0 # for log-log relationship
#' lenwt_b <- 3.0 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 2 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA,
#'                          EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodPattern <- 2 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 1 # Ratio of females to males at age zero
#' FinalSex_Pmax <- 1 # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_A50 <- 35 # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- 60 # Logistic sex change relationship parameters (95% of max probability)
#' mat_A50 <- c(20, 20) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(30, 30) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA,                          EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_A50 <- c(15, 15) # females, males - Logistic age selectivity relationship parameters
#' sel_A95 <- c(25, 25) # females, males - Logistic age selectivity relationship parameters
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- c(15, 15) # females, males - Logistic age fish retention at age parameters (inflection point)
#' ret_A95 <- c(25, 25) # females, males - Logistic age fish retention at age parameters (95% of maximum retention)
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 0.07 # natural mortality  (year-1)
#' FMort <- 0.07 # estimate of fishing mortality, e.g. from catch curve analysis
#' Res=CalcYPRAndSPRForFMort(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                           lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
#'                           ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
#'                           mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
#'                           DiscMort, Steepness, SRrel_Type, NatMort, FMort)
#' @export
CalcYPRAndSPRForFMort <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                  lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
                                  ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                                  mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
                                  DiscMort, Steepness, SRrel_Type, NatMort, FMort) {

  # number of model time steps
  nTimeSteps <- 1 + (MaxModelAge / TimeStep)

  # Ages. As species is relatively short-lived, specify small time steps (e.g. 0.1 years)
  Ages <- seq(0,MaxModelAge,TimeStep)

  # calculate length at age - von Bertalanffy growth curve
  if (is.na(EstLenAtAge[1,1])) {
    FemLenAtAge <- Linf[1] * (1 - exp(-vbK[1] * (Ages - tzero[1])))
    MalLenAtAge <- Linf[2] * (1 - exp(-vbK[2] * (Ages - tzero[2])))
  }
  if (!is.na(EstLenAtAge[1,1])) {
    FemLenAtAge <- EstLenAtAge[,1]
    MalLenAtAge <- EstLenAtAge[,2]
  }

  # set negative lengths to zero
  FemLenAtAge[which(FemLenAtAge<0)] = 0
  MalLenAtAge[which(MalLenAtAge<0)] = 0

  # calculate weight (g) at age
  if (is.na(EstWtAtAge[1,1])) {
    if (WLrel_Type == 1) { # power relationship
      FemWtAtAge <- (lenwt_a * FemLenAtAge ^ lenwt_b) / 1000 # weight at length, kg
      MalWtAtAge <- (lenwt_a * MalLenAtAge ^ lenwt_b) / 1000 # weight at length, kg
    }
    if (WLrel_Type == 2) { # log-log relationship
      FemWtAtAge <- exp(ln_lenwt_a + lenwt_b * log(FemLenAtAge)) / 1000 # weight at length, kg
      MalWtAtAge <- exp(ln_lenwt_a + lenwt_b * log(MalLenAtAge)) / 1000 # weight at length, kg
    }
  }
  if (!is.na(EstWtAtAge[1,1])) {
    FemWtAtAge <- EstWtAtAge[,1]
    MalWtAtAge <- EstWtAtAge[,2]
  }


  # for sex-changing species, calculate proportion female at age
  if (ReprodPattern == 1) { # protogynous hermaphroditism (female to male sex change)
    PropFemAtAge = NA
  }
  if (ReprodPattern == 2) { # protogynous hermaphroditism (female to male sex change)
    PropFemAtAge = 1 - (FinalSex_Pmax / (1 + exp(-log(19) * (Ages-FinalSex_A50) / (FinalSex_A95 - FinalSex_A50))))
  }
  if (ReprodPattern == 3) { # protandrous hermaphroditism (male to female sex change)
    PropFemAtAge = FinalSex_Pmax / (1 + exp(-log(19) * (Ages-FinalSex_A50) / (FinalSex_A95 - FinalSex_A50)))
  }

  # calculate proportion mature at age
  if (is.na(EstMatAtAge[1,1])) {
    FemPropMatAtAge <- 1 / (1 + exp(-log(19) * (Ages-mat_A50[1]) / (mat_A95[1] - mat_A50[1])))
    MalPropMatAtAge <- 1 / (1 + exp(-log(19) * (Ages-mat_A50[2]) / (mat_A95[2] - mat_A50[2])))
  }
  if (!is.na(EstMatAtAge[1,1])) {
    FemPropMatAtAge <- EstMatAtAge[,1]
    MalPropMatAtAge <- EstMatAtAge[,2]
  }

  # Calculate gear selectivity at age
  FemGearSelAtAge <- 1/(1+exp(-log(19) * (Ages - sel_A50[1]) / (sel_A95[1] - sel_A50[1])))
  MalGearSelAtAge <- 1/(1+exp(-log(19) * (Ages - sel_A50[2]) / (sel_A95[2] - sel_A50[2])))

  # Calculate fish retention at age
  FemRetProbAtAge <- ret_Pmax[1] / (1+exp(-log(19) * (Ages - ret_A50[1]) / (ret_A95[1] - ret_A50[1])))
  MalRetProbAtAge <- ret_Pmax[2] / (1+exp(-log(19) * (Ages - ret_A50[2]) / (ret_A95[2] - ret_A50[2])))

  # Calculate selectivity of landings at age
  FemSelLandAtAge <- FemGearSelAtAge * FemRetProbAtAge
  MalSelLandAtAge <- MalGearSelAtAge * MalRetProbAtAge

  # Calculate selectivity of discards at age
  FemSelDiscAtAge <- FemGearSelAtAge * (1 - FemRetProbAtAge)
  MalSelDiscAtAge <- MalGearSelAtAge * (1 - MalRetProbAtAge)

  # calculate female and male mortality at age associated with discarding of undersize fish
  FemDiscFAtAge <- FMort * DiscMort * FemSelDiscAtAge
  MalDiscFAtAge <- FMort * DiscMort * MalSelDiscAtAge

  # calculate female and male mortality at age associated with landings
  FemLandFAtAge <- FMort * FemSelLandAtAge
  MalLandFAtAge <- FMort * MalSelLandAtAge

  # calculate (total) female and male fishing mortality at age
  if (DiscMort == 0) {
    FemFAtAge <- FMort * FemRetProbAtAge
    MalFAtAge <- FMort * MalRetProbAtAge
  } else {
    FemFAtAge <- FMort * (FemSelLandAtAge + (DiscMort * FemSelDiscAtAge))
    MalFAtAge <- FMort * (MalSelLandAtAge + (DiscMort * MalSelDiscAtAge))
  }

  # calculate female and male total mortality at age
  FemZAtAge <- FemFAtAge + NatMort
  MalZAtAge <- MalFAtAge + NatMort

  # calculate relative unfished female and male survival at age
  UnfishFemSurvAtAge <- rep(0,nTimeSteps)
  UnfishMalSurvAtAge <- rep(0,nTimeSteps)
  tempUnfishFemSurvAtAge <- rep(0,nTimeSteps)
  tempUnfishMalSurvAtAge <- rep(0,nTimeSteps)
  for (j in 1:nTimeSteps) {
    if (j==1) { # age=0
      UnfishFemSurvAtAge[j] <- InitRatioFem
      UnfishMalSurvAtAge[j] <- (1 - InitRatioFem)
    }
    else if (j < nTimeSteps) {
      UnfishFemSurvAtAge[j] <- UnfishFemSurvAtAge[j-1] * exp(-NatMort * TimeStep)
      UnfishMalSurvAtAge[j] <- UnfishMalSurvAtAge[j-1] * exp(-NatMort * TimeStep)
    }
    else if (j==nTimeSteps) { # maximum model age, plus group
      UnfishFemSurvAtAge[j] <- (UnfishFemSurvAtAge[j-1] * exp(-NatMort * TimeStep)) /
        (1 - exp(-NatMort * TimeStep))
      UnfishMalSurvAtAge[j] <- (UnfishMalSurvAtAge[j-1] * exp(-NatMort * TimeStep)) /
        (1 - exp(-NatMort * TimeStep))
    }

    if (ReprodPattern > 1) { # hermaphroditic species
      # assuming that, regardless of fishing pattern, the sex ratio at age is determined by
      # by the logistic curve describing probability of sex change at age
      tempUnfishFemSurvAtAge[j] = PropFemAtAge[j] * (UnfishFemSurvAtAge[j] + UnfishMalSurvAtAge[j])
      tempUnfishMalSurvAtAge[j] = (1 - PropFemAtAge[j]) * (UnfishFemSurvAtAge[j] + UnfishMalSurvAtAge[j])
      UnfishFemSurvAtAge[j] = tempUnfishFemSurvAtAge[j]
      UnfishMalSurvAtAge[j] = tempUnfishMalSurvAtAge[j]
    }
  }

  # calculate female and male unfished spawning biomass at age
  UnfishFemBiomAtAge <- UnfishFemSurvAtAge * FemPropMatAtAge * FemWtAtAge
  UnfishMalBiomAtAge <- UnfishMalSurvAtAge * MalPropMatAtAge * MalWtAtAge

  # Unfished spawning biomass in kg
  UnfishFemSpawnBiom <- sum(UnfishFemBiomAtAge)
  UnfishMalSpawnBiom <- sum(UnfishMalBiomAtAge)
  UnfishCombSexSpawnBiom <- UnfishFemSpawnBiom + UnfishMalSpawnBiom

  # calculate female and male survival at age for fished population
  FishedFemSurvAtAge <- rep(0,nTimeSteps)
  FishedMalSurvAtAge <- rep(0,nTimeSteps)
  tempFishedFemSurvAtAge <- rep(0,nTimeSteps)
  tempFishedMalSurvAtAge <- rep(0,nTimeSteps)
  for (j in 1:nTimeSteps) {
    if (j==1) { # age=0
      FishedFemSurvAtAge[j] <- InitRatioFem
      FishedMalSurvAtAge[j] <- (1 - InitRatioFem)
    }
    else if (j < nTimeSteps) {
      FishedFemSurvAtAge[j] <- FishedFemSurvAtAge[j-1] * exp(-FemZAtAge[j-1] * TimeStep)
      FishedMalSurvAtAge[j] <- FishedMalSurvAtAge[j-1] * exp(-MalZAtAge[j-1] * TimeStep)
    }
    else if (j==nTimeSteps) { # maximum model age, plus group
      FishedFemSurvAtAge[j] <- FishedFemSurvAtAge[j-1] * exp(-FemZAtAge[j-1] * TimeStep) /
        (1 - exp(-(FMort+NatMort)*TimeStep))
      FishedMalSurvAtAge[j] <- FishedMalSurvAtAge[j-1] * exp(-MalZAtAge[j-1] * TimeStep) /
        (1 - exp(-(FMort+NatMort)*TimeStep))
    }

    if (ReprodPattern > 1) { # hermaphroditic species
      # assuming that, regardless of fishing pattern, the sex ratio at age is determined by
      # by the logistic curve describing probability of sex change at age
      tempFishedFemSurvAtAge[j] = PropFemAtAge[j] * (FishedFemSurvAtAge[j] + FishedMalSurvAtAge[j])
      tempFishedMalSurvAtAge[j] = (1 - PropFemAtAge[j]) * (FishedFemSurvAtAge[j] + FishedMalSurvAtAge[j])
      FishedFemSurvAtAge[j] = tempFishedFemSurvAtAge[j]
      FishedMalSurvAtAge[j] = tempFishedMalSurvAtAge[j]
    }
  }

  # calculate female and male mature biomass at age for fished population
  FishedFemBiomAtAge <- FishedFemSurvAtAge * FemPropMatAtAge * FemWtAtAge
  FishedMalBiomAtAge <- FishedMalSurvAtAge * MalPropMatAtAge * MalWtAtAge

  # calculate female and male catch at age (in numbers) - Baranov catch equation
  FemCatchAtAgeNum <- FishedFemSurvAtAge * (FemFAtAge/FemZAtAge) *
    (1 - exp(-(FemZAtAge * TimeStep)))
  MalCatchAtAgeNum <- FishedMalSurvAtAge * (MalFAtAge/MalZAtAge) *
    (1 - exp(-(MalZAtAge * TimeStep)))
  #sum(FemCatchAtAgeNum)

  # calculate female and male catch at age (in biomass)
  FemCatchAtAge <- FemCatchAtAgeNum * FemWtAtAge
  MalCatchAtAge <- MalCatchAtAgeNum * MalWtAtAge

  # calculate yield per recruit in kg
  YPR <- sum(FemCatchAtAge) + sum(MalCatchAtAge)

  # calculate female and male spawning biomass per recruit in kg
  FishFemSpawnBiom <- sum(FishedFemBiomAtAge)
  FishMalSpawnBiom <- sum(FishedMalBiomAtAge)
  FishCombSexSpawnBiom <- FishFemSpawnBiom + FishMalSpawnBiom

  # calculate spawning potential ratio (SPR)
  Fem_SPR <- FishFemSpawnBiom / UnfishFemSpawnBiom
  Mal_SPR <- FishMalSpawnBiom / UnfishMalSpawnBiom
  CombSex_SPR <- (FishFemSpawnBiom + FishMalSpawnBiom)  / (UnfishFemSpawnBiom + UnfishMalSpawnBiom)

  if (ReprodPattern == 1) { # gonochoristic species
    UnfishSpawnBiom = UnfishFemSpawnBiom
    FishSpawnBiom = FishFemSpawnBiom
  }
  if (ReprodPattern > 1) { # hermaphroditic species
    UnfishSpawnBiom = UnfishCombSexSpawnBiom
    FishSpawnBiom = FishCombSexSpawnBiom
  }

  if (SRrel_Type == 1) { # Beverton-Holt

    # Define Beverton and Holt stock-recruitment relationship
    # to account for impacts of fishing on recruitment, i.e. through its impact on spawning biomass
    # (which is typically ignored in standard per recruit models)
    BH_SRRa <- (UnfishSpawnBiom / 1.0) * ((1-Steepness) / (4*Steepness))
    BH_SRRb <- (Steepness - 0.2) / (0.8 * Steepness * 1.0)

    # Check that the specified initial recruitment can be recovered,
    # given the S_R parameters and unfished female biomass
    # (Check <- UnfishSpawnBiom/(BH_SRRa + BH_SRRb * UnfishSpawnBiom))

    # calculate equilibrium recruitment
    BH_Equil_Rec <- (FishSpawnBiom-BH_SRRa) / (BH_SRRb*FishSpawnBiom)

    # Equations from Norm Hall - Beverton and Holt and Ricker relationship,
    # when both parameterised using steepness
  }
  if (SRrel_Type == 2) { # Ricker
    BH_SRRa = NA
    BH_SRRb = NA
    BH_Equil_Rec = NA
  }

  # Alternative method for deriving equilbrium recruitment
  if (SRrel_Type == 1) { # Beverton-Holt
    Equil_Rec = ((4 * Steepness * FishSpawnBiom) - (1 - Steepness) * UnfishSpawnBiom) /
      (5 * (Steepness - 0.2) * FishSpawnBiom)
  }

  if (SRrel_Type == 2) { # Ricker
    # Ricker
    Equil_Rec <- (UnfishSpawnBiom * (1- ((4 * log(UnfishSpawnBiom/FishSpawnBiom)) /
                                           (5 * log(5 * Steepness))))) / FishSpawnBiom
  }

  # calculate equilibrium catch
  Equil_Catch <- Equil_Rec * YPR

  # calculate equilibrium female spawning biomass
  Equil_FemSpBiom <- Equil_Rec * FishFemSpawnBiom
  Equil_MalSpBiom <- Equil_Rec * FishMalSpawnBiom

  # calculate equilibrium model SPR
  Equilmod_FemRelBiom <- Equil_FemSpBiom / UnfishFemSpawnBiom
  Equilmod_MalRelBiom <- Equil_MalSpBiom / UnfishMalSpawnBiom
  Equilmod_CombSexRelBiom <- (Equil_FemSpBiom + Equil_MalSpBiom) / (UnfishFemSpawnBiom + UnfishMalSpawnBiom)

  Results = list(Ages = Ages,
                 YPR = YPR,
                 Fem_SPR = Fem_SPR,
                 Mal_SPR = Mal_SPR,
                 CombSex_SPR = CombSex_SPR,
                 UnfishFemSurvAtAge = UnfishFemSurvAtAge,
                 UnfishMalSurvAtAge = UnfishMalSurvAtAge,
                 UnfishFemBiomAtAge = UnfishFemBiomAtAge,
                 UnfishMalBiomAtAge = UnfishMalBiomAtAge,
                 FishedFemSurvAtAge = FishedFemSurvAtAge,
                 FishedMalSurvAtAge = FishedMalSurvAtAge,
                 FishedFemBiomAtAge = FishedFemBiomAtAge,
                 FishedMalBiomAtAge = FishedMalBiomAtAge,
                 FemCatchAtAgeNum = FemCatchAtAgeNum,
                 MalCatchAtAgeNum = MalCatchAtAgeNum,
                 FemCatchAtAge = FemCatchAtAge,
                 MalCatchAtAge = MalCatchAtAge,
                 BH_SRRa = BH_SRRa,
                 BH_SRRb = BH_SRRb,
                 BH_Equil_Rec = BH_Equil_Rec,
                 Equil_Rec = Equil_Rec,
                 Equil_Catch = Equil_Catch,
                 Equil_FemSpBiom = Equil_FemSpBiom,
                 Equil_MalSpBiom = Equil_MalSpBiom,
                 Equilmod_FemRelBiom = Equilmod_FemRelBiom,
                 Equilmod_MalRelBiom = Equilmod_MalRelBiom,
                 Equilmod_CombSexRelBiom = Equilmod_CombSexRelBiom,
                 FemLenAtAge = FemLenAtAge,
                 MalLenAtAge = MalLenAtAge,
                 FemWtAtAge = FemWtAtAge,
                 MalWtAtAge = MalWtAtAge,
                 PropFemAtAge = PropFemAtAge,
                 FemPropMatAtAge = FemPropMatAtAge,
                 MalPropMatAtAge = MalPropMatAtAge,
                 FemGearSelAtAge = FemGearSelAtAge,
                 MalGearSelAtAge = MalGearSelAtAge,
                 FemRetProbAtAge = FemRetProbAtAge,
                 MalRetProbAtAge = MalRetProbAtAge,
                 FemSelLandAtAge = FemSelLandAtAge,
                 MalSelLandAtAge = MalSelLandAtAge,
                 FemSelDiscAtAge = FemSelDiscAtAge,
                 MalSelDiscAtAge = MalSelDiscAtAge,
                 FemDiscFAtAge = FemDiscFAtAge,
                 MalDiscFAtAge = MalDiscFAtAge,
                 FemLandFAtAge = FemLandFAtAge,
                 MalLandFAtAge = MalLandFAtAge,
                 FemFAtAge = FemFAtAge,
                 MalFAtAge = MalFAtAge,
                 FemZAtAge = FemZAtAge,
                 MalZAtAge = MalZAtAge)

  return(Results)
}


#' Get outputs from per recruit analysis across a range of fishing mortality values
#'
#' This function provides outputs associated with per recruit analysis, and an
#' extended form of this analysis with a Beverton-Holt stock recruitment relationship to account
#' for potential impacts of fishing on recruitment. Outputs are provided for a range of
#' fishing mortality values, including the current, estimated value.
#'
#' Function requires values for the initial sex ratio (i.e. of juvenile recruits),
#' maximum age, model timestep (e.g. 1 = annual, 1/12 = monthly), von Bertalanffy growth
#' parameters, weight-length parameters, i.e using a power relationship, logistic length at
#' maturity parameters, logistic selectivity parameters, steepness parameter of the stock-recruitment
#' relationship, natural mortality (NatMort), fishing mortality (FMort.
#' (Numbers and a Numeric Vectors)
#' @param params maximum age considered in model (MaxModelAge), model timestep, e.g. 1 = annual, 1/12 = monthly (TimeStep),
#' von Bertalanffy growth parameters (Linf, vbK, tzero), Female and male estimated lengths-at-age (set to NA if
#' using von Bertalanffy growth equation, or set the von Bertalanffy growth parameters to NA and input as data frame
#' (columns for sex) from any growth model (EstLenAtAge), weight-length parameters, using a power or log-log relationship
#' (lenwt_a, ln_lenwt_a, lenwt_b), weight-length relationship type (WLrel_Type), Female and male estimated
#' weights-at-ages  (set to NA if weight-length relationship parameters, or set the weight-length parameters
#' to NA and input these as data frame (columns for sex), (EstWtAtAge), reproductive pattern, 1=gonochoristic,
#' 2=protogynous, 3=protandrous (ReprodPattern), initial sex ratio of juvenile recruits (InitRatioFem),
#' Maximum proportion of final sex for hermaphroditic species (set to NA for gonochoristic species) (FinalSex_Pmax),
#' logistic sex change parameters for hermaphroditic species (set to NA for gonochoristic species) (FinalSex_A50, FinalSex_A95),
#' logistic length at maturity parameters (mat_L50, mat_L95), (set to NA if using age at maturity parameters, or set the
#' age at maturity parameters to NA and input these as data fra,e (columns for sex) (EstMatAtAge), logistic gear selectivity
#' parameters (sel_A50, sel_A95),logistic retention curve parameters (ret_Pmax, ret_A50, ret_A95), discard mortality rate
#' (DiscMort), steepness parameter of the stock-recruitment relationship (Steepness), stock-recruitment relationship type,
#' 1 = Beverton-Holt, 2=Ricker (SRrel_Type), natural mortality (NatMort), current fishing mortality (Current_F).
#' @return Age at each timestep (Ages), female and male length at age (FemLenAtAge, MalLenAtAge),
#' female and male weight at age (FemWtAtAge, FemWtAtAge), proportion female at age (PropFemAtAge), female and male
#' proportion mature at age (FemPropMatAtAge, MalPropMatAtAge), female and male selectivity at age (FemSelAtAge, MalSelAtAge),
#' yield-per-recruit (YPR), unfished female, male and combined sex spawning potential ratio (Fem_SPR, Mal_SPR,
#' CombSex_SPR), unfished female and male per recruit survival at age (UnfishFemSurvAtAge, UnfishMalSurvAtAge),
#' unfished female and male mature biomass at age (UnfishFemBiomAtAge, UnfishMalBiomAtAge), fished female and
#' male spawning potential ratio (Fem_SPR, Mal_SPR), fished female and male per recruit survival at age
#' (UnfishFemSurvAtAge, UnfishMalSurvAtAge), fished female and male mature biomass at age (UnfishFemBiomAtAge,
#' UnfishMalBiomAtAge), female and male per recruit catches at age in numbers (FemCatchAtAgeNum, MalCatchAtAgeNum)
#' and biomass (FemCatchAtAge, MalCatchAtAge), and for the extended model, equilbrium recruitment (Equil_Rec),
#' equilbrium catch (Equil_Catch), equilbrium female and male spawning biomass (Equil_FemSpBiom, Equil_MalSpBiom) and
#' relative female, male and combined sex spawning biomass (Equilmod_FemRelBiom, Equilmod_MalRelBiom, Equilmod_CombSexRelBiom),
#' female and male gear selectivity at age (FemSelAtAge, MalSelAtAge), female and male retention at age probabilities,
#' (FemRetProbAtAge, MalRetProbAtAge), selectivity of female and male fish landings (FemSelLandAtAge, MalSelLandAtAge),
#' female and male selectivity of discards (FemSelDiscAtAge, MalSelDiscAtAge), female and male fishing mortality
#' associated with discarding (FemDiscFAtAge, MalDiscFAtAge), female and male fishing mortality associated with
#' fish landings (FemLandFAtAge, MalLandFAtAge), female and male total fishing mortality at age (FemFAtAge, MalFAtAge),
#' female and male total mortality at age (FemZAtAge, MalZAtAge)range of fishing mortality values applied for which per quantities are calculated
#' (FishMort = seq(0,2,0.01)), equilbrium recruitment vs FMort (Equil_Rec), equilbrium catch vs FMort (Equil_Catch),
#' equilbrium female spawning biomass vs FMort (Equil_FemSpBiom), equilbrium spawning potential ratio vs FMort
#' (Equilmod_SPR), maximum yield per recruit (maxypr), maximum equilbrium catch (maxeqCatch), fishing mortality
#' associated with maxypr (Fmax), fishing mortality associated with maxeqCatch (FmaxeqCatch), YPR vs FishMort
#' (YPRResults), Equil_Catch vs FishMort (EquilCatchResults), Equilmod_FemRelBiom vs FishMort (Fem_SPRResults),
#' Equilmod_MalRelBiom vs FishMort (Mal_SPRResults), Equilmod_CombSexRelBiom vs FishMort (CombSex_SPRResults),
#' Equil_Rec vs FishMort (EquilRecResults)
#' @examples
#' # Example 1. Non-hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 20 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be approprriate to use a smaller time step)
#' Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.5, 0.5) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstFemLenAtAge=NA,
#'                           EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- NA # for log-log relationship
#' lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA,
#'                          EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#'
#' mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA,
#'                           EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- c(2.5, 2.5) # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- c(3.5, 3.5) # females, males - Logistic age fish retention at age parameters
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 0.2 # natural mortality  (year-1)
#' Current_F <- 0.07 # estimate of fishing mortality, e.g. from catch curve analysis
#' Res = GetPerRecruitResults(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                            lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
#'                            ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
#'                            mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
#'                            DiscMort, Steepness, SRrel_Type, NatMort, Current_F)
#' # Example 2: hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 100 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be approprriate to use a smaller time step)
#' Linf <- c(1000, 1000) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.1, 0.1) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstMalLenAtAge=NA,
#'                           EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- -11.0 # for log-log relationship
#' lenwt_b <- 3.0 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 2 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA,
#'                          EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodPattern <- 2 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 1 # Ratio of females to males at age zero
#' FinalSex_Pmax <- 1 # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_A50 <- 35 # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- 60 # Logistic sex change relationship parameters (95% of max probability)
#' mat_A50 <- c(20, 20) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(30, 30) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA,                          EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_A50 <- c(15, 15) # females, males - Logistic age selectivity relationship parameters
#' sel_A95 <- c(25, 25) # females, males - Logistic age selectivity relationship parameters
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- c(15, 15) # females, males - Logistic age fish retention at age parameters (inflection point)
#' ret_A95 <- c(25, 25) # females, males - Logistic age fish retention at age parameters (95% of maximum retention)
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 0.07 # natural mortality  (year-1)
#' Current_F <- 0.07 # estimate of fishing mortality, e.g. from catch curve analysis
#' Res = GetPerRecruitResults(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                            lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
#'                            ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
#'                            mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
#'                            DiscMort, Steepness, SRrel_Type, NatMort, Current_F)
#' @export
GetPerRecruitResults <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                 lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
                                 ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                                 mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
                                 DiscMort, Steepness, SRrel_Type, NatMort, Current_F) {

  FishMort <- seq(0,2,0.01)
  nFVals <- length(FishMort) # fishing mortality
  YPRResults <- rep(0,nFVals)
  Fem_SPRResults <- rep(0,nFVals)
  Mal_SPRResults <- rep(0,nFVals)
  CombSex_SPRResults <- rep(0,nFVals)

  EquilRecResults <- rep(0,nFVals)
  EquilCatchResults <- rep(0,nFVals)
  Equil_FemSpBiomResults <- rep(0,nFVals)
  Equil_MalSpBiomResults <- rep(0,nFVals)

  Equilmod_FemRelBiomResults <- rep(0,nFVals)
  Equilmod_MalRelBiomResults <- rep(0,nFVals)
  Equilmod_CombSexRelBiomResults <- rep(0,nFVals)


  for (k in 1:nFVals) {
    FMort = FishMort[k]
    Res = CalcYPRAndSPRForFMort(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
                                ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                                mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
                                DiscMort, Steepness, SRrel_Type, NatMort, FMort)
    # per recruit results
    YPRResults[k] <- Res$YPR
    Fem_SPRResults[k] = Res$Fem_SPR
    Mal_SPRResults[k] = Res$Mal_SPR
    CombSex_SPRResults[k] = Res$CombSex_SPR
    # extended model results
    EquilRecResults[k] = Res$Equil_Rec
    EquilCatchResults[k] = Res$Equil_Catch
    Equil_FemSpBiomResults[k] = Res$Equil_FemSpBiom
    Equil_MalSpBiomResults[k] = Res$Equil_MalSpBiom
    Equilmod_FemRelBiomResults[k] = Res$Equilmod_FemRelBiom
    Equilmod_MalRelBiomResults[k] = Res$Equilmod_MalRelBiom
    Equilmod_CombSexRelBiomResults[k] = Res$Equilmod_CombSexRelBiom
  }
  maxypr <- max(YPRResults) # maximum yield per recruit
  maxeqCatch <- max(EquilCatchResults) # maximum equilbrium catch
  Fmax <- FishMort[which(YPRResults==maxypr)] # fishing mortality at ypr maximum
  FmaxeqCatch <- FishMort[which(EquilCatchResults==maxeqCatch)] # fishing mortality at maxeqCatch

  # calc FMSY
  F_MSY=FmaxeqCatch
  x=which(FishMort==F_MSY)
  # calc fem biomass ratio, at BMSY, 0.5BMSY and 1.2BMSY

  if (ReprodPattern == 1) { # gonochoristic species
    BMSY_Thresh=Equilmod_FemRelBiomResults[x]
  }
  if (ReprodPattern > 1) { # hermaphroditic species
    BMSY_Thresh=Equilmod_CombSexRelBiomResults[x]
  }
  BMSY_Lim=0.5*BMSY_Thresh
  BMSY_Targ=1.2*BMSY_Thresh

  # get results for current F
  FMort = Current_F
  Res2 = CalcYPRAndSPRForFMort(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                               lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
                               ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                               mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
                               DiscMort, Steepness, SRrel_Type, NatMort, FMort)



  Results = list(Ages = Res2$Ages,
                 YPR = Res2$YPR,
                 Fem_SPR = Res2$Fem_SPR,
                 Mal_SPR = Res2$Mal_SPR,
                 CombSex_SPR = Res2$CombSex_SPR,
                 FemLenAtAge = Res2$FemLenAtAge,
                 MalLenAtAge = Res2$MalLenAtAge,
                 FemWtAtAge = Res2$FemWtAtAge,
                 MalWtAtAge = Res2$MalWtAtAge,
                 PropFemAtAge = Res2$PropFemAtAge,
                 FemPropMatAtAge = Res2$FemPropMatAtAge,
                 MalPropMatAtAge = Res2$MalPropMatAtAge,
                 UnfishFemSurvAtAge = Res2$UnfishFemSurvAtAge,
                 UnfishMalSurvAtAge = Res2$UnfishMalSurvAtAge,
                 UnfishFemBiomAtAge = Res2$UnfishFemBiomAtAge,
                 UnfishMalBiomAtAge = Res2$UnfishMalBiomAtAge,
                 FishedFemSurvAtAge = Res2$FishedFemSurvAtAge,
                 FishedMalSurvAtAge = Res2$FishedMalSurvAtAge,
                 FishedFemBiomAtAge = Res2$FishedFemBiomAtAge,
                 FishedMalBiomAtAge = Res2$FishedMalBiomAtAge,
                 FemCatchAtAgeNum = Res2$FemCatchAtAgeNum,
                 MalCatchAtAgeNum = Res2$MalCatchAtAgeNum,
                 FemCatchAtAge = Res2$FemCatchAtAge,
                 MalCatchAtAge = Res2$MalCatchAtAge,
                 Equil_Rec = Res2$Equil_Rec,
                 Equil_Catch = Res2$Equil_Catch,
                 Equil_FemSpBiom = Res2$Equil_FemSpBiom,
                 Equil_MalSpBiom = Res2$Equil_MalSpBiom,
                 Equilmod_FemRelBiom = Res2$Equilmod_FemRelBiom,
                 Equilmod_MalRelBiom = Res2$Equilmod_MalRelBiom,
                 Equilmod_CombSexRelBiom = Res2$Equilmod_CombSexRelBiom,
                 BH_SRRa = Res2$BH_SRRa,
                 BH_SRRb = Res2$BH_SRRb,
                 BH_Equil_Rec = Res2$BH_Equil_Rec,
                 FemGearSelAtAge = Res2$FemGearSelAtAge,
                 MalGearSelAtAge = Res2$MalGearSelAtAge,
                 FemRetProbAtAge = Res2$FemRetProbAtAge,
                 MalRetProbAtAge = Res2$MalRetProbAtAge,
                 FemSelLandAtAge = Res2$FemSelLandAtAge,
                 MalSelLandAtAge = Res2$MalSelLandAtAge,
                 FemSelDiscAtAge = Res2$FemSelDiscAtAge,
                 MalSelDiscAtAge = Res2$MalSelDiscAtAge,
                 FemDiscFAtAge = Res2$FemDiscFAtAge,
                 MalDiscFAtAge = Res2$MalDiscFAtAge,
                 FemLandFAtAge = Res2$FemLandFAtAge,
                 MalLandFAtAge = Res2$MalLandFAtAge,
                 FemFAtAge = Res2$FemFAtAge,
                 MalFAtAge = Res2$MalFAtAge,
                 FemZAtAge = Res2$FemZAtAge,
                 MalZAtAge = Res2$MalZAtAge,
                 FishMort = FishMort,
                 maxypr = maxypr,
                 F_MSY = F_MSY,
                 Fmax = Fmax,
                 BMSY_Targ=BMSY_Targ,
                 BMSY_Thresh=BMSY_Thresh,
                 BMSY_Lim=BMSY_Lim,
                 FmaxeqCatch = FmaxeqCatch,
                 YPRResults = YPRResults,
                 EquilCatchResults = EquilCatchResults,
                 Fem_SPRResults = Fem_SPRResults,
                 Mal_SPRResults = Mal_SPRResults,
                 CombSex_SPRResults = CombSex_SPRResults,
                 Equilmod_FemRelBiomResults = Equilmod_FemRelBiomResults,
                 Equilmod_MalRelBiomResults = Equilmod_MalRelBiomResults,
                 Equilmod_CombSexRelBiomResults = Equilmod_CombSexRelBiomResults,
                 EquilRecResults = EquilRecResults)
  return(Results)

}

# function to calculate x axis maximum and interval
Get_xaxis_scale <- function(x_data) {

  xx = floor(log10(max(x_data)))
  temp_xmin = floor(min(0.8 * x_data))
  xmin = round(temp_xmin, -xx)
  temp_xmax = ceiling(1.2*max(x_data))
  temp_xmax2 = round(temp_xmax, -xx)
  xint = (temp_xmax2-xmin)/4
  xmax = xmin + (4 * xint)
  results = list(xmin = xmin,
                 xmax = xmax,
                 xint = xint)


  return(results)
}


# function to calculate y axis maximum and interval
Get_yaxis_scale <- function(y_data) {

  # get maximum y value
  x = max(1.35*y_data)

  # get number of digits
  xx = floor(log10(max(y_data)))

  # get rounded maximum value
  xxx = round(x, -xx)

  # get interval, with first number of 1, 2 or 5
  xxxx = xxx / 3

  # get y interval
  yint = round(xxxx, -floor(log10(max(xxxx))))

  ymax = yint * 4

  results = list(ymax = ymax,
                 yint = yint)

  return(results)
}

#' Get plots associated with per recruit analysis and extended analysis
#' with a stock-recruitment relationship
#'
#' This function provides a range of plots associated with per recruit analysis, and an
#' extended form of this analysis with a Beverton-Holt stock recruitment relationship to account
#' for potential impacts of fishing on recruitment. Outputs include information on biology
#' inputs and various per recruit outputs.
#'
#' Function requires values for the initial sex ratio (i.e. of juvenile recruits),
#' maximum age, model timestep (e.g. 1 = annual, 1/12 = monthly), von Bertalanffy growth
#' parameters, weight-length parameters, i.e using a power relationship, logistic length at
#' maturity parameters, logistic selectivity parameters, steepness parameter of the stock-recruitment
#' relationship, natural mortality (NatMort), current fishing mortality (Current_F)
#'
#' (Numbers and a Numeric Vectors)
#' @param maximum age considered in model (MaxModelAge), model timestep, e.g. 1 = annual, 1/12 = monthly (TimeStep),
#' von Bertalanffy growth parameters (Linf, vbK, tzero), Female and male estimated lengths-at-age (set to NA if
#' using von Bertalanffy growth equation, or set the von Bertalanffy growth parameters to NA and input as data frame
#' (columns for sex) from any growth model (EstLenAtAge), weight-length parameters, using a power or log-log relationship
#' (lenwt_a, ln_lenwt_a, lenwt_b), weight-length relationship type (WLrel_Type), Female and male estimated
#' weights-at-ages  (set to NA if weight-length relationship parameters, or set the weight-length parameters
#' to NA and input these as data frame (columns for sex), (EstWtAtAge), reproductive pattern, 1=gonochoristic,
#' 2=protogynous, 3=protandrous (ReprodPattern), initial sex ratio of juvenile recruits (InitRatioFem),
#' Maximum proportion of final sex for hermaphroditic species (set to NA for gonochoristic species) (FinalSex_Pmax),
#' logistic sex change parameters for hermaphroditic species (set to NA for gonochoristic species) (FinalSex_A50, FinalSex_A95),
#' logistic length at maturity parameters (mat_L50, mat_L95), (set to NA if using age at maturity parameters, or set the
#' age at maturity parameters to NA and input these as data fra,e (columns for sex) (EstMatAtAge), logistic gear selectivity
#' parameters (sel_A50, sel_A95),logistic retention curve parameters (ret_Pmax, ret_A50, ret_A95), discard mortality rate
#' (DiscMort), steepness parameter of the stock-recruitment relationship (Steepness), stock-recruitment relationship type,
#' 1 = Beverton-Holt, 2=Ricker (SRrel_Type), natural mortality (NatMort), option for plotting reference points on final plot,
#'1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points (RefPointPlotOpt), current fishing mortality (Current_F).
#' @return
#' plots inlcuding sex specific growth curves, weight at age curves, selectivity vs maturity at age curves,
#' total and natural mortality at age curves, fished vs unfished survival at age, fished vs unfished
#' female mature biomass at age, catch at age, at current estimated fishing mortality, yield per recruit
#' and equilibrium catch vs fishing mortality values, i.e. FishMort = seq(0,3,0.01), spawning potential
#' ratio and relative  spawning biomass vs FishMort, equilbrium recruitment vs FishMort.
#' @examples
#' # Example 1. Non-hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 20 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be approprriate to use a smaller time step)
#' Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.5, 0.5) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstMalLenAtAge=NA,
#'                           EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- NA # for log-log relationship
#' lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA,
#'                          EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#'
#' mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA,
#'                           EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- c(2.5, 2.5) # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- c(3.5, 3.5) # females, males - Logistic age fish retention at age parameters
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 0.2 # natural mortality  (year-1)
#' Current_F <- 0.07 # estimate of fishing mortality, e.g. from catch curve analysis
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' PlotPerRecruitResults(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                       lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
#'                       ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
#'                       mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
#'                       DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F)
#' # Example 2: hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 100 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be approprriate to use a smaller time step)
#' Linf <- c(1000, 1000) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.1, 0.1) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstMalLenAtAge=NA,
#'                           EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- -11.0 # for log-log relationship
#' lenwt_b <- 3.0 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 2 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA,
#'                          EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodPattern <- 2 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 1 # Ratio of females to males at age zero
#' FinalSex_Pmax <- 1 # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_A50 <- 35 # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- 60 # Logistic sex change relationship parameters (95% of max probability)
#' mat_A50 <- c(20, 20) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(30, 30) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA,                          EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_A50 <- c(15, 15) # females, males - Logistic age selectivity relationship parameters
#' sel_A95 <- c(25, 25) # females, males - Logistic age selectivity relationship parameters
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- c(15, 15) # females, males - Logistic age fish retention at age parameters (inflection point)
#' ret_A95 <- c(25, 25) # females, males - Logistic age fish retention at age parameters (95% of maximum retention)
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 0.07 # natural mortality  (year-1)
#' Current_F <- 0.07 # estimate of fishing mortality, e.g. from catch curve analysis
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' PlotPerRecruitResults(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                       lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
#'                       ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
#'                       mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
#'                       DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F)
#' @export
PlotPerRecruitResults <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                  lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
                                  ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                                  mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
                                  DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F) {

  Res = GetPerRecruitResults(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
                             ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                             mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
                             DiscMort, Steepness, SRrel_Type, NatMort, Current_F)


  #Plot 1:
  par(mfrow = c(2,2), mar=c(3.5,4,2,2),
      oma=c(1,1,1,1),tck=-0.03,mgp = c(3, 0.5, 0))

  # plot growth curve
  y1=max(Res$FemLenAtAge)
  y2=max(Res$MalLenAtAge)
  if (y1 > y2) {
    ylims = Get_yaxis_scale(Res$FemLenAtAge)
  } else {
    ylims = Get_yaxis_scale(Res$MalLenAtAge)
  }
  ymax = ylims$ymax; yint = ylims$yint
  plot(Res$Ages,Res$FemLenAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,MaxModelAge),
       col="red",yaxt="n",xaxt="n",ylab="",xlab="")
  lines(Res$Ages, Res$MalLenAtAge,col="blue")
  axis(1,at=seq(0,MaxModelAge,1), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,MaxModelAge,1), labels=seq(0,MaxModelAge,1),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Total length (mm"),plain(")"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('bottomright', col=c("red","blue"),legend=c("Female","Male"),bty='n', cex=0.8,lwd=1.75)

  # plot weight at age
  y1=max(Res$FemWtAtAge)
  y2=max(Res$MalWtAtAge)
  if (y1 > y2) {
    ylims = Get_yaxis_scale(Res$FemWtAtAge)
  } else {
    ylims = Get_yaxis_scale(Res$MalWtAtAge)
  }
  ymax = ylims$ymax; yint = ylims$yint
  plot(Res$Ages,Res$FemWtAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,MaxModelAge),col="red",yaxt="n",xaxt="n",
       ylab="",xlab="")
  lines(Res$Ages,Res$MalWtAtAge,col="blue")
  axis(1,at=seq(0,MaxModelAge,1), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,MaxModelAge,1), labels=seq(0,MaxModelAge,1),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Total weight (kg"),plain(")"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('bottomright', col=c("red","blue"),legend=c("Female","Male"),bty='n', cex=0.8,lwd=1.75)

  # plot female maturity and selectivity at age
  plot(Res$Ages,Res$FemPropMatAtAge,"l", pch=16,frame.plot=F,ylim=c(0,1),xlim=c(0,MaxModelAge),col="red",yaxt="n",xaxt="n",
       ylab="",xlab="",cex=0.8)
  if (DiscMort == 0) {
    lines(Res$Ages, Res$FemRetProbAtAge, "l", col="red",lty="dotted", cex=0.8)
    legend('bottomright', col=c("red","red"),legend=c("Fem. mature","Fem. reten."),
           lty=c("solid","dotted"),bty='n', cex=0.8,lwd=1.75)
  } else {
    lines(Res$Ages,Res$FemSelLandAtAge, "l", col="red",lty="dotted",cex=0.8)
    lines(Res$Ages,Res$FemSelDiscAtAge, "l", col="brown",lty="dotted",cex=0.8)
    legend('bottomright', col=c("red","red","brown"),legend=c("Fem. mature","Fem. land.","Fem. disc."),
           lty=c("solid","dotted","dotted"),bty='n', cex=0.8,lwd=1.75)
  }
  axis(1,at=seq(0,MaxModelAge,1), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,1,0.5), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,MaxModelAge,1), labels=seq(0,MaxModelAge,1),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,1,0.5), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Proportion"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  if (ReprodPattern > 1) {
    lines(Res$Ages,Res$PropFemAtAge, "l", col="black",lty="solid",cex=0.8)
    legend('topright', col="black",legend="Prop. Fem.",
           lty="solid",bty='n', cex=0.8,lwd=1.75)
  }

  # plot male maturity and selectivity at age
  plot(Res$Ages, Res$MalPropMatAtAge,"l", pch=16, frame.plot=F,ylim=c(0,1),xlim=c(0,MaxModelAge),col="blue",yaxt="n",xaxt="n",
       ylab="",xlab="", cex=0.8)
  if (DiscMort == 0) {
    lines(Res$Ages, Res$MalRetProbAtAge, "l", col="blue",lty="dotted", cex=0.8)
    legend('bottomright', col=c("blue","blue"),legend=c("Male mature","Male reten."),
           lty=c("solid","dotted"),bty='n', cex=0.8,lwd=1.75)
  } else {
    lines(Res$Ages, Res$FemSelLandAtAge, "l", col="blue",lty="dotted", cex=0.8)
    lines(Res$Ages, Res$FemSelDiscAtAge, "l", col="purple",lty="dotted", cex=0.8)
    legend('bottomright', col=c("blue","blue","purple"),legend=c("Male mature","Male land.", "Male disc."),
           lty=c("solid","dotted","dotted"),bty='n', cex=0.8,lwd=1.75)
  }
  axis(1,at=seq(0,MaxModelAge,1), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,1,0.5), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,MaxModelAge,1), labels=seq(0,MaxModelAge,1),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,1,0.5), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Proportion"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  if (ReprodPattern > 1) {
    lines(Res$Ages,1-Res$PropFemAtAge, "l", col="black",lty="solid",cex=0.8)
    legend('topright', col="black",legend="Prop. Male",
           lty="solid",bty='n', cex=0.8,lwd=1.75)
  }

  #plot 2:
  par(mfrow = c(3,2), mar=c(3.5,4,2,2),
      oma=c(1,1,1,1),tck=-0.03,mgp = c(3, 0.5, 0))

  # plot female mortality at age
  ylims = Get_yaxis_scale(Res$FemFAtAge)
  ymax = ylims$ymax; yint = ylims$yint
  plot(Res$Ages, Res$FemFAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,MaxModelAge),col="red",yaxt="n",xaxt="n",
       ylab="",xlab="")
  # lines(Res$Ages,Res$FemZAtAge,lty="dotted",col="red")
  lines(Res$Ages,Res$FemDiscFAtAge,lty="dotted",col="brown")
  lines(Res$Ages,Res$FemLandFAtAge,lty="dotted",col="purple")
  lines(Res$Ages,rep(NatMort,length(Res$Ages)),lty="dashed")
  axis(1,at=seq(0,MaxModelAge,1), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,MaxModelAge,1), labels=seq(0,MaxModelAge,1),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Mortality") ~ (year^{-1}))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col=c("red","brown","purple","black"),lty=c("solid","dotted","dotted","dashed"),
         legend=c("Fem. F","Fem. DiscF","Fem. LandF","Fem. M"),bty='n', cex=1.0,lwd=1.75)

  # plot male mortality at age
  ylims = Get_yaxis_scale(Res$MalFAtAge)
  ymax = ylims$ymax; yint = ylims$yint
  plot(Res$Ages, Res$MalFAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,MaxModelAge),col="blue",yaxt="n",xaxt="n",
       ylab="",xlab="")
  # lines(Res$Ages, Res$MalZAtAge,lty="dotted",col="blue")
  lines(Res$Ages,Res$MalDiscFAtAge,lty="dotted",col="brown")
  lines(Res$Ages,Res$MalLandFAtAge,lty="dotted",col="purple")
  lines(Res$Ages,rep(NatMort,length(Res$Ages)),lty="dashed")
  axis(1,at=seq(0,MaxModelAge,1), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,MaxModelAge,1), labels=seq(0,MaxModelAge,1),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Mortality") ~ (year^{-1}))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col=c("blue","brown","purple","black"),lty=c("solid","dotted","dotted","dashed"),
         legend=c("Mal. F","Mal. DiscF","Mal.LandF","Mal. M"),bty='n', cex=1.0,lwd=1.75)

  # plot fished and unfished female survival in terms of numbers given specified current fully-selected fishing mortality
  plot(Res$Ages, Res$UnfishFemSurvAtAge,"l",frame.plot=F,ylim=c(0,InitRatioFem),xlim=c(0,MaxModelAge),col="red",yaxt="n",xaxt="n",
       ylab="",xlab="")
  lines(Res$Ages, Res$FishedFemSurvAtAge,col="red",lty="dotted")
  axis(1,at=seq(0,MaxModelAge,1), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,InitRatioFem,0.1), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,MaxModelAge,1), labels=seq(0,MaxModelAge,1),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,InitRatioFem,0.1), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Rel. survival"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col="red",lty=c("solid","dotted"),
         legend=c("Fem. unfish","Fem. fish"),bty='n', cex=1.0,lwd=1.75)

  # plot fished and unfished male survival in terms of numbers given specified current fully-selected fishing mortality
  plot(Res$Ages, Res$UnfishMalSurvAtAge,"l",frame.plot=F,ylim=c(0,InitRatioFem),xlim=c(0,MaxModelAge),col="blue",yaxt="n",xaxt="n",
       ylab="",xlab="")
  lines(Res$Ages, Res$FishedMalSurvAtAge,col="blue",lty="dotted")
  axis(1,at=seq(0,MaxModelAge,1), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,InitRatioFem,0.1), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,MaxModelAge,1), labels=seq(0,MaxModelAge,1),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,InitRatioFem,0.1), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Rel. survival"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col="blue",lty=c("solid","dotted"),
         legend=c("Mal. unfish","Mal. fish"),bty='n', cex=1.0,lwd=1.75)

  # plot fished and unfished mature female biomass at age given specified current fully-selected fishing mortality
  ylims = Get_yaxis_scale(Res$UnfishFemBiomAtAge)
  ymax = ylims$ymax; yint = ylims$yint
  plot(Res$Ages, Res$UnfishFemBiomAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,MaxModelAge),
       col="red",yaxt="n",xaxt="n",ylab="",xlab="")
  lines(Res$Ages, Res$FishedFemBiomAtAge,col="red",lty="dotted")
  axis(1,at=seq(0,MaxModelAge,1), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,MaxModelAge,1), labels=seq(0,MaxModelAge,1),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Biom. at age (kg"),plain(")"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col=c("red","red"),lty=c("solid","dotted"),legend=c("Fem. unfish","Mal. unfish"),bty='n', cex=1.0,lwd=1.75)

  # plot fished and unfished mature male biomass at age given specified current fully-selected fishing mortality
  ylims = Get_yaxis_scale(Res$UnfishMalBiomAtAge)
  ymax = ylims$ymax; yint = ylims$yint
  plot(Res$Ages, Res$UnfishMalBiomAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,MaxModelAge),
       col="blue",yaxt="n",xaxt="n",ylab="",xlab="")
  lines(Res$Ages, Res$FishedMalBiomAtAge,col="blue",lty="dotted")
  axis(1,at=seq(0,MaxModelAge,1), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,MaxModelAge,1), labels=seq(0,MaxModelAge,1),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Biom. at age (kg"),plain(")"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col=c("blue","blue"),lty=c("solid","dotted"),legend=c("Fem. unfish","Mal. unfish"),bty='n', cex=1.0,lwd=1.75)

  #Plot 3:

  # plot female and male catch at age, given specified current fully-selected fishing mortality
  par(mfrow = c(3,2), mar=c(3.5,4,2,2),
      oma=c(1,1,1,1),tck=-0.03,mgp = c(3, 0.5, 0))

  FemCatchNumAtAgeProp <- Res$FemCatchAtAgeNum / sum(Res$FemCatchAtAgeNum)
  MalCatchNumAtAgeProp <- Res$MalCatchAtAgeNum / sum(Res$MalCatchAtAgeNum)
  ylims = Get_yaxis_scale(FemCatchNumAtAgeProp)
  ymax = ylims$ymax; yint = ylims$yint
  plot(Res$Ages, FemCatchNumAtAgeProp,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,MaxModelAge),
       col="red",yaxt="n",xaxt="n",ylab="",xlab="")
  lines(Res$Ages, MalCatchNumAtAgeProp,col="blue","l")
  axis(1,at=seq(0,MaxModelAge,1), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,MaxModelAge,1), labels=seq(0,MaxModelAge,1),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext("Catch at age Prop.",las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col=c("red","blue"),lty="solid",legend=c("females","males"),bty='n', cex=0.8,lwd=1.75)


  # plot yield per recruit (per recruit analysis) and equilibrium catch (equilibrium age-structured model)
  # given specified current fully-selected fishing mortality
  ylims = Get_yaxis_scale(Res$YPRResults)
  ymax = ylims$ymax; yint = ylims$yint

  plot(Res$FishMort, Res$YPRResults,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,max(Res$FishMort)),
       col="black",yaxt="n",xaxt="n",ylab="",xlab="")
  points(Current_F, Res$YPR,cex=1.2,col="black",pch=16)
  lines(Res$FishMort, Res$EquilCatchResults,col="blue")
  points(Current_F, Res$Equil_Catch, cex=1.2,col="blue",pch=16)
  axis(1,at=seq(0,max(Res$FishMort),0.5), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,max(Res$FishMort),0.5), labels=seq(0,max(Res$FishMort),0.5),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("YPR / Eq.Catch (kg"),plain(")"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col=c("black","blue"),lty=c("solid","solid"),
         legend=c("YPR","Eq.catch"),bty='n', cex=0.8,lwd=1.75)

  # plot spawning potential ratio (per recruit analysis) and relative spawning biomass (equilbrium age-structured model)
  # given specified current fully-selected fishing mortality, for female
  plot(Res$FishMort, Res$Fem_SPRResults,"l",frame.plot=F,ylim=c(0,1.0),xlim=c(0,max(Res$FishMort)),
       col="red",yaxt="n",xaxt="n",ylab="",xlab="", lty="dotted")
  lines(Res$FishMort, Res$Equilmod_FemRelBiomResults,col="red",lty="dotted")
  points(Current_F, Res$Fem_SPR,cex=1.2,col="red",pch=16)
  points(Current_F, Res$Equilmod_FemRelBiom,cex=1.2,col="red",pch=1)
  axis(1,at=seq(0,max(Res$FishMort),0.5), cex.axis=0.8, lwd=1.75,lab=F) # x axis
  axis(2,at=seq(0,1,0.2), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,max(Res$FishMort),0.5), labels=seq(0,max(Res$FishMort),0.5),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add x labels
  axis(2,at=seq(0,1,0.2), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Biom. ratio"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col=c("red","red"),lty=c("dotted","solid"),
         legend=c("Fem SPR","Fem. Rel biom"),bty='n', cex=0.8,lwd=1.75)

  # plot spawning potential ratio (per recruit analysis) and relative spawning biomass (equilbrium age-structured model)
  # given specified current fully-selected fishing mortality, for each sex
  plot(Res$FishMort, Res$Mal_SPRResults,"l",frame.plot=F,ylim=c(0,1.0),xlim=c(0,max(Res$FishMort)),
       col="blue",yaxt="n",xaxt="n",ylab="",xlab="", lty="dotted")
  lines(Res$FishMort, Res$Equilmod_MalRelBiomResults,col="blue",lty="solid")
  points(Current_F, Res$Mal_SPR,cex=1.2,col="blue",pch=16)
  points(Current_F, Res$Equilmod_MalRelBiom,cex=1.2,col="blue",pch=1)
  axis(1,at=seq(0,max(Res$FishMort),0.5), cex.axis=0.8, lwd=1.75,lab=F) # x axis
  axis(2,at=seq(0,1,0.2), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,max(Res$FishMort),0.5), labels=seq(0,max(Res$FishMort),0.5),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add x labels
  axis(2,at=seq(0,1,0.2), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Biom. ratio"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col=c("blue","blue"),lty=c("dotted","solid"),
         legend=c("Mal SPR","Mal. Rel biom"),bty='n', cex=0.8,lwd=1.75)

  # plot spawning potential ratio (per recruit analysis) and relative spawning biomass (equilbrium age-structured model)
  # given specified current fully-selected fishing mortality, for combined sexes
  plot(Res$FishMort, Res$CombSex_SPRResults,"l",frame.plot=F,ylim=c(0,1.0),xlim=c(0,max(Res$FishMort)),
       col="black",yaxt="n",xaxt="n",ylab="",xlab="", lty="dotted")
  lines(Res$FishMort, Res$Equilmod_CombSexRelBiomResults,col="black",lty="solid")
  points(Current_F, Res$CombSex_SPR,cex=1.2,col="black",pch=16)
  points(Current_F, Res$Equilmod_CombSexRelBiom,cex=1.2,col="black",pch=1)
  axis(1,at=seq(0,max(Res$FishMort),0.5), cex.axis=0.8, lwd=1.75,lab=F) # x axis
  axis(2,at=seq(0,1,0.2), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,max(Res$FishMort),0.5), labels=seq(0,max(Res$FishMort),0.5),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add x labels
  axis(2,at=seq(0,1,0.2), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Biom. ratio"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col=c("black","black"),lty=c("dotted","solid"),
         legend=c("CombSex SPR","CombSex Rel biom"),bty='n', cex=0.8,lwd=1.75)

  # plot equilibrium recruitment vs F
  # ylims = Get_yaxis_scale(Res$EquilRecResults)
  # ymax = ylims$ymax; yint = ylims$yint
  ymax = 1.0
  yint = 0.2

  plot(Res$FishMort, Res$EquilRecResults,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,max(Res$FishMort)),
       col="red",yaxt="n",xaxt="n",ylab="",xlab="")
  points(Current_F, Res$Equil_Rec, cex=1.2,col="red",pch=16)
  axis(1,at=seq(0,max(Res$FishMort),0.5), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,max(Res$FishMort),0.5), labels=seq(0,max(Res$FishMort),0.5),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Equil. Recruitment"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)

}


#' Deterministic plot of yield vs F from per recruit analysis and extended analysis
#' with a stock-recruitment relationship
#'
#' This function provides a deterministic plot of yield vs F from per recruit analysis and extended analysis
#' with a stock-recruitment relationship
#'
#' Function requires values for the initial sex ratio (i.e. of juvenile recruits),
#' maximum age, model timestep (e.g. 1 = annual, 1/12 = monthly), von Bertalanffy growth
#' parameters, weight-length parameters, i.e using a power relationship, logistic length at
#' maturity parameters, logistic selectivity parameters, steepness parameter of the stock-recruitment
#' relationship, natural mortality (NatMort), current fishing mortality (Current_F)
#'
#' (Numbers and a Numeric Vectors)
#' @param maximum age considered in model (MaxModelAge), model timestep, e.g. 1 = annual, 1/12 = monthly (TimeStep),
#' von Bertalanffy growth parameters (Linf, vbK, tzero), Female and male estimated lengths-at-age (set to NA if
#' using von Bertalanffy growth equation, or set the von Bertalanffy growth parameters to NA and input as data frame
#' (columns for sex) from any growth model (EstLenAtAge), weight-length parameters, using a power or log-log relationship
#' (lenwt_a, ln_lenwt_a, lenwt_b), weight-length relationship type (WLrel_Type), Female and male estimated
#' weights-at-ages  (set to NA if weight-length relationship parameters, or set the weight-length parameters
#' to NA and input these as data frame (columns for sex), (EstWtAtAge), reproductive pattern, 1=gonochoristic,
#' 2=protogynous, 3=protandrous (ReprodPattern), initial sex ratio of juvenile recruits (InitRatioFem),
#' Maximum proportion of final sex for hermaphroditic species (set to NA for gonochoristic species) (FinalSex_Pmax),
#' logistic sex change parameters for hermaphroditic species (set to NA for gonochoristic species) (FinalSex_A50, FinalSex_A95),
#' logistic length at maturity parameters (mat_L50, mat_L95), (set to NA if using age at maturity parameters, or set the
#' age at maturity parameters to NA and input these as data fra,e (columns for sex) (EstMatAtAge), logistic gear selectivity
#' parameters (sel_A50, sel_A95),logistic retention curve parameters (ret_Pmax, ret_A50, ret_A95), discard mortality rate
#' (DiscMort), steepness parameter of the stock-recruitment relationship (Steepness), stock-recruitment relationship type,
#' 1 = Beverton-Holt, 2=Ricker (SRrel_Type), natural mortality (NatMort), option for plotting reference points on final plot,
#'1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points (RefPointPlotOpt), current fishing mortality (Current_F).
#' @return
#' deterministic plot of yield vs F from per recruit analysis and extended analysis
#' with a stock-recruitment relationship.
#' @examples
#' # Example 1. Non-hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 20 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be approprriate to use a smaller time step)
#' Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.5, 0.5) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstMalLenAtAge=NA,
#'                           EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- NA # for log-log relationship
#' lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA,
#'                          EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#'
#' mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA,
#'                           EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- c(2.5, 2.5) # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- c(3.5, 3.5) # females, males - Logistic age fish retention at age parameters
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 0.2 # natural mortality  (year-1)
#' Current_F <- 0.07 # estimate of fishing mortality, e.g. from catch curve analysis
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' PlotPerRecruit_Yield_no_err(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                       lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
#'                       ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
#'                       mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
#'                       DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F)
#' # Example 2: hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 100 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be approprriate to use a smaller time step)
#' Linf <- c(1000, 1000) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.1, 0.1) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstMalLenAtAge=NA,
#'                           EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- -11.0 # for log-log relationship
#' lenwt_b <- 3.0 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 2 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA,
#'                          EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodPattern <- 2 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 1 # Ratio of females to males at age zero
#' FinalSex_Pmax <- 1 # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_A50 <- 35 # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- 60 # Logistic sex change relationship parameters (95% of max probability)
#' mat_A50 <- c(20, 20) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(30, 30) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA,                          EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_A50 <- c(15, 15) # females, males - Logistic age selectivity relationship parameters
#' sel_A95 <- c(25, 25) # females, males - Logistic age selectivity relationship parameters
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- c(15, 15) # females, males - Logistic age fish retention at age parameters (inflection point)
#' ret_A95 <- c(25, 25) # females, males - Logistic age fish retention at age parameters (95% of maximum retention)
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 0.07 # natural mortality  (year-1)
#' Current_F <- 0.07 # estimate of fishing mortality, e.g. from catch curve analysis
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' PlotPerRecruit_Yield_no_err(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                       lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
#'                       ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
#'                       mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
#'                       DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F)
#' @export
PlotPerRecruit_Yield_no_err <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                        lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
                                        ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                                        mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
                                        DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F) {

  Res = GetPerRecruitResults(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
                             ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                             mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
                             DiscMort, Steepness, SRrel_Type, NatMort, Current_F)

  # F vs YPR and equil catch
  ymax = max(Res$YPRResults)*1.5; yint = 0.2
  # plot(Res$FishMort, Res$YPRResults, "l", frame.plot=F, ylim=c(0, ymax), xlim=c(0, max(Res$FishMort)),
  #      col="black", yaxt="n", xaxt="n", ylab="", xlab="")
  # points(Current_F, Res$YPR, cex=1.2, col="black", pch=16)
  # lines(Res$FishMort, Res$EquilCatchResults, col="black", lty="dotted")
  plot(Res$FishMort, Res$EquilCatchResults, "l", frame.plot=F, ylim=c(0, ymax), xlim=c(0, max(Res$FishMort)),
       col="black", yaxt="n", xaxt="n", ylab="", xlab="")
  points(Current_F, Res$Equil_Catch, cex=1.2, col="black", pch=16)
  axis(1, at=seq(0, max(Res$FishMort), 0.5), cex.axis=1, lwd=1.75, lab=F)
  axis(2, at=seq(0, ymax, yint), cex.axis=1, lwd=1.75, lab=F)
  axis(1, at=seq(0, max(Res$FishMort), 0.5), labels=seq(0, max(Res$FishMort), 0.5),
       cex.axis=1, line=0.5, las=1, lwd=1.5, tick=F)
  axis(2, at=seq(0, ymax, yint), cex.axis=1, line=0.5, las=1, lwd=1.5, tick=F)
  mtext(expression(paste(plain("Equilibrium catch (kg)"))), las=3, side=2, line=3, cex=1, lwd=1.75)
  mtext(expression(paste(italic("F") ~ (year^{-1}))), las=1, side=1, line=3, cex=1, lwd=1.75)
  # legend("topright", col=c("black", "black"), lty=c("solid", "dotted"), legend=c("YPR", "Eq. catch"),
         # bty="n", cex=0.8, lwd=1.75, inset = 0.05)
  legend("topleft", col="black", pch = 16, lty=0, legend="Estimate",
         bty="n", cex=0.8, inset = 0.05)
}

#' Deterministic plot of SPR and relative equilibrium biomass vs F from per recruit analysis and extended analysis
#' with a stock-recruitment relationship
#'
#' This function provides a deterministic plot of SPR and relative equilibrium biomass  vs F from per recruit analysis and extended analysis
#' with a stock-recruitment relationship
#'
#' Function requires values for the initial sex ratio (i.e. of juvenile recruits),
#' maximum age, model timestep (e.g. 1 = annual, 1/12 = monthly), von Bertalanffy growth
#' parameters, weight-length parameters, i.e using a power relationship, logistic length at
#' maturity parameters, logistic selectivity parameters, steepness parameter of the stock-recruitment
#' relationship, natural mortality (NatMort), current fishing mortality (Current_F)
#'
#' @param maximum age considered in model (MaxModelAge), model timestep, e.g. 1 = annual, 1/12 = monthly (TimeStep),
#' von Bertalanffy growth parameters (Linf, vbK, tzero), Female and male estimated lengths-at-age (set to NA if
#' using von Bertalanffy growth equation, or set the von Bertalanffy growth parameters to NA and input as data frame
#' (columns for sex) from any growth model (EstLenAtAge), weight-length parameters, using a power or log-log relationship
#' (lenwt_a, ln_lenwt_a, lenwt_b), weight-length relationship type (WLrel_Type), Female and male estimated
#' weights-at-ages  (set to NA if weight-length relationship parameters, or set the weight-length parameters
#' to NA and input these as data frame (columns for sex), (EstWtAtAge), reproductive pattern, 1=gonochoristic,
#' 2=protogynous, 3=protandrous (ReprodPattern), initial sex ratio of juvenile recruits (InitRatioFem),
#' Maximum proportion of final sex for hermaphroditic species (set to NA for gonochoristic species) (FinalSex_Pmax),
#' logistic sex change parameters for hermaphroditic species (set to NA for gonochoristic species) (FinalSex_A50, FinalSex_A95),
#' logistic length at maturity parameters (mat_L50, mat_L95), (set to NA if using age at maturity parameters, or set the
#' age at maturity parameters to NA and input these as data fra,e (columns for sex) (EstMatAtAge), logistic gear selectivity
#' parameters (sel_A50, sel_A95),logistic retention curve parameters (ret_Pmax, ret_A50, ret_A95), discard mortality rate
#' (DiscMort), steepness parameter of the stock-recruitment relationship (Steepness), stock-recruitment relationship type,
#' 1 = Beverton-Holt, 2=Ricker (SRrel_Type), natural mortality (NatMort), option for plotting reference points on final plot,
#'1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points (RefPointPlotOpt), current fishing mortality (Current_F).
#' @return
#' deterministic plot of SPR and relative equilibrium biomass  vs F from per recruit analysis and extended analysis
#' with a stock-recruitment relationship.
#' @examples
#' # Example 1. Non-hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 20 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be approprriate to use a smaller time step)
#' Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.5, 0.5) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstMalLenAtAge=NA,
#'                           EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- NA # for log-log relationship
#' lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA,
#'                          EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#'
#' mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA,
#'                           EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- c(2.5, 2.5) # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- c(3.5, 3.5) # females, males - Logistic age fish retention at age parameters
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 0.2 # natural mortality  (year-1)
#' Current_F <- 0.2 # estimate of fishing mortality, e.g. from catch curve analysis
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' PlotPerRecruit_Biom_no_err(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                       lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
#'                       ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
#'                       mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
#'                       DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F)
#' # Example 2: hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 100 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be approprriate to use a smaller time step)
#' Linf <- c(1000, 1000) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.1, 0.1) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstMalLenAtAge=NA,
#'                           EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- -11.0 # for log-log relationship
#' lenwt_b <- 3.0 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 2 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA,
#'                          EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodPattern <- 2 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 1 # Ratio of females to males at age zero
#' FinalSex_Pmax <- 1 # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_A50 <- 35 # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- 60 # Logistic sex change relationship parameters (95% of max probability)
#' mat_A50 <- c(20, 20) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(30, 30) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA,                          EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_A50 <- c(15, 15) # females, males - Logistic age selectivity relationship parameters
#' sel_A95 <- c(25, 25) # females, males - Logistic age selectivity relationship parameters
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- c(15, 15) # females, males - Logistic age fish retention at age parameters (inflection point)
#' ret_A95 <- c(25, 25) # females, males - Logistic age fish retention at age parameters (95% of maximum retention)
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 0.07 # natural mortality  (year-1)
#' Current_F <- 0.07 # estimate of fishing mortality, e.g. from catch curve analysis
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' PlotPerRecruit_Biom_no_err(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                       lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
#'                       ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
#'                       mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
#'                       DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F)
#' @export
PlotPerRecruit_Biom_no_err <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                       lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
                                       ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                                       mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
                                       DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F) {

  Res = GetPerRecruitResults(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
                             ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                             mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
                             DiscMort, Steepness, SRrel_Type, NatMort, Current_F)

  # F vs SPR and Brel
  # plot(Res$FishMort, Res$Fem_SPRResults, "l", frame.plot=F, ylim=c(0, 1), xlim=c(0, max(Res$FishMort)),
  #      col="black", yaxt="n", xaxt="n", ylab="", xlab="")
  # lines(Res$FishMort, Res$Equilmod_FemRelBiomResults, col="black", lty="dotted")
  # points(Current_F, Res$Fem_SPR, cex=1.2, col="black", pch=16)

  plot(Res$FishMort, Res$Equilmod_FemRelBiomResults, "l", frame.plot=F, ylim=c(0, 1), xlim=c(0, max(Res$FishMort)),
       col="black", yaxt="n", xaxt="n", ylab="", xlab="")
  points(Current_F, Res$Equilmod_FemRelBiom, cex=1.2, col="black", pch=16)
  axis(1, at=seq(0, max(Res$FishMort), 0.5), cex.axis=1, lwd=1.75, lab=F)
  axis(2, at=seq(0, 1, 0.2), cex.axis=1, lwd=1.75, lab=F)
  axis(1, at=seq(0, max(Res$FishMort), 0.5), labels = seq(0, max(Res$FishMort), 0.5),
       cex.axis=1, line=0.5, las=1, lwd=1.5, tick=F)
  axis(2, at=seq(0, 1, 0.2), cex.axis=1, line=0.5, las=1, lwd=1.5, tick=F)
  mtext(expression(paste(plain("Relative spawning biomass"))), las=3, side=2, line=3, cex=1, lwd=1.75)
  mtext(expression(paste(italic("F") ~ (year^{-1}))), las=1, side=1, line=3, cex=1, lwd=1.75)
  # legend("topright", col=c("black", "black"), lty=c("solid", "dotted"), legend=c("SPR", "Rel. biomass"),
  #        bty="n", cex=0.8, lwd=1.75, inset = 0.05)
  legend("topleft", col="black", pch = 16, lty=0, legend="Estimate",
  bty="n", cex=0.8, inset = 0.05)

  # RefPointPlotOpt <- 0 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
  if (RefPointPlotOpt == 1) {
    lines(abline(h=0.4, col="green"))
    lines(abline(h=0.3, col="orange"))
    lines(abline(h=0.2, col="red"))
    legend('topright', col=c("green","orange","red"),lty=c("solid","solid","solid"),
           legend=c("0.4","0.3","0.2"),bty='n', cex=0.8, lwd=1.75)
  }
  if (RefPointPlotOpt == 2) {
    lines(abline(h=Res$BMSY_Targ, col="green"))
    lines(abline(h=Res$BMSY_Thresh, col="orange"))
    lines(abline(h=Res$BMSY_Lim, col="red"))
    legend('topright', col=c("green","orange","red"),lty=c("solid","solid","solid"),
           legend=c("1.2BMSY","BMSY","0.5BMSY"),bty='n', cex=0.8,lwd=1.75)
  }
}


#' Get outputs from per recruit analysis across a range of fishing mortality values, allowing for error in
#' M, h and F
#'
#' This function provides outputs associated with per recruit analysis, and an
#' extended form of this analysis with a Beverton-Holt stock recruitment relationship to account
#' for potential impacts of fishing on recruitment. Outputs are provided for a range of
#' fishing mortality values, including the current, estimated value, with specified error in M, h and F.
#'
#' Function requires values for the initial sex ratio (i.e. of juvenile recruits),
#' maximum age, model timestep (e.g. 1 = annual, 1/12 = monthly), von Bertalanffy growth
#' parameters, weight-length parameters, i.e using a power relationship, logistic length at
#' maturity parameters, logistic selectivity parameters, steepness parameter of the stock-recruitment
#' relationship with associated error, natural mortality with associated error (NatMort), fishing mortality
#' with associated error (FMort.
#' (Numbers and a Numeric Vectors)
#' @param params maximum age considered in model (MaxModelAge), model timestep, e.g. 1 = annual, 1/12 = monthly (TimeStep),
#' von Bertalanffy growth parameters (Linf, vbK, tzero), Female and male estimated lengths-at-age (set to NA if
#' using von Bertalanffy growth equation, or set the von Bertalanffy growth parameters to NA and input as data frame
#' (columns for sex) from any growth model (EstLenAtAge), weight-length parameters, using a power or log-log relationship
#' (lenwt_a, ln_lenwt_a, lenwt_b), weight-length relationship type (WLrel_Type), Female and male estimated
#' weights-at-ages  (set to NA if weight-length relationship parameters, or set the weight-length parameters
#' to NA and input these as data frame (columns for sex), (EstWtAtAge), reproductive pattern, 1=gonochoristic,
#' 2=protogynous, 3=protandrous (ReprodPattern), initial sex ratio of juvenile recruits (InitRatioFem),
#' Maximum proportion of final sex for hermaphroditic species (set to NA for gonochoristic species) (FinalSex_Pmax),
#' logistic sex change parameters for hermaphroditic species (set to NA for gonochoristic species) (FinalSex_A50, FinalSex_A95),
#' logistic length at maturity parameters (mat_L50, mat_L95), (set to NA if using age at maturity parameters, or set the
#' age at maturity parameters to NA and input these as data fra,e (columns for sex) (EstMatAtAge), logistic gear selectivity
#' parameters (sel_A50, sel_A95),logistic retention curve parameters (ret_Pmax, ret_A50, ret_A95), discard mortality rate
#' (DiscMort), steepness parameter of the stock-recruitment relationship (Steepness) and associated standard deviation
#' (Steepness_sd), stock-recruitment relationship type, 1 = Beverton-Holt, 2=Ricker (SRrel_Type), natural mortality
#' (NatMort) and associated standard deviation (NatMort_sd), current fishing mortality (Current_F) and associated
#' standard deviation (Current_F_sd).
#' @return range of fishing mortality values for calculating per recruit quantities (PerRec_FValues)
#' estimated female SPR and equilibrium relative biomass values for trials with random M, h and F,
#' (Fem_SPR_Vals, Equil_RelFemSpBiom_Vals), estimated relative biomass at maximum sustainable yield
#' by trial (BMSY_Vals), estimated female spawning potential ratio values by F for each trial (Sim_FemSPR)
#' estimated female relative biomass values by F for each trial (Sim_Equil_RelFemSpBiom), median,
#' lower 95 and upper 95 percentile values for SPR for current F (EstFemSPR, EstLow95FemSPR, EstUp95FemSPR), median,
#' lower 95 and upper 95 percentile values for equilibrium relative female spawning biomass, median
#' BMSY ratio (EstBMSYratio), summary matrix containing key outputs from analysis (ResSummary_with_err)
#' @examples
#' # Example. Non-hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 20 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be approprriate to use a smaller time step)
#' Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.5, 0.5) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstFemLenAtAge=NA,
#'                           EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- NA # for log-log relationship
#' lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA,
#'                          EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#'
#' mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA,
#'                           EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- c(2.5, 2.5) # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- c(3.5, 3.5) # females, males - Logistic age fish retention at age parameters
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' Steepness_sd <- 0.025
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort <- 0.2 # natural mortality  (year-1)
#' NatMort_sd <- 0.025
#' Current_F <- 0.1 # estimate of fishing mortality, e.g. from catch curve analysis
#' Current_F_sd <- 0.005
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' nReps = 50
#' GetPerRecruitResults_with_err(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                               lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodPattern,
#'                               InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95, mat_A50, mat_A95,
#'                               EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95, DiscMort,
#'                               Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd, Current_F,
#'                               Current_F_sd, nReps)
#' @export
GetPerRecruitResults_with_err <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                          lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodPattern,
                                          InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95, mat_A50, mat_A95,
                                          EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95, DiscMort,
                                          Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd, Current_F,
                                          Current_F_sd, nReps) {



  FValues = rnorm(nReps, Current_F, Current_F_sd)
  hValues = rnorm(nReps, Steepness, Steepness_sd)
  MValues = rnorm(nReps, NatMort, NatMort_sd)


  for (i in 1:nReps) {
    FMort = FValues[i]
    Steepness = hValues[i]
    NatMort = MValues[i]
    PREst = GetPerRecruitResults(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                 lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
                                 ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                                 mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
                                 DiscMort, Steepness, SRrel_Type, NatMort, FMort)

    if (i == 1) {
      Fem_SPR_Vals = rep(0, nReps)
      Equil_RelFemSpBiom_Vals = rep(0, nReps)
      BMSY_Vals = rep(0,nReps)
      Sim_FemSPR = data.frame(matrix(nrow=nReps, ncol=length(PREst$FishMort)))
      colnames(Sim_FemSPR) = PREst$FishMort
      Sim_Equil_RelFemSpBiom = Sim_FemSPR
    }

    Fem_SPR_Vals[i] = PREst$Fem_SPR
    Equil_RelFemSpBiom_Vals[i] = PREst$Equilmod_FemRelBiom
    BMSY_Vals[i] = PREst$BMSY_Thresh
    Sim_FemSPR[i,] = PREst$Fem_SPRResults
    Sim_Equil_RelFemSpBiom[i,] = PREst$Equilmod_FemRelBiomResults

    cat("i",i,'\n')
  }

  # Save key outputs (median and 95% confidence levels)
  EstFemSPR = round(median(Fem_SPR_Vals),3)
  EstLow95FemSPR = as.numeric(round(quantile(Fem_SPR_Vals, 0.025),3))
  EstUp95FemSPR = as.numeric(round(quantile(Fem_SPR_Vals, 0.975),3))
  EstEquilRelFemSpBiom = round(median(Equil_RelFemSpBiom_Vals),3)
  Low95EquilRelFemSpBiom = as.numeric(round(quantile(Equil_RelFemSpBiom_Vals, 0.025),3))
  Upp95EquilRelFemSpBiom = as.numeric(round(quantile(Equil_RelFemSpBiom_Vals, 0.975),3))
  EstBMSYratio = round(median(BMSY_Vals),3)
  ResSummary_with_err <- data.frame(EstFemSPR, EstLow95FemSPR, EstUp95FemSPR, EstEquilRelFemSpBiom, Low95EquilRelFemSpBiom, Upp95EquilRelFemSpBiom, EstBMSYratio)
  colnames(ResSummary_with_err)=c("SPR", "LowSPR", "UppSPR", "EquilSB", "LowEquilSB", "UppEquilSB", "BMSYratio")

  Results = list(PerRec_FValues = PREst$FishMort,
                 Fem_SPR_Vals=Fem_SPR_Vals,
                 Equil_RelFemSpBiom_Vals=Equil_RelFemSpBiom_Vals,
                 BMSY_Vals=BMSY_Vals,
                 Sim_FemSPR=Sim_FemSPR,
                 Sim_Equil_RelFemSpBiom=Sim_Equil_RelFemSpBiom,
                 EstFemSPR=EstFemSPR,
                 EstLow95FemSPR=EstLow95FemSPR,
                 EstUp95FemSPR=EstUp95FemSPR,
                 EstEquilRelFemSpBiom=EstEquilRelFemSpBiom,
                 Low95EquilRelFemSpBiom=Low95EquilRelFemSpBiom,
                 Upp95EquilRelFemSpBiom=Upp95EquilRelFemSpBiom,
                 EstBMSYratio=EstBMSYratio,
                 ResSummary_with_err=ResSummary_with_err)

  return(Results)

}


#' Plot of SPR and relative equilibrium biomass vs F from per recruit analysis and extended analysis
#' with a stock-recruitment relationship, with error for M, h and F
#'
#' This function provides a plot of SPR and relative equilibrium biomass vs F from per recruit analysis
#' and extended analysis with a stock-recruitment relationship, with specified error for M, h and F
#'
#' Function requires values for the initial sex ratio (i.e. of juvenile recruits),
#' maximum age, model timestep (e.g. 1 = annual, 1/12 = monthly), von Bertalanffy growth
#' parameters, weight-length parameters, i.e using a power relationship, logistic length at
#' maturity parameters, logistic selectivity parameters, steepness parameter of the stock-recruitment
#' relationship, natural mortality (NatMort), current fishing mortality (Current_F)
#' @param params maximum age considered in model (MaxModelAge), model timestep, e.g. 1 = annual, 1/12 = monthly (TimeStep),
#' von Bertalanffy growth parameters (Linf, vbK, tzero), Female and male estimated lengths-at-age (set to NA if
#' using von Bertalanffy growth equation, or set the von Bertalanffy growth parameters to NA and input as data frame
#' (columns for sex) from any growth model (EstLenAtAge), weight-length parameters, using a power or log-log relationship
#' (lenwt_a, ln_lenwt_a, lenwt_b), weight-length relationship type (WLrel_Type), Female and male estimated
#' weights-at-ages  (set to NA if weight-length relationship parameters, or set the weight-length parameters
#' to NA and input these as data frame (columns for sex), (EstWtAtAge), reproductive pattern, 1=gonochoristic,
#' 2=protogynous, 3=protandrous (ReprodPattern), initial sex ratio of juvenile recruits (InitRatioFem),
#' Maximum proportion of final sex for hermaphroditic species (set to NA for gonochoristic species) (FinalSex_Pmax),
#' logistic sex change parameters for hermaphroditic species (set to NA for gonochoristic species) (FinalSex_A50, FinalSex_A95),
#' logistic length at maturity parameters (mat_L50, mat_L95), (set to NA if using age at maturity parameters, or set the
#' age at maturity parameters to NA and input these as data fra,e (columns for sex) (EstMatAtAge), logistic gear selectivity
#' parameters (sel_A50, sel_A95),logistic retention curve parameters (ret_Pmax, ret_A50, ret_A95), discard mortality rate
#' (DiscMort), steepness parameter of the stock-recruitment relationship (Steepness) and associated standard deviation
#' (Steepness_sd), stock-recruitment relationship type, 1 = Beverton-Holt, 2=Ricker (SRrel_Type), natural mortality
#' (NatMort) and associated standard deviation (NatMort_sd), current fishing mortality (Current_F) and associated
#' standard deviation (Current_F_sd).
#' @return plot of SPR and relative equilibrium biomass vs F from per recruit analysis
#' and extended analysis with a stock-recruitment relationship, with specified error for M, h and F
#' @examples
#' # Example. Non-hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 20 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be approprriate to use a smaller time step)
#' Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.5, 0.5) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstFemLenAtAge=NA,
#'                           EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- NA # for log-log relationship
#' lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA,
#'                          EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#'
#' mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA,
#'                           EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- c(2.5, 2.5) # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- c(3.5, 3.5) # females, males - Logistic age fish retention at age parameters
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' Steepness_sd <- 0.025
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort <- 0.2 # natural mortality  (year-1)
#' NatMort_sd <- 0.025
#' Current_F <- 0.1 # estimate of fishing mortality, e.g. from catch curve analysis
#' Current_F_sd <- 0.005
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' nReps = 50
#' PlotPerRecruit_Biom_with_err(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                              lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodPattern,
#'                              InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95, mat_A50, mat_A95,
#'                              EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95, DiscMort,
#'                              Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd, Current_F,
#'                              Current_F_sd, RefPointPlotOpt, nReps)
#' @export
PlotPerRecruit_Biom_with_err <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                         lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodPattern,
                                         InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95, mat_A50, mat_A95,
                                         EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95, DiscMort,
                                         Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd, Current_F,
                                         Current_F_sd, RefPointPlotOpt, nReps) {

  # get BMSY reference points
  res = GetPerRecruitResults(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge,
                             ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                             mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95,
                             DiscMort, Steepness, SRrel_Type, NatMort, Current_F)

  Res=GetPerRecruitResults_with_err(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                    lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodPattern,
                                    InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95, mat_A50, mat_A95,
                                    EstMatAtAge, sel_A50, sel_A95, ret_Pmax, ret_A50, ret_A95, DiscMort,
                                    Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd, Current_F,
                                    Current_F_sd, nReps)

  # Plot per recruit outputs with uncertainty
  EqB_med = apply(Res$Sim_Equil_RelFemSpBiom,2,quantile, probs=c(0.5))
  EqB_lw = apply(Res$Sim_Equil_RelFemSpBiom,2,quantile, probs=c(0.025))
  EqB_hi = apply(Res$Sim_Equil_RelFemSpBiom,2,quantile, probs=c(0.975))
  plot(Res$PerRec_FValues, EqB_med, "l", frame.plot=F, ylim=c(0,1.0), xlim=c(0,max(Res$PerRec_FValues)),
       col="black", yaxt="n", xaxt="n", ylab="", xlab="")
  if (RefPointPlotOpt == 1) {
    lines(abline(h = 0.4, col = "green"))
    lines(abline(h = 0.3, col = "orange"))
    lines(abline(h = 0.2, col = "red"))
    legend("topright", col = c("green", "orange", "red"), lty = c("solid", "solid", "solid"),
           legend = c("0.4", "0.3", "0.2"), bty = "n", cex = 0.8, lwd = 1.75)
  }
  if (RefPointPlotOpt == 2) {
    lines(abline(h = res$BMSY_Targ, col = "green"))
    lines(abline(h = res$BMSY_Thresh, col = "orange"))
    lines(abline(h = res$BMSY_Lim, col = "red"))
    legend("topright", col = c("green", "orange", "red"), lty = c("solid", "solid", "solid"),
           legend = c("1.2BMSY", "BMSY",  "0.5BMSY"), bty = "n", cex = 0.8, lwd = 1.75)
  }
  points(Current_F, Res$EstEquilRelFemSpBiom, cex=1.2, col="black", pch=16)
  arrows(Current_F, Res$Low95EquilRelFemSpBiom, Current_F, Res$Upp95EquilRelFemSpBiom, length=0.05, angle=90, code=3)
  polygon(c(res$FishMort,rev(res$FishMort)),c(EqB_lw,rev(EqB_hi)), col=grey(0.5,0.25),
          border=grey(0.5,0.25))
  axis(1, at=seq(0, max(res$FishMort), 0.5), cex.axis=1, lwd=1, lab=F)
  axis(2, at=seq(0, 1, 0.2), cex.axis=1, lwd=1, lab=F)
  axis(1, at=seq(0, max(res$FishMort), 0.5), labels = seq(0, max(res$FishMort), 0.5),
       cex.axis=1, line=0.5, las=1, lwd=1, tick=F)
  axis(2, at=seq(0, 1, 0.2), cex.axis=1, line=0.5, las=1, lwd=1, tick=F)
  mtext(expression(paste(plain("Relative spawning biomass"))), las=3, side=2, line=3, cex=1, lwd=1.75)
  mtext(expression(paste(italic("F") ~ (year^{-1}))), las=1, side=1, line=3, cex=1, lwd=1.75)
  legend("topleft", col="black", pch = 16, legend="Estimate",
         bty="n", cex=0.8, lty=0, inset = 0.05)

}


