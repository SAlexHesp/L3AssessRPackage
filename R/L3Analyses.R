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

# Alex Hesp, August 2022
# Department of Primary Industries and Regional Development
# Catch curve and per recruit analysis package

# **************************
# general plotting functions
# **************************

 #' Get maximum and minimum x values and x interval
 #'
 #' Used for setting plot defaults
 #'
 #' @keywords internal
 #'
 #' @param x_data x axis data for plot
 #'
 #' @return list object with minimum and maximum x axis values and x axis interval
 Get_xaxis_scale <- function(x_data) {

   # modified from code provided by
   # https://peltiertech.com/calculate-nice-axis-scales-in-excel-vba/

   xmax_data = max(x_data)
   xmin_data = min(x_data)

   xpow = log10(xmax_data-xmin_data)
   xint = 10 ^ (xpow-round(xpow,0))

   if (xint>=0 & xint<2.5) {
     xint = 0.2
   }
   if (xint>=2.5 & xint<5) {
     xint = 0.5
   }
   if (xint>=5 & xint<7.5) {
     xint = 1
   }
   if (xint>=7.5) {
     xint = 2
   }

   xint = xint * 10^round(xpow,0) # major ticks
   xmin = xint * round(xmin_data/xint,0)
   xmax = xint * (round(xmax_data / xint,0) + 1)

   results = list(xmin = xmin,
                  xmax = xmax,
                  xint = xint)

 }

 #' Get maximum and minimum y values and y interval
 #'
 #' Used for setting plot defaults
 #'
 #' @keywords internal
 #'
 #' @param y_data y axis data for plot
 #'
 #' @return list object with minimum and maximum y axis values and y axis interval
 Get_yaxis_scale <- function(y_data) {

   # modified from code provided by
   # https://peltiertech.com/calculate-nice-axis-scales-in-excel-vba/

   ymax_data = max(1.1 * y_data)
   ymin_data = min(y_data)

   ypow = log10(ymax_data-ymin_data)
   yint = 10 ^ (ypow-round(ypow,0))

   if (yint>=0 & yint<2.5) {
     yint = 0.2
   }
   if (yint>=2.5 & yint<5) {
     yint = 0.5
   }
   if (yint>=5 & yint<7.5) {
     yint = 1
   }
   if (yint>=7.5) {
     yint = 2
   }

   yint = yint * 10^round(ypow,0) # major ticks
   ymin = yint * round(ymin_data/yint,0)
   ymax = yint * (round(ymax_data / yint,0) + 1)

   results = list(ymin = ymin,
                  ymax = ymax,
                  yint = yint)

 }


# **************************
# growth functions
# **************************

 #' Schnute growth function
 #'
 #' Calculates estimated length at age from Schnute growth curve, given the age, growth parameters and
 #' reference ages
 #'
 #' @keywords internal
 #'
 #' @param Age specified age
 #' @param t1 first reference age
 #' @param t2 second reference age
 #' @param y1 length at t1
 #' @param y2 length at t2
 #' @param a growth curve parameter
 #' @param a growth curve parameter
 #'
 #' @return Estimated length at a specified age
 SchnuteGrowthfunction <- function (Age, t1, t2, y1, y2, a, b) {
   # Schnute's versatile growth model

   # Reference:
   #    Schnute, J. T. and Richards, L. J.  (1990).  A unified aproach to the
   #    analysis of fish growth, maturity and survivorship data.  Can. J.
   #    Fish. Aquat. Sci. 47: 24-40.

   # If the age lies below the theoretical age at which the length is zero,
   # the estimated length is assumed to be zero.

   # robustify routine, prevent y2 being less than y1
   if (y2 < y1+2) {
     y2 = y1+2
     # cat("SchnuteGrowthfunction: y1", y1,"y2",y2,'\n')
   }

   # ' Determine which equation is to be used
   if (a == 0) {
     if (b == 0) {
       # Eqn(18)
       y = y1 * exp(log(y2 / y1) * (Age - t1) / (t2 - t1))

       if (y < 0) {
         y = 0
       } # y < 0
     } else {          #(b == 0)
       # Eqn(17)
       # First, let's work out tzero
       #cat("Age", Age, "y1",y1,"y2",y2,"b",b,'\n')
       tzero = t1 - (y1 ^ b) * (t2 - t1) / (y2 ^ b - y1 ^ b)

       if (Age < tzero) {
         y = 0
       } else {
         v = (y1 ^ b + (y2 ^ b - y1 ^ b) * (Age - t1) / (t2 - t1))
         y = v ^ (1 / b)
       }
     }
   }# a == 0

   if (a != 0) {
     if (b == 0) {
       # Eqn(16)
       y = y1 * exp(log(y2 / y1) * (1 - exp(-a * (Age - t1))) / (1 - exp(-a * (t2 - t1))))
       if (y < 0) {
         y = 0 }
     } else {
       # Eqn(15)
       # First. let's work out tzero
       if (1 + (y1 ^ b) * (1 - exp(-a * (t2 - t1))) / (y2 ^ b - y1 ^ b) <= 0) {
         tzero = t1 - log(1E-4) / a
       } else {
         tzero = t1 - log(1 + (y1 ^ b) * (1 - exp(-a * (t2 - t1))) / (y2 ^ b - y1 ^ b)) / a
       }
       if (is.nan(tzero)) {
         cat("SchnuteGrowthfunction: Problem calculating tzero",'\n')
       }
       # cat("Age",Age, "tzero",tzero, '\n')
       if (Age < tzero) {
         y = 0
       } else {
         v = (y1 ^ b + (y2 ^ b - y1 ^ b)
              * (1 - exp(-a * (Age - t1))) / (1 - exp(-a * (t2 - t1))))
         y = v ^ (1 / b)
       } # else
     } # b == 0
   }  # a != 0

   return(y)
 } # end function


 #' Inverse Schnute growth function
 #'
 #' Calculate age given length, from Schnute growth function
 #'
 #' @keywords internal
 #'
 #' @param FishLen specified length
 #' @param t1 first reference age
 #' @param t2 second reference age
 #' @param y1 length at t1
 #' @param y2 length at t2
 #' @param a growth curve parameter
 #' @param a growth curve parameter
 #'
 #' @return Age at specified length
 InverseSchnuteGrowthfunction <- function (FishLen, t1, t2, y1, y2, a, b) {

   # return age from Length, given Schnute growth parameters
   # (can be used when both a and b are not equal to zero)
   Linf = ((exp(a*t2)*y2^b-exp(a*t1)*y1^b)/(exp(a*t2)-exp(a*t1)))^(1/b)


   if (FishLen >= Linf) {
     FishLen=Linf
     Age=MaxAge
   } else {
     Age=log(1-(FishLen^b-y1^b)/(y2^b-y1^b)*(1-exp(-a*(t2-t1))))/-a+t1
   }


   return(Age)
 }

 #' Calculate expected length after year of growth
 #'
 #' Calculate expected length after year of growth, given initial length (from von Bertalanffy or Schnute function)
 #'
 #' @keywords internal
 #'
 #' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
 #' @param TimeStep model timestep (fraction of year)
 #' @param GrowthParamsForSex c(Linf,vbK) von Bertalanffy, or c(t1,t2,y1,y2,a,b) Schnute
 #' @param RefnceAgesForSex c(t1,t2)
 #' @param midpt mid points of length classes
 #'
 #' @return Expected lengths (ExpLen)
 CalcLengthAfterGrowthForTimetep <- function(GrowthCurveType, TimeStep, GrowthParamsForSex, RefnceAgesForSex, midpt) {

   # calculate mean length after growth for current timestep, given growth curve type,
   # mean length of the current length class and initial length

   if (GrowthCurveType == 1) { # von Bert
     Linf = GrowthParams[1]
     vbK = GrowthParams[2]
     ExpLen = midpt + (Linf - midpt) * (1 - exp(-vbK*TimeStep))
   }
   if (GrowthCurveType == 2) { # Schnute

     t1 = RefnceAgesForSex[1]
     t2 = RefnceAgesForSex[2]
     y1 = GrowthParamsForSex[1]
     y2 = GrowthParamsForSex[2]
     a = GrowthParamsForSex[3]
     b = GrowthParamsForSex[4]

     nObs <- length(midpt)
     StartAge = rep(NA,nObs)
     ExpLen = rep(NA,nObs)
     for (i in 1:nObs) {
       FishLen = midpt[i]
       StartAge[i] = InverseSchnuteGrowthfunction(FishLen, t1, t2, y1, y2, a, b)
       Age <- StartAge[i] + TimeStep
       ExpLen[i] = SchnuteGrowthfunction(Age, t1, t2, y1, y2, a, b)
       #cat("Age",Age,"MeanLen",MeanLen[i],'\n')
     }
   }

   return(ExpLen)
 }

#******************************************
# Catch curve analyses - size-based methods
#******************************************

#' Calculate gillnet selectivity
#'
#' Calculate gillnet selectivity at length, by mesh and overall (assuming equal fishing intensity
#' among meshes) using the method of Kirkwood and Walker (1986), given model parameter values
#' (theta1 and theta2)
#'
#' @param theta1 parameter of Kirkwood and Walker (1986) model
#' @param theta2 parameter of Kirkwood and Walker (1986) model
#' @param MeshSize_mm vector of gillnet mesh sizes, in mm
#' @param nMeshes number of gillnet meshes
#' @param nLenCl number of length classes
#' @param midpt mid points of length classes
#'
#' @return Selectivity at length for each mesh (SelAtLengthForMesh), overall length-based
#' selectivity (SelAtLength)
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

#' Calculate logistic length-based selectivity (or retention)
#'
#' Calculate logistic length-based selectivity or retention (asymptotic curve)
#'
#' @keywords internal
#'
#' @param L50 length at which 50 percent of fish are selected into the fishery
#' @param L95 length at which 95 percent of fish are selected into the fishery
#' @param nLenCl number of length classes
#' @param midpt mid points of length classes
#'
#' @return Selectivity at length (SelAtLength)
CalcLogisticSelOrReten <- function(L50, L95, midpt) {

  nLenCl = length(midpt)
  SelAtLength = rep(0,nLenCl)
  SelAtLength = 1 / (1 + exp(-log(19) * (midpt - L50) / (L95-L50)))

  return(SelAtLength)

}

#' Calculate size distribution of 1+ year old recruits
#'
#' Calculate size distribution of 1+ year old recruits, given mean length
#' at age 1 and specified CV for lengths at age, assuming a normal distribution
#'
#' @keywords internal
#'
#' @param MeanSizeAtAge mean length at age from growth curve
#' @param CVSizeAtAge coefficient of variation (CV), common across all mean lengths at age
#' @param lbnd lower bounds of length classes
#' @param ubnd upper bounds of length classes
#' @param midpt mid points of length classes
#' @param nLenCl number of length classes
#'
#' @return Selectivity at length (SelAtLength)
CalcSizeDistOfRecruits <- function(MeanSizeAtAge, CVSizeAtAge, lbnd, ubnd, midpt, nLenCl) {

  # Mean length and SD of 0+ recruits, at 1 year of age
  if (is.vector(MeanSizeAtAge)) {
    MeanLenRec <- MeanSizeAtAge[1]
    # robustify for positive tzero values
    if (MeanLenRec < 5) MeanLenRec = 5 # i.e. assuming mean size at least 20 mm at first time step
    SDAgeOneRecruits = MeanLenRec * CVSizeAtAge

    RecLenDist = rep(0,nLenCl)
    RecLenDist = pnorm(ubnd, mean=MeanLenRec, sd=SDAgeOneRecruits, lower.tail = T) -
      pnorm(lbnd, mean=MeanLenRec, sd=SDAgeOneRecruits,lower.tail = T)

    RecLenDist = RecLenDist / sum(RecLenDist)
  }
  if (is.matrix(MeanSizeAtAge) | is.data.frame(MeanSizeAtAge)) {
    RecLenDist <- data.frame(matrix(nrow = 2, ncol = nLenCl))
    colnames(RecLenDist) <- midpt
    RecLenDist = as.matrix(RecLenDist)
    MeanLenRec = rep(0,2)
    SDAgeOneRecruits = rep(0,2)

    for (i in 1:2) {
      MeanLenRec[i] <- MeanSizeAtAge[i,1]
      # robustify for positive tzero values
      if (MeanLenRec[i] < 5) MeanLenRec[i] = 5 # i.e. assuming mean size at least 50 mm at 1 year of age
      if (length(CVSizeAtAge)==1) { # age and length-based catch curve, estimated cv value
        SDAgeOneRecruits[i] = MeanLenRec[i] * CVSizeAtAge[i]
      }
      if (length(CVSizeAtAge)==2) { # length-based catch curve, sex-specific inputted cv values
        SDAgeOneRecruits[i] = MeanLenRec[i] * CVSizeAtAge[i]
      }

      RecLenDist[i,] = pnorm(ubnd, mean=MeanLenRec[i], sd=SDAgeOneRecruits[i], lower.tail = T) -
        pnorm(lbnd, mean=MeanLenRec[i], sd=SDAgeOneRecruits[i],lower.tail = T)
      RecLenDist[i,] = RecLenDist[i,] / sum(RecLenDist[i,])
    }
  }

  return(RecLenDist)

}

#' Visualise length at age data when applying a length transition matrix
#'
#' This function simulates growth trajectories for individual fish for visualising length
#' at age data produced when applying a length transition matrix, specified
#' mean sizes at ages and associated variation (common CV for mean lengths at age). Note
#' that this diagnostic analysis does not account for any fishing mortality (or selectivity)
#' effects on variation in data.
#'
#' @param nFish number of growth trajectories for individual fish to simulate
#' @param TimeStep model timestep (fraction of year)
#' @param MinAge minimum age
#' @param MinAge maximum age
#' @param CVSizeAtAge coefficient of variation (CV), common across all mean lengths at age
#' @param lbnd mid points of length classes
#' @param midpt mid points of length classes
#' @param ubnd mid points of length classes
#' @param nLenCl number of length classes
#' @param MainLabel plot label
#' @param xaxis_lab y axis label
#' @param yaxis_lab x axis label
#' @param xmax maximum x axis value
#' @param xint x axis interval
#' @param ymax maximum y axis value
#' @param yint y axis interval
#'
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
#' Growth_params = c(Linf, vbK, tzero) # currently only set up for von Bertalanffy growth equation
#' CVSizeAtAge = 0.05
#' nFish = 100
#' TimeStep = 1
#' VisualiseGrowthApplyingLTM(nFish, TimeStep, MaxAge, Growth_params, CVSizeAtAge, lbnd, midpt, ubnd,
#'                            nLenCl, MainLabel=NA, xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA, ymax=NA, yint=NA)
#' @export
VisualiseGrowthApplyingLTM <- function (nFish, TimeStep, MaxAge, Growth_params, CVSizeAtAge, lbnd, midpt, ubnd,
                                        nLenCl, MainLabel, xaxis_lab, yaxis_lab, xmax, xint, ymax, yint) {

  # a bit of code to check growth transition matrix is likely to
  # give sensible results - not part of catch curve model

  # note, this doesn't account for mortality effects
  Ages = seq(TimeStep,MaxAge,TimeStep)
  if (is.na(xaxis_lab)) xaxis_lab = "Age (y)"
  if (is.na(yaxis_lab)) yaxis_lab = "Length (mm)"
  xlims = Get_xaxis_scale(Ages)
  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint
  ylims = Get_yaxis_scale(midpt)
  if (is.na(ymax)) ymax = ylims$ymax
  if (is.na(yint)) yint = ylims$yint

  Linf = Growth_params[1]
  vbK = Growth_params[2]
  tzero = Growth_params[3]
  MeanSizeAtAge = Linf * (1 - exp (-vbK * (Ages - tzero))) # e.g. as calculated from a von Bertalanffy growth curve
  MeanEndingLength = midpt + (Linf - midpt) * (1 - exp(-vbK*TimeStep)) # e.g. as calculated from a von Bertalanffy growth curve
  TimestepGrowthSizeInc = MeanEndingLength-midpt # amount of annual growth with respect to initial length

  LTM = CalcLTM_cpp(TimestepGrowthSizeInc, CVSizeAtAge, lbnd, midpt, ubnd, nLenCl)

  FishLenAtAge = data.frame(matrix(nrow=nFish, ncol=length(Ages)))
  colnames(FishLenAtAge) = Ages
  FishLenAtAge = as.matrix(FishLenAtAge)

  for (j in 1:nFish) { # fish
    k = 0
    for (kk in Ages) { # age
     k=k+1
      if (k ==1) {
        rndlen = rnorm(1,MeanSizeAtAge[1],MeanSizeAtAge[1]*CVSizeAtAge)
        FishLenAtAge[j,k] = rndlen
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

        FishLenAtAge[j,k] = midpt[rndlc]
      } # else
    } # age

  } # fish

  plot(Ages, MeanSizeAtAge,  "l", cex.main=1, frame.plot=F, xlab = xaxis_lab,
       ylab = yaxis_lab, xlim=c(0,xmax), ylim=c(0,ymax), main=MainLabel)
  for (j in 1:nFish) { # fish
    lines(Ages, FishLenAtAge[j,], cex=0.6, col=j)
  }
  lines(Ages, MeanSizeAtAge, lwd=2, cex=0.6, col="black")
}


#' Get results for a fitted length-based catch curve
#'
#' This function fits a length-based catch curve with length-based selectivity, based on
#' either 1) specified length-based selectivity inputted as a vector, or 2) an estimated asymptotic logistic
#' selectivity curve, and specified growth curve (either von Bertalanffy or Schnute growth model). The model is
#' to a sample of fish length frequency data, by minimising the negative log-likelihood associated
#' with the parameters and data, using nlminb. It provides various statistical outputs in include
#' convergence statistics, parameter estimates and associated 95 percent confidence limits and associated
#' variance-covariance matrix, calculated using the MASS package.
#'
#' @param params vector of model parameters in log space (params) to be estimated
#' @param GrowthCurveType 1 = von Bertalanffy, 2 = Schnute
#' @param GrowthParams c(Linf, vbK, CVSizeAtAge) single sex von Bertalanffy, or data.frame(Linf=Linf, vbK=vbK, CVSizeAtAge=CVSizeAtAge),
#' both sexes von Bertalanffy, or c(y1, y2, a, b) single sex Schnute, or data.frame(y1=y1, y2=y2, a=a, b=b), both sexes Schnute
#' @param RefnceAges reference ages for Schnute function (set to NA if growth based on another function)
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param SelectivityType 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' @param ObsRetCatchFreqAtLen observed frequency of fish at length in sample
#' @param lbnd lower bounds of length classes
#' @param ubnd upper bounds of length classes
#' @param midpt mid points of length classes
#' @param SelectivityVec selectivity at length
#' @param PropReleased proportion of fish that are released
#' @param ObsDiscCatchFreqAtLen proportion of fish that are released
#' @param DiscMort Proportion of fish that die following to capture and release
#' @param CVSizeAtAge coefficient of variation (CV), common across all mean lengths at age
#' @param MaxAge maximum age considered in model
#' @param NatMort natural mortality
#' @param TimeStep model timestep (e.g. 1 = annual, 1/12 = monthly)
#'
#' @return negative log-likelihood (nll), nlminb convergence diagnostic (convergence)
#' sample size (SampleSize), growth parameter estimates with lower and upper 95 percent
#' confidence limits (ParamEst), point estimates for estimated parameters (params)
#' and variance-covariance matrix (vcov.Params), gear selectivity at length (SelAtLength), retention at length (RetAtLength),
#' selectivity of landings at length (SelLandAtLength), growth curve (MeanSizeAtAge), midpoint of each length class (midpt),
#' mean length after 1 year from growth curve, given initial length (MeanEndingLength), mean change in length after 1 year,
#' from initial length - note, assuming normal a distribution allows for possibility of negative growth
#' if above asymptotic length (TimestepGrowthSizeInc), length distribution of 1+ year old recruits (RecLenDist),
#' expected retained catches, at length (RetCatchAtLen), expected discarded catch at length (DiscCatchAtLen),
#' expected total (discard plus retained) catches at length (TotCatchAtLen), proportion of catch in each length classes for sexes combined
#' and separate (ExpRetCatchPropInLenClass, ExpRetCatchPropInLenClass_Fem, ExpRetCatchPropInLenClass_Mal),
#' expected retained catch proportions given integer age (ExpRetCatchPropLengthGivenIntAge, ExpRetCatchPropLengthGivenIntAge_Fem,
#' ExpRetCatchPropLengthGivenIntAge_Mal), observed catch in each length class for sexes combined (ObsRetCatchFreqAtLen), CV for modelled lengthsat age, around mean
#' length at MaxAge, for sexes combined or separate (CV_LenAtMaxAge, FemCV_LenAtMaxAge, MalCV_LenAtMaxAge), fishing mortality
#' at length, for sexes combined or separate (FAtLen, FAtLen_Fem, FAtLen_Mal), total mortality at length, for sexes combined
#' or separate (ZAtLen, ZAtLen_Fem, ZAtLen_Mal), fishing mortality associated with landings at length, for sexes combined or
#' separate (FAtLen, FAtLen_Fem, FAtLen_Mal), fishing mortality associated with capture and retention at length, for sexes combined or
#' separate (FAtLenReten, FAtLenReten_Fem, FAtLenReten_Mal), fishing mortality associated with capture and discarding at
#' length, for sexes combined or separate (FAtLenDisc, FAtLenDisc_Fem, FAtLenDisc_Mal)
#'
#' @examples
#' # Simulate data
#' SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
#' set.seed(123)
#' MaxAge = 30
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = 0.2
#' MaxLen = 1500
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 500) # L50, L95 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95 for retention
#' DiscMort = 0 # proportion of fish that die due to natural mortality
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.08
#' GrowthParams = c(Linf, vbK, CVSizeAtAge)
#' RefnceAges = NA
#' # # 2 sexes, von Bertalanffy
#' # GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' # Linf = c(700,850)
#' # vbK = c(0.3,0.2)
#' # CVSizeAtAge = c(0.08,0.08)
#' # GrowthParams = data.frame(Linf=Linf, vbK=vbK, CVSizeAtAge=CVSizeAtAge)
#' # RefnceAges = NA
#' # 1 sex, Schnute
#' # GrowthCurveType = 2 # 1 = von Bertalanffy, 2 = Schnute
#' # t1 = 1 # growth - Schnute
#' # t2 = 10 # growth - Schnute
#' # y1 = 400 # growth - Schnute
#' # y2 = 1000 # growth - Schnute
#' # a = 0.1 # growth - Schnute
#' # b = 2.0 # growth - Schnute
#' # GrowthParams = c(y1, y2, a, b)
#' # RefnceAges = c(t1,t2)
#' # CVSizeAtAge = 0.08
#' # # 2 sexes, Schnute
#' # GrowthCurveType = 2 # 1 = von Bertalanffy, 2 = Schnute
#' # t1 = c(1,1) # growth - Schnute
#' # t2 = c(10,10) # growth - Schnute
#' # y1 = c(435.3,400.3) # growth - Schnute
#' # y2 = c(1089.1,1034.5) # growth - Schnute
#' # a = c(0.044,0.136) # growth - Schnute
#' # b = c(2.748,1.971) # growth - Schnute
#' # CVSizeAtAge = c(0.08, 0.08)
#' # GrowthParams = data.frame(y1=y1, y2=y2, a=a, b=b)
#' # RefnceAges = data.frame(t1=t1,t2=t2)
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen
#' ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsRetCatchFreqAtLen)
#' PropReleased = NA # proportion of fish released, vector including mean and sd (option probably now obselete)
#' midpt=Res$midpt
#' lbnd=Res$lbnd
#' ubnd=Res$ubnd
#' InitFishMort = 0.25 # specify starting parameters
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
#' InitL50 = 600
#' InitDelta = 100
#' params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
#' Res=GetLengthBasedCatchCurveResults(params, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
#'                                           lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
#' #
#' # # Example with selectivity specified as a vector
#' # # Simulate data
#' # SampleSize=5000
#' # set.seed(123)
#' # MaxAge = 30
#' # TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' # NatMort = 4.22/MaxAge
#' # FishMort = 0.2
#' # MaxLen = 1500
#' # LenInc = 20
#' # MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
#' # SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' # lbnd = seq(0,MaxLen - LenInc, LenInc)
#' # midpt = lbnd + (LenInc/2)
#' # SelectivityVec = 1 / (1 + exp(-log(19)*(midpt-600)/(700-600)))
#' # SelParams = c(NA, NA) # L50, L95 for gear selectivity
#' # RetenParams = c(NA, NA) # L50, L95 for retention
#' # DiscMort = 0 # proportion of fish that die due to natural mortality
#' #  # single sex, von Bertalanffy
#' # GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' # Linf = 800
#' # vbK = 0.2
#' # CVSizeAtAge = 0.08
#' # GrowthParams = c(Linf, vbK, CVSizeAtAge)
#' # RefnceAges = NA
#' # Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#' #                        SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen
#' # ObsDiscCatchFreqAtLen = NA # vector including mean and sd
#' # midpt=Res$midpt
#' # lbnd=Res$lbnd
#' # ubnd=Res$ubnd
#' # InitFishMort = 0.25 # specify starting parameters
#' # InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
#' # params = c(InitFishMort_logit)
#' # Res=GetLengthBasedCatchCurveResults(params, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
#' #                                      lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
#' @export
GetLengthBasedCatchCurveResults <- function (params, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                             lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
{
  nlmb <- nlminb(params, CalcObjFunc_LengthBasedCatchCurve,
                 gradient = NULL, hessian = TRUE)
  (hess.out = optimHess(nlmb$par, CalcObjFunc_LengthBasedCatchCurve))
  (vcov.Params = solve(hess.out))
  (ses = sqrt(diag(vcov.Params)))
  # logit space
  temp = c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1])
  # inverse logit transformed value
  EstFMort = 1/(1+exp(-temp))

  if (SelectivityType == 1) {
    ParamEst = t(data.frame(FMort = round(EstFMort, 2)))
  }
  if (SelectivityType == 2) { # include estimates for logistic gear selectivity parameters

    EstL50_sel = exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
    EstDelta_sel = exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
    if (length(params==3)) {
      ParamEst = t(data.frame(FMort = round(EstFMort, 3), L50_sel = round(EstL50_sel, 3), Delta_sel = round(EstDelta_sel, 3)))
    }
    if (length(params==5)) {
      EstL50_ret = exp(c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
      EstDelta_ret = exp(c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5]))
      ParamEst = t(data.frame(FMort = round(EstFMort, 3), L50_sel = round(EstL50_sel, 3), Delta_sel = round(EstDelta_sel, 3),
                              L50_ret = round(EstL50_ret, 3), Delta_ret = round(EstDelta_ret, 3)))
    }
  }
  colnames(ParamEst) = c("Estimate", "lw_95%CL", "up_95%CL")

  # store some diagnostic outputs from model
  params = nlmb$par
  CatchCurveType=1 #1=length-based, 2=age and length based
  res = AgeAndLengthBasedCatchCurvesCalcs(params, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge, CatchCurveType,
                                          MLL, SelectivityType, lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)


  SampleSize = sum(ObsRetCatchFreqAtLen)
  nll = nlmb$objective
  convergence = nlmb$convergence
  SelAtLength = res$SelAtLength
  RetAtLength = res$RetAtLength
  SelLandAtLength = res$SelLandAtLength
  GrowthModelType = res$GrowthModelType
  ExpRetCatchPropInLenClass = res$ExpRetCatchPropInLenClass
  ExpRetCatchPropInLenClass_Fem = res$ExpRetCatchPropInLenClass_Fem
  ExpRetCatchPropInLenClass_Mal = res$ExpRetCatchPropInLenClass_Mal
  ExpDiscCatchPropInLenClass = res$ExpDiscCatchPropInLenClass
  ExpDiscCatchPropInLenClass_Fem = res$ExpDiscCatchPropInLenClass_Fem
  ExpDiscCatchPropInLenClass_Mal = res$ExpDiscCatchPropInLenClass_Mal
  ExpRetCatchPropLengthGivenIntAge = res$ExpRetCatchPropLengthGivenIntAge
  ExpRetCatchPropLengthGivenIntAge_Fem = res$ExpRetCatchPropLengthGivenIntAge_Fem
  ExpRetCatchPropLengthGivenIntAge_Mal = res$ExpRetCatchPropLengthGivenIntAge_Mal
  RetCatchAtLen = res$RetCatchAtLen # retained catch
  TotCatchAtLen = res$TotCatchAtLen # total catch
  DiscCatchAtLen = res$DiscCatchAtLen # discarded catch
  FAtLen = res$FAtLen
  FAtLen_Fem = res$FAtLen_Fem
  FAtLen_Mal = res$FAtLen_Mal
  ZAtLen = res$ZAtLen
  ZAtLen_Fem = res$ZAtLen_Fem
  ZAtLen_Mal = res$ZAtLen_Mal
  FAtLenReten = res$FAtLenReten
  FAtLenReten_Fem = res$FAtLenReten_Fem
  FAtLenReten_Mal = res$FAtLenReten_Mal
  FAtLenDisc = res$FAtLenDisc
  FAtLenDisc_Fem = res$FAtLenDisc_Fem
  FAtLenDisc_Mal = res$FAtLenDisc_Mal

  # get expected proportion of fish that are released
  EstPropReleased = sum(res$DiscCatchAtLen) / sum(res$TotCatchAtLen)

  # get expected mean length of released fish
  ExpMeanLenReleased = sum((res$DiscCatchAtLen * res$midpt)) / sum(res$DiscCatchAtLen)

  # calculate approximate cv for lengths at max age (growth diagnostic)
  # single / combined sex
  MinAge = floor(TimeStep)
  nAgeCl = length(MinAge:MaxAge)

  if (GrowthModelType ==  1 | GrowthModelType == 3) { # single sex
    if (TimeStep==1) {
      Probs = as.vector(unlist(res$ExpRetCatchPropLengthGivenIntAge[nAgeCl,]))
      RandFreq=rep(midpt,(as.vector(rmultinom(1,1000,Probs))))
      MeanRandFreq=mean(RandFreq)
      sdRandFreq=sd(RandFreq)
      CV_LenAtMaxAge=sdRandFreq/MeanRandFreq
    } else {
      CV_LenAtMaxAge = NA
    }
    FemCV_LenAtMaxAge = NA
    MalCV_LenAtMaxAge = NA
  }

  # females
  if (GrowthModelType ==  2 | GrowthModelType == 4) { # separate sex
    if (TimeStep==1) {
      FemProbs = as.vector(unlist(res$ExpRetCatchPropLengthGivenIntAge_Fem[nAgeCl,]))
      FemRandFreq=rep(midpt,(as.vector(rmultinom(1,1000,FemProbs))))
      FemMeanRandFreq=mean(FemRandFreq)
      FemsdRandFreq=sd(FemRandFreq)
      FemCV_LenAtMaxAge=FemsdRandFreq/FemMeanRandFreq

      # males
      MalProbs = as.vector(unlist(res$ExpRetCatchPropLengthGivenIntAge_Mal[nAgeCl,]))
      MalRandFreq=rep(midpt,(as.vector(rmultinom(1,1000,MalProbs))))
      MalMeanRandFreq=mean(MalRandFreq)
      MalsdRandFreq=sd(MalRandFreq)
      MalCV_LenAtMaxAge=MalsdRandFreq/MalMeanRandFreq
    } else {
      FemCV_LenAtMaxAge = NA
      MalCV_LenAtMaxAge = NA
    }
    CV_LenAtMaxAge = NA
  }


  Results = list(nll = nll,
                 convergence = convergence,
                 SampleSize = SampleSize,
                 ParamEst = ParamEst,
                 params = nlmb$par,
                 vcov.Params = vcov.Params,
                 SelAtLength=res$SelAtLength,
                 RetAtLength = res$RetAtLength,
                 SelLandAtLength = res$SelLandAtLength,
                 MeanSizeAtAge=res$MeanSizeAtAge,
                 midpt=res$midpt,
                 MeanEndingLength=res$MeanEndingLength,
                 TimestepGrowthSizeInc=res$TimestepGrowthSizeInc,
                 RecLenDist=res$RecLenDist,
                 ExpRetCatchPropInLenClass=res$ExpRetCatchPropInLenClass,
                 ExpRetCatchPropInLenClass_Fem=res$ExpRetCatchPropInLenClass_Fem,
                 ExpRetCatchPropInLenClass_Mal=res$ExpRetCatchPropInLenClass_Mal,
                 ExpDiscCatchPropInLenClass = ExpDiscCatchPropInLenClass,
                 ExpDiscCatchPropInLenClass_Fem = ExpDiscCatchPropInLenClass_Fem,
                 ExpDiscCatchPropInLenClass_Mal = ExpDiscCatchPropInLenClass_Mal,
                 ObsRetCatchFreqAtLen=res$ObsRetCatchFreqAtLen,
                 ExpRetCatchPropLengthGivenIntAge=ExpRetCatchPropLengthGivenIntAge,
                 ExpRetCatchPropLengthGivenIntAge_Fem=ExpRetCatchPropLengthGivenIntAge_Fem,
                 ExpRetCatchPropLengthGivenIntAge_Mal=ExpRetCatchPropLengthGivenIntAge_Mal,
                 CV_LenAtMaxAge=CV_LenAtMaxAge,
                 FemCV_LenAtMaxAge=FemCV_LenAtMaxAge,
                 MalCV_LenAtMaxAge=MalCV_LenAtMaxAge,
                 RetCatchAtLen = RetCatchAtLen,
                 DiscCatchAtLen = DiscCatchAtLen,
                 TotCatchAtLen = TotCatchAtLen,
                 EstPropReleased = EstPropReleased,
                 ExpMeanLenReleased = ExpMeanLenReleased,
                 FAtLen = FAtLen,
                 FAtLen_Fem = FAtLen_Fem,
                 FAtLen_Mal = FAtLen_Mal,
                 ZAtLen = ZAtLen,
                 ZAtLen_Fem = ZAtLen_Fem,
                 ZAtLen_Mal = ZAtLen_Mal,
                 FAtLenReten = FAtLenReten,
                 FAtLenReten_Fem = FAtLenReten_Fem,
                 FAtLenReten_Mal = FAtLenReten_Mal,
                 FAtLenDisc = FAtLenDisc,
                 FAtLenDisc_Fem = FAtLenDisc_Fem,
                 FAtLenDisc_Mal = FAtLenDisc_Mal)

  return(Results)
}


#' Get NLL for length based catch curve
#'
#' @keywords internal
#'
#' @param params estimated parameters, including fishing mortality and logistic selectivity parameters when estimated
#'
#' @return negative log-likelihood (NLL)
CalcObjFunc_LengthBasedCatchCurve <- function(params) {
  # get NLL for length based catch curve, for optimisation

  CatchCurveType=1 #1=length-based, 2=age and length based
  Res = AgeAndLengthBasedCatchCurvesCalcs(params, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge, CatchCurveType,
                                          MLL, SelectivityType, lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)

  ExpRetCatchPropInLenClass = Res$ExpRetCatchPropInLenClass
  ExpDiscCatchPropInLenClass = Res$ExpDiscCatchPropInLenClass

  NLL_PropReleased = 0
  if (!is.na(PropReleased[1])) {
    EstPropReleased = sum(Res$DiscCatchAtLen) / sum(Res$TotCatchAtLen)
    NLL_PropReleased = -dnorm(EstPropReleased, PropReleased[1], PropReleased[2], log=TRUE)
  }

  NLL_RetCatch = CalcNLLMargLengthComposition(ObsRetCatchFreqAtLen, ExpRetCatchPropInLenClass)

  NLL_DiscCatch = 0
  if (!is.na(ObsDiscCatchFreqAtLen[1])) {
    NLL_DiscCatch = CalcNLLMargLengthComposition(ObsDiscCatchFreqAtLen, ExpDiscCatchPropInLenClass)
  }

  NLL = NLL_RetCatch + NLL_PropReleased + NLL_DiscCatch + Res$L50_Pen + Res$L95_Pen

  cat("NLL", NLL, " NLL_DiscCatch ", NLL_DiscCatch, " NLL_PropReleased ", NLL_PropReleased,
      " FMort", 1/(1+exp(-params[1])), " L50_sel ", exp(params[2]), " delta_sel ", exp(params[3]),
      " L50_ret ", exp(params[4]), " delta_ret ", exp(params[5]),
      " L50_Pen ", Res$L50_Pen, " L95_Pen " ,Res$L95_Pen, '\n')

  return(NLL)

}


#' Get NLL for age and length-based catch curve
#'
#' @keywords internal
#' @param estimated parameters, including fishing mortality, growth and logistic selectivity parameters
#' @return negative log-likelihood (NLL)
CalcObjFunc_AgeAndLengthBasedCatchCurve <- function(params) {
  # get NLL for length based catch curve, for optimisation

  CatchCurveType=2
  GrowthParams = NA
  RefnceAges = NA
  CVSizeAtAge = NA
  GrowthCurveType = 1 # von Bertalanffy
  Res = AgeAndLengthBasedCatchCurvesCalcs(params, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge, CatchCurveType,
                                          MLL, SelectivityType, lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)


  L50_Pen = Res$L50_Pen
  L95_Pen = Res$L95_Pen
  GrowthModelType = Res$GrowthModelType # 1 = von Bert single sex, 2 = von Bert, 2 sexes,
  # 3 = Schnute - single sex, 4 = Schunte 2 sexes

  nLenCl = length(midpt)
  MinAge = floor(TimeStep)
  nAgeCl = length(MinAge:MaxAge)

  if (GrowthModelType == 1 | GrowthModelType == 3) {
    # get NLL for marginal length composition
    ExpRetCatchPropInLenClass = Res$ExpRetCatchPropInLenClass
    Length_NLL = CalcNLLMargLengthComposition(ObsRetCatchFreqAtLen, ExpRetCatchPropInLenClass)

    # get NLL for age at length observations
    ExpRetCatchPropAtIntAge = as.matrix(Res$ExpRetCatchPropAtIntAge)
    #cat("ExpRetCatchPropAtIntAge",ExpRetCatchPropAtIntAge,'\n')
    ExpRetCatchPropLengthGivenIntAge = as.matrix(Res$ExpRetCatchPropLengthGivenIntAge)
    #cat("ExpRetCatchPropLengthGivenIntAge",ExpRetCatchPropLengthGivenIntAge,'\n')
    ExpRetCatchPropIntAgeGivenLength = CalcExpCatchPropIntAgeGivenLength_cpp(nLenCl, nAgeCl, ExpRetCatchPropLengthGivenIntAge, ExpRetCatchPropAtIntAge)

    if (TimeStep == 1) {
      ObsCatchFreqAtLengthAndIntAge = ObsRetCatchFreqAtLengthAndAge
    } else {
      ObsCatchFreqAtLengthAndIntAge = ConvertObsDataFromDecAgesToIntegerAges(TimeStep, MaxAge, nLenCl, ObsRetCatchFreqAtLengthAndAge)
      #cat("ObsCatchFreqAtLengthAndIntAge",ObsCatchFreqAtLengthAndIntAge,'\n')
    }
    CondAgeAtLengthNLL = CalcNLLCondAgeAtLength_cpp(nLenCl, nAgeCl, ObsCatchFreqAtLengthAndIntAge, ExpRetCatchPropIntAgeGivenLength)
  }
  if (GrowthModelType == 2 | GrowthModelType == 4) {

    # get NLL for marginal length composition
    # females
    ExpRetCatchPropInLenClass_Fem = Res$ExpRetCatchPropInLenClass_Fem
    ObsRetCatchFreqAtLen_Fem = unlist(ObsRetCatchFreqAtLen[1,])
    Length_NLL_Fem = CalcNLLMargLengthComposition(ObsRetCatchFreqAtLen_Fem, ExpRetCatchPropInLenClass_Fem)
    # males
    ExpRetCatchPropInLenClass_Mal = Res$ExpRetCatchPropInLenClass_Mal
    ObsRetCatchFreqAtLen_Mal = unlist(ObsRetCatchFreqAtLen[2,])
    Length_NLL_Mal = CalcNLLMargLengthComposition(ObsRetCatchFreqAtLen_Mal, ExpRetCatchPropInLenClass_Mal)
    Length_NLL = Length_NLL_Fem + Length_NLL_Mal

    # get NLL for age at length observations
    # females
    ExpRetCatchPropAtIntAge_Fem = as.matrix(Res$ExpRetCatchPropAtIntAge_Fem)
    ExpRetCatchPropLengthGivenIntAge_Fem = as.matrix(Res$ExpRetCatchPropLengthGivenIntAge_Fem)
    ObsRetCatchFreqAtLengthAndAge_Fem = unlist(ObsRetCatchFreqAtLengthAndAge[,,1])
    if (TimeStep == 1) {
      ObsCatchFreqAtLengthAndIntAge_Fem = ObsRetCatchFreqAtLengthAndAge_Fem
    } else {
      ObsCatchFreqAtLengthAndIntAge_Fem = ConvertObsDataFromDecAgesToIntegerAges(TimeStep, MaxAge, nLenCl, ObsRetCatchFreqAtLengthAndAge_Fem)
    }
    ExpRetCatchPropIntAgeGivenLength_Fem = CalcExpCatchPropIntAgeGivenLength_cpp(nLenCl, nAgeCl, ExpRetCatchPropLengthGivenIntAge_Fem, ExpRetCatchPropAtIntAge_Fem)
    CondAgeAtLengthNLL_Fem = CalcNLLCondAgeAtLength_cpp(nLenCl, nAgeCl, ObsCatchFreqAtLengthAndIntAge_Fem, ExpRetCatchPropIntAgeGivenLength_Fem)
    # males
    ExpRetCatchPropAtIntAge_Mal = as.matrix(Res$ExpRetCatchPropAtIntAge_Mal)
    ExpRetCatchPropLengthGivenIntAge_Mal = as.matrix(Res$ExpRetCatchPropLengthGivenIntAge_Mal)
    ObsRetCatchFreqAtLengthAndAge_Mal = unlist(ObsRetCatchFreqAtLengthAndAge[,,2])
    if (TimeStep == 1) {
      ObsCatchFreqAtLengthAndIntAge_Mal = ObsRetCatchFreqAtLengthAndAge_Mal
    } else {
      ObsCatchFreqAtLengthAndIntAge_Mal = ConvertObsDataFromDecAgesToIntegerAges(TimeStep, MaxAge, nLenCl, ObsRetCatchFreqAtLengthAndAge_Mal)
    }
    ExpRetCatchPropIntAgeGivenLength_Mal = CalcExpCatchPropIntAgeGivenLength_cpp(nLenCl, nAgeCl, ExpRetCatchPropLengthGivenIntAge_Mal, ExpRetCatchPropAtIntAge_Mal)
    CondAgeAtLengthNLL_Mal = CalcNLLCondAgeAtLength_cpp(nLenCl, nAgeCl, ObsCatchFreqAtLengthAndIntAge_Mal, ExpRetCatchPropIntAgeGivenLength_Mal)

    CondAgeAtLengthNLL = CondAgeAtLengthNLL_Fem + CondAgeAtLengthNLL_Mal
  }

  NLL = Length_NLL + CondAgeAtLengthNLL + L50_Pen + L95_Pen
  cat("NLL", NLL, "Length_NLL", Length_NLL, "CondAgeAtLengthNLL", CondAgeAtLengthNLL, "L50_Pen", L50_Pen,
      "L95_Pen",L95_Pen, "\n")
  cat("F ", 1/(1+exp(-params[1])), " other params ",exp(params[2:length(params)]), "\n")
  cat("", "\n")

  return(NLL)

}


#' Get NLL for marginal length composition data
#'
#' @keywords internal
#'
#' @param ObsRetCatchFreqAtLen observed frequencies of fish in length classes
#' @param ExpCatchPropInLenClass observed catch proportions in length classes
#'
#' @return negative log-likelihood (Length_NLL)
CalcNLLMargLengthComposition <- function(ObsRetCatchFreqAtLen, ExpRetCatchPropInLenClass) {

  Length_NLL = -sum(ObsRetCatchFreqAtLen * log(ExpRetCatchPropInLenClass + 1E-4))

  return(Length_NLL)
}

#' Get growth parameter values for catch curve
#'
#' Get growth parameter values for catch curve given specified catch curve, growth and selectivity model options
#'
#' @keywords internal
#' @param params estimated model parameters (varies, depending on growth curve type, selectivity type and catch curve type)
#' @param GrowthParams c(Linf, vbK, CVSizeAtAge) single sex von Bertalanffy, or data.frame(Linf=Linf, vbK=vbK, CVSizeAtAge=CVSizeAtAge),
#' both sexes von Bertalanffy, or c(y1, y2, a, b) single sex Schnute, or data.frame(y1=y1, y2=y2, a=a, b=b), both sexes Schnute
#' @param RefnceAges reference ages for Schnute function (set to NA if growth based on another function)
#' @param CatchCurveType 1=length based, 2=age and length based
#' @param SelectivityType 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#'
#' @return values of growth parameters given specified catch curve, growth and selectivity model option
GetGrowthAndSelectivityParams <- function(params, GrowthParams, RefnceAges, CatchCurveType, SelectivityType) {

  L50 = NA
  L95 = NA
  L50_ret = NA
  L95_ret = NA
  Linf = NA
  vbK = NA
  t1 = NA
  t2 = NA
  y1 = NA
  y2 = NA
  a = NA
  b = NA
  CVSizeAtAge = NA
  GrowthModelType = NA

  # get growth and selectivity parameters for length based catch curve model
  if (CatchCurveType == 1) { # length based catch curve
    if (SelectivityType == 2) { # logistic selectivity
      L50 = exp(params[2])
      L95 = L50 + exp(params[3])
      if (length(params)==5) { #
        L50_ret = exp(params[4])
        L95_ret = L50_ret + exp(params[5])
      }
    }
    if (GrowthCurveType == 1) { # von Bertalanffy
      if(is.vector(GrowthParams)) { # single or combined sex
        Linf = GrowthParams[1]
        vbK = GrowthParams[2]
        GrowthModelType = 1
      } else { #  separate sexes
        Linf = GrowthParams[,1]
        vbK = GrowthParams[,2]
        GrowthModelType = 2
      }
    } else { # Schnute
      if(is.vector(GrowthParams)) { # single or combined sex
        t1 = c(RefnceAges[1],RefnceAges[1])
        t2 = c(RefnceAges[2],RefnceAges[2])
        y1 = c(GrowthParams[1],GrowthParams[1])
        y2 = c(GrowthParams[2],GrowthParams[2])
        a = c(GrowthParams[3],GrowthParams[3])
        b = c(GrowthParams[4],GrowthParams[4])
        GrowthModelType = 3
      } else { #  separate sexes
        t1=RefnceAges[,1]; t2=RefnceAges[,2]
        y1=GrowthParams[,1]; y2=GrowthParams[,2]
        a=GrowthParams[,3]; b=GrowthParams[,4]
        GrowthModelType = 4
      }
    }
  }

  # get growth and selectivity parameters for length and age based catch curve model
  if (CatchCurveType == 2) { # age and length based catch curve
    if (GrowthCurveType == 1) { # von Bertalanffy
      if (SelectivityType == 1 & length(params)==4) { # selectivity vector input, single sex input
        Linf = exp(params[2])
        vbK = exp(params[3])
        CVSizeAtAge = exp(params[4])
        GrowthModelType = 1
      }
      if (SelectivityType == 2 & length(params)==6) { # logistic selectivity, single sex input
        L50 = exp(params[2])
        L95 = L50 + exp(params[3])
        Linf = exp(params[4])
        vbK = exp(params[5])
        CVSizeAtAge = exp(params[6])
        GrowthModelType = 1
      }
      if (SelectivityType == 1 & length(params)==6) { # selectivity vector input, separate sex input
        Linf = exp(c(params[2],params[3]))
        vbK = exp(c(params[4],params[5]))
        CVSizeAtAge = exp(params[6])
        GrowthModelType = 2
      }
      if (SelectivityType == 2 & length(params)==8) { # logistic selectivity, separate sex input
        L50 = exp(params[2])
        L95 = L50 + exp(params[3])
        Linf = exp(c(params[4],params[5]))
        vbK = exp(c(params[6],params[7]))
        CVSizeAtAge = exp(params[8])
        GrowthModelType = 2
      }
    } # von Bertalanffy

    if (GrowthCurveType == 2) { # Schnute
      # not yet implemeted!!!
      GrowthModelType = 3
      # GrowthModelType = 4
    }
  } # age and length based catch curve

  # calc Linf for Schnute growth curve, for penalty function
  if (GrowthModelType == 3) { # Schnute
    Linf = ((exp(a*t2)*y2^b-exp(a*t1)*y1^b)/(exp(a*t2)-exp(a*t1)))^(1/b)
  }
  if (GrowthModelType == 4) { # Schnute
    Linf1 = ((exp(a[1]*t2[1])*y2[1]^b[1]-exp(a[1]*t1[1])*y1[1]^b[1])/(exp(a[1]*t2[1])-exp(a[1]*t1[1])))^(1/b[1])
    Linf2 = ((exp(a[2]*t2[2])*y2[2]^b[2]-exp(a[2]*t1[2])*y1[2]^b[2])/(exp(a[2]*t2[2])-exp(a[2]*t1[2])))^(1/b[2])
    Linf = max(Linf1,Linf2)
  }

  result = list(L50 = L50,
                L95 = L95,
                L50_ret = L50_ret,
                L95_ret = L95_ret,
                Linf = Linf,
                vbK = vbK,
                t1 = t1,
                t2 = t2,
                y1 = y1,
                y2 = y2,
                a = a,
                b = b,
                CVSizeAtAge = CVSizeAtAge,
                GrowthModelType = GrowthModelType)

  return(result)
}


#' Get inputs for length transition matrices
#'
#' Get key inputs required to calculate length transition matrices, for models with either one or two sexes,
#' and using von Bertalanffy or Schnute growth model
#'
#' @keywords internal
#'
#' @param MaxAge maximum age considered in model
#' @param TimeStep model timestep (e.g. 1 = annual, 1/12 = monthly)
#' @param nLenCl number of length classes
#' @param midpt mid points of length classes
#' @param GrowthModelType 1=von Bertalanffy single sex, 2=von Bertalanffy both sexes, 3=Schnute single sex,
#' 4=Schnute both sexes
#' @param Linf asymptotic length (von Bertalanffy growth curve)
#' @param vbK growth coefficient (von Bertalanffy growth curve)
#' @param t1 first reference age (Schnute growth curve)
#' @param t2 second reference age (Schnute growth curve)
#' @param y1 length at first reference age (Schnute growth curve)
#' @param y2 length at second reference age (Schnute growth curve)
#' @param a Schnute growth curve parameter
#' @param b Schnute growth curve parameter
#' @param GrowthParams growth parameters (used for GrowthModelType 2 or 4)
#' @param RefnceAges Schnute reference ages (used for GrowthModelType 4)
#'
#' @return mean size at age from growth curve (MeanSizeAtAge), estimated length after one year for each length class
#' mid point, given growth curve (MeanEndingLength) and associated change in length after one year of growth (TimestepGrowthSizeInc)
GetInputsForLengthTransitionMatrices <- function(MaxAge, TimeStep, nLenCl, midpt, GrowthModelType, Linf, vbK,
                                                 t1, t2, y1, y2, a, b, GrowthParams, RefnceAges) {

  # growth calcs for single sex catch curve (or using combined growth curve for both sexes)
  Ages = seq(TimeStep,MaxAge,TimeStep)

  if (GrowthModelType == 1) { # von Bertalanffy, single sex
    MeanSizeAtAge = Linf * (1 - exp (-vbK * Ages))
    MeanEndingLength = midpt + (Linf - midpt) * (1 - exp(-vbK*TimeStep))
    TimestepGrowthSizeInc = MeanEndingLength-midpt # amount of annual growth with respect to initial length
  }
  if (GrowthModelType == 2) { # von Bertalanffy, separate sex
    MeanSizeAtAge <- data.frame(matrix(nrow = 2, ncol = length(Ages)))
    colnames(MeanSizeAtAge) <- Ages
    MeanEndingLength <- data.frame(matrix(nrow = 2, ncol = nLenCl))
    colnames(MeanEndingLength) <- midpt
    TimestepGrowthSizeInc = MeanEndingLength
    RecLenDist = MeanEndingLength
    for (i in 1:2) {
      MeanSizeAtAge[i,] = Linf[i] * (1 - exp (-vbK[i] * Ages))
      MeanEndingLength[i,] = midpt + (Linf[i] - midpt) * (1 - exp(-vbK[i]*TimeStep))
      TimestepGrowthSizeInc[i,] = MeanEndingLength[i,] - midpt # amount of annual growth with respect to initial length
    }
  }
  if (GrowthModelType == 3) { # Schnute, single sex
    MeanSizeAtAge = rep(0,length(Ages))
    k = 0
    for (t in Ages) {
      k = k + 1

      t1 = RefnceAges[1]
      t2 = RefnceAges[2]
      y1 = GrowthParams[1]
      y2 = GrowthParams[2]
      a = GrowthParams[3]
      b = GrowthParams[4]
      MeanSizeAtAge[k] = SchnuteGrowthfunction(t, t1, t2, y1, y2, a, b)
    }
    RefnceAgesForSex = RefnceAges
    GrowthParamsForSex = GrowthParams
    MeanEndingLength = CalcLengthAfterGrowthForTimetep(GrowthCurveType=2, TimeStep, GrowthParamsForSex, RefnceAgesForSex, midpt)
    TimestepGrowthSizeInc = MeanEndingLength-midpt # amount of annual growth with respect to initial length
  }
  if (GrowthModelType == 4) { # Schnute, separate sex
    MeanSizeAtAge <- data.frame(matrix(nrow = 2, ncol = length(Ages)))
    colnames(MeanSizeAtAge) <- Ages
    MeanEndingLength <- data.frame(matrix(nrow = 2, ncol = nLenCl))
    colnames(MeanEndingLength) <- midpt
    TimestepGrowthSizeInc = MeanEndingLength
    RecLenDist = MeanEndingLength
    for (i in 1:2) {
      t1=RefnceAges[i,1]
      t2=RefnceAges[i,2]
      y1=GrowthParams[i,1]
      y2=GrowthParams[i,2]
      a=GrowthParams[i,3]
      b=GrowthParams[i,4]
      k = 0
      for (t in Ages) {
        k = k + 1
        MeanSizeAtAge[i,k] =  SchnuteGrowthfunction(t, t1, t2, y1, y2, a, b)
      }
      GrowthParamsForSex = c(y1, y2, a, b)
      RefnceAgesForSex = c(t1, t2)
      MeanEndingLength[i,] = CalcLengthAfterGrowthForTimetep(GrowthCurveType=2, TimeStep, GrowthParamsForSex, RefnceAgesForSex, midpt)
      TimestepGrowthSizeInc[i,] = MeanEndingLength[i,] - midpt
    }
  }

  result = list(Ages = Ages,
                MeanSizeAtAge = MeanSizeAtAge,
                MeanEndingLength = MeanEndingLength,
                TimestepGrowthSizeInc = TimestepGrowthSizeInc)

  return(result)
}


#' Calculations for length and age and lengths-based catch curves
#'
#' This function undertakes a range of calculations required by the length and age and lengths-based catch curves
#'
#' @keywords internal
#'
#' @param params vector of model parameters in log space (params) to be estimated
#' @param GrowthCurveType 1 = von Bertalanffy, 2 = Schnute
#' @param GrowthParams c(Linf, vbK, CVSizeAtAge) single sex von Bertalanffy, or data.frame(Linf=Linf, vbK=vbK, CVSizeAtAge=CVSizeAtAge),
#' both sexes von Bertalanffy, or c(y1, y2, a, b) single sex Schnute, or data.frame(y1=y1, y2=y2, a=a, b=b), both sexes Schnute
#' @param RefnceAges reference ages for Schnute function (set to NA if growth based on another function)
#' @param CVSizeAtAge coefficient of variation (CV), common across all mean lengths at age
#' @param GrowthCurveType 1 = von Bertalanffy, 2 = Schnute
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param SelectivityType 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' @param DiscMort proportion of fish that die following catch and release
#' @param lbnd lower bounds of length classes
#' @param ubnd upper bounds of length classes
#' @param midpt mid points of length classes
#' @param SelectivityVec selectivity at length
#' @param MaxAge maximum age considered in model
#' @param NatMort natural mortality
#' @param TimeStep model timestep (e.g. 1 = annual, 1/12 = monthly)
#'
#' @return gear selectivity at length (SelAtLength), selectivity of landings at length (SelLandAtLength),
#'  mean size at age from growth curve (MeanSizeAtAge), retention probability at length (RetAtLength),
#' mid points of length classes (midpt), mean length after 1 year of growth given initial length of midpt
#' (MeanEndingLength), associated change in length after 1 year of growth (TimestepGrowthSizeInc), length distribution
#' of 1+ year old recruits (RecLenDist), retained catch proportions at each length for females, males and sexes combined
#' (ExpRetCatchPropInLenClass, ExpRetCatchPropInLenClass_Fem, ExpRetCatchPropInLenClass_Mal), expected retained catches at integer
#' age classes for females, males and sexes combined (RetCatchAtIntAge, RetCatchAtIntAge_Fem, RetCatchAtIntAge_Mal), retained
#' catch proportions at integer age classes for females, males and sexes combined (ExpRetCatchPropAtIntAge, ExpRetCatchPropAtIntAge_Fem, ExpRetCatchPropAtIntAge_Mal),
#' relative catches at length and age (Catch), expected catch proportions at length given integer age class for females, males and sexes combined
#' (ExpRetCatchPropLengthGivenIntAge, ExpRetCatchPropLengthGivenIntAge_Fem, ExpRetCatchPropLengthGivenIntAge_Mal), penalty values
#' for L50 and L95 selectivity parameters (L50_Pen, L95_Pen), growth model type, 1=single sex von Bertalanffy, 2=2 sex
#' von Bertalanffy, 3=1 sex Schnute, 4=2 sex Schnute (GrowthModelType), expected total (retained plus discarded) catches at
#' integer age and length (TotCatchAtIntAgeLen, TotCatchAtIntAgeLen_Fem, TotCatchAtIntAgeLen_Mal), retained catches at
#' length (RetCatchAtIntAgeLen, RetCatchAtIntAgeLen_Fem, RetCatchAtIntAgeLen_Mal), discarded catches at length (DiscCatchAtIntAgeLen,
#' DiscCatchAtIntAgeLen_Fem, DiscCatchAtIntAgeLen_Mal), total catches (discarded plus retained) at length (TotCatchAtLen,
#' TotCatchAtLen_Fem, TotCatchAtLen_Mal), retained catches at length (RetCatchAtLen, RetCatchAtLen_Fem, RetCatchAtLen_Mal),
#' discarded catches at length (DiscCatchAtLen, DiscCatchAtLen_Fem, RetCatchAtLen_Mal), fishing mortality
#' at length, for sexes combined or separate (FAtLen, FAtLen_Fem, FAtLen_Mal), total mortality at length, for sexes combined
#' or separate (ZAtLen, ZAtLen_Fem, ZAtLen_Mal), fishing mortality associated with landings at length, for sexes combined or
#' separate (FAtLen, FAtLen_Fem, FAtLen_Mal), fishing mortality associated with capture and retention at length, for sexes combined or
#' separate (FAtLenReten, FAtLenReten_Fem, FAtLenReten_Mal), fishing mortality associated with capture and discarding at
#' length, for sexes combined or separate (FAtLenDisc, FAtLenDisc_Fem, FAtLenDisc_Mal)
#' @export
#'
AgeAndLengthBasedCatchCurvesCalcs <- function (params, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge, CatchCurveType,
                                               MLL, SelectivityType, lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
{

  # get parameters for specified growth curve and catch curve type
  res = GetGrowthAndSelectivityParams(params, GrowthParams, RefnceAges, CatchCurveType, SelectivityType)

  L50 = res$L50; L95 = res$L95
  L50_ret = res$L50_ret; L95_ret = res$L95_ret

  Linf = res$Linf; vbK = res$vbK
  t1 = res$t1; t2 = res$t2; y1 = res$y1; y2 = res$y2; a = res$a; b = res$b
  GrowthModelType = res$GrowthModelType # 1 = single sex von Bertalanffy, 2 = 2 sex von Bertalanffy,
  # 3 = 1 sex Schnute, 4 = 2 sex Schnute

  if (CatchCurveType == 2) { # age and length based catch curve
    CVSizeAtAge = res$CVSizeAtAge
  }

  # ensure L95 < Linf
  L95_Pen = 0
  L50_Pen = 0
  if (SelectivityType == 2) { # logistic selectivity
    if (L95 > max(Linf)) {
      L95_Pen = (L95 - max(Linf))^2
      L95 = max(Linf)
    }
    # calculate L50 penalty
    if (L50 > (L95 - 2)) {
      L50_Pen = ((L95 - 2) - L50)^2
      L50 = L95 - 2
    }
  }

  # selectivity
  if (SelectivityType == 1) { # inputted as vector
    SelAtLength = SelectivityVec
  }
  if (SelectivityType == 2) { # logistic gear selectivity
    SelAtLength = CalcLogisticSelOrReten(L50, L95, midpt)
  }

  # retention
  if (!is.na(L50_ret)) { # calculate retention curve
    RetAtLength = CalcLogisticSelOrReten(L50_ret, L95_ret, midpt)
  } else {
    if (is.na(MLL)) { # specifying retention as 1, i.e. all fish caught are retained
      RetAtLength = rep(1,length(midpt))
    } else { # knife edge retention at MLL
      RetAtLength = rep(1E-20,length(midpt))
      RetAtLength[which(midpt>=MLL)]=1
    }
  }

  # get key inputs for length transition matrices
  nLenCl = length(midpt)
  Res = GetInputsForLengthTransitionMatrices(MaxAge, TimeStep, nLenCl, midpt, GrowthModelType, Linf, vbK, t1, t2, y1, y2, a, b, GrowthParams, RefnceAges)
  MeanSizeAtAge = Res$MeanSizeAtAge
  # cat("MeanSizeAtAge",MeanSizeAtAge,'\n')
  TimestepGrowthSizeInc = Res$TimestepGrowthSizeInc
  MeanEndingLength = Res$MeanEndingLength

  # get size distribution for juvenile recruits
  CV_BothSexes = c(CVSizeAtAge, CVSizeAtAge)
  RecLenDist = CalcSizeDistOfRecruits(MeanSizeAtAge, CV_BothSexes, lbnd, ubnd, midpt, nLenCl) # length distribution of 1+ recruits

  nTimeSteps = length(seq(TimeStep,MaxAge,TimeStep))
  Ages = seq(TimeStep, MaxAge, TimeStep)
  # cat("Ages",Ages,'\n')
  AgeCl = floor(round(Ages,6))
  # cat("AgeCl",AgeCl,'\n')
  nAgeCl = length(unique(AgeCl))
  MinAge = min(AgeCl)

  if (GrowthModelType ==  1 | GrowthModelType == 3) { # single sex
    # Calculate length transition matrix
    LTM = CalcLTM_cpp(TimestepGrowthSizeInc, CVSizeAtAge, lbnd, midpt, ubnd, nLenCl) # length-transition matrix

    # calculate catches
    InitRecNumber = 1.0
    CatchCurveResults = CalcCatches_AgeAndLengthBasedCatchCurves_cpp(params, NatMort, RecLenDist, InitRecNumber, MaxAge, TimeStep, nTimeSteps,
                                                                     nLenCl, midpt, RetAtLength, SelAtLength, DiscMort, LTM)

    # fishing mortality and total mortality at length
    FAtLen = CatchCurveResults$FAtLen
    ZAtLen = CatchCurveResults$ZAtLen
    FAtLenReten = CatchCurveResults$FAtLenReten
    FAtLenDisc = CatchCurveResults$FAtLenDisc
    FAtLen_Fem = NA; FAtLen_Mal = NA
    ZAtLen_Fem = NA; ZAtLen_Mal = NA
    FAtLenReten_Fem = NA; FAtLenReten_Mal = NA
    FAtLenDisc_Fem = NA; FAtLenDisc_Mal = NA

    # selectivity of landings (accounting for MLL, if exists)
    SelLandAtLength = CatchCurveResults$SelLandAtLength

    # catches at length
    TotCatchAtLen = CatchCurveResults$TotCatchAtLen # total catches (released + retained)
    RetCatchAtLen = CatchCurveResults$RetCatchAtLen
    DiscCatchAtLen = CatchCurveResults$DiscCatchAtLen
    ExpRetCatchPropInLenClass = RetCatchAtLen / sum(RetCatchAtLen)
    ExpDiscCatchPropInLenClass = DiscCatchAtLen / sum(DiscCatchAtLen)

    # catches at length and decimal age
    TotCatchAtDecAgeLen = CatchCurveResults$TotCatchAtDecAgeLen # total catches (released + retained)
    RetCatchAtDecAgeLen = CatchCurveResults$RetCatchAtDecAgeLen
    DiscCatchAtDecAgeLen = CatchCurveResults$DiscCatchAtDecAgeLen

    # catches at integer ages
    TotCatchAtIntAgeLen <- data.frame(matrix(0, nrow = nAgeCl, ncol = nLenCl))
    colnames(TotCatchAtIntAgeLen) <- midpt
    TotCatchAtIntAgeLen = as.matrix(TotCatchAtIntAgeLen)
    RetCatchAtIntAgeLen = TotCatchAtIntAgeLen
    DiscCatchAtIntAgeLen = TotCatchAtIntAgeLen

    for (t in 1:nTimeSteps) {
      if (AgeCl[1]<1) {
        x = AgeCl[t]+1
      } else {
        x = AgeCl[t]
      }

      TotCatchAtIntAgeLen[x,] = TotCatchAtIntAgeLen[x,] + CatchCurveResults$TotCatchAtDecAgeLen[t,]
      RetCatchAtIntAgeLen[x,] = RetCatchAtIntAgeLen[x,] + CatchCurveResults$RetCatchAtDecAgeLen[t,]
      DiscCatchAtIntAgeLen[x,] = DiscCatchAtIntAgeLen[x,] + CatchCurveResults$DiscCatchAtDecAgeLen[t,]
    }

    # expected prop at integer age
    RetCatchAtIntAge = as.vector(unlist(rowSums(RetCatchAtIntAgeLen)))
    # cat("RetCatchAtIntAge", RetCatchAtIntAge,'\n')
    # plot(MinAge:MaxAge,RetCatchAtIntAge)
    ExpRetCatchPropAtIntAge = RetCatchAtIntAge / sum(RetCatchAtIntAge)
    ExpRetCatchPropLengthGivenIntAge = CalcExpRetCatchPropLengthGivenIntAge(MinAge, MaxAge, midpt, RetCatchAtIntAge, RetCatchAtIntAgeLen)

    TotCatchAtLen_Fem = NA; TotCatchAtLen_Mal = NA
    RetCatchAtLen_Fem = NA; RetCatchAtLen_Mal = NA
    DiscCatchAtLen_Fem = NA; DiscCatchAtLen_Mal = NA
    TotCatchAtIntAgeLen_Fem = NA; TotCatchAtIntAgeLen_Mal = NA
    RetCatchAtIntAgeLen_Fem = NA; RetCatchAtIntAgeLen_Mal = NA
    DiscCatchAtIntAgeLen_Fem = NA; DiscCatchAtIntAgeLen_Mal = NA
    ExpRetCatchPropLengthGivenIntAge_Fem = NA; ExpRetCatchPropLengthGivenIntAge_Mal = NA
    ExpRetCatchPropInLenClass_Fem = NA; ExpRetCatchPropInLenClass_Mal = NA
    ExpDiscCatchPropInLenClass_Fem = NA; ExpDiscCatchPropInLenClass_Mal = NA
    RetCatchAtIntAge_Fem = NA; RetCatchAtIntAge_Mal = NA
    ExpRetCatchPropAtIntAge_Fem = NA; ExpRetCatchPropAtIntAge_Mal = NA

  } # single sex (or combined growth) catch curve


  if (GrowthModelType ==  2 | GrowthModelType == 4) { # separate sex

    TimestepGrowthSizeInc_Fem = as.vector(unlist(TimestepGrowthSizeInc[1,]))
    if (CatchCurveType == 1) CVSizeAtAge_Fem = CVSizeAtAge[1]
    if (CatchCurveType == 2) CVSizeAtAge_Fem = CVSizeAtAge
    LTM_Fem = CalcLTM_cpp(TimestepGrowthSizeInc_Fem, CVSizeAtAge_Fem, lbnd, midpt, ubnd, nLenCl) # length-transition matrix - females

    TimestepGrowthSizeInc_Mal = as.vector(unlist(TimestepGrowthSizeInc[2,]))
    if (CatchCurveType == 1) CVSizeAtAge_Mal = CVSizeAtAge[2]
    if (CatchCurveType == 2) CVSizeAtAge_Mal = CVSizeAtAge
    LTM_Mal = CalcLTM_cpp(TimestepGrowthSizeInc_Mal, CVSizeAtAge_Mal, lbnd, midpt, ubnd, nLenCl) # length-transition matrix - males

    RecLenDist_Fem = as.vector(unlist(RecLenDist[1,]))
    InitRecNumber = 0.5
    CatchCurveResults_Fem = CalcCatches_AgeAndLengthBasedCatchCurves_cpp(params, NatMort, RecLenDist_Fem, InitRecNumber, MaxAge, TimeStep, nTimeSteps,
                                                                         nLenCl, midpt, RetAtLength, SelAtLength, DiscMort, LTM_Fem)

    RecLenDist_Mal = as.vector(unlist(RecLenDist[2,]))
    InitRecNumber = 0.5
    CatchCurveResults_Mal = CalcCatches_AgeAndLengthBasedCatchCurves_cpp(params, NatMort, RecLenDist_Mal, InitRecNumber, MaxAge, TimeStep, nTimeSteps,
                                                                         nLenCl, midpt, RetAtLength, SelAtLength, DiscMort, LTM_Mal)


    # fishing mortality and total mortality at length
    FAtLen_Fem = CatchCurveResults_Fem$FAtLen
    ZAtLen_Fem = CatchCurveResults_Fem$ZAtLen
    FAtLenReten_Fem = CatchCurveResults_Fem$FAtLenReten
    FAtLenDisc_Fem = CatchCurveResults_Fem$FAtLenDisc
    FAtLen_Mal = CatchCurveResults_Mal$FAtLen
    ZAtLen_Mal = CatchCurveResults_Mal$ZAtLen
    FAtLenReten_Mal = CatchCurveResults_Mal$FAtLenReten
    FAtLenDisc_Mal = CatchCurveResults_Mal$FAtLenDisc
    FAtLen = NA
    ZAtLen = NA
    FAtLenReten = NA
    FAtLenDisc = NA

    # selectivity of landings (accounting for MLL, if exists)
    SelLandAtLength = CatchCurveResults_Fem$SelLandAtLength

    # catches at length
    TotCatchAtLen_Fem = CatchCurveResults_Fem$TotCatchAtLen # total catches (released + retained)
    TotCatchAtLen_Mal = CatchCurveResults_Mal$TotCatchAtLen
    TotCatchAtLen = TotCatchAtLen_Fem + TotCatchAtLen_Mal

    # expected prop at length
    RetCatchAtLen_Fem = CatchCurveResults_Fem$RetCatchAtLen
    RetCatchAtLen_Mal = CatchCurveResults_Mal$RetCatchAtLen
    RetCatchAtLen = RetCatchAtLen_Fem + RetCatchAtLen_Mal
    ExpRetCatchPropInLenClass_Fem = RetCatchAtLen_Fem / (sum(RetCatchAtLen_Fem) + sum(RetCatchAtLen_Mal)) # prop (with respect to all fish)
    ExpRetCatchPropInLenClass_Mal = RetCatchAtLen_Mal / (sum(RetCatchAtLen_Fem) + sum(RetCatchAtLen_Mal)) # prop (with respect to all fish)
    ExpRetCatchPropInLenClass = RetCatchAtLen / sum(RetCatchAtLen)

    # discarded catch at length
    DiscCatchAtLen_Fem = CatchCurveResults_Fem$DiscCatchAtLen
    DiscCatchAtLen_Mal = CatchCurveResults_Mal$DiscCatchAtLen
    DiscCatchAtLen = DiscCatchAtLen_Fem + DiscCatchAtLen_Mal
    ExpDiscCatchPropInLenClass_Fem = DiscCatchAtLen_Fem / (sum(DiscCatchAtLen_Fem) + sum(DiscCatchAtLen_Mal)) # prop (with respect to all fish)
    ExpDiscCatchPropInLenClass_Mal = DiscCatchAtLen_Mal / (sum(DiscCatchAtLen_Fem) + sum(DiscCatchAtLen_Mal)) # prop (with respect to all fish)
    ExpDiscCatchPropInLenClass = DiscCatchAtLen / sum(DiscCatchAtLen)

    # convert to integer ages
    TotCatchAtIntAgeLen_Fem <- data.frame(matrix(0, nrow = nAgeCl, ncol = nLenCl))
    colnames(TotCatchAtIntAgeLen_Fem) <- midpt
    TotCatchAtIntAgeLen_Fem = as.matrix(TotCatchAtIntAgeLen_Fem)
    TotCatchAtIntAgeLen_Mal = TotCatchAtIntAgeLen_Fem
    RetCatchAtIntAgeLen_Fem = TotCatchAtIntAgeLen_Fem
    RetCatchAtIntAgeLen_Mal = TotCatchAtIntAgeLen_Fem
    DiscCatchAtIntAgeLen_Fem = TotCatchAtIntAgeLen_Fem
    DiscCatchAtIntAgeLen_Mal = TotCatchAtIntAgeLen_Fem
    RetCatchAtIntAge_Fem = rep(0,nAgeCl)
    RetCatchAtIntAge_Mal = rep(0,nAgeCl)

    for (t in 1:nTimeSteps) {
      if (AgeCl[1]<1) {
        x = AgeCl[t]+1
      } else {
        x = AgeCl[t]
      }
      TotCatchAtIntAgeLen_Fem[x,] = TotCatchAtIntAgeLen_Fem[x,] + CatchCurveResults_Fem$TotCatchAtDecAgeLen[t,]
      TotCatchAtIntAgeLen_Mal[x,] = TotCatchAtIntAgeLen_Mal[x,] + CatchCurveResults_Mal$TotCatchAtDecAgeLen[t,]
      RetCatchAtIntAgeLen_Fem[x,] = RetCatchAtIntAgeLen_Fem[x,] + CatchCurveResults_Fem$RetCatchAtDecAgeLen[t,]
      RetCatchAtIntAgeLen_Mal[x,] = RetCatchAtIntAgeLen_Mal[x,] + CatchCurveResults_Mal$RetCatchAtDecAgeLen[t,]
      DiscCatchAtIntAgeLen_Fem[x,] = DiscCatchAtIntAgeLen_Fem[x,] + CatchCurveResults_Fem$DiscCatchAtDecAgeLen[t,]
      DiscCatchAtIntAgeLen_Mal[x,] = DiscCatchAtIntAgeLen_Mal[x,] + CatchCurveResults_Mal$DiscCatchAtDecAgeLen[t,]
      RetCatchAtIntAge_Fem[x] = RetCatchAtIntAge_Fem[x] + sum(RetCatchAtIntAgeLen_Fem[x,])
      RetCatchAtIntAge_Mal[x] = RetCatchAtIntAge_Mal[x] + sum(RetCatchAtIntAgeLen_Mal[x,])
    }
    TotCatchAtIntAgeLen = TotCatchAtIntAgeLen_Fem + TotCatchAtIntAgeLen_Mal
    RetCatchAtIntAgeLen = RetCatchAtIntAgeLen_Fem + RetCatchAtIntAgeLen_Mal
    DiscCatchAtIntAgeLen = DiscCatchAtIntAgeLen_Fem + DiscCatchAtIntAgeLen_Mal

    # expected retained catch prop at
    ExpRetCatchPropAtIntAge_Fem = RetCatchAtIntAge_Fem / sum(RetCatchAtIntAge_Fem)
    ExpRetCatchPropAtIntAge_Mal = RetCatchAtIntAge_Mal / sum(RetCatchAtIntAge_Mal)
    RetCatchAtIntAge = RetCatchAtIntAge_Fem + RetCatchAtIntAge_Mal
    ExpRetCatchPropAtIntAge = RetCatchAtIntAge / sum(RetCatchAtIntAge)

    # expected prop at age, given length
    ExpRetCatchPropLengthGivenIntAge_Fem = CalcExpRetCatchPropLengthGivenIntAge(MinAge, MaxAge, midpt, RetCatchAtIntAge_Fem, RetCatchAtIntAgeLen_Fem)
    ExpRetCatchPropLengthGivenIntAge_Mal = CalcExpRetCatchPropLengthGivenIntAge(MinAge, MaxAge, midpt, RetCatchAtIntAge_Mal, RetCatchAtIntAgeLen_Mal)
    ExpRetCatchPropLengthGivenIntAge=NA
  }


  Results = list(SelAtLength = SelAtLength,
                 RetAtLength = RetAtLength,
                 SelLandAtLength = SelLandAtLength,
                 MeanSizeAtAge = MeanSizeAtAge,
                 midpt = midpt,
                 MeanEndingLength = MeanEndingLength,
                 TimestepGrowthSizeInc = TimestepGrowthSizeInc,
                 RecLenDist = RecLenDist,
                 ExpRetCatchPropInLenClass = ExpRetCatchPropInLenClass,
                 ExpRetCatchPropInLenClass_Fem = ExpRetCatchPropInLenClass_Fem,
                 ExpRetCatchPropInLenClass_Mal = ExpRetCatchPropInLenClass_Mal,
                 ExpDiscCatchPropInLenClass = ExpDiscCatchPropInLenClass,
                 ExpDiscCatchPropInLenClass_Fem = ExpDiscCatchPropInLenClass_Fem,
                 ExpDiscCatchPropInLenClass_Mal = ExpDiscCatchPropInLenClass_Mal,
                 RetCatchAtIntAge = RetCatchAtIntAge,
                 RetCatchAtIntAge_Fem = RetCatchAtIntAge_Fem,
                 RetCatchAtIntAge_Mal = RetCatchAtIntAge_Mal,
                 ExpRetCatchPropAtIntAge = ExpRetCatchPropAtIntAge,
                 ExpRetCatchPropAtIntAge_Fem = ExpRetCatchPropAtIntAge_Fem,
                 ExpRetCatchPropAtIntAge_Mal = ExpRetCatchPropAtIntAge_Mal,
                 ExpRetCatchPropLengthGivenIntAge = ExpRetCatchPropLengthGivenIntAge,
                 ExpRetCatchPropLengthGivenIntAge_Fem = ExpRetCatchPropLengthGivenIntAge_Fem,
                 ExpRetCatchPropLengthGivenIntAge_Mal = ExpRetCatchPropLengthGivenIntAge_Mal,
                 L50_Pen = L50_Pen,
                 L95_Pen = L95_Pen,
                 GrowthModelType = GrowthModelType,
                 TotCatchAtIntAgeLen = TotCatchAtIntAgeLen,
                 TotCatchAtIntAgeLen_Fem = TotCatchAtIntAgeLen_Fem,
                 TotCatchAtIntAgeLen_Mal = TotCatchAtIntAgeLen_Mal,
                 RetCatchAtIntAgeLen = RetCatchAtIntAgeLen,
                 RetCatchAtIntAgeLen_Fem = RetCatchAtIntAgeLen_Fem,
                 RetCatchAtIntAgeLen_Mal = RetCatchAtIntAgeLen_Mal,
                 DiscCatchAtIntAgeLen = DiscCatchAtIntAgeLen,
                 DiscCatchAtIntAgeLen_Fem = DiscCatchAtIntAgeLen_Fem,
                 DiscCatchAtIntAgeLen_Mal = DiscCatchAtIntAgeLen_Mal,
                 TotCatchAtLen = TotCatchAtLen,
                 TotCatchAtLen_Fem = TotCatchAtLen_Fem,
                 TotCatchAtLen_Mal = TotCatchAtLen_Mal,
                 RetCatchAtLen = RetCatchAtLen,
                 RetCatchAtLen_Fem = RetCatchAtLen_Fem,
                 RetCatchAtLen_Mal = RetCatchAtLen_Mal,
                 DiscCatchAtLen = DiscCatchAtLen,
                 DiscCatchAtLen_Fem = DiscCatchAtLen_Fem,
                 RetCatchAtLen_Mal = RetCatchAtLen_Mal,
                 FAtLen = FAtLen,
                 FAtLen_Fem = FAtLen_Fem,
                 FAtLen_Mal = FAtLen_Mal,
                 ZAtLen = ZAtLen,
                 ZAtLen_Fem = ZAtLen_Fem,
                 ZAtLen_Mal = ZAtLen_Mal,
                 FAtLenReten = FAtLenReten,
                 FAtLenReten_Fem = FAtLenReten_Fem,
                 FAtLenReten_Mal = FAtLenReten_Mal,
                 FAtLenDisc = FAtLenDisc,
                 FAtLenDisc_Fem = FAtLenDisc_Fem,
                 FAtLenDisc_Mal = FAtLenDisc_Mal)

  return(Results)
}


#' Get results for a fitted age and length-based catch curve
#'
#' This function fits a length-based catch curve with length-based selectivity, based on
#' either 1) specified length-based selectivity inputted as a vector or 2) estimated asymptotic logistic
#' selectivity curve, and estimated von Bertalalnffy growth curve. The model is fitted to a
#' sample of fish length and age data, by minimising the overall negative log-likelihood,
#' including the NLL associated with the marginal length composition and a conditional age at length NLL,
#' given the parameters (selectivity, growth and mortality) and data, using nlminb.
#' It provides various statistical outputs in include convergence statistics, parameter estimates
#' and associated 95 percent confidence limits and associated variance-covariance matrix, calculated using
#' the MASS package.
#'
#' @param params vector of model parameters (in log space) to be estimated
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param SelectivityType 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' @param ObsRetCatchFreqAtLen observed frequencies in length classes
#' @param ObsRetCatchFreqAtLengthAndAge observed frequencies in length and age classes
#' @param lbnd lower bounds of length classes
#' @param ubnd upper bounds of length classes
#' @param midpt mid points of length classes
#' @param SelectivityVec selectivity at length
#' @param DiscMort proportion of fish that die following capture and release
#' @param MaxAge maximum age considered in model
#' @param NatMort natural mortality
#' @param TimeStep model timestep (e.g. 1 = annual, 1/12 = monthly)
#'
#' @return negative log-likelihood (nll), nlminb convergence diagnostic (convergence)
#' sample size (SampleSize), growth parameter estimates with lower and upper 95 percent
#' confidence limits (ParamEst), point estimates for estimated parameters (params)
#' and variance-covariance matrix (vcov.Params), gear selectivity at length (SelAtLength), retention
#' at length (RetAtLength), selectivity of landings at length (SelLandAtLength),
#' growth curve (MeanSizeAtAge), midpoint of each length class (midpt), mean length after 1 year from
#' growth curve, given initial length (MeanEndingLength), mean change in length after 1 year,
#' from initial length - note, assuming normal a distribution allows for possibility of negative growth
#' if above asyptotic length (TimestepGrowthSizeInc), length distribution of 1+ year old recruits (RecLenDist),
#' expected retained catches at length and integer age class (RetCatchAtIntAgeLen), expected catches at length (ExpCatchAtLen),
#' catch proportions at length (ExpCatchPropInLenClass), expected retained catches at integer age classes (RetCatchAtIntAge),
#' catch proportions at integer age class of females, males and sexes combined (ExpRetCatchPropAtIntAge, ExpRetCatchPropAtIntAge_Fem,
#' ExpRetCatchPropAtIntAge_Mal), catch proportions at length and integer age class for females, males and sexes combined
#' (ExpRetCatchPropLengthGivenIntAge, ExpRetCatchPropLengthGivenIntAge_Fem, ExpRetCatchPropLengthGivenIntAge_Mal),
#' observed catch frequency at length data (ObsRetCatchFreqAtLen), observed catch frequency at length and age
#' (ObsRetCatchFreqAtLengthAndAge), approximate CV at length corresponding to maximum age (CV_LenAtMaxAge, FemCV_LenAtMaxAge,
#' MalCV_LenAtMaxAge)
#'
#' @examples
#' # Simulate data
#' set.seed(123)
#' SampleSize=5000
#' MaxAge = 26
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' MinAge = floor(TimeStep)
#' nAgeCl = length(MinAge:MaxAge)
#' nTimeSteps = length(seq(TimeStep,MaxAge,TimeStep))
#' NatMort = 4.22/MaxAge
#' FishMort = 0.2
#' MaxLen = 1100
#' LenInc = 10
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified is knife-edged at MLL
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0
#' # # single sex, von Bertalanffy
#' # GrowthCurveType = 1 # 1 = von Bertalanffy
#' # Linf = 800
#' # vbK = 0.2
#' # CVSizeAtAge = 0.08
#' # GrowthParams = c(Linf, vbK, CVSizeAtAge)
#' # RefnceAges = NA
#' # 2 sexes, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy
#' Linf = c(700,850)
#' vbK = c(0.25,0.2)
#' CVSizeAtAge = c(0.05,0.05)
#' RefnceAges = NA
#' GrowthParams = data.frame(Linf=Linf, vbK=vbK, CVSizeAtAge=CVSizeAtAge)
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' lbnd=Res$lbnd
#' midpt=Res$midpt
#' ubnd=Res$ubnd
#' # # get data - 1 sex (or combined sexes)
#' # ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen # 1 sex
#' # ObsRetCatchFreqAtLengthAndAge = as.matrix(Res$ObsCatchFreqAtLengthAndDecAge) # 1 sex
#' # get data - 2 sexes
#' ObsRetCatchFreqAtLen <- data.frame(matrix(nrow = 2, ncol = length(midpt))) # 2 sex
#' colnames(ObsRetCatchFreqAtLen) <- midpt
#' ObsRetCatchFreqAtLen[1,] = Res$ObsRetCatchFreqAtLen_Fem
#' ObsRetCatchFreqAtLen[2,] = Res$ObsRetCatchFreqAtLen_Mal
#' ObsRetCatchFreqAtLengthAndAge = array(c(unlist(Res$ObsCatchFreqAtLengthAndDecAge_Fem), unlist(Res$ObsCatchFreqAtLengthAndDecAge_Mal)),
#'                                    c(nTimeSteps, length(midpt), 2), dimnames=list(rownames(Res$ObsCatchFreqAtLengthAndDecAge_Fem),
#'                                                                                   colnames(Res$ObsCatchFreqAtLengthAndDecAge_Fem)))
#' # # get params - 1 sex
#' # InitFishMort = 0.3 # specify starting parameters
#' # InitL50 = 320
#' # InitDelta = 50 # L95-L50
#' # InitLinf = 800
#' # InitvbK = 0.2
#' # InitCVSizeAtAge = 0.05
#' # get params - 2 sexes
#' InitFishMort = 0.3 # specify starting parameters
#' InitL50 = 320
#' InitDelta = 50 # L95-L50
#' InitLinf = c(800,800)
#' InitvbK = c(0.25,0.25)
#' InitCVSizeAtAge = 0.05
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' params = c(InitFishMort_logit, log(c(InitL50, InitDelta, InitLinf, InitvbK, InitCVSizeAtAge)))
#' FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
#' # # Example with specified selectivity vector
#' # # Simulate data
#' # set.seed(123)
#' # SampleSize=5000
#' # MaxAge = 26
#' # TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' # MinAge = floor(TimeStep)
#' # nAgeCl = length(MinAge:MaxAge)
#' # nTimeSteps = length(seq(TimeStep,MaxAge,TimeStep))
#' # NatMort = 4.22/MaxAge
#' # FishMort = 0.2
#' # MaxLen = 1100
#' # LenInc = 10
#' # MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA, otherwise retention is knife-edged at MLL
#' # SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' # lbnd = seq(0,MaxLen - LenInc, LenInc)
#' # midpt = lbnd + (LenInc/2)
#' # SelectivityVec = 1 / (1 + exp(-log(19)*(midpt-400)/(500-400)))
#' # SelParams = c(NA, NA) # L50, L95 for gear selectivity
#' # RetenParams = c(NA, NA) # L50, L95 for retention
#' # DiscMort = 0
#' # # single sex, von Bertalanffy
#' # GrowthCurveType = 1 # 1 = von Bertalanffy
#' # Linf = 800
#' # vbK = 0.2
#' # CVSizeAtAge = 0.08
#' # GrowthParams = c(Linf, vbK, CVSizeAtAge)
#' # RefnceAges = NA
#' # Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#' #                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # lbnd=Res$lbnd
#' # midpt=Res$midpt
#' # ubnd=Res$ubnd
#' # # get data - 1 sex (or combined sexes)
#' # ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen # 1 sex
#' # ObsRetCatchFreqAtLengthAndAge = as.matrix(Res$ObsCatchFreqAtLengthAndDecAge) # 1 sex
#' # # get params - 1 sex
#' # InitFishMort = 0.3 # specify starting parameters
#' # InitLinf = 800
#' # InitvbK = 0.2
#' # InitCVSizeAtAge = 0.05
#' # InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' # params = c(InitFishMort_logit, log(c(InitLinf, InitvbK, InitCVSizeAtAge)))
#' # FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#' #                                                lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
#' @export
GetAgeAndLengthBasedCatchCurveResults <- function (params, MLL, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                                   lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
{

  nlmb <- nlminb(params, CalcObjFunc_AgeAndLengthBasedCatchCurve,
                 gradient = NULL, hessian = TRUE, control=list(eval.max=1000, iter.max=1000))

  (hess.out = optimHess(nlmb$par, CalcObjFunc_AgeAndLengthBasedCatchCurve))
  (vcov.Params = solve(hess.out))
  (ses = sqrt(diag(vcov.Params)))

  # logit space
  temp = c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1])
  # inverse logit transformed value
  EstFMort = 1/(1+exp(-temp));


  if (SelectivityType == 1 & length(params)==4) { # selectivity vector input, single sex input
    EstLinf = exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
    Estk = exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
    EstCV = exp(c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) *  ses[4]))
    ParamEst = t(data.frame(FMort = round(EstFMort, 3), Linf = round(EstLinf, 3),
                            vbK = round(Estk, 3), CV = round(EstCV, 3)))
  }
  if (SelectivityType == 2 & length(params)==6) { # logistic selectivity, single sex input
    EstL50 = exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
    EstDelta = exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
    EstLinf = exp(c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
    Estk = exp(c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5]))
    EstCV = exp(c(nlmb$par[6], nlmb$par[6] + c(-1.96, 1.96) * ses[6]))
    ParamEst = t(data.frame(FMort = round(EstFMort, 3), SelL50 = round(EstL50, 3), Delta = round(EstDelta, 3),
                            Linf = round(EstLinf, 3), vbK = round(Estk, 3), CV = round(EstCV, 3)))
  }
  if (SelectivityType == 1 & length(params)==6) {   # selectivity vector input, separate sex input
    EstLinf_F = exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
    EstLinf_M = exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
    Estk_F = exp(c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
    Estk_M = exp(c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5]))
    EstCV = exp(c(nlmb$par[6], nlmb$par[6] + c(-1.96, 1.96) *  ses[6]))
    ParamEst = t(data.frame(FMort = round(EstFMort, 3), Linf_F = round(EstLinf_F, 3), Linf_M = round(EstLinf_M, 3),
                            vbK_F = round(Estk_F, 3), vbK_M = round(Estk_M, 3), CV = round(EstCV, 3)))
  }
  if (SelectivityType == 2 & length(params)==8) { # selectivity vector input, separate sex input # logistic selectivity, separate sex input
    EstL50 = exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
    EstDelta = exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
    EstLinf_F = exp(c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
    EstLinf_M = exp(c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5]))
    Estk_F = exp(c(nlmb$par[6], nlmb$par[6] + c(-1.96, 1.96) * ses[6]))
    Estk_M = exp(c(nlmb$par[7], nlmb$par[7] + c(-1.96, 1.96) * ses[7]))
    EstCV = exp(c(nlmb$par[8], nlmb$par[8] + c(-1.96, 1.96) * ses[8]))
    ParamEst = t(data.frame(FMort = round(EstFMort, 3), SelL50 = round(EstL50, 3), Delta = round(EstDelta, 3),
                            Linf_F = round(EstLinf_F, 3), Linf_M = round(EstLinf_M, 3),
                            vbK_F = round(Estk_F, 3), vbK_M = round(Estk_M, 3), CV = round(EstCV, 3)))
  }
  colnames(ParamEst) = c("Estimate", "lw_95%CL", "up_95%CL")

  # store some diagnostic outputs from model
  CatchCurveType=2
  GrowthCurveType=1 # von Bertalanffy
  params = nlmb$par
  GrowthParams = NA
  RefnceAges= NA
  CVSizeAtAge = NA
  res = AgeAndLengthBasedCatchCurvesCalcs(params, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge, CatchCurveType,
                                          MLL, SelectivityType, lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
  SampleSize = sum(ObsRetCatchFreqAtLen)
  nll = nlmb$objective
  convergence = nlmb$convergence

  ExpRetCatchPropLengthGivenIntAge = res$ExpRetCatchPropLengthGivenIntAge
  ExpRetCatchPropLengthGivenIntAge_Fem = res$ExpRetCatchPropLengthGivenIntAge_Fem
  ExpRetCatchPropLengthGivenIntAge_Mal = res$ExpRetCatchPropLengthGivenIntAge_Mal
  GrowthModelType = res$GrowthModelType

  # calculate approximate cv for lengths at max age (growth diagnostic)
  # single / combined sex
  if (GrowthModelType ==  1 | GrowthModelType == 3) { # single sex
    if (TimeStep==1) {
      Probs = as.vector(unlist(res$ExpRetCatchPropLengthGivenIntAge[MaxAge,]))
      RandFreq=rep(midpt,(as.vector(rmultinom(1,1000,Probs))))
      MeanRandFreq=mean(RandFreq)
      sdRandFreq=sd(RandFreq)
      CV_LenAtMaxAge=sdRandFreq/MeanRandFreq
    } else {
      CV_LenAtMaxAge = NA
    }
    FemCV_LenAtMaxAge = NA
    MalCV_LenAtMaxAge = NA
  }

  # females
  if (GrowthModelType ==  2 | GrowthModelType == 4) { # separate sex
    if (TimeStep==1) {
      FemProbs = as.vector(unlist(res$ExpRetCatchPropLengthGivenIntAge_Fem[MaxAge,]))
      FemRandFreq=rep(midpt,(as.vector(rmultinom(1,1000,FemProbs))))
      FemMeanRandFreq=mean(FemRandFreq)
      FemsdRandFreq=sd(FemRandFreq)
      FemCV_LenAtMaxAge=FemsdRandFreq/FemMeanRandFreq

      # males
      MalProbs = as.vector(unlist(res$ExpRetCatchPropLengthGivenIntAge_Mal[MaxAge,]))
      MalRandFreq=rep(midpt,(as.vector(rmultinom(1,1000,MalProbs))))
      MalMeanRandFreq=mean(MalRandFreq)
      MalsdRandFreq=sd(MalRandFreq)
      MalCV_LenAtMaxAge=MalsdRandFreq/MalMeanRandFreq
    } else {
      FemCV_LenAtMaxAge = NA
      MalCV_LenAtMaxAge = NA
    }
    CV_LenAtMaxAge = NA
  }


  Results = list(nll = nll,
                 convergence = convergence,
                 SampleSize = SampleSize,
                 ParamEst = ParamEst,
                 params = nlmb$par,
                 vcov.Params = vcov.Params,
                 SelAtLength=res$SelAtLength,
                 RetAtLength = res$RetAtLength,
                 SelLandAtLength = res$SelLandAtLength,
                 MeanSizeAtAge=res$MeanSizeAtAge,
                 midpt=res$midpt,
                 MeanEndingLength=res$MeanEndingLength,
                 TimestepGrowthSizeInc=res$TimestepGrowthSizeInc,
                 RecLenDist=res$RecLenDist,
                 RetCatchAtIntAgeLen=res$RetCatchAtIntAgeLen,
                 ExpCatchPropInLenClass=res$ExpCatchPropInLenClass,
                 RetCatchAtIntAge=res$RetCatchAtIntAge,
                 ExpRetCatchPropAtIntAge=res$ExpRetCatchPropAtIntAge,
                 ExpRetCatchPropAtIntAge_Fem=res$ExpRetCatchPropAtIntAge_Fem,
                 ExpRetCatchPropAtIntAge_Mal=res$ExpRetCatchPropAtIntAge_Mal,
                 ExpRetCatchPropLengthGivenIntAge=ExpRetCatchPropLengthGivenIntAge,
                 ExpRetCatchPropLengthGivenIntAge_Fem=ExpRetCatchPropLengthGivenIntAge_Fem,
                 ExpRetCatchPropLengthGivenIntAge_Mal=ExpRetCatchPropLengthGivenIntAge_Mal,
                 ObsRetCatchFreqAtLen=ObsRetCatchFreqAtLen,
                 ObsRetCatchFreqAtLengthAndAge=ObsRetCatchFreqAtLengthAndAge,
                 CV_LenAtMaxAge=CV_LenAtMaxAge,
                 FemCV_LenAtMaxAge=FemCV_LenAtMaxAge,
                 MalCV_LenAtMaxAge=MalCV_LenAtMaxAge)

  return(Results)
}


#' Convert observed length and age data from decimal ages to integer age classes
#'
#' Convert observed length and age data from decimal ages to integer age classes, used in age and
#' length based catch curve analyses
#'
#' @keywords internal
#'
#' @param TimeStep model timestep
#' @param MaxAge maximum age considered by model
#' @param nLenCl number of length classes
#' @param ObsRetCatchFreqAtLengthAndAge observed decimal age and length frequency data
#'
#' @return ObsCatchFreqAtLengthAndIntAge
ConvertObsDataFromDecAgesToIntegerAges <- function(TimeStep, MaxAge, nLenCl, ObsRetCatchFreqAtLengthAndAge) {

  Ages = seq(TimeStep,MaxAge,TimeStep)
  nTimeSteps = length(Ages)
  AgeCl = floor(round(Ages,6))
  nAgeCl = length(unique(AgeCl))

  ObsCatchFreqAtLengthAndIntAge <- data.frame(matrix(0, nrow = nAgeCl, ncol = nLenCl))
  colnames(ObsCatchFreqAtLengthAndIntAge) <- midpt
  ObsCatchFreqAtLengthAndIntAge = as.matrix(ObsCatchFreqAtLengthAndIntAge)

  for (t in 1:nTimeSteps) {
    if (AgeCl[1]<1) {
      x = AgeCl[t]+1
    } else {
      x = AgeCl[t]
    }

    ObsCatchFreqAtLengthAndIntAge[x,] = ObsCatchFreqAtLengthAndIntAge[x,] + ObsRetCatchFreqAtLengthAndAge[t,]
  }

  return(ObsCatchFreqAtLengthAndIntAge)
}

#' Calculate expected retained catch proportions in length classes, within integer age classes
#'
#' Calculate expected retained catch proportions in length classes, within integer age classes, used in age and
#' length based catch curve analyses
#'
#' @keywords internal
#'
#' @param MinAge minimum age
#' @param MaxAge maximum age considered by model
#' @param RetCatchAtIntAge expected retained catches at integer age classes
#' @param RetCatchAtIntAgeLen expected retained catches at integer age classes and length classes
#'
#' @return ExpRetCatchPropLengthGivenIntAge
CalcExpRetCatchPropLengthGivenIntAge <- function(MinAge, MaxAge, midpt, RetCatchAtIntAge, RetCatchAtIntAgeLen) {

  nAgeCl = length(MinAge:MaxAge)
  nLenCl = length(midpt)
  ExpRetCatchPropLengthGivenIntAge <- data.frame(matrix(0, nrow = nAgeCl, ncol = nLenCl))
  colnames(ExpRetCatchPropLengthGivenIntAge) <- midpt
  ExpRetCatchPropLengthGivenIntAge = as.matrix(ExpRetCatchPropLengthGivenIntAge)

  x = 0
  for (t in seq(MinAge,MaxAge,1)) {
    x = x + 1
    if (RetCatchAtIntAge[x]==0) {
      ExpRetCatchPropLengthGivenIntAge[x,] = 0
    } else {
      ExpRetCatchPropLengthGivenIntAge[x,] = RetCatchAtIntAgeLen[x,] / RetCatchAtIntAge[x]
    }
  }

  return(ExpRetCatchPropLengthGivenIntAge)
}



#' Simulate some age frequency data
#'
#' Simulate age frequency data with specified selectivity and mortality parameters
#'
#' @param SampleSize required sample size
#' @param MinAge minimum age
#' @param MaxAge maximum age
#' @param SelA50 age at 50 percent selectivity
#' @param SelA95 age at 95 percent selectivity
#' @param NatMort natural mortality
#' @param FMort fully-selected fishing mortality
#'
#' @return Age classes (Ages), selectivity at age (SelAtAge), expected proportions at age (PropAtAge), observed
#' catch frequency at age (CatchSample)
#'
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
#' SampleSize = 1000 # required number of fish for age sample
#' Res=SimAgeFreqData(SampleSize, MinAge, MaxAge, SelA50, SelA95, NatMort, FMort)
#' ObsAgeFreq = unlist(as.vector(Res$CatchSample))
#' @export
SimAgeFreqData <- function(SampleSize, MinAge, MaxAge, SelA50, SelA95, NatMort, FMort) {

  Ages = MinAge:MaxAge
  SelAtAge = 1 / (1 + exp(-log(19) * (Ages - SelA50) / (SelA95 - SelA50)))
  FAtAge = SelAtAge * FMort
  ZAtAge = NatMort + FAtAge
  N = numeric(length(Ages))

  k=1
  N[k] = 1
  for (i in seq(MinAge+1,MaxAge,1)) {
    k=k+1
    if (i < MaxAge) {
      N[k] = N[k-1] * exp(-ZAtAge[k-1])
    } else {
      N[k] = N[k-1] * exp(-ZAtAge[k-1] / (1 - exp(-ZAtAge[k])))
    }
  }
  CatchAtAge = N * (FAtAge / ZAtAge) * (1 - exp(-ZAtAge)) # catch at age
  PropAtAge = CatchAtAge / sum(CatchAtAge)
  CatchSample = rmultinom(n=1, size=SampleSize, prob=PropAtAge)

  results = list(Ages=Ages,
                 SelAtAge=SelAtAge,
                 PropAtAge=PropAtAge,
                 CatchSample=CatchSample)
}


#' Simulate age and length frequency data
#'
#' Function simulates age and length frequency data with specified selectivity, mortality and growth parameters
#'
#' @param SampleSize required same size, single value generates same number of retained and discarded fish,
#' # can also specify vector, e.g. SampleSize=c(10000, 2000), produces 10000 retained fish, 2000 discarded fish
#' @param MaxAge maximum age
#' @param TimeStep model timestep (e.g. 1 = annual, 1/12 = monthly)
#' @param NatMort natural mortality
#' @param FishMort fully-selected fishing mortality
#' @param MaxLen maximum length
#' @param LenInc length class interval
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param SelectivityType 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' @param L50 length at 50 percent selectivity
#' @param L95 length at 95 percent selectivity
#' @param SelectivityVec selectivity at length, set to NA if selectivity parameters inputted
#' @param DiscMort Proportion of fish that die following to capture and release
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
#' @param GrowthParams c(Linf, vbK) single sex von Bertalanffy, or data.frame(Linf=Linf, vbK=vbK) separate sex von Bertalanffy,
#' c(y1, y2, a, b) single sex Schnute, or data.frame(y1=y1, y2=y2, a=a, b=b) separate sex Schnute
#' @param RefnceAges Schnute reference ages, either c(t1,t2) single sex, c(t1=t1, t2=t2) separate sex, or set to NA if using von Bertalanffy model
#' @param CVSizeAtAge coefficient of variation for mean lengths at age
#'
#' @return Observed catch frequency at age, ObsCatchFreqAtAge (Vector), Observed catch proportions at
#' age, ObsRelCatchAtAge (Vectors), Observed catch frequency at length, ObsRetCatchFreqAtLen (Vector),
#' observed catch proportions at length, ObsRelCatchAtLen (Vectors), Observed catch frequency at length
#' and age, ObsRetCatchFreqAtLengthAndAge (Matrix), age classes, ObsAgeCl (Vector), length class
#' mid points, ObslenClMidPt (Vector)
#' @examples
#' set.seed(123)
#' SampleSize=10000 # generate 10000 retained fish and 10000 discarded fish,
#' # can also specify vector, e.g. SampleSize=c(10000, 2000), produces 10000 retained fish, 2000 discarded fish
#' MaxAge = 20
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 0.2
#' FishMort = 0.2
#' MaxLen = 1000
#' LenInc = 20
#' midpt = seq(0,MaxLen - LenInc, LenInc) + (LenInc/2)
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified is knife-edged at MLL
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 500) # L50, L95 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95 for retention
#' DiscMort = 0.5
#' GrowthCurveType = 1 # 1 = von Bert, 2 = Schnute
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.08
#' GrowthParams = c(Linf, vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' @export
SimLenAndAgeFreqData <- function(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                                 SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge) {

  # Simulate some observed catch at length data and age and length data
  DecAges = seq(TimeStep,MaxAge,TimeStep)
  nTimeSteps = length(DecAges)
  lbnd = seq(0,MaxLen - LenInc, LenInc)
  ubnd = lbnd + LenInc
  midpt = lbnd + (LenInc/2)
  nLenCl = length(midpt)
  FishMort_logit = log(FishMort/(1-FishMort)) # logit transform (so F is always between 0 and 1)

  # get parameters for specified growth curve
  # For single sex/combined sex data simulation model, to reduce coding, applying 2-sex model with same parameters
  # for each sex, and then combining data later
  if (GrowthCurveType == 1) { # von Bertalanffy
    if (is.vector(GrowthParams)) { # params inputted as combined sex
      Linf = c(GrowthParams[1],GrowthParams[1])
      vbK = c(GrowthParams[2],GrowthParams[2])
      CVSizeAtAge = c(CVSizeAtAge, CVSizeAtAge)
    }
    if (is.data.frame(GrowthParams)) { # params inputted as separate sex
      Linf = GrowthParams[,1]
      vbK = GrowthParams[,2]
    }
    GrowthModelType = 2 # von Bert, sep. sex
  } # GrowthCurveType = 1

  if (GrowthCurveType == 2) { # Schnute
    if (is.vector(GrowthParams)) { # combined sex
      t1 = c(RefnceAges[1],RefnceAges[1])
      t2 = c(RefnceAges[2],RefnceAges[2])
      y1 = c(GrowthParams[1],GrowthParams[1])
      y2 = c(GrowthParams[2],GrowthParams[2])
      a = c(GrowthParams[3],GrowthParams[3])
      b = c(GrowthParams[4],GrowthParams[4])
      CVSizeAtAge = c(CVSizeAtAge, CVSizeAtAge)
      GrowthParamsForSex = c(y1, y2, a, b)
      RefnceAgesForSex = c(t1, t2)
      GrowthParams = data.frame(y1=y1, y2=y2, a=a, b=b)
      RefnceAges = data.frame(t1=t1,t2=t2)
    }
    if (is.data.frame(GrowthParams)) { # separate sex
      t1=RefnceAges[,1]; t2=RefnceAges[,2]
      y1=GrowthParams[,1]; y2=GrowthParams[,2]
      a=GrowthParams[,3]; b=GrowthParams[,4]
    }
    GrowthModelType = 4 # Schnute, sep. sex
  } # GrowthCurveType = 2

  # get key inputs for length transition matrices
  Res = GetInputsForLengthTransitionMatrices(MaxAge, TimeStep, nLenCl, midpt, GrowthModelType,
                                             Linf, vbK, t1, t2, y1, y2, a, b, GrowthParams, RefnceAges)

  MeanSizeAtAge = Res$MeanSizeAtAge
  TimestepGrowthSizeInc = Res$TimestepGrowthSizeInc
  MeanEndingLength = Res$MeanEndingLength

  # selectivity
  if (SelectivityType == 1) { # inputted as vector
    SelAtLength = SelectivityVec
  }
  if (SelectivityType == 2) { # logistic gear selectivity
    L50=SelParams[1]
    L95=L50 + SelParams[2]
    SelAtLength = CalcLogisticSelOrReten(L50, L95, midpt)
  }
  # cat("SelAtLength",SelAtLength,'\n')

  # retention
  if (!is.na(RetenParams[1])) { # calculate retention curve
    L50=RetenParams[1]
    L95=L50+RetenParams[2]
    RetAtLength = CalcLogisticSelOrReten(L50, L95, midpt)
  } else {
    if (is.na(MLL)) { # specifying retention as 1, i.e. all fish caught are retained
      RetAtLength = rep(1,length(midpt))
    } else { # knife edge retention at MLL
      RetAtLength = rep(1E-20,length(midpt))
      RetAtLength[which(midpt>=MLL)]=1
    }
  }
  # cat("RetAtLength",RetAtLength,'\n')

  # size distribution of 1+ recruits
  RecLenDist = CalcSizeDistOfRecruits(MeanSizeAtAge, CVSizeAtAge, lbnd, ubnd, midpt, nLenCl)

  TimestepGrowthSizeInc_Fem = as.vector(unlist(TimestepGrowthSizeInc[1,]))
  CVSizeAtAge_Fem = CVSizeAtAge[1]
  LTM_Fem = CalcLTM_cpp(TimestepGrowthSizeInc_Fem, CVSizeAtAge_Fem, lbnd, midpt, ubnd, nLenCl) # length-transition matrix - females

  TimestepGrowthSizeInc_Mal = as.vector(unlist(TimestepGrowthSizeInc[2,]))
  CVSizeAtAge_Mal = CVSizeAtAge[2]
  LTM_Mal = CalcLTM_cpp(TimestepGrowthSizeInc_Mal, CVSizeAtAge_Mal, lbnd, midpt, ubnd, nLenCl) # length-transition matrix - males

  RecLenDist_Fem = as.vector(unlist(RecLenDist[1,]))
  InitRecNumber = 0.5
  CatchCurveResults_Fem = CalcCatches_AgeAndLengthBasedCatchCurves_cpp(FishMort_logit, NatMort, RecLenDist_Fem, InitRecNumber, MaxAge, TimeStep, nTimeSteps,
                                                                       nLenCl, midpt, RetAtLength, SelAtLength, DiscMort, LTM_Fem)

  RecLenDist_Mal = as.vector(unlist(RecLenDist[2,]))
  InitRecNumber = 0.5
  CatchCurveResults_Mal = CalcCatches_AgeAndLengthBasedCatchCurves_cpp(FishMort_logit, NatMort, RecLenDist_Mal, InitRecNumber, MaxAge, TimeStep, nTimeSteps,
                                                                       nLenCl, midpt, RetAtLength, SelAtLength, DiscMort, LTM_Mal)

  # retained catch at length
  RetCatchAtLen_Fem = CatchCurveResults_Fem$RetCatchAtLen
  RetCatchAtLen_Mal = CatchCurveResults_Mal$RetCatchAtLen
  RetCatchAtLen = RetCatchAtLen_Fem + RetCatchAtLen_Mal

  # discarded catch at length
  DiscCatchAtLen_Fem = CatchCurveResults_Fem$DiscCatchAtLen
  DiscCatchAtLen_Mal = CatchCurveResults_Mal$DiscCatchAtLen
  DiscCatchAtLen = DiscCatchAtLen_Fem + DiscCatchAtLen_Mal

  # total catch at length (discarded pus retained catches)
  TotCatchAtLen_Fem = CatchCurveResults_Fem$TotCatchAtLen
  TotCatchAtLen_Mal = CatchCurveResults_Mal$TotCatchAtLen
  TotCatchAtLen = TotCatchAtLen_Fem + TotCatchAtLen_Mal

  # retained catch at decimal and and length
  SelLandAtLength = CatchCurveResults_Fem$SelLandAtLength # selectivity of landings
  RetCatchAtDecAgeLen_Fem = as.matrix(unlist(CatchCurveResults_Fem$RetCatchAtDecAgeLen)) # retained catches, decimal ages and lengths
  RetCatchAtDecAgeLen_Mal = as.matrix(unlist(CatchCurveResults_Mal$RetCatchAtDecAgeLen))
  RetCatchAtDecAgeLen = RetCatchAtDecAgeLen_Fem + RetCatchAtDecAgeLen_Mal # single sex or combined sexes

  # discarded catch at decimal and and length
  DiscCatchAtDecAgeLen_Fem = as.matrix(unlist(CatchCurveResults_Fem$DiscCatchAtDecAgeLen)) # retained catches, decimal ages and lengths
  DiscCatchAtDecAgeLen_Mal = as.matrix(unlist(CatchCurveResults_Mal$DiscCatchAtDecAgeLen))
  DiscCatchAtDecAgeLen = DiscCatchAtDecAgeLen_Fem + DiscCatchAtDecAgeLen_Mal # single sex or combined sexes

  # total catch at decimal and and length (discards plus retained catches)
  TotCatchAtDecAgeLen_Fem = as.matrix(unlist(CatchCurveResults_Fem$TotCatchAtDecAgeLen)) # retained catches, decimal ages and lengths
  TotCatchAtDecAgeLen_Mal = as.matrix(unlist(CatchCurveResults_Mal$TotCatchAtDecAgeLen))
  TotCatchAtDecAgeLen = TotCatchAtDecAgeLen_Fem + TotCatchAtDecAgeLen_Mal # single sex or combined sexes

  # retained catch proportions at length
  ExpRetCatchPropAtLen_Fem = RetCatchAtLen_Fem / sum(RetCatchAtLen_Fem)
  ExpRetCatchPropAtLen_Mal = RetCatchAtLen_Mal / sum(RetCatchAtLen_Mal)
  ExpRetCatchPropAtLen = RetCatchAtLen / sum(RetCatchAtLen)

  # discarded catch proportions at length
  ExpDiscCatchPropAtLen_Fem = DiscCatchAtLen_Fem / sum(DiscCatchAtLen_Fem)
  ExpDiscCatchPropAtLen_Mal = DiscCatchAtLen_Mal / sum(DiscCatchAtLen_Mal)
  ExpDiscCatchPropAtLen = DiscCatchAtLen / sum(DiscCatchAtLen)

  # generate observed retained catch frequencies at length
  SampleSize_Fem = SampleSize[1] * (sum(RetCatchAtLen_Fem) / (sum(RetCatchAtLen_Fem) + sum(RetCatchAtLen_Mal)))
  SampleSize_Fem = round(SampleSize_Fem,0)
  ObsRetCatchFreqAtLen_Fem = as.vector(rmultinom(1, SampleSize_Fem, ExpRetCatchPropAtLen_Fem))
  ObsRelCatchAtLen_Fem = ObsRetCatchFreqAtLen_Fem/sum(ObsRetCatchFreqAtLen_Fem)

  SampleSize_Mal = SampleSize[1] * (sum(RetCatchAtLen_Mal) / (sum(RetCatchAtLen_Fem) + sum(RetCatchAtLen_Mal)))
  SampleSize_Mal = round(SampleSize_Mal,0)
  ObsRetCatchFreqAtLen_Mal = as.vector(rmultinom(1, SampleSize_Mal, ExpRetCatchPropAtLen_Mal))
  ObsRelCatchAtLen_Mal = ObsRetCatchFreqAtLen_Mal/sum(ObsRetCatchFreqAtLen_Mal)
  ObsRetCatchFreqAtLen = ObsRetCatchFreqAtLen_Fem + ObsRetCatchFreqAtLen_Mal # combined sexes
  ObsRelCatchAtLen = ObsRetCatchFreqAtLen/sum(ObsRetCatchFreqAtLen)

  # is.na(MLL)
  # is.na(SelParams)
  # is.na(RetenParams)

  # generate observed discarded catch frequencies at length
  if (is.na(MLL) & is.na(RetenParams[1])) {
    ObsDiscCatchFreqAtLen = NA
    ObsDiscCatchFreqAtLen_Fem = NA
    ObsDiscCatchFreqAtLen_Mal = NA
  } else {

    if (length(SampleSize)==2) {
      DiscSampleSize = SampleSize[2] # specifying sample sizes for both retained and discarded fish
    } else {
      DiscSampleSize = SampleSize # set discarded fish sample size equal to retained fish sample size
    }
    # cat("DiscSampleSize",DiscSampleSize,'\n')
    # cat("DiscCatchAtLen_Fem",DiscCatchAtLen_Fem,'\n')
    # cat("ExpDiscCatchPropAtLen_Fem",ExpDiscCatchPropAtLen_Fem,'\n')

    DiscSampleSize_Fem = DiscSampleSize * (sum(DiscCatchAtLen_Fem) / (sum(DiscCatchAtLen_Fem) + sum(DiscCatchAtLen_Mal)))
    DiscSampleSize_Fem = round(DiscSampleSize_Fem,0)
    ObsDiscCatchFreqAtLen_Fem = as.vector(rmultinom(1, DiscSampleSize_Fem, ExpDiscCatchPropAtLen_Fem))
    ObsRelDiscCatchAtLen_Fem = ObsDiscCatchFreqAtLen_Fem/sum(ObsDiscCatchFreqAtLen_Fem)

    DiscSampleSize_Mal = DiscSampleSize * (sum(DiscCatchAtLen_Mal) / (sum(DiscCatchAtLen_Fem) + sum(DiscCatchAtLen_Mal)))
    DiscSampleSize_Mal = round(DiscSampleSize_Mal,0)
    ObsDiscCatchFreqAtLen_Mal = as.vector(rmultinom(1, DiscSampleSize_Mal, ExpDiscCatchPropAtLen_Mal))
    ObsRelDiscCatchAtLen_Mal = ObsDiscCatchFreqAtLen_Mal/sum(ObsDiscCatchFreqAtLen_Mal)
    ObsDiscCatchFreqAtLen = ObsDiscCatchFreqAtLen_Fem + ObsDiscCatchFreqAtLen_Mal # combined sexes
    ObsRelDiscCatchAtLen = ObsDiscCatchFreqAtLen/sum(ObsDiscCatchFreqAtLen)
  }

  # retained catches at decimal age
  RetCatchAtDecAge_Fem = rowSums(RetCatchAtDecAgeLen_Fem)
  RetCatchAtDecAge_Mal = rowSums(RetCatchAtDecAgeLen_Mal)
  RetCatchAtDecAge = RetCatchAtDecAge_Fem + RetCatchAtDecAge_Mal

  # Discarded catches at decimal age
  DiscCatchAtDecAge_Fem = rowSums(DiscCatchAtDecAgeLen_Fem)
  DiscCatchAtDecAge_Mal = rowSums(DiscCatchAtDecAgeLen_Mal)
  DiscCatchAtDecAge = DiscCatchAtDecAge_Fem + DiscCatchAtDecAge_Mal

  # Discarded catches at decimal age
  TotCatchAtDecAge_Fem = rowSums(TotCatchAtDecAgeLen_Fem)
  TotCatchAtDecAge_Mal = rowSums(TotCatchAtDecAgeLen_Mal)
  TotCatchAtDecAge = TotCatchAtDecAge_Fem + TotCatchAtDecAge_Mal

  # catch proportions at decimal age
  ExpCatchPropAtDecAge_Fem = RetCatchAtDecAge_Fem / sum(RetCatchAtDecAge_Fem)
  ExpCatchPropAtDecAge_Mal = RetCatchAtDecAge_Mal / sum(RetCatchAtDecAge_Mal)
  ExpCatchPropAtDecAge = RetCatchAtDecAge / sum(RetCatchAtDecAge)

  ExpRetCatchPropLengthGivenDecAge_Fem <- data.frame(matrix(nrow = nTimeSteps, ncol = nLenCl))
  colnames(ExpRetCatchPropLengthGivenDecAge_Fem) <- midpt
  ExpRetCatchPropLengthGivenDecAge_Fem = as.matrix(ExpRetCatchPropLengthGivenDecAge_Fem)
  ExpRetCatchPropLengthGivenDecAge_Mal = ExpRetCatchPropLengthGivenDecAge_Fem

  for (j in 1:nTimeSteps) {
    if (sum(RetCatchAtDecAgeLen_Fem[j,])==0) {
      ExpRetCatchPropLengthGivenDecAge_Fem[j,] = 0
    } else {
      ExpRetCatchPropLengthGivenDecAge_Fem[j,] = RetCatchAtDecAgeLen_Fem[j,] / sum(RetCatchAtDecAgeLen_Fem[j,])
    }
    if (sum(RetCatchAtDecAgeLen_Mal[j,])==0) {
      ExpRetCatchPropLengthGivenDecAge_Mal[j,] = 0
    } else {
      ExpRetCatchPropLengthGivenDecAge_Mal[j,] = RetCatchAtDecAgeLen_Mal[j,] / sum(RetCatchAtDecAgeLen_Mal[j,])
    }
  }

  ExpCatchPropDecAgeGivenLength_Fem = ExpRetCatchPropLengthGivenDecAge_Fem
  ExpCatchPropDecAgeGivenLength_Mal = ExpRetCatchPropLengthGivenDecAge_Fem
  ObsCatchPropDecAgeAtLength_Fem = ExpRetCatchPropLengthGivenDecAge_Fem
  ObsCatchPropDecAgeAtLength_Mal = ExpRetCatchPropLengthGivenDecAge_Fem
  FractDenom_Fem = rep(0,nLenCl); FractDenom_Mal = rep(0,nLenCl)
  for (i in 1:nLenCl) {
    for (j in 1:nTimeSteps) {
      FractDenom_Fem[i] = FractDenom_Fem[i] + (ExpRetCatchPropLengthGivenDecAge_Fem[j,i] * ExpCatchPropAtDecAge_Fem[j])
      if (FractDenom_Fem[i] == 0) {
        FractDenom_Fem[i] = 1E-20
      }
      FractDenom_Mal[i] = FractDenom_Mal[i] + (ExpRetCatchPropLengthGivenDecAge_Mal[j,i] * ExpCatchPropAtDecAge_Mal[j])
      if (FractDenom_Mal[i] == 0) {
        FractDenom_Mal[i]= 1E-20
      }
    }
    for (j in 1:nTimeSteps) {
      ExpCatchPropDecAgeGivenLength_Fem[j,i] = (ExpRetCatchPropLengthGivenDecAge_Fem[j,i] * ExpCatchPropAtDecAge_Fem[j]) / FractDenom_Fem[i]
      ExpCatchPropDecAgeGivenLength_Mal[j,i] = (ExpRetCatchPropLengthGivenDecAge_Mal[j,i] * ExpCatchPropAtDecAge_Mal[j]) / FractDenom_Mal[i]
    }
  }

  # observed catch frequencies at length and decimal age
  ObsCatchFreqAtLengthAndDecAge_Fem = ExpRetCatchPropLengthGivenDecAge_Fem
  ObsCatchFreqAtLengthAndDecAge_Mal = ExpRetCatchPropLengthGivenDecAge_Fem
  ObsCatchFreqAtLengthAndDecAge = ExpRetCatchPropLengthGivenDecAge_Fem
  for (i in 1:nLenCl) {
    ObsCatchFreqAtLengthAndDecAge_Fem[,i] = rmultinom(1, ObsRetCatchFreqAtLen_Fem[i], ExpCatchPropDecAgeGivenLength_Fem[,i])
    ObsCatchFreqAtLengthAndDecAge_Mal[,i] = rmultinom(1, ObsRetCatchFreqAtLen_Mal[i], ExpCatchPropDecAgeGivenLength_Mal[,i])
    ObsCatchFreqAtLengthAndDecAge[,i] = ObsCatchFreqAtLengthAndDecAge_Fem[,i] + ObsCatchFreqAtLengthAndDecAge_Mal[,i]
  }

  # generate data for individual fish (age classes and mid points of length classes)
  ObsDecAge_Fem = rep(NA, SampleSize_Fem); ObsDecAge_Mal = rep(NA, SampleSize_Mal); ObsDecAge = rep(NA, SampleSize[1])
  ObsAgeCl_Fem = rep(NA, SampleSize_Fem); ObsAgeCl_Mal = rep(NA, SampleSize_Mal); ObsAgeCl = rep(NA, SampleSize[1])
  ObsLenClMidPt_Fem = rep(NA, SampleSize_Fem); ObsLenClMidPt_Mal = rep(NA, SampleSize_Mal); ObsLenClMidPt = rep(NA, SampleSize[1])

  # females
  strt=1; fnsh=0
  for (j in 1:nTimeSteps) {
    for (i in 1:nLenCl) {
      x=ObsCatchFreqAtLengthAndDecAge_Fem[j,i] # number of fish in current length and age class
      if(x>0) {
        fnsh=strt+x-1
        ObsDecAge_Fem[strt:fnsh]=j*TimeStep
        ObsLenClMidPt_Fem[strt:fnsh]=midpt[i]
        strt=strt+x
      }
    }
  }
  # males
  strt=1; fnsh=0
  for (j in 1:nTimeSteps) {
    for (i in 1:nLenCl) {
      x=ObsCatchFreqAtLengthAndDecAge_Mal[j,i] # number of fish in current length and age class
      if(x>0) {
        fnsh=strt+x-1
        ObsDecAge_Mal[strt:fnsh]=j*TimeStep
        ObsLenClMidPt_Mal[strt:fnsh]=midpt[i]
        strt=strt+x
      }
    }
  }

  ObsAgeCl_Fem = floor(round(ObsDecAge_Fem,6)); ObsAgeCl_Mal = floor(round(ObsDecAge_Mal,6))
  ObsAgeCl = c(ObsAgeCl_Fem, ObsAgeCl_Mal) # decimal ages combined sexes
  ObsDecAge = c(ObsDecAge_Fem, ObsDecAge_Mal) # decimal ages combined sexes
  ObsLenClMidPt = c(ObsLenClMidPt_Fem, ObsLenClMidPt_Mal)

  # get frequencies at integer ages
  MinAge = floor(TimeStep)
  nAgeCl = length(MinAge:MaxAge)
  ObsCatchFreqAtLengthAndIntAge_Fem <- data.frame(matrix(0,nrow = nAgeCl, ncol = nLenCl))
  ObsCatchFreqAtLengthAndIntAge_Fem = as.matrix(ObsCatchFreqAtLengthAndIntAge_Fem)
  colnames(ObsCatchFreqAtLengthAndIntAge_Fem)=midpt

  ObsCatchFreqAtLengthAndIntAge_Mal=ObsCatchFreqAtLengthAndIntAge_Fem
  ObsCatchFreqAtLengthAndIntAge=ObsCatchFreqAtLengthAndIntAge_Fem
  nTimeSteps = length(seq(TimeStep,MaxAge,TimeStep))

  for (i in 1:length(ObsAgeCl_Fem)) {
    tempAge = ObsAgeCl_Fem[i]
    tempLen = ObsLenClMidPt_Fem[i]
    j = tempAge - MinAge + 1
    k = which(midpt==tempLen)
    ObsCatchFreqAtLengthAndIntAge_Fem[j,k] = ObsCatchFreqAtLengthAndIntAge_Fem[j,k] + 1
  }
  for (i in 1:length(ObsAgeCl_Mal)) {
    tempAge = ObsAgeCl_Mal[i]
    tempLen = ObsLenClMidPt_Mal[i]
    j = tempAge - MinAge + 1
    k = which(midpt==tempLen)
    ObsCatchFreqAtLengthAndIntAge_Mal[j,k] = ObsCatchFreqAtLengthAndIntAge_Mal[j,k] + 1
  }
  ObsCatchFreqAtLengthAndIntAge = ObsCatchFreqAtLengthAndIntAge_Fem + ObsCatchFreqAtLengthAndIntAge_Mal

  if (is.vector(GrowthParams)) { # params inputted as combined sex
    ObsRetCatchFreqAtLen_Fem = NA
    ObsRetCatchFreqAtLen_Mal = NA
    RetCatchAtDecAge_Fem = NA
    RetCatchAtDecAge_Mal = NA
    DiscCatchAtDecAge_Fem = NA
    DiscCatchAtDecAge_Mal = NA
    ExpCatchPropAtDecAge_Fem = NA
    ExpCatchPropAtDecAge_Mal = NA
    ExpRetCatchPropLengthGivenDecAge_Fem = NA
    ExpRetCatchPropLengthGivenDecAge_Mal = NA
    ExpCatchPropDecAgeGivenLength_Fem = NA
    ExpCatchPropDecAgeGivenLength_Mal = NA
    ObsCatchFreqAtLengthAndDecAge_Fem = NA
    ObsCatchFreqAtLengthAndDecAge_Mal = NA
    ObsLenClMidPt_Fem = NA
    ObsLenClMidPt_Mal = NA
  }
  if (is.data.frame(GrowthParams)) { # params inputted as separate sex
    RetCatchAtDecAge = NA
    DiscCatchAtDecAge = NA
    ObsCatchFreqAtLengthAndDecAge = NA
    ObsAge = NA
    ObsLenClMidPt = NA
  }

  # get expected proportion and associated sd of fish that are released
  PropReleased = sum(DiscCatchAtLen) / sum(TotCatchAtLen)
  PropReleased_sd = sqrt(PropReleased * (1 - PropReleased))/SampleSize

  # get expected mean length and associated sd of released fish
  MeanLenReleased = sum((DiscCatchAtLen * midpt)) / sum(DiscCatchAtLen)
  MeanLenReleased_sd = sqrt((sum((Res$DiscCatchAtLen * (Res$midpt^2))) / sum(Res$DiscCatchAtLen)) - MeanLenReleased^2)

  Results = list(DecAges = DecAges,
                 RecLenDist = RecLenDist,
                 MeanSizeAtAge = MeanSizeAtAge,
                 TimestepGrowthSizeInc = TimestepGrowthSizeInc,
                 MeanEndingLength = MeanEndingLength,
                 lbnd = lbnd,
                 midpt = midpt,
                 ubnd = ubnd,
                 RetAtLength = RetAtLength,
                 SelLandAtLength = SelLandAtLength,
                 DiscCatchAtLen = DiscCatchAtLen,
                 RetCatchAtLen = RetCatchAtLen,
                 TotCatchAtLen = TotCatchAtLen,
                 SampleSize_Fem = SampleSize_Fem,
                 SampleSize_Mal = SampleSize_Mal,
                 ObsRetCatchFreqAtLen = ObsRetCatchFreqAtLen,
                 ObsRetCatchFreqAtLen_Fem = ObsRetCatchFreqAtLen_Fem,
                 ObsRetCatchFreqAtLen_Mal = ObsRetCatchFreqAtLen_Mal,
                 ObsDiscCatchFreqAtLen = ObsDiscCatchFreqAtLen,
                 ObsDiscCatchFreqAtLen_Fem = ObsDiscCatchFreqAtLen_Fem,
                 ObsDiscCatchFreqAtLen_Mal = ObsDiscCatchFreqAtLen_Mal,
                 RetCatchAtDecAge_Fem = RetCatchAtDecAge_Fem,
                 RetCatchAtDecAge_Mal = RetCatchAtDecAge_Mal,
                 RetCatchAtDecAge = RetCatchAtDecAge,
                 DiscCatchAtDecAge_Fem = DiscCatchAtDecAge_Fem,
                 DiscCatchAtDecAge_Mal = DiscCatchAtDecAge_Mal,
                 DiscCatchAtDecAge = DiscCatchAtDecAge,
                 TotCatchAtDecAge_Fem = TotCatchAtDecAge_Fem,
                 TotCatchAtDecAge_Mal = TotCatchAtDecAge_Mal,
                 TotCatchAtDecAge = TotCatchAtDecAge,
                 ExpCatchPropAtDecAge = ExpCatchPropAtDecAge,
                 ExpCatchPropAtDecAge_Fem = ExpCatchPropAtDecAge_Fem,
                 ExpCatchPropAtDecAge_Mal = ExpCatchPropAtDecAge_Mal,
                 ExpRetCatchPropLengthGivenDecAge_Fem = ExpRetCatchPropLengthGivenDecAge_Fem,
                 ExpRetCatchPropLengthGivenDecAge_Mal = ExpRetCatchPropLengthGivenDecAge_Mal,
                 ExpCatchPropDecAgeGivenLength_Fem = ExpCatchPropDecAgeGivenLength_Fem,
                 ExpCatchPropDecAgeGivenLength_Mal = ExpCatchPropDecAgeGivenLength_Mal,
                 ObsCatchFreqAtLengthAndDecAge_Fem = ObsCatchFreqAtLengthAndDecAge_Fem,
                 ObsCatchFreqAtLengthAndDecAge_Mal = ObsCatchFreqAtLengthAndDecAge_Mal,
                 ObsCatchFreqAtLengthAndDecAge = ObsCatchFreqAtLengthAndDecAge,
                 ObsCatchFreqAtLengthAndIntAge_Fem = ObsCatchFreqAtLengthAndIntAge_Fem,
                 ObsCatchFreqAtLengthAndIntAge_Mal = ObsCatchFreqAtLengthAndIntAge_Mal,
                 ObsCatchFreqAtLengthAndIntAge = ObsCatchFreqAtLengthAndIntAge,
                 ObsDecAge = ObsDecAge,
                 ObsDecAge_Fem = ObsDecAge_Fem,
                 ObsDecAge_Mal = ObsDecAge_Mal,
                 ObsAgeCl = ObsAgeCl,
                 ObsAgeCl_Fem = ObsAgeCl_Fem,
                 ObsAgeCl_Mal = ObsAgeCl_Mal,
                 ObsLenClMidPt = ObsLenClMidPt,
                 ObsLenClMidPt_Fem = ObsLenClMidPt_Fem,
                 ObsLenClMidPt_Mal = ObsLenClMidPt_Mal,
                 PropReleased = PropReleased,
                 PropReleased_sd = PropReleased_sd,
                 MeanLenReleased = MeanLenReleased,
                 MeanLenReleased_sd = MeanLenReleased_sd)
  return(Results)

}


#' Produce plot of fitted length-based catch curve to length composition data for retained catches
#'
#' @param params vector of model parameters (in log space) to be estimated (if FittedRes=NA)
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param SelectivityType 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' @param ObsRetCatchFreqAtLen observed frequencies in length classes
#' @param lbnd lower bounds of length classes
#' @param ubnd upper bounds of length classes
#' @param midpt mid points of length classes
#' @param SelectivityVec selectivity at length
#' @param PropReleased proportion of fish that are released
#' @param ObsDiscCatchFreqAtLen Mean length of fish that are released
#' @param DiscMort selectivity at length
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
#' @param GrowthParams c(Linf, vbK) single sex von Bertalanffy, or data.frame(Linf=Linf, vbK=vbK) separate sex von Bertalanffy,
#' c(y1, y2, a, b) single sex Schnute, or data.frame(y1=y1, y2=y2, a=a, b=b) separate sex Schnute
#' @param RefnceAges Schnute reference ages, either c(t1,t2) single sex, c(t1=t1, t2=t2) separate sex, or set to NA if using von Bertalanffy model
#' @param MaxAge maximum age considered in model
#' @param NatMort natural mortality
#' @param TimeStep model timestep (e.g. 1 = annual, 1/12 = monthly)
#' @param MainLabel plot label
#' @param xaxis_lab y axis label
#' @param yaxis_lab x axis label
#' @param xmax maximum x axis value
#' @param xint x axis interval
#' @param ymax maximum y axis value
#' @param yint y axis interval
#' @param PlotCLs logical (TRUE=plot 95 percent confidence limits for line)
#' @param FittedRes saved results from GetLengthBasedCatchCurveResults function (model will be refitted if set to NA)
#' @param nReps number of random parameter sets from parametric resampling to generate outputs with error
#'
#' @return plot of fitted length-based catch curve to length composition data for retained fish
#'
#' @examples
#' # Simulate data
#' SampleSize=5000
#' set.seed(123)
#' MaxAge = 40
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = 0.15
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified is knife-edged at MLL
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 500) # L50, L95 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95 for retention
#' DiscMort = 0.5
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = 900
#' vbK = 0.11
#' CVSizeAtAge = 0.08
#' GrowthParams = c(Linf, vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen
#' ObsDiscCatchFreqAtLen = Res$ObsDiscCatchFreqAtLen
#' PropReleased=NA
#' midpt=Res$midpt
#' lbnd=Res$lbnd
#' ubnd=Res$ubnd
#' InitFishMort = 0.25 # specify starting parameters
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
#' InitL50 = 300
#' InitDelta = 100
#' params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
#' FittedRes=GetLengthBasedCatchCurveResults(params, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
#'                                           lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge,
#'                                           MaxAge, NatMort, TimeStep)
#' PlotLengthBasedCatchCurve_RetCatch(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
#'                                  SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
#'                                  xaxis_lab=NA, yaxis_lab=NA, xmax=1500, xint=500,
#'                                  ymax=0.15, yint=0.05, PlotCLs=TRUE, FittedRes, nReps=200)
#' @export
PlotLengthBasedCatchCurve_RetCatch <- function(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                             SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                             MaxAge, NatMort, TimeStep, MainLabel, xaxis_lab, yaxis_lab, xmax, xint,
                                             ymax, yint, PlotCLs, FittedRes, nReps) {

  .pardefault <- par(no.readonly = TRUE) # store current par settings

  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    res =  FittedRes
  } else {
    res=GetLengthBasedCatchCurveResults(params, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                        lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
  }

  params = res$params
  vcov.params = res$vcov.Params
  ExpCatchAtLen = res$ExpRetCatchPropInLenClass
  ObsRelCatchAtLen = ObsRetCatchFreqAtLen/sum(ObsRetCatchFreqAtLen)
  SelAtLength = res$SelAtLength
  RetAtLength = res$RetAtLength
  SelLandAtLength = res$SelLandAtLength

  set.seed(123)
  sims = data.frame(MASS::mvrnorm(n = nReps, params, vcov.params))
  EstPropAtLen.sim = data.frame(matrix(nrow = nReps, ncol = length(midpt)))

  for (j in 1:nReps) {
    params = unlist(sims[j,])
    CatchCurveType=1 #1=length-based, 2=age and length based
    Res=AgeAndLengthBasedCatchCurvesCalcs(params, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge, CatchCurveType,
                                          MLL, SelectivityType, lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)

    EstPropAtLen.sim[j,] = unlist(Res$ExpRetCatchPropInLenClass)
    cat("j",j,'\n')
  }

  EstProp.sim = apply(EstPropAtLen.sim, 2, median)
  EstProp.sim_low = apply(EstPropAtLen.sim, 2, quantile, probs = 0.025)
  EstProp.sim_up = apply(EstPropAtLen.sim, 2, quantile, probs = 0.975)

  if (is.na(xaxis_lab)) xaxis_lab = "Length (mm)"
  if (is.na(yaxis_lab)) yaxis_lab = "Prop. (retained catch)"
  xlims = Get_xaxis_scale(ubnd)
  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint
  ylims = Get_yaxis_scale(ObsRelCatchAtLen)
  if (is.na(ymax)) ymax = ylims$ymax
  if (is.na(yint)) yint = ylims$yint

  plot(midpt, ObsRelCatchAtLen, "p", main=MainLabel, cex.main=1.2, pch=16, cex=0.8, xaxt = "n", yaxt = "n",
       xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax))
  if (is.data.frame(GrowthParams)) {
    lines(midpt, Res$ExpRetCatchPropInLenClass_Fem, lty="dotted", col="dark grey")
    lines(midpt, Res$ExpRetCatchPropInLenClass_Mal, lty="dotted", col="blue")
    legend("bottomleft", legend=c("Female","Male"), y.intersp = 1.0, inset=c(0.05,0.05),
           cex = 0.8, bty="n", lty="dotted", col=c("dark grey","blue"))
  }
  if (PlotCLs == TRUE) {
    sm1 = spline(Res$midpt, EstProp.sim_low, n=100, method="natural")
    sm2 = spline(Res$midpt, EstProp.sim_up, n=100, method="natural")
    sm1$y[which(sm1$y<0)]=0
    sm2$y[which(sm2$y<0)]=0
    if (!is.na(MLL)) {
      sm1$y[which(sm1$x<MLL)]=0
      sm2$y[which(sm2$x<MLL)]=0
    }
    x = c(sm1$x, rev(sm2$x)) # using shading for 95% CLs
    y = c(sm1$y, rev(sm2$y))
    polygon(x,y, col="pink",border=NA)
  }
  points(midpt, ObsRelCatchAtLen, col="black", pch=16, cex=0.8)
  points(midpt, ExpCatchAtLen, col="red", pch=1, cex=0.8)
  axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
  axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
  axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
  axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
  # inverse logit transformed value
  params = res$params # from point estimate, not from resampled values
  Fval = round(1/(1+exp(-params[1])),2)
  Fest = bquote("F =" ~ .(Fval) ~ y^-1)
  if (SelectivityType==1) {
    legend("topright", pch=-1, legend=as.expression(Fest),
           lty="solid",col="black", bty='n', cex=0.8,lwd=-1, y.intersp=1.2, adj=0)
  }
  if (SelectivityType==2) {
    L50est=paste("L50 =",round(exp(params[2]),0),"mm")
    L95est=paste("L95 - L50 =",round(exp(params[3]),0),"mm")
    legend("topright", pch=-1, legend=c(as.expression(Fest), L50est, L95est),
           lty="solid",col="black", bty='n', cex=0.8,lwd=-1, y.intersp=1.2)
  }
  legend("topleft", legend=c("Observed","Estimated"), y.intersp = 1.0, inset=c(0.13,0),
         lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16,1), col=c("black","red"))

  par(.pardefault)

}


#' Produce plot of fitted length-based catch curve to length composition data for discarded catches
#'
#' @param params vector of model parameters (in log space) to be estimated (if FittedRes=NA)
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param SelectivityType 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' @param ObsRetCatchFreqAtLen observed frequencies in length classes
#' @param lbnd lower bounds of length classes
#' @param ubnd upper bounds of length classes
#' @param midpt mid points of length classes
#' @param SelectivityVec selectivity at length
#' @param PropReleased proportion of fish that are released
#' @param ObsDiscCatchFreqAtLen Mean length of fish that are released
#' @param DiscMort selectivity at length
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
#' @param GrowthParams c(Linf, vbK) single sex von Bertalanffy, or data.frame(Linf=Linf, vbK=vbK) separate sex von Bertalanffy,
#' c(y1, y2, a, b) single sex Schnute, or data.frame(y1=y1, y2=y2, a=a, b=b) separate sex Schnute
#' @param RefnceAges Schnute reference ages, either c(t1,t2) single sex, c(t1=t1, t2=t2) separate sex, or set to NA if using von Bertalanffy model
#' @param MaxAge maximum age considered in model
#' @param NatMort natural mortality
#' @param TimeStep model timestep (e.g. 1 = annual, 1/12 = monthly)
#' @param MainLabel plot label
#' @param xaxis_lab y axis label
#' @param yaxis_lab x axis label
#' @param xmax maximum x axis value
#' @param xint x axis interval
#' @param ymax maximum y axis value
#' @param yint y axis interval
#' @param PlotCLs logical (TRUE=plot 95 percent confidence limits for line)
#' @param FittedRes saved results from GetLengthBasedCatchCurveResults function (model will be refitted if set to NA)
#' @param nReps number of random parameter sets from parametric resampling to generate outputs with error
#'
#' @return plot of fitted length-based catch curve to length composition data for discarded fish (when length data available for discards)
#'
#' @examples
#' # Simulate data
#' SampleSize=5000
#' set.seed(123)
#' MaxAge = 40
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = 0.15
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified is knife-edged at MLL
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 500) # L50, L95 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95 for retention
#' DiscMort = 0.5
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = 900
#' vbK = 0.11
#' CVSizeAtAge = 0.08
#' GrowthParams = c(Linf, vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen
#' ObsDiscCatchFreqAtLen = Res$ObsDiscCatchFreqAtLen
#' PropReleased=NA
#' midpt=Res$midpt
#' lbnd=Res$lbnd
#' ubnd=Res$ubnd
#' InitFishMort = 0.25 # specify starting parameters
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
#' InitL50 = 300
#' InitDelta = 100
#' params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
#' FittedRes=GetLengthBasedCatchCurveResults(params, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
#'                                           lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge,
#'                                           MaxAge, NatMort, TimeStep)
#' PlotLengthBasedCatchCurve_DiscCatch(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
#'                                  SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
#'                                  xaxis_lab=NA, yaxis_lab=NA, xmax=1500, xint=500,
#'                                  ymax=0.15, yint=0.05, PlotCLs=TRUE, FittedRes, nReps=200)
#' @export
PlotLengthBasedCatchCurve_DiscCatch <- function(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                               SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                               MaxAge, NatMort, TimeStep, MainLabel, xaxis_lab, yaxis_lab, xmax, xint,
                                               ymax, yint, PlotCLs, FittedRes, nReps) {

  .pardefault <- par(no.readonly = TRUE) # store current par settings

  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    res =  FittedRes
  } else {
    res=GetLengthBasedCatchCurveResults(params, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                        lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
  }

  params = res$params
  vcov.params = res$vcov.Params
  ExpCatchAtLen = res$ExpDiscCatchPropInLenClass
  ObsRelCatchAtLen = ObsDiscCatchFreqAtLen/sum(ObsDiscCatchFreqAtLen)
  SelAtLength = res$SelAtLength
  RetAtLength = res$RetAtLength
  SelLandAtLength = res$SelLandAtLength

  set.seed(123)
  sims = data.frame(MASS::mvrnorm(n = nReps, params, vcov.params))
  EstPropAtLen.sim = data.frame(matrix(nrow = nReps, ncol = length(midpt)))

  for (j in 1:nReps) {
    params = unlist(sims[j,])
    CatchCurveType=1 #1=length-based, 2=age and length based
    Res=AgeAndLengthBasedCatchCurvesCalcs(params, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge, CatchCurveType,
                                          MLL, SelectivityType, lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)

    EstPropAtLen.sim[j,] = unlist(Res$ExpDiscCatchPropInLenClass)
    cat("j",j,'\n')
  }

  EstProp.sim = apply(EstPropAtLen.sim, 2, median)
  EstProp.sim_low = apply(EstPropAtLen.sim, 2, quantile, probs = 0.025)
  EstProp.sim_up = apply(EstPropAtLen.sim, 2, quantile, probs = 0.975)

  if (is.na(xaxis_lab)) xaxis_lab = "Length (mm)"
  if (is.na(yaxis_lab)) yaxis_lab = "Prop. (discarded catch)"
  xlims = Get_xaxis_scale(ubnd)
  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint
  ylims = Get_yaxis_scale(ObsRelCatchAtLen)
  if (is.na(ymax)) ymax = ylims$ymax
  if (is.na(yint)) yint = ylims$yint

  plot(midpt, ObsRelCatchAtLen, "p", main=MainLabel, cex.main=1.2, pch=16, cex=0.8, xaxt = "n", yaxt = "n",
       xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax))
  if (is.data.frame(GrowthParams)) {

    lines(midpt, Res$ExpDiscCatchPropInLenClass_Fem, lty="dotted", col="dark grey")
    lines(midpt, Res$ExpDiscCatchPropInLenClass_Mal, lty="dotted", col="blue")
    legend("bottomright", legend=c("Female","Male"), y.intersp = 1.0, inset=c(0.05,0.05),
           cex = 0.8, bty="n", lty="dotted", col=c("dark grey","blue"))
  }
  if (PlotCLs == TRUE) {
    sm1 = spline(Res$midpt, EstProp.sim_low, n=100, method="natural")
    sm2 = spline(Res$midpt, EstProp.sim_up, n=100, method="natural")
    sm1$y[which(sm1$y<0)]=0
    sm2$y[which(sm2$y<0)]=0
    if (!is.na(MLL)) {
      if ((FittedRes$ParamEst[2,1] + FittedRes$ParamEst[3,1])<MLL) {
        sm1$y[which(sm1$x>MLL)]=0
        sm2$y[which(sm2$x>MLL)]=0
      }
    }
    x = c(sm1$x, rev(sm2$x)) # using shading for 95% CLs
    y = c(sm1$y, rev(sm2$y))
    polygon(x,y, col="pink",border=NA)
  }
  points(midpt, ObsRelCatchAtLen, col="black", pch=16, cex=0.8)
  points(midpt, ExpCatchAtLen, col="red", pch=1, cex=0.8)
  axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
  axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
  axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
  axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
  # inverse logit transformed value
  params = res$params # from point estimate, not from resampled values
  Fval = round(1/(1+exp(-params[1])),2)
  Fest = bquote("F =" ~ .(Fval) ~ y^-1)
  if (SelectivityType==1) {
    legend("topright", pch=-1, legend=as.expression(Fest),
           lty="solid",col="black", bty='n', cex=0.8,lwd=-1, y.intersp=1.2, adj=0)
  }
  if (SelectivityType==2) {
    L50est=paste("L50 =",round(exp(params[2]),0),"mm")
    L95est=paste("L95 - L50 =",round(exp(params[3]),0),"mm")
    legend("topright", pch=-1, legend=c(as.expression(Fest), L50est, L95est),
           lty="solid",col="black", bty='n', cex=0.8,lwd=-1, y.intersp=1.2)
  }
  legend("topleft", legend=c("Observed","Estimated"), y.intersp = 1.0, inset=c(0.13,0),
         lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16,1), col=c("black","red"))

  par(.pardefault)

}


#' Show estimated selectivity curve from length-based catch curve model
#'
#' @param params vector of model parameters (in log space) to be estimated (if FittedRes=NA)
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param SelectivityType 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' @param ObsRetCatchFreqAtLen observed frequencies in length classes
#' @param lbnd lower bounds of length classes
#' @param ubnd upper bounds of length classes
#' @param midpt mid points of length classes
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param SelectivityVec selectivity at length
#' @param PropReleased proportion of fish that are released
#' @param ObsDiscCatchFreqAtLen mean lenght of released fish
#' @param DiscMort Proportion of fish that die following to capture and release
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
#' @param GrowthParams c(Linf, vbK) single sex von Bertalanffy, or data.frame(Linf=Linf, vbK=vbK) separate sex von Bertalanffy,
#' c(y1, y2, a, b) single sex Schnute, or data.frame(y1=y1, y2=y2, a=a, b=b) separate sex Schnute
#' @param RefnceAges Schnute reference ages, either c(t1,t2) single sex, c(t1=t1, t2=t2) separate sex, or set to NA if using von Bertalanffy model
#' @param MaxAge maximum age considered in model
#' @param NatMort natural mortality
#' @param TimeStep model timestep (e.g. 1 = annual, 1/12 = monthly)
#' @param MainLabel plot label
#' @param xaxis_lab y axis label
#' @param yaxis_lab x axis label
#' @param xmax maximum x axis value
#' @param xint x axis interval
#' @param ymax maximum y axis value
#' @param yint y axis interval
#' @param PlotCLs logical (TRUE=plot 95 percent confidence limits for line)
#' @param FittedRes saved results from GetLengthBasedCatchCurveResults function (model will be refitted if set to NA)
#' @param nReps number of random parameter sets from parametric resampling to generate outputs with error
#'
#' @return Plot of selectivity curve estimated by a length-based catch curve
#'
#' @examples
#' # Simulate data
#' SampleSize=5000
#' set.seed(123)
#' MaxAge = 40
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = 0.15
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified is knife-edged at MLL
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 500) # L50, L95 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95 for retention
#' DiscMort = 0.5
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = 900
#' vbK = 0.11
#' CVSizeAtAge = 0.08
#' GrowthParams = c(Linf, vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen
#' ObsDiscCatchFreqAtLen = Res$ObsDiscCatchFreqAtLen
#' PropReleased=NA
#' midpt=Res$midpt
#' lbnd=Res$lbnd
#' ubnd=Res$ubnd
#' InitFishMort = 0.25 # specify starting parameters
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
#' InitL50 = 300
#' InitDelta = 100
#' params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
#' FittedRes=GetLengthBasedCatchCurveResults(params, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
#'                                           lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge,
#'                                           MaxAge, NatMort, TimeStep)
#' PlotLengthBasedCatchCurve_Selectivity(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityVec,
#'                                       PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
#'                                       MaxAge, NatMort, TimeStep, MainLabel=NA, xaxis_lab=NA, yaxis_lab=NA, xmax=1200, xint=200,
#'                                       ymax=1.0, yint=0.2, PlotCLs=TRUE, FittedRes, nReps=200)
#' @export
PlotLengthBasedCatchCurve_Selectivity <- function(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityVec,
                                                  PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                                  MaxAge, NatMort, TimeStep, MainLabel, xaxis_lab, yaxis_lab, xmax, xint,
                                                  ymax, yint, PlotCLs, FittedRes, nReps) {

  .pardefault <- par(no.readonly = TRUE) # store current par settings

  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    res =  FittedRes
  } else {
    res=GetLengthBasedCatchCurveResults(params, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                        lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)

  }

  FishLen = seq(0,max(ubnd),1)
  nFishLen = length(FishLen)

  if (PlotCLs == TRUE) {
    if (SelectivityType == 2) {
      params = res$params
      vcov.params = res$vcov.Params
      set.seed(123)
      sims = data.frame(MASS::mvrnorm(n = nReps, params, vcov.params))
      SelAtLength.sim = data.frame(matrix(nrow = nReps, ncol = nFishLen))
      for (j in 1:nReps) {
        ParamVals = exp(unlist(sims[j,]))
        if (SelectivityType == 2) {
          # L95Val = (log(19)/ParamVals[3]) + ParamVals[2]

          SelAtLength.sim[j,] = 1 / (1 + exp(-log(19) * (FishLen - ParamVals[2]) /
                                               ParamVals[3]))
        }
        cat("j",j,"ParamVals",ParamVals,  '\n')
      }
      EstProp.sim = as.vector(apply(SelAtLength.sim, 2, median, na.rm=T))
      EstProp.sim_low = as.vector(apply(SelAtLength.sim, 2, quantile, probs = 0.025, na.rm=T))
      EstProp.sim_up = as.vector(apply(SelAtLength.sim, 2, quantile, probs = 0.975, na.rm=T))

    }
  }


  if (SelectivityType == 1) { # input selectivity as vector
    EstProp.sim = SelAtLength
    EstProp.sim_low = SelAtLength
    EstProp.sim_up = SelAtLength
  }

  if (is.na(xaxis_lab)) xaxis_lab = "Length (mm)"
  if (is.na(yaxis_lab)) yaxis_lab = "Selectivity"
  xlims = Get_xaxis_scale(midpt)
  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint
  if (is.na(ymax)) ymax = 1.2
  if (is.na(yint)) yint = 0.2

  plot(res$mid, res$SelAtLength, "p", main=MainLabel, cex.main=1.2, pch=1, cex=0.6,
       xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax),
       ylim=c(0,ymax), col="red")

  if (PlotCLs == TRUE) {
    sm1 = spline(FishLen, EstProp.sim_low, n=100, method="natural")
    sm2 = spline(FishLen, EstProp.sim_up, n=100, method="natural")
    x = c(sm1$x, rev(sm2$x)) # using shading for 95% CLs
    y = c(sm1$y, rev(sm2$y))
    polygon(x,y, col="pink",border=NA)

  }

  points(res$mid, res$SelAtLength, col="red", cex=0.6)
  axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
  axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
  axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
  axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)

  if (SelectivityType==2) { # logistic selectivity
    L50est=paste("L50 =",round(exp(params[2]),0),"mm")
    L95est=paste("L95 - L50 =",round(exp(params[3]),0),"mm")
    if(is.na(MLL)) {
      legend("topleft", pch=-1, legend=c(L50est, L95est), lty="solid",col="black",
             bty='n', cex=0.6,lwd=-1, y.intersp=1.0)
    }
  }
  if (is.numeric(MLL)) {
    abline(v=MLL)
    legend("topleft", pch=-1, legend=c(L50est, L95est, paste("MLL =",MLL)), lty="solid",col="black",
           bty='n', cex=0.6,lwd=-1, y.intersp=1.0)
  }

  par(.pardefault)
}


#' Show estimated fishing mortality at length relationships from length-based catch curve model
#'
#' This function provides plots of estimated fishing mortality at length relationships from a length-based
#' catch curve model, associated with fish retention and discard mortality, fish retention, fish discard mortality,
#' and 'fishing mortality' used to calculate to fish numbers that are released.
#'
#' @param params vector of model parameters (in log space) to be estimated (if FittedRes=NA)
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param SelectivityType 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' @param ObsRetCatchFreqAtLen observed frequencies in length classes
#' @param lbnd lower bounds of length classes
#' @param ubnd upper bounds of length classes
#' @param midpt mid points of length classes
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param SelectivityVec selectivity at length
#' @param PropReleased proportion of fish that are released
#' @param ObsDiscCatchFreqAtLen mean length of released fish
#' @param DiscMort Proportion of fish that die following to capture and release
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
#' @param GrowthParams c(Linf, vbK) single sex von Bertalanffy, or data.frame(Linf=Linf, vbK=vbK) separate sex von Bertalanffy,
#' c(y1, y2, a, b) single sex Schnute, or data.frame(y1=y1, y2=y2, a=a, b=b) separate sex Schnute
#' @param RefnceAges Schnute reference ages, either c(t1,t2) single sex, c(t1=t1, t2=t2) separate sex, or set to NA if using von Bertalanffy model
#' @param MaxAge maximum age considered in model
#' @param NatMort natural mortality
#' @param TimeStep model timestep (e.g. 1 = annual, 1/12 = monthly)
#' @param xmax maximum x axis value
#' @param xint x axis interval
#' @param ymax maximum y axis value
#' @param yint y axis interval
#'
#' @return Plots of fishing mortality at length relationships estimated by a length-based catch curve model
#'
#' @examples
#' # Simulate data
#' SampleSize=5000
#' set.seed(123)
#' MaxAge = 40
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = 0.15
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified is knife-edged at MLL
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 500) # L50, L95 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95 for retention
#' DiscMort = 0.5
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = 900
#' vbK = 0.11
#' CVSizeAtAge = 0.08
#' GrowthParams = c(Linf, vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen
#' ObsDiscCatchFreqAtLen = Res$ObsDiscCatchFreqAtLen
#' PropReleased=NA
#' midpt=Res$midpt
#' lbnd=Res$lbnd
#' ubnd=Res$ubnd
#' InitFishMort = 0.25 # specify starting parameters
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
#' InitL50 = 300
#' InitDelta = 100
#' params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
#' FittedRes=GetLengthBasedCatchCurveResults(params, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
#'                                           lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge,
#'                                           MaxAge, NatMort, TimeStep)
#' PlotLengthBasedCatchCurve_Mortality(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityVec,
#'                                     PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
#'                                     MaxAge, NatMort, TimeStep, xmax=NA, xint=NA, ymax=NA, yint=NA, FittedRes)
#' @export
PlotLengthBasedCatchCurve_Mortality <- function(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityVec,
                                                PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                                MaxAge, NatMort, TimeStep, xmax, xint, ymax, yint, FittedRes) {

  .pardefault <- par(no.readonly = TRUE) # store current par settings

  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    res =  FittedRes
  } else {
    res=GetLengthBasedCatchCurveResults(params, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                        lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)

  }


  CatchCurveType=1 #1=length-based, 2=age and length based
  params = res$params
  res = AgeAndLengthBasedCatchCurvesCalcs(params, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge, CatchCurveType,
                                          MLL, SelectivityType, lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
  GrowthModelType = res$GrowthModelType

  xlims = Get_xaxis_scale(midpt)
  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint
  if (is.na(ymax)) ymax = 0.5
  if (is.na(yint)) yint = 0.1

  par(mfrow=c(2,2), mar=c(5,4,2,2), oma=c(2,2,2,2))
  yaxis_lab = expression(paste(italic("F") ~ (year^{-1})))
  xaxis_lab = "Length(mm)"

  # F at age
  if (GrowthModelType ==  1 | GrowthModelType == 3) { # single sex
    plot(midpt, res$FAtLen, "o", main="F (retention + discard mortality)", cex.main=1.0, pch=1, cex=0.6,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=NA, frame=F, xlim=c(0,xmax),
         ylim=c(0,ymax), col="red")
    axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
    axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
    axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    mtext(expression(paste(italic("F") ~ (year^{-1}))),las=3,side=2,line=2.5,cex=1.0,lwd=1.75)
  }
  if (GrowthModelType ==  2 | GrowthModelType == 4) { # separate sex
    plot(midpt, res$FAtLen_Fem, "o", main="F (retention + discard mortality)", cex.main=1.0, pch=1, cex=0.6,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=NA, frame=F, xlim=c(0,xmax),
         ylim=c(0,ymax), col="red")
    axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
    axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
    axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    mtext(expression(paste(italic("F") ~ (year^{-1}))),las=3,side=2,line=2.5,cex=1.0,lwd=1.75)
    lines(midpt, res$FAtLen_Mal, col="blue")
    legend('topright', col=c("red","blue"),lty="solid",legend=c("females","males"),bty='n', cex=0.8,lwd=1.75)
  }

  # F at age, associated with retention
  if (GrowthModelType ==  1 | GrowthModelType == 3) { # single sex
    plot(midpt, res$FAtLenReten, "o", main="F (retention)", cex.main=1.0, pch=1, cex=0.6,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=NA, frame=F, xlim=c(0,xmax),
         ylim=c(0,ymax), col="blue")
    axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
    axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
    axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    mtext(expression(paste(italic("F") ~ (year^{-1}))),las=3,side=2,line=3,cex=1.0,lwd=1.75)
  }
  if (GrowthModelType ==  2 | GrowthModelType == 4) { # separate sex
    plot(midpt, res$FAtLenReten_Fem, "o", main="F (retention)", cex.main=1.0, pch=1, cex=0.6,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=NA, frame=F, xlim=c(0,xmax),
         ylim=c(0,ymax), col="blue")
    axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
    axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
    axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    mtext(expression(paste(italic("F") ~ (year^{-1}))),las=3,side=2,line=3,cex=1.0,lwd=1.75)
    lines(midpt, res$FAtLenReten_Mal, col="blue")
    legend('topright', col=c("red","blue"),lty="solid",legend=c("females","males"),bty='n', cex=0.8,lwd=1.75)
  }

  # F at age, associated with fish capture and post-release mortality
  if (GrowthModelType ==  1 | GrowthModelType == 3) { # single sex
    plot(midpt, res$FAtLenDisc*DiscMort, "o", main="F (discard mortality)", cex.main=1.0, pch=1, cex=0.6,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=NA, frame=F, xlim=c(0,xmax),
         ylim=c(0,ymax), col="brown")
    axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
    axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
    axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    mtext(expression(paste(italic("F") ~ (year^{-1}))),las=3,side=2,line=2.5,cex=1.0,lwd=1.75)
  }
  if (GrowthModelType ==  2 | GrowthModelType == 4) { # separate sex
    plot(midpt, res$FAtLenDisc_Fem*DiscMort, "o", main="F (discard mortality)", cex.main=1.0, pch=1, cex=0.6,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=NA, frame=F, xlim=c(0,xmax),
         ylim=c(0,ymax), col="brown")
    axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
    axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
    axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    mtext(expression(paste(italic("F") ~ (year^{-1}))),las=3,side=2,line=2.5,cex=1.0,lwd=1.75)
    lines(midpt, res$FAtLenDisc_Mal*DiscMort, col="blue")
    legend('topright', col=c("red","blue"),lty="solid",legend=c("females","males"),bty='n', cex=0.8,lwd=1.75)
  }


  # F at age, associated with fish capture and release
  if (GrowthModelType ==  1 | GrowthModelType == 3) { # single sex
    plot(midpt, res$FAtLenDisc, "o", main="F (capture + release)", cex.main=1.0, pch=1, cex=0.6,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=NA, frame=F, xlim=c(0,xmax),
         ylim=c(0,ymax), col="red")
    axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
    axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
    axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    mtext(expression(paste(italic("F") ~ (year^{-1}))),las=3,side=2,line=2.5,cex=1.0,lwd=1.75)
  }
  if (GrowthModelType ==  2 | GrowthModelType == 4) { # separate sex
    plot(midpt, res$FAtLenDisc_Fem, "o", main="F (capture + release)", cex.main=1.0, pch=1, cex=0.6,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=NA, frame=F, xlim=c(0,xmax),
         ylim=c(0,ymax), col="red")
    axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
    axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
    axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    mtext(expression(paste(italic("F") ~ (year^{-1}))),las=3,side=2,line=2.5,cex=1.0,lwd=1.75)
    lines(midpt, res$FAtLenDisc_Mal, col="blue")
    legend('topright', col=c("red","blue"),lty="solid",legend=c("females","males"),bty='n', cex=0.8,lwd=1.75)
  }
  par(.pardefault)
}



#' Show fit of age and length catch curve model to marginal length composition
#'
#' This function provides a plot of the fit of an age and length-based catch curve model, with length-
#' based selectivity to marginal length composition data. The model is fitted to a sample of fish
#' length and age data, by minimising the overall negative log-likelihood, including the NLL associated
#' with the marginal length composition and a conditional age at length NLL,
#' given the parameters (selectivity, growth and mortality) and data, using nlminb.
#' It provides various statistical outputs in include convergence statistics, parameter estimates
#' and associated 95 percent confidence limits and associated variance-covariance matrix, calculated using
#' the MASS package.
#'
#' @param params vector of model parameters (in log space) to be estimated (if FittedRes=NA)
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param SelectivityType 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' @param ObsRetCatchFreqAtLen observed frequencies in length classes
#' @param ObsRetCatchFreqAtLengthAndAge observed frequencies in length and age classes
#' @param lbnd lower bounds of length classes
#' @param ubnd upper bounds of length classes
#' @param midpt mid points of length classes
#' @param SelectivityVec selectivity at length
#' @param DiscMort Proportion of fish that die following to capture and release
#' @param MaxAge maximum age considered in model
#' @param NatMort natural mortality
#' @param TimeStep model timestep (e.g. 1 = annual, 1/12 = monthly)
#' @param MainLabel plot label
#' @param xaxis_lab y axis label
#' @param yaxis_lab x axis label
#' @param xmax maximum x axis value
#' @param xint x axis interval
#' @param ymax maximum y axis value
#' @param yint y axis interval
#' @param PlotCLs logical (TRUE=plot 95 percent confidence limits for line)
#' @param FittedRes saved results from GetLengthBasedCatchCurveResults function (model will be refitted if set to NA)
#' @param nReps number of random parameter sets from parametric resampling to generate outputs with error
#'
#' @return Plot of fitted age and length-based catch curve to marginal length composition
#'
#' @examples
#' # Simulate data
#' set.seed(123)
#' SampleSize=5000
#' MaxAge = 26
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' MinAge = floor(TimeStep)
#' nAgeCl = length(MinAge:MaxAge)
#' nTimeSteps = length(seq(TimeStep,MaxAge,TimeStep))
#' NatMort = 4.22/MaxAge
#' FishMort = 0.2
#' MaxLen = 1100
#' LenInc = 10
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified is knife-edged at MLL
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0
#' # # single sex, von Bertalanffy
#' # GrowthCurveType = 1 # 1 = von Bertalanffy
#' # Linf = 800
#' # vbK = 0.2
#' # CVSizeAtAge = 0.08
#' # GrowthParams = c(Linf, vbK, CVSizeAtAge)
#' # RefnceAges = NA
#' # 2 sexes, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy
#' Linf = c(700,850)
#' vbK = c(0.25,0.2)
#' CVSizeAtAge = c(0.05,0.05)
#' RefnceAges = NA
#' GrowthParams = data.frame(Linf=Linf, vbK=vbK, CVSizeAtAge=CVSizeAtAge)
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' lbnd=Res$lbnd
#' midpt=Res$midpt
#' ubnd=Res$ubnd
#' # # get data - 1 sex (or combined sexes)
#' # ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen # 1 sex
#' # ObsRetCatchFreqAtLengthAndAge = as.matrix(Res$ObsCatchFreqAtLengthAndDecAge) # 1 sex
#' # get data - 2 sexes
#' ObsRetCatchFreqAtLen <- data.frame(matrix(nrow = 2, ncol = length(midpt))) # 2 sex
#' colnames(ObsRetCatchFreqAtLen) <- midpt
#' ObsRetCatchFreqAtLen[1,] = Res$ObsRetCatchFreqAtLen_Fem
#' ObsRetCatchFreqAtLen[2,] = Res$ObsRetCatchFreqAtLen_Mal
#' ObsRetCatchFreqAtLengthAndAge = array(c(unlist(Res$ObsCatchFreqAtLengthAndDecAge_Fem), unlist(Res$ObsCatchFreqAtLengthAndDecAge_Mal)),
#'                                    c(nTimeSteps, length(midpt), 2), dimnames=list(rownames(Res$ObsCatchFreqAtLengthAndDecAge_Fem),
#'                                                                                   colnames(Res$ObsCatchFreqAtLengthAndDecAge_Fem)))
#' # # get params - 1 sex
#' # InitFishMort = 0.3 # specify starting parameters
#' # InitL50 = 320
#' # InitDelta = 50 # L95-L50
#' # InitLinf = 800
#' # InitvbK = 0.2
#' # InitCVSizeAtAge = 0.05
#' # get params - 2 sexes
#' InitFishMort = 0.3 # specify starting parameters
#' InitL50 = 320
#' InitDelta = 50 # L95-L50
#' InitLinf = c(800,800)
#' InitvbK = c(0.25,0.25)
#' InitCVSizeAtAge = 0.05
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' params = c(InitFishMort_logit, log(c(InitL50, InitDelta, InitLinf, InitvbK, InitCVSizeAtAge)))
#' FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
#' # Plot. Note, can skip above step and set FittedRes=NA (plot function will be slower)
#'PlotAgeLengthCatchCurve_MargLength(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                   lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep, MainLabel=NA,
#'                                   xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA,
#'                                   ymax=0.10, yint=0.02, PlotCLs=TRUE, FittedRes, nReps=200)
#' @export
PlotAgeLengthCatchCurve_MargLength <- function(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                               lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep, MainLabel,
                                               xaxis_lab, yaxis_lab, xmax, xint,
                                               ymax, yint, PlotCLs, FittedRes, nReps) {

  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    res =  FittedRes
  } else {

    res=GetAgeAndLengthBasedCatchCurveResults(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                              lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
  }

  params = res$params
  vcov.params = res$vcov.Params
  ExpCatchAtLen = res$ExpCatchAtLen
  if (is.vector(ObsRetCatchFreqAtLen)) {
    ObsRelCatchAtLen = ObsRetCatchFreqAtLen/sum(ObsRetCatchFreqAtLen)
  }
  if (is.data.frame(ObsRetCatchFreqAtLen)) {
    ObsRelCatchAtLen = data.frame(matrix(nrow = 2, ncol = length(midpt)))
    colnames(ObsRelCatchAtLen) = midpt
    ObsRelCatchAtLen[1,] = ObsRetCatchFreqAtLen[1,]/(sum(ObsRetCatchFreqAtLen[1,])+sum(ObsRetCatchFreqAtLen[2,]))
    ObsRelCatchAtLen[2,] = ObsRetCatchFreqAtLen[2,]/(sum(ObsRetCatchFreqAtLen[1,])+sum(ObsRetCatchFreqAtLen[2,]))
  }
  set.seed(123)
  sims = data.frame(MASS::mvrnorm(n = nReps, params, vcov.params))
  EstPropAtLen.sim = data.frame(matrix(nrow = nReps, ncol = length(midpt)))
  EstPropAtLen.sim_Fem = EstPropAtLen.sim
  EstPropAtLen.sim_Mal = EstPropAtLen.sim
  GrowthParams = NA
  RefnceAges = NA

  for (j in 1:nReps) {
    params = unlist(sims[j,])
    CatchCurveType=2 #1=length-based, 2=age and length based
    GrowthCurveType=1 # von Bertalanffy
    Res=AgeAndLengthBasedCatchCurvesCalcs(params, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge, CatchCurveType,
                                          MLL, SelectivityType, lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
    EstPropAtLen.sim[j,] = unlist(Res$ExpRetCatchPropInLenClass)
    EstPropAtLen.sim_Fem[j,] = unlist(Res$ExpRetCatchPropInLenClass_Fem)
    EstPropAtLen.sim_Mal[j,] = unlist(Res$ExpRetCatchPropInLenClass_Mal)
    cat("j",j,'\n')
  }

  # combined sexes
  if (is.vector(ObsRetCatchFreqAtLen)) {
    EstProp.sim = apply(EstPropAtLen.sim, 2, median)
    EstProp.sim_low = apply(EstPropAtLen.sim, 2, quantile, probs = 0.025)
    EstProp.sim_up = apply(EstPropAtLen.sim, 2, quantile, probs = 0.975)
  }
  # separate sexes - females
  if (is.data.frame(ObsRetCatchFreqAtLen)) {
    EstPropF.sim = apply(EstPropAtLen.sim_Fem, 2, median)
    EstPropF.sim_low = apply(EstPropAtLen.sim_Fem, 2, quantile, probs = 0.025)
    EstPropF.sim_up = apply(EstPropAtLen.sim_Fem, 2, quantile, probs = 0.975)
    # separate sexes - males
    EstPropM.sim = apply(EstPropAtLen.sim_Mal, 2, median)
    EstPropM.sim_low = apply(EstPropAtLen.sim_Mal, 2, quantile, probs = 0.025)
    EstPropM.sim_up = apply(EstPropAtLen.sim_Mal, 2, quantile, probs = 0.975)
  }

  if (is.na(xaxis_lab)) xaxis_lab = "Length (mm)"
  xlims = Get_xaxis_scale(ubnd)
  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint
  ylims = Get_yaxis_scale(ObsRelCatchAtLen)
  if (is.na(ymax)) ymax = ylims$ymax
  if (is.na(yint)) yint = ylims$yint

  # combined sexes
  if (is.vector(ObsRetCatchFreqAtLen)) {
    if (is.na(yaxis_lab)) yaxis_lab = "Proportion"
    plot(midpt, ObsRelCatchAtLen, "p", main=MainLabel, cex.main=1.2, pch=16, cex=0.8,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax))
    points(midpt, ObsRelCatchAtLen, col="black", pch=16, cex=0.8)
    if (PlotCLs == TRUE) {
      # x = c(Res$midpt,rev(Res$midpt)) # using shading for 95% CLs
      # y = c(EstProp.sim_low, rev(EstProp.sim_up))
      # polygon(x,y,col="pink",border=NA)
      sm1 = spline(Res$midpt, EstProp.sim_low, n=100, method="natural")
      sm2 = spline(Res$midpt, EstProp.sim_up, n=100, method="natural")
      x = c(sm1$x, rev(sm2$x)) # using shading for 95% CLs
      y = c(sm1$y, rev(sm2$y))
      polygon(x,y,col="pink",border=NA)
    }
    points(midpt, ObsRelCatchAtLen, col="black", pch=16, cex=0.8)
    points(midpt, EstProp.sim, col="red", pch=1, cex=0.8)
    axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
    axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
    axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
  }

  # separate sexes
  if (is.data.frame(ObsRetCatchFreqAtLen)) {
    if (is.na(yaxis_lab)) yaxis_lab1 = "Proportion - Females"
    plot(midpt, ObsRelCatchAtLen[1,], "p", main=MainLabel, cex.main=1.2, pch=16, cex=0.8,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab1,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax))
    points(midpt, ObsRelCatchAtLen[1,], col="black", pch=16, cex=0.8)
    if (PlotCLs == TRUE) {
      # x = c(Res$midpt,rev(Res$midpt)) # using shading for 95% CLs
      # y = c(EstPropF.sim_low, rev(EstPropF.sim_up))
      # polygon(x,y,col="pink",border=NA)
      sm1 = spline(Res$midpt, EstPropF.sim_low, n=100, method="natural")
      sm2 = spline(Res$midpt, EstPropF.sim_up, n=100, method="natural")
      x = c(sm1$x, rev(sm2$x)) # using shading for 95% CLs
      y = c(sm1$y, rev(sm2$y))
      polygon(x,y,col="pink",border=NA)
    }
    points(midpt, ObsRelCatchAtLen[1,], col="black", pch=16, cex=0.8)
    points(midpt, EstPropF.sim, col="red", pch=1, cex=0.8)
    axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
    axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
    axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)

    if (is.na(yaxis_lab)) yaxis_lab2 = "Proportion - Males"
    plot(midpt, ObsRelCatchAtLen[2,], "p", main=MainLabel, cex.main=1.2, pch=16, cex=0.8,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab2,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax))
    points(midpt, ObsRelCatchAtLen[2,], col="black", pch=16, cex=0.8)
    if (PlotCLs == TRUE) {
      # x = c(Res$midpt,rev(Res$midpt)) # using shading for 95% CLs
      # y = c(EstPropM.sim_low, rev(EstPropM.sim_up))
      # polygon(x,y,col="light blue",border=NA)
      sm1 = spline(Res$midpt, EstPropM.sim_low, n=100, method="natural")
      sm2 = spline(Res$midpt, EstPropM.sim_up, n=100, method="natural")
      x = c(sm1$x, rev(sm2$x)) # using shading for 95% CLs
      y = c(sm1$y, rev(sm2$y))
      polygon(x,y,col="light blue",border=NA)
    }
    points(midpt, ObsRelCatchAtLen[2,], col="black", pch=16, cex=0.8)
    points(midpt, EstPropM.sim, col="blue", pch=1, cex=0.8)
    axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
    axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
    axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
  }
}


#' Show estimated growth curve from age and length catch curve model
#'
#' This function provides a plot of growth curve estimated from an age and length-based catch curve model, with length-
#' based selectivity. The model is fitted to a sample of fish length and age data, by minimising the overall negative log-likelihood,
#' including the NLL associated with the marginal length composition and a conditional age at length NLL,
#' given the parameters (selectivity, growth and mortality) and data, using nlminb.
#' It provides various statistical outputs in include convergence statistics, parameter estimates
#' and associated 95 percent confidence limits and associated variance-covariance matrix, calculated using
#' the MASS package.
#'
#' @param params vector of model parameters (in log space) to be estimated (if FittedRes=NA)
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param SelectivityType 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' @param ObsRetCatchFreqAtLen observed frequencies in length classes
#' @param ObsRetCatchFreqAtLengthAndAge observed frequencies in length and age classes
#' @param lbnd lower bounds of length classes
#' @param ubnd upper bounds of length classes
#' @param midpt mid points of length classes
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param SelectivityVec selectivity at length
#' @param DiscMort proportion of fish that die following capture and release
#' @param MaxAge maximum age considered in model
#' @param NatMort natural mortality
#' @param TimeStep model timestep (e.g. 1 = annual, 1/12 = monthly)
#' @param MainLabel plot label
#' @param xaxis_lab y axis label
#' @param yaxis_lab x axis label
#' @param xmax maximum x axis value
#' @param xint x axis interval
#' @param ymax maximum y axis value
#' @param yint y axis interval
#' @param PlotCLs logical (TRUE=plot 95 percent confidence limits for line)
#' @param FittedRes saved results from GetLengthBasedCatchCurveResults function (model will be refitted if set to NA)
#' @param nReps number of random parameter sets from parametric resampling to generate outputs with error
#'
#' @return Plot showing growth curve estimated by age and length-based catch curve
#'
#' @examples
#' # Simulate data
#' set.seed(123)
#' SampleSize=5000
#' MaxAge = 26
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' MinAge = floor(TimeStep)
#' nAgeCl = length(MinAge:MaxAge)
#' nTimeSteps = length(seq(TimeStep,MaxAge,TimeStep))
#' NatMort = 4.22/MaxAge
#' FishMort = 0.2
#' MaxLen = 1100
#' LenInc = 10
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified is knife-edged at MLL
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0
#' # # single sex, von Bertalanffy
#' # GrowthCurveType = 1 # 1 = von Bertalanffy
#' # Linf = 800
#' # vbK = 0.2
#' # CVSizeAtAge = 0.08
#' # GrowthParams = c(Linf, vbK, CVSizeAtAge)
#' # RefnceAges = NA
#' # 2 sexes, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy
#' Linf = c(700,850)
#' vbK = c(0.25,0.2)
#' CVSizeAtAge = c(0.05,0.05)
#' RefnceAges = NA
#' GrowthParams = data.frame(Linf=Linf, vbK=vbK, CVSizeAtAge=CVSizeAtAge)
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' lbnd=Res$lbnd
#' midpt=Res$midpt
#' ubnd=Res$ubnd
#' # # get data - 1 sex (or combined sexes)
#' # ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen # 1 sex
#' # ObsRetCatchFreqAtLengthAndAge = as.matrix(Res$ObsCatchFreqAtLengthAndDecAge) # 1 sex
#' # get data - 2 sexes
#' ObsRetCatchFreqAtLen <- data.frame(matrix(nrow = 2, ncol = length(midpt))) # 2 sex
#' colnames(ObsRetCatchFreqAtLen) <- midpt
#' ObsRetCatchFreqAtLen[1,] = Res$ObsRetCatchFreqAtLen_Fem
#' ObsRetCatchFreqAtLen[2,] = Res$ObsRetCatchFreqAtLen_Mal
#' ObsRetCatchFreqAtLengthAndAge = array(c(unlist(Res$ObsCatchFreqAtLengthAndDecAge_Fem), unlist(Res$ObsCatchFreqAtLengthAndDecAge_Mal)),
#'                                    c(nTimeSteps, length(midpt), 2), dimnames=list(rownames(Res$ObsCatchFreqAtLengthAndDecAge_Fem),
#'                                                                                   colnames(Res$ObsCatchFreqAtLengthAndDecAge_Fem)))
#' # # get params - 1 sex
#' # InitFishMort = 0.3 # specify starting parameters
#' # InitL50 = 320
#' # InitDelta = 50 # L95-L50
#' # InitLinf = 800
#' # InitvbK = 0.2
#' # InitCVSizeAtAge = 0.05
#' # get params - 2 sexes
#' InitFishMort = 0.3 # specify starting parameters
#' InitL50 = 320
#' InitDelta = 50 # L95-L50
#' InitLinf = c(800,800)
#' InitvbK = c(0.25,0.25)
#' InitCVSizeAtAge = 0.05
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' params = c(InitFishMort_logit, log(c(InitL50, InitDelta, InitLinf, InitvbK, InitCVSizeAtAge)))
#' FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
#' # plot
#' PlotAgeLengthCatchCurve_Growth(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep, MainLabel=NA,
#'                                xaxis_lab=NA, yaxis_lab=NA, xmax=40, xint=10,
#'                                ymax=1000, yint=200, PlotCLs=TRUE, FittedRes, nReps=200)
#' @export
PlotAgeLengthCatchCurve_Growth <- function(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                           lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep, MainLabel,
                                           xaxis_lab, yaxis_lab, xmax, xint,
                                           ymax, yint, PlotCLs, FittedRes, nReps) {


  # generate data for individual fish (age classes and mid points of length classes)
  if (is.vector(ObsRetCatchFreqAtLen)) {
    SampleSize = sum(ObsRetCatchFreqAtLen)
    ObsAge = rep(NA, SampleSize)
    ObslenClMidPt = rep(NA, SampleSize)
  }
  if (is.data.frame(ObsRetCatchFreqAtLen)) {
    SampleSize_F = sum(ObsRetCatchFreqAtLen[1,])
    SampleSize_M = sum(ObsRetCatchFreqAtLen[2,])
    ObsAge_F = rep(NA, SampleSize_F)
    ObsAge_M = rep(NA, SampleSize_M)
    ObslenClMidPt_F = rep(NA, SampleSize_F)
    ObslenClMidPt_M = rep(NA, SampleSize_M)
  }
  # dim(ObsRetCatchFreqAtLengthAndAge)
  strt=1; fnsh=0
  strtF=1; fnshF=0
  strtM=1; fnshM=0
  nLenCl = length(midpt)
  MinAge = floor(TimeStep)
  nAgeCl = length(MinAge:MaxAge)
  DecAges = seq(TimeStep,MaxAge,TimeStep)
  nTimeSteps = length(DecAges)

  for (i in 1:nTimeSteps) {
    for (j in 1:nLenCl) {
      if (is.vector(ObsRetCatchFreqAtLen)) { # single sex
        x=ObsRetCatchFreqAtLengthAndAge[i,j] # number of fish in current length and age class
        if(x>0) {
          fnsh=strt+x-1
          ObsAge[strt:fnsh]=i*TimeStep
          ObslenClMidPt[strt:fnsh]=midpt[j]
          strt=strt+x
        }
      }
      if (is.data.frame(ObsRetCatchFreqAtLen)) { # 2 sexes
        # females
        x=ObsRetCatchFreqAtLengthAndAge[i,j,1] # number of females in current length and age class
        if(x>0) {
          fnshF=strtF+x-1
          ObsAge_F[strtF:fnshF]=i*TimeStep
          ObslenClMidPt_F[strtF:fnshF]=midpt[j]
          strtF=strtF+x
        }
        # males
        x=ObsRetCatchFreqAtLengthAndAge[i,j,2] # number of males in current length and age class
        if(x>0) {
          fnshM=strtM+x-1
          ObsAge_M[strtM:fnshM]=i*TimeStep
          ObslenClMidPt_M[strtM:fnshM]=midpt[j]
          strtM=strtM+x
        }
      }
    }
  }

  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    res =  FittedRes
  } else {
    res=GetAgeAndLengthBasedCatchCurveResults(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                              lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
  }
  params = res$params
  vcov.params = res$vcov.Params
  set.seed(123)
  sims = data.frame(MASS::mvrnorm(n = nReps, params, vcov.params))
  EstLenAtAge.sim = data.frame(matrix(nrow = nReps, ncol = nTimeSteps))
  EstLenAtAgeF.sim = EstLenAtAge.sim
  EstLenAtAgeM.sim = EstLenAtAge.sim

  for (j in 1:nReps) {
    ParamVals = exp(unlist(sims[j,]))
    if (is.vector(ObsRetCatchFreqAtLen)) { # combined sex
      if (SelectivityType == 1) { # input vector
        EstLenAtAge.sim[j,] = ParamVals[2] * (1 - exp(-ParamVals[3]*(DecAges)))
      }
      if (SelectivityType == 2) { # estimated
        EstLenAtAge.sim[j,] = ParamVals[4] * (1 - exp(-ParamVals[5]*(DecAges)))
      }
    }
    if (is.data.frame(ObsRetCatchFreqAtLen)) { # separate sexes
      if (SelectivityType == 1) { # input vector
        EstLenAtAgeF.sim[j,] = ParamVals[2] * (1 - exp(-ParamVals[4]*(DecAges)))
        EstLenAtAgeM.sim[j,] = ParamVals[3] * (1 - exp(-ParamVals[5]*(DecAges)))
      }
      if (SelectivityType == 2) { # estimated
        EstLenAtAgeF.sim[j,] = ParamVals[4] * (1 - exp(-ParamVals[6]*(DecAges)))
        EstLenAtAgeM.sim[j,] = ParamVals[5] * (1 - exp(-ParamVals[7]*(DecAges)))
      }
    }

    cat("j",j,'\n')
  }
  if (is.vector(ObsRetCatchFreqAtLen)) { # combined sex
    EstProp.sim = as.vector(apply(EstLenAtAge.sim, 2, median))
    EstProp.sim_low = as.vector(apply(EstLenAtAge.sim, 2, quantile, probs = 0.025))
    EstProp.sim_up = as.vector(apply(EstLenAtAge.sim, 2, quantile, probs = 0.975))
  }
  if (is.data.frame(ObsRetCatchFreqAtLen)) { # separate sexes
    EstPropF.sim = as.vector(apply(EstLenAtAgeF.sim, 2, median))
    EstPropF.sim_low = as.vector(apply(EstLenAtAgeF.sim, 2, quantile, probs = 0.025))
    EstPropF.sim_up = as.vector(apply(EstLenAtAgeF.sim, 2, quantile, probs = 0.975))
    EstPropM.sim = as.vector(apply(EstLenAtAgeM.sim, 2, median))
    EstPropM.sim_low = as.vector(apply(EstLenAtAgeM.sim, 2, quantile, probs = 0.025))
    EstPropM.sim_up = as.vector(apply(EstLenAtAgeM.sim, 2, quantile, probs = 0.975))
  }

  # AgeClasses = MinAge:MaxAge
  if (is.na(xaxis_lab)) xaxis_lab = "Age (y)"
  xlims = Get_xaxis_scale(DecAges)
  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint

  if (is.vector(ObsRetCatchFreqAtLen)) { # combined sex
    if (is.na(yaxis_lab)) yaxis_lab = "Length (mm)"
    ylims = Get_yaxis_scale(midpt)
    if (is.na(ymax)) ymax = ylims$ymax
    if (is.na(yint)) yint = ylims$yint
    plot(ObsAge, ObslenClMidPt, "p", main=MainLabel, cex.main=1.2, pch=16, cex=0.6,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax))
    points(ObsAge, ObslenClMidPt, col="black", cex=0.6)
    points(DecAges,EstProp.sim, col="red", cex=0.6)
    if (PlotCLs == TRUE) {
      x = c(DecAges,rev(DecAges)) # using shading for 95% CLs
      y = c(EstProp.sim_low, rev(EstProp.sim_up))
      polygon(x,y,col="pink",border=NA)
    }
    points(ObsAge, ObslenClMidPt, col="black", cex=0.6)
    points(DecAges, EstProp.sim, col="red", pch=1, cex=0.6)
    axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
    axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
    axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
  }
  if (is.data.frame(ObsRetCatchFreqAtLen)) { # combined sex
    # females
    ylims = Get_yaxis_scale(midpt)
    if (is.na(ymax)) ymax = ylims$ymax
    if (is.na(yint)) yint = ylims$yint
    if (is.na(yaxis_lab)) yaxis_lab1 = "Length females"
    plot(ObsAge_F, ObslenClMidPt_F, "p", main=MainLabel, cex.main=1.2, pch=16, cex=0.6,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2), ylab=list(yaxis_lab1,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax))
    points(ObsAge_F, ObslenClMidPt_F, col="black", cex=0.6)
    points(DecAges,EstPropF.sim, col="red", cex=0.6)
    if (PlotCLs == TRUE) {
      x = c(DecAges,rev(DecAges)) # using shading for 95% CLs
      y = c(EstPropF.sim_low, rev(EstPropF.sim_up))
      polygon(x,y,col="pink",border=NA)
    }
    points(ObsAge_F, ObslenClMidPt_F, col="black", cex=0.6)
    points(DecAges, EstPropF.sim, col="red", pch=1, cex=0.6)
    axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
    axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
    axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    # males
    if (is.na(yaxis_lab)) yaxis_lab2 = "Length males"
    ylims = Get_yaxis_scale(midpt)
    if (is.na(ymax)) ymax = ylims$ymax
    if (is.na(yint)) yint = ylims$yint
    plot(ObsAge_M, ObslenClMidPt_M, "p", main=MainLabel, cex.main=1.0, pch=16, cex=0.6,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2), ylab=list(yaxis_lab2,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax))
    points(ObsAge_M, ObslenClMidPt_M, col="black", cex=0.6)
    points(DecAges,EstPropM.sim, col="blue", cex=0.6)
    if (PlotCLs == TRUE) {
      x = c(DecAges,rev(DecAges)) # using shading for 95% CLs
      y = c(EstPropM.sim_low, rev(EstPropM.sim_up))
      polygon(x,y,col="light blue",border=NA)
    }
    points(ObsAge_M, ObslenClMidPt_M, col="black", cex=0.6)
    points(DecAges, EstPropM.sim, col="blue", pch=1, cex=0.6)
    axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
    axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
    axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
    axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
  }
}


#' Show estimated selectivity curve from age and length catch curve model
#'
#' This function provides a plot of the selectivity curve estimated from an age and length-based
#' catch curve model. The model is fitted to a sample of fish length and age data, by
#' minimising the overall negative log-likelihood, including the NLL associated with the
#' marginal length composition and a conditional age at length NLL,
#' given the parameters (selectivity, growth and mortality) and data, using nlminb.
#' It provides various statistical outputs in include convergence statistics, parameter estimates
#' and associated 95 percent confidence limits and associated variance-covariance matrix, calculated using
#' the MASS package.
#'
#' @param params vector of model parameters (in log space) to be estimated (if FittedRes=NA)
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param SelectivityType 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' @param ObsRetCatchFreqAtLen observed frequencies in length classes
#' @param ObsRetCatchFreqAtLengthAndAge observed frequencies in length and age classes
#' @param lbnd lower bounds of length classes
#' @param ubnd upper bounds of length classes
#' @param midpt mid points of length classes
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param SelectivityVec selectivity at length
#' @param DiscMort Proportion of fish that die following to capture and release
#' @param MaxAge maximum age considered in model
#' @param NatMort natural mortality
#' @param TimeStep model timestep (e.g. 1 = annual, 1/12 = monthly)
#' @param MainLabel plot label
#' @param xaxis_lab y axis label
#' @param yaxis_lab x axis label
#' @param xmax maximum x axis value
#' @param xint x axis interval
#' @param ymax maximum y axis value
#' @param yint y axis interval
#' @param PlotCLs logical (TRUE=plot 95 percent confidence limits for line)
#' @param FittedRes saved results from GetLengthBasedCatchCurveResults function (model will be refitted if set to NA)
#' @param nReps number of random parameter sets from parametric resampling to generate outputs with error
#'
#' @return Plot of selectivity curve estimated by an age and length-based catch curve
#'
#' @examples
#' # Simulate data
#' set.seed(123)
#' SampleSize=5000
#' MaxAge = 26
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' MinAge = floor(TimeStep)
#' nAgeCl = length(MinAge:MaxAge)
#' nTimeSteps = length(seq(TimeStep,MaxAge,TimeStep))
#' NatMort = 4.22/MaxAge
#' FishMort = 0.2
#' MaxLen = 1100
#' LenInc = 10
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified is knife-edged at MLL
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0
#' # # single sex, von Bertalanffy
#' # GrowthCurveType = 1 # 1 = von Bertalanffy
#' # Linf = 800
#' # vbK = 0.2
#' # CVSizeAtAge = 0.08
#' # GrowthParams = c(Linf, vbK, CVSizeAtAge)
#' # RefnceAges = NA
#' # 2 sexes, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy
#' Linf = c(700,850)
#' vbK = c(0.25,0.2)
#' CVSizeAtAge = c(0.05,0.05)
#' RefnceAges = NA
#' GrowthParams = data.frame(Linf=Linf, vbK=vbK, CVSizeAtAge=CVSizeAtAge)
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' lbnd=Res$lbnd
#' midpt=Res$midpt
#' ubnd=Res$ubnd
#' # # get data - 1 sex (or combined sexes)
#' # ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen # 1 sex
#' # ObsRetCatchFreqAtLengthAndAge = as.matrix(Res$ObsCatchFreqAtLengthAndDecAge) # 1 sex
#' # get data - 2 sexes
#' ObsRetCatchFreqAtLen <- data.frame(matrix(nrow = 2, ncol = length(midpt))) # 2 sex
#' colnames(ObsRetCatchFreqAtLen) <- midpt
#' ObsRetCatchFreqAtLen[1,] = Res$ObsRetCatchFreqAtLen_Fem
#' ObsRetCatchFreqAtLen[2,] = Res$ObsRetCatchFreqAtLen_Mal
#' ObsRetCatchFreqAtLengthAndAge = array(c(unlist(Res$ObsCatchFreqAtLengthAndDecAge_Fem), unlist(Res$ObsCatchFreqAtLengthAndDecAge_Mal)),
#'                                    c(nTimeSteps, length(midpt), 2), dimnames=list(rownames(Res$ObsCatchFreqAtLengthAndDecAge_Fem),
#'                                                                                   colnames(Res$ObsCatchFreqAtLengthAndDecAge_Fem)))
#' # # get params - 1 sex
#' # InitFishMort = 0.3 # specify starting parameters
#' # InitL50 = 320
#' # InitDelta = 50 # L95-L50
#' # InitLinf = 800
#' # InitvbK = 0.2
#' # InitCVSizeAtAge = 0.05
#' # get params - 2 sexes
#' InitFishMort = 0.3 # specify starting parameters
#' InitL50 = 320
#' InitDelta = 50 # L95-L50
#' InitLinf = c(800,800)
#' InitvbK = c(0.25,0.25)
#' InitCVSizeAtAge = 0.05
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' params = c(InitFishMort_logit, log(c(InitL50, InitDelta, InitLinf, InitvbK, InitCVSizeAtAge)))
#' FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
#' # Note, can skip above step and set FittedRes=NA (plot function will be slower)
#' # plot
#' PlotAgeLengthCatchCurve_Selectivity(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                     lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep, MainLabel=NA,
#'                                     xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA,
#'                                     ymax=NA, yint=NA, PlotCLs=TRUE, FittedRes, nReps=200)
#' @export
PlotAgeLengthCatchCurve_Selectivity <- function(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                                lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep, MainLabel,
                                                xaxis_lab, yaxis_lab, xmax, xint,
                                                ymax, yint, PlotCLs, FittedRes, nReps) {

  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    res =  FittedRes
  } else {
    res=GetAgeAndLengthBasedCatchCurveResults(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                              lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
  }

  if (SelectivityType == 1) {
    SelAtLength = res$SelAtLength
  }
  if (SelectivityType == 2) {
    params = res$params
    vcov.params = res$vcov.Params
    set.seed(123)
    sims = data.frame(MASS::mvrnorm(n = nReps, params, vcov.params))
    SelAtLength.sim = data.frame(matrix(nrow = nReps, ncol = length(midpt)))

    for (j in 1:nReps) {
      ParamVals = exp(unlist(sims[j,]))

      if (SelectivityType == 2) {
        SelAtLength.sim[j,] = 1 / (1 + exp(-log(19) * (midpt - ParamVals[2]) /
                                             (ParamVals[3])))
      }
      cat("j",j,'\n')
    }
    EstProp.sim = as.vector(apply(SelAtLength.sim, 2, median))
    EstProp.sim_low = as.vector(apply(SelAtLength.sim, 2, quantile, probs = 0.025))
    EstProp.sim_up = as.vector(apply(SelAtLength.sim, 2, quantile, probs = 0.975))

  }

  if (SelectivityType == 1) { # input selectivity as vector
    EstProp.sim = SelAtLength
    EstProp.sim_low = SelAtLength
    EstProp.sim_up = SelAtLength
  }

  if (is.na(xaxis_lab)) xaxis_lab = "Length (mm)"
  if (is.na(yaxis_lab)) yaxis_lab = "Selectivity"
  xlims = Get_xaxis_scale(midpt)
  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint
  if (is.na(ymax)) ymax = 1.2
  if (is.na(yint)) yint = 0.2

  plot(midpt, EstProp.sim, "p", main=MainLabel, cex.main=1.2, pch=1, cex=0.6,
       xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax),
       ylim=c(0,ymax), col="red")
  if (PlotCLs == TRUE) {
    # x = c(midpt,rev(midpt)) # using shading for 95% CLs
    # y = c(EstProp.sim_low, rev(EstProp.sim_up))
    # polygon(x,y,col="pink",border=NA)
    sm1 = spline(midpt, EstProp.sim_low, n=100, method="natural")
    sm2 = spline(midpt, EstProp.sim_up, n=100, method="natural")
    sm1$y[which(sm1$y<0)]=0
    sm2$y[which(sm1$y>1)]=1
    x = c(sm1$x, rev(sm2$x)) # using shading for 95% CLs
    y = c(sm1$y, rev(sm2$y))
    polygon(x,y, col="pink",border=NA)

  }
  sm1 = spline(midpt, EstProp.sim, n=100, method="natural")
  lines(sm1$x,sm1$y, col="red", cex=0.6)
  points(midpt, EstProp.sim, col="red", pch=1, cex=0.6)
  axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
  axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
  axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
  axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)

  if (SelectivityType==2) { # logistic selectivity
    L50est=paste("L50 =",round(exp(params[2]),0),"mm")
    L95est=paste("L95 =",round(exp(params[3]),0),"mm")
    legend("topleft", pch=-1, legend=c(L50est, L95est), lty="solid",col="black",
           bty='n', cex=0.6,lwd=-1, y.intersp=1.0)
  }
}

#' Show estimated proportions at age in each length class, from age and length catch curve model
#'
#' This function provides a plot of the estimated proportions at age in each length class,
#' as estimated from an age and length-based  catch curve model with length-based selectivity. The model is
#' fitted to a sample of fish length and age data, by minimising the overall negative log-likelihood,
#' including the NLL associated with the marginal length composition and a conditional age at length NLL,
#' given the parameters (selectivity, growth and mortality) and data, using nlminb.
#' It provides various statistical outputs in include convergence statistics, parameter estimates
#' and associated 95 percent  confidence limits and associated variance-covariance matrix, calculated using
#' the MASS package.
#'
#' @param params vector of model parameters (in log space) to be estimated (if FittedRes=NA)
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param SelectivityType 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' @param ObsRetCatchFreqAtLen observed frequencies in length classes
#' @param ObsRetCatchFreqAtLengthAndAge observed frequencies in length and age classes
#' @param lbnd lower bounds of length classes
#' @param ubnd upper bounds of length classes
#' @param midpt mid points of length classes
#' @param SelectivityVec selectivity at length
#' @param DiscMort proportion of fish that die following capture and release
#' @param MaxAge maximum age considered in model
#' @param NatMort natural mortality
#' @param TimeStep model timestep (e.g. 1 = annual, 1/12 = monthly)
#' @param MainLabel plot label
#' @param xaxis_lab y axis label
#' @param yaxis_lab x axis label
#' @param xmax maximum x axis value
#' @param xint x axis interval
#' @param ymax maximum y axis value
#' @param yint y axis interval
#' @param FittedRes saved results from GetLengthBasedCatchCurveResults function (model will be refitted if set to NA)
#'
#' @return Plot of fitted age and length-based catch curve to age at length data
#'
#' @examples
#' # Simulate data
#' set.seed(123)
#' SampleSize=5000
#' MaxAge = 26
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' MinAge = floor(TimeStep)
#' nAgeCl = length(MinAge:MaxAge)
#' nTimeSteps = length(seq(TimeStep,MaxAge,TimeStep))
#' NatMort = 4.22/MaxAge
#' FishMort = 0.2
#' MaxLen = 1100
#' LenInc = 10
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified is knife-edged at MLL
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0
#' # # single sex, von Bertalanffy
#' # GrowthCurveType = 1 # 1 = von Bertalanffy
#' # Linf = 800
#' # vbK = 0.2
#' # CVSizeAtAge = 0.08
#' # GrowthParams = c(Linf, vbK, CVSizeAtAge)
#' # RefnceAges = NA
#' # 2 sexes, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy
#' Linf = c(700,850)
#' vbK = c(0.25,0.2)
#' CVSizeAtAge = c(0.05,0.05)
#' RefnceAges = NA
#' GrowthParams = data.frame(Linf=Linf, vbK=vbK, CVSizeAtAge=CVSizeAtAge)
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' lbnd=Res$lbnd
#' midpt=Res$midpt
#' ubnd=Res$ubnd
#' # # get data - 1 sex (or combined sexes)
#' # ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen # 1 sex
#' # ObsRetCatchFreqAtLengthAndAge = as.matrix(Res$ObsCatchFreqAtLengthAndDecAge) # 1 sex
#' # get data - 2 sexes
#' ObsRetCatchFreqAtLen <- data.frame(matrix(nrow = 2, ncol = length(midpt))) # 2 sex
#' colnames(ObsRetCatchFreqAtLen) <- midpt
#' ObsRetCatchFreqAtLen[1,] = Res$ObsRetCatchFreqAtLen_Fem
#' ObsRetCatchFreqAtLen[2,] = Res$ObsRetCatchFreqAtLen_Mal
#' ObsRetCatchFreqAtLengthAndAge = array(c(unlist(Res$ObsCatchFreqAtLengthAndDecAge_Fem), unlist(Res$ObsCatchFreqAtLengthAndDecAge_Mal)),
#'                                    c(nTimeSteps, length(midpt), 2), dimnames=list(rownames(Res$ObsCatchFreqAtLengthAndDecAge_Fem),
#'                                                                                   colnames(Res$ObsCatchFreqAtLengthAndDecAge_Fem)))
#' # # get params - 1 sex
#' # InitFishMort = 0.3 # specify starting parameters
#' # InitL50 = 320
#' # InitDelta = 50 # L95-L50
#' # InitLinf = 800
#' # InitvbK = 0.2
#' # InitCVSizeAtAge = 0.05
#' # get params - 2 sexes
#' InitFishMort = 0.3 # specify starting parameters
#' InitL50 = 320
#' InitDelta = 50 # L95-L50
#' InitLinf = c(800,800)
#' InitvbK = c(0.25,0.25)
#' InitCVSizeAtAge = 0.05
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' params = c(InitFishMort_logit, log(c(InitL50, InitDelta, InitLinf, InitvbK, InitCVSizeAtAge)))
#' FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
#' # Plot. Note, can skip above step and set FittedRes=NA (plot function will be slower)
#' PlotAgeLengthCatchCurve_Cond_AL(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep, MainLabel=NA,
#'                                 xaxis_lab=NA, yaxis_lab=NA, xmax=40, xint=10,
#'                                 ymax=NA, yint=NA, FittedRes)
#' @export
PlotAgeLengthCatchCurve_Cond_AL <- function(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                            lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep, MainLabel,
                                            xaxis_lab, yaxis_lab, xmax, xint, ymax, yint, FittedRes) {

  .pardefault <- par(no.readonly = TRUE) # store current par settings

  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    res =  FittedRes
  } else {
    res=GetAgeAndLengthBasedCatchCurveResults(params, MLL, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                              lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
  }

  if (is.na(xaxis_lab)) xaxis_lab = "AgeClass (y)"
  if (is.na(yaxis_lab)) yaxis_lab = "Proportion"
  MinAge = floor(TimeStep)
  AgeClasses = MinAge:MaxAge
  nAgeCl = length(MinAge:MaxAge)
  xlims = Get_xaxis_scale(AgeClasses)
  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint
  if (is.na(ymax)) ymax = 1.0
  if (is.na(yint)) yint = 0.2
  nLenCl = length(midpt)

  # predicted conditional proportion of age given length
  # single sex
  if (is.vector(ObsRetCatchFreqAtLen)) {
    ExpRetCatchPropAtIntAge=res$ExpRetCatchPropAtIntAge
    ExpRetCatchPropLengthGivenIntAge=res$ExpRetCatchPropLengthGivenIntAge
    ExpRetCatchPropIntAgeGivenLength <- data.frame(matrix(nrow = nAgeCl, ncol = nLenCl))
    colnames(ExpRetCatchPropIntAgeGivenLength) <- midpt
    ObsCatchPropAgeAtLength = ExpRetCatchPropIntAgeGivenLength
    ExpRetCatchPropIntAgeGivenLength = CalcExpCatchPropIntAgeGivenLength_cpp(nLenCl, nAgeCl, ExpRetCatchPropLengthGivenIntAge, ExpRetCatchPropAtIntAge)
  }

  # 2 sexes
  if (is.data.frame(ObsRetCatchFreqAtLen)) {
    ExpRetCatchPropLengthGivenIntAge_Fem=res$ExpRetCatchPropLengthGivenIntAge_Fem
    ExpRetCatchPropLengthGivenIntAge_Mal=res$ExpRetCatchPropLengthGivenIntAge_Mal
    ExpRetCatchPropAtIntAge_Fem=res$ExpRetCatchPropAtIntAge_Fem
    ExpRetCatchPropAtIntAge_Mal=res$ExpRetCatchPropAtIntAge_Mal
    FractDenom_Fem = rep(0,nLenCl)
    FractDenom_Mal = rep(0,nLenCl)
    ExpRetCatchPropIntAgeGivenLength_Fem <- data.frame(matrix(nrow = nAgeCl, ncol = nLenCl))
    colnames(ExpRetCatchPropIntAgeGivenLength_Fem) <- midpt
    ExpRetCatchPropIntAgeGivenLength_Mal = ExpRetCatchPropIntAgeGivenLength_Fem
    ObsCatchPropAgeAtLength_Fem = ExpRetCatchPropIntAgeGivenLength_Fem
    ObsCatchPropAgeAtLength_Mal = ExpRetCatchPropIntAgeGivenLength_Fem
    ExpRetCatchPropIntAgeGivenLength_Fem = CalcExpCatchPropIntAgeGivenLength_cpp(nLenCl, nAgeCl, ExpRetCatchPropLengthGivenIntAge_Fem, ExpRetCatchPropAtIntAge_Fem)
    ExpRetCatchPropIntAgeGivenLength_Mal = CalcExpCatchPropIntAgeGivenLength_cpp(nLenCl, nAgeCl, ExpRetCatchPropLengthGivenIntAge_Mal, ExpRetCatchPropAtIntAge_Mal)
  }

  # plot observed and expected age proportions for each length class
  # single sex
  if (is.vector(ObsRetCatchFreqAtLen)) {
    if (TimeStep == 1) {
      ObsCatchFreqAtLengthAndIntAge = ObsRetCatchFreqAtLengthAndAge
    } else {
      ObsCatchFreqAtLengthAndIntAge = ConvertObsDataFromDecAgesToIntegerAges(TimeStep, MaxAge, nLenCl, ObsRetCatchFreqAtLengthAndAge)
    }
    par(mfcol=c(3,3), mar=c(3.5,3.5,1,1), oma=c(1,1,1,0), tck=-0.03)
    k=0
    for (i in 1:nLenCl) {
      k=k+1
      if (k==10) k=1

      ObsCatchPropAgeAtLength[,i] = ObsCatchFreqAtLengthAndIntAge[,i] / sum(ObsCatchFreqAtLengthAndIntAge[,i])
      plot(1:nAgeCl, ObsCatchPropAgeAtLength[,i], "p", main='', cex.main=1.0,
           pch=16, cex=0.8, cex.main=0.8, xaxt = "n", yaxt = "n", xlab=NA,
           ylab=NA, frame=F, xlim=c(0,xmax), ylim=c(0,ymax), col="black")
      lines(1:nAgeCl, ExpRetCatchPropIntAgeGivenLength[,i], col="red")
      axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
      axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)

      if (k==1 | k==2 | k==3) {
        axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
        mtext(yaxis_lab,las=3,side=2,line=3, cex=0.8)
      }
      if (k==3 | k==6 | k==9) {
        axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
        mtext(xaxis_lab,las=1,side=1,line=3,adj=0.5, cex=0.8)
      }
      if (is.na(MainLabel)) {
        MainLab=paste(lbnd[i],"-",ubnd[i],"mm")
      } else {
        MainLab=MainLabel
      }
      mtext(MainLab,las=1,side=3,line=0,adj=0.5, cex=0.6)
      legend("topright", legend=c("Obs","Exp"), lty="solid", lwd=c(-1,1),
             pch=c(16,-1), col=c("black","red"), bty='n', cex=0.8)
    }
  }

  # single sex
  if (is.data.frame(ObsRetCatchFreqAtLen)) {
    if (TimeStep == 1) {
      # dim(ObsRetCatchFreqAtLengthAndAge)
      ObsCatchFreqAtLengthAndIntAge_Fem = ObsRetCatchFreqAtLengthAndAge[,,1]
      ObsCatchFreqAtLengthAndIntAge_Mal = ObsRetCatchFreqAtLengthAndAge[,,2]
    } else {
      ObsRetCatchFreqAtLengthAndAge_Fem = as.matrix(ObsRetCatchFreqAtLengthAndAge[,,1])
      ObsRetCatchFreqAtLengthAndAge_Mal = as.matrix(ObsRetCatchFreqAtLengthAndAge[,,2])
      ObsCatchFreqAtLengthAndIntAge_Fem = ConvertObsDataFromDecAgesToIntegerAges(TimeStep, MaxAge, nLenCl, ObsRetCatchFreqAtLengthAndAge_Fem)
      ObsCatchFreqAtLengthAndIntAge_Mal = ConvertObsDataFromDecAgesToIntegerAges(TimeStep, MaxAge, nLenCl, ObsRetCatchFreqAtLengthAndAge_Mal)
    }

    # females
    par(mfcol=c(3,3), mar=c(3.5,3.5,1,1), oma=c(1,1,1,0), tck=-0.03)
    k=0
    for (i in 1:nLenCl) {
      k=k+1
      if (k==10) k=1
      ObsCatchPropAgeAtLength_Fem[,i] = ObsCatchFreqAtLengthAndIntAge_Fem[,i] / sum(ObsCatchFreqAtLengthAndIntAge_Fem[,i])
      plot(AgeClasses, ObsCatchPropAgeAtLength_Fem[,i], "p", main='', cex.main=1.0,
           pch=16, cex=0.8, cex.main=0.8, xaxt = "n", yaxt = "n", xlab=NA,
           ylab=NA, frame=F, xlim=c(0,xmax), ylim=c(0,ymax), col="black")
      lines(AgeClasses, ExpRetCatchPropIntAgeGivenLength_Fem[,i], col="red")
      axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
      axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)

      if (k==1 | k==2 | k==3) {
        axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
        mtext(yaxis_lab,las=3,side=2,line=3, cex=0.8)
      }
      if (k==3 | k==6 | k==9) {
        axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
        mtext(xaxis_lab,las=1,side=1,line=3,adj=0.5, cex=0.8)
      }
      if (is.na(MainLabel)) {
        MainLab=paste(lbnd[i],"-",ubnd[i],"mm")
      } else {
        MainLab=MainLabel
      }
      mtext(MainLab,las=1,side=3,line=0,adj=0.5, cex=0.6)
      legend("topright", legend=c("Fem_Obs","Fem_Exp"), lty="solid", lwd=c(-1,1),
             pch=c(16,-1), col=c("black","red"), bty='n', cex=0.8)
    }
    # males
    par(mfcol=c(3,3), mar=c(3.5,3.5,1,1), oma=c(1,1,1,0), tck=-0.03)
    k=0
    for (i in 1:nLenCl) {
      k=k+1
      if (k==10) k=1
      ObsCatchPropAgeAtLength_Mal[,i] = ObsCatchFreqAtLengthAndIntAge_Mal[,i] / sum(ObsCatchFreqAtLengthAndIntAge_Mal[,i])
      plot(AgeClasses, ObsCatchPropAgeAtLength_Mal[,i], "p", main='', cex.main=1.0,
           pch=16, cex=0.8, cex.main=0.8, xaxt = "n", yaxt = "n", xlab=NA,
           ylab=NA, frame=F, xlim=c(0,xmax), ylim=c(0,ymax), col="black")
      lines(AgeClasses, ExpRetCatchPropIntAgeGivenLength_Mal[,i], col="blue")
      axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
      axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)

      if (k==1 | k==2 | k==3) {
        axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
        mtext(yaxis_lab,las=3,side=2,line=3, cex=0.8)
      }
      if (k==3 | k==6 | k==9) {
        axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
        mtext(xaxis_lab,las=1,side=1,line=3,adj=0.5, cex=0.8)
      }
      if (is.na(MainLabel)) {
        MainLab=paste(lbnd[i],"-",ubnd[i],"mm")
      } else {
        MainLab=MainLabel
      }
      mtext(MainLab,las=1,side=3,line=0,adj=0.5, cex=0.6)
      legend("topright", legend=c("Mal_Obs","Mal_Exp"), lty="solid", lwd=c(-1,1),
             pch=c(16,-1), col=c("black","blue"), bty='n', cex=0.8)
    }
  }

  # reset default par options
  par(.pardefault)
}


#' Calculate values for length-converted catch curve with von Bertalanffy growth
#'
#' Return outputs associated with length-converted catch curve (Pauly, 1990), required to fit the model
#'
#' @keywords internal
#'
#' @param GrowthParams c(Linf,vbK) von Bertalanffy growth parameters
#' @param ObsRetCatchFreqAtLen observed frequencies in length classes
#' @param MinFreq minimum frequency to include
#' @param lbnd lower bounds of length classes
#' @param ubnd upper bounds of length classes
#' @param midpt mid points of length classes
#'
#' @return length class with peak frequency (PeakLencl), last length class included in analysis (LastLenCl),
#' observed data to which catch curve is fitted (ObsCatchFreqAtLen2), relatives age at midpoints of length classes
#' (Age_midptlencl), natural logarithms of changes in numbers with respect to time taken to grow through
#' each length class (ln_n_dt)
#'
Calcs_PaulyLenConvertCatchCurve <- function(GrowthParams, ObsRetCatchFreqAtLen,
                                            MinFreq, lbnd, midpt, ubnd) {

  Linf = GrowthParams[1]
  vbK = GrowthParams[2]

  # determine length class at peak frequency. If there is more than one length class
  # with the same "peak" frequency, choose the larger one
  PeakFrequency <- max(ObsRetCatchFreqAtLen)
  PeakLencl <- which(ObsRetCatchFreqAtLen==PeakFrequency)
  if (length(PeakLencl)>1) PeakLencl = max(PeakLencl)

  # Select all length classes with > specified number of required observations
  # where length category is below Linf. (Need to remove categories for which
  # ln(n*dl/dt) becomes negative)
  nLenCl = length(lbnd)
  k = PeakLencl
  while (ObsRetCatchFreqAtLen[k] >= MinFreq & ubnd[k] < Linf)  {
    LastLenCl <- k
    k = k + 1
    if (k > nLenCl)
      stop("problem with while loop")
  }

  # number of length categories for analysis
  nLenCl_CC <- length(ObsRetCatchFreqAtLen[PeakLencl:LastLenCl])

  # data to which catch curve is fitted
  ObsCatchFreqAtLen2 <- rep(NA,nLenCl)
  ObsCatchFreqAtLen2[PeakLencl:LastLenCl] <- ObsRetCatchFreqAtLen[PeakLencl:LastLenCl]

  # calc Age_midptlencl, delta_t and Obs_ln_n_dt
  Age_lbndlencl <- rep(0,nLenCl_CC)
  Age_ubndlencl <- rep(0,nLenCl_CC)
  Age_midptlencl <- rep(0,nLenCl_CC)
  DeltaT_yrs <- rep(0,nLenCl_CC)
  Obs_ln_n_dt <- rep(0,nLenCl_CC)
  i=0
  for (j in PeakLencl:LastLenCl) {
    i=i+1
    Age_lbndlencl[i] = (-1 / vbK * log(1 - lbnd[j] / Linf))
    Age_ubndlencl[i] = (-1 / vbK * log(1 - ubnd[j] / Linf))
    Age_midptlencl[i] = (-1 / vbK * log(1 - midpt[j] / Linf))
    DeltaT_yrs[i] = Age_ubndlencl[i] - Age_lbndlencl[i]
    Obs_ln_n_dt[i] = log(ObsRetCatchFreqAtLen[j] / DeltaT_yrs[i])
  }

  results = list(PeakLencl = PeakLencl,
                 LastLenCl = LastLenCl,
                 ObsCatchFreqAtLen2 = ObsCatchFreqAtLen2,
                 Age_midptlencl = Age_midptlencl,
                 Obs_ln_n_dt = Obs_ln_n_dt)

  return(results)
}

#' Calculate values for length-converted catch curve with Schnute growth
#'
#' Return outputs associated with length-converted catch curve (see Pauly, 1990), with growth described using the
#' Schnute (1981) growth model, required to fit the model
#'
#' @keywords internal
#'
#' @param GrowthParams c(y1,y2,a,b) Schnute growth curve parameters
#' @param RefnceAges c(t1,t2) Schnute growth curve reference ages
#' @param ObsRetCatchFreqAtLen observed frequencies in length classes
#' @param MinFreq minimum frequency to include
#' @param lbnd lower bounds of length classes
#' @param midpt mid points of length classes
#' @param ubnd upper bounds of length classes
#'
#' @return length class with peak frequency (PeakLencl), last length class included in analysis (LastLenCl),
#' observed data to which catch curve is fitted (ObsCatchFreqAtLen2), relatives age at midpoints of length classes
#' (Age_midptlencl), natural logarithms of changes in numbers with respect to time taken to grow through
#' each length class (ln_n_dt)
#'
Calcs_LenCovertCatchCurve_Schnute <- function(GrowthParams, RefnceAges, ObsRetCatchFreqAtLen,
                                                  MinFreq, lbnd, midpt, ubnd) {


  t1 = RefnceAges[1]
  t2 = RefnceAges[2]
  y1 = GrowthParams[1]
  y2 = GrowthParams[2]
  a = GrowthParams[3]
  b = GrowthParams[4]

  # determine length class at peak frequency. If there is more than one length class
  # with the same "peak" frequency, choose the larger one
  PeakFrequency <- max(ObsRetCatchFreqAtLen)
  PeakLencl <- which(ObsRetCatchFreqAtLen==PeakFrequency)
  if (length(PeakLencl)>1) PeakLencl = max(PeakLencl)

  # Select all length classes with > specified number of required observations
  # where length category is below Linf. (Need to remove categories for which
  # ln(n*dl/dt) becomes negative)
  Linf = ((exp(a*t2)*y2^b-exp(a*t1)*y1^b)/(exp(a*t2)-exp(a*t1)))^(1/b)
  nLenCl = length(lbnd)
  k = PeakLencl
  LastLenCl <- k
  for (i in ObsRetCatchFreqAtLen[k:length(ObsRetCatchFreqAtLen)]) {
    if (ObsRetCatchFreqAtLen[k] >= MinFreq) {
      if (ubnd[k]<Linf) {
        LastLenCl <- k
        k = k + 1
      }}
  }

  # number of length categories for analysis
  nLenCl_CC <- length(ObsRetCatchFreqAtLen[PeakLencl:LastLenCl])

  # data to which catch curve is fitted
  ObsCatchFreqAtLen2 <- rep(NA,nLenCl)
  ObsCatchFreqAtLen2[PeakLencl:LastLenCl] <- ObsRetCatchFreqAtLen[PeakLencl:LastLenCl]

  # calc Age_midptlencl, delta_t and ln_n_dt
  Age_lbndlencl <- rep(0,nLenCl_CC)
  Age_ubndlencl <- rep(0,nLenCl_CC)
  Age_midptlencl <- rep(0,nLenCl_CC)
  Age_midptlencl <- rep(0,nLenCl_CC)
  DeltaT_yrs <- rep(0,nLenCl_CC)
  Obs_ln_n_dt <- rep(0,nLenCl_CC)
  i=0
  for (j in PeakLencl:LastLenCl) {
    i=i+1

    Age_lbndlencl[i] = InverseSchnuteGrowthfunction(lbnd[j], t1, t2, y1, y2, a, b)
    Age_ubndlencl[i] = InverseSchnuteGrowthfunction(ubnd[j], t1, t2, y1, y2, a, b)
    Age_midptlencl[i] = InverseSchnuteGrowthfunction(midpt[j], t1, t2, y1, y2, a, b)
    DeltaT_yrs[i] = max(0,Age_ubndlencl[i] - Age_lbndlencl[i])
    Obs_ln_n_dt[i] = log(ObsRetCatchFreqAtLen[j] / DeltaT_yrs[i])
  }


  results = list(PeakLencl = PeakLencl,
                 LastLenCl = LastLenCl,
                 ObsCatchFreqAtLen2=ObsCatchFreqAtLen2,
                 Age_midptlencl=Age_midptlencl,
                 Obs_ln_n_dt=Obs_ln_n_dt)

  return(results)
}

#' Fit length-converted catch model and return results
#'
#' Fit specified length-converted catch curve model with growth described either by a von
#' Bertalanffy growth model or Schunute growth model (with a and b not equal to zero)
#'
#' @param ModelType 1=von Bertalanffy growth, 2=Schnute growth
#' @param GrowthParams c(Linf, vbK) von Bertalanffy growth or c(y1, y2, a, b) Schnute growth
#' @param RefnceAges Schnute growth curve reference ages
#' @param ObsRetCatchFreqAtLen observed frequencies in length classes
#' @param MinFreq minimum frequency to include
#' @param lbnd lower bounds of length classes
#' @param midpt mid points of length classes
#' @param ubnd upper bounds of length classes
#'
#' @return sample size for observed data (SampleSize), Estimates of parameters and associated 95 percent
#' confidence limits (ParamEst), length class with greatest frequency of fish (PeakLencl), last length
#' class considered in analysis (LastLenCl), observed frequency data to which catch curve is fitted (ObsCatchFreqAtLen2),
#' relative ages at mid points of length classes (Age_midptlencl), natural logarithms of observed and expected numbers in length
#' classes divided by time interval taken to grow through length class (Obs_ln_n_dt, Est_ln_n_dt), lower and upper
#' 95 percent confidence limits for (Est_ln_n_dtlow, Est_ln_n_dtup)
#'
#' @examples
#' # Simulate data
#' SampleSize=1000
#' set.seed(123)
#' MaxAge = 30
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = 0.2
#' MaxLen = 1500
#' LenInc = 50
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified is knife-edged at MLL
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 500) # L50, L95 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95 for retention
#' DiscMort = 0
#' # # von Bertalanffy
#' # GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' # Linf = 800
#' # vbK = 0.2
#' # CVSizeAtAge = 0.08
#' # GrowthParams = c(Linf, vbK)
#' # RefnceAges = NA
#' # Schnute
#' GrowthCurveType = 2 # 1 = von Bertalanffy, 2 = Schnute
#' t1 = 1 # growth - Schnute
#' t2 = 10 # growth - Schnute
#' y1 = 400 # growth - Schnute
#' y2 = 1000 # growth - Schnute
#' a = 0.1 # growth - Schnute
#' b = 2.0 # growth - Schnute
#' CVSizeAtAge = 0.08
#' GrowthParams = c(y1, y2, a, b)
#' RefnceAges = c(t1,t2)
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#'
#' ObsRetCatchFreqAtLen = as.vector(Res$ObsRetCatchFreqAtLen)
#' MinFreq = 20 # set minimum frequency for larger lengths for analysis
#' # note, this needs to be high enough so that data for ln(n/dt) vs relative age for essentially straight
#' # line - if not, Z will be biased!!!
#' midpt=Res$midpt
#' lbnd=Res$lbnd
#' ubnd=Res$ubnd
#' ModelType = 2 # 1 = von Bertalanffy growth curve (Pauly), 2 = length-converted catch curve - Schnute growth curve)
#' res=GetLenConvCatchCurveResults(ModelType, GrowthParams, RefnceAges, ObsRetCatchFreqAtLen,
#'                                 MinFreq, lbnd, midpt, ubnd)
#' @export
GetLenConvCatchCurveResults <- function(ModelType, GrowthParams, RefnceAges, ObsRetCatchFreqAtLen,
                                        MinFreq, lbnd, midpt, ubnd) {

  # Pauly's length-converted catch curve
  if (ModelType == 1) {
    res=Calcs_PaulyLenConvertCatchCurve(GrowthParams, ObsRetCatchFreqAtLen,
                                        MinFreq, lbnd, midpt, ubnd)
  }
  if (ModelType == 2) { # length conv catch curve using Schnute function
    res=Calcs_LenCovertCatchCurve_Schnute(GrowthParams, RefnceAges, ObsRetCatchFreqAtLen,
                                          MinFreq, lbnd, midpt, ubnd)
  }

  # run the linear catch curve model.
  Age_midptlencl=res$Age_midptlencl
  Obs_ln_n_dt=res$Obs_ln_n_dt
  mod1 <- lm(Obs_ln_n_dt~Age_midptlencl)

  # get linear model parameters
  yint <- as.numeric(mod1$coefficients[1])
  ZMort <- - as.numeric(mod1$coefficients[2])

  # get standard errors for the parameters
  yint_se <- as.numeric(summary(mod1)$coef[3])
  ZMort_se <- as.numeric(summary(mod1)$coef[4])

  # get approx 95% confidence limits for parameters
  Est_ZMort = c(ZMort, ZMort + c(-1.96, 1.96) * ZMort_se)
  Est_yIntercept = c(yint, yint + c(-1.96, 1.96) * yint_se)
  ParamEst = t(data.frame(Est_ZMort=round(Est_ZMort,2), Est_yIntercept=round(Est_yIntercept,2)))

  # store results in data frame
  colnames(ParamEst) = c("Estimate","lw_95%CL","up_95%CL")

  yint_low <- yint - (1.96 * yint_se)
  yint_up <- yint + (1.96 * yint_se)
  ZMort_low <- ZMort - (1.96 * ZMort_se)
  ZMort_up <- ZMort + (1.96 * ZMort_se)

  # get 95% confidence limits for fitted line
  Age_midptlencl=res$Age_midptlencl
  newdat <- data.frame(Age_midptlencl=Age_midptlencl)
  conf_int = predict(mod1, newdata = newdat, interval = 'confidence')
  Est_ln_n_dt = conf_int[,1]
  Est_ln_n_dtlow = conf_int[,2]
  Est_ln_n_dtup = conf_int[,3]

  # get additional info from fitted model
  Age_midptlencl = res$Age_midptlencl
  ObsCatchFreqAtLen2 = res$ObsCatchFreqAtLen2
  PeakLencl = res$PeakLencl
  LastLenCl = res$LastLenCl
  Obs_ln_n_dt = res$Obs_ln_n_dt
  nLenCl_CC <- length(res$PeakLencl:res$LastLenCl)

  # store sample size
  SampleSize = sum(ObsRetCatchFreqAtLen)

  results = list(SampleSize = SampleSize,
                 ParamEst = ParamEst,
                 PeakLencl = PeakLencl,
                 LastLenCl = LastLenCl,
                 ObsCatchFreqAtLen2 = ObsCatchFreqAtLen2,
                 Age_midptlencl = Age_midptlencl,
                 Obs_ln_n_dt = Obs_ln_n_dt,
                 Est_ln_n_dt = Est_ln_n_dt,
                 Est_ln_n_dtlow = Est_ln_n_dtlow,
                 Est_ln_n_dtup = Est_ln_n_dtup)

  return(results)
}

#' Plot results for specified length-converted catch model
#'
#' Plot results for specified length-converted catch curve model, including either Pauly's method,
#' or a model allowing for Schnute growth with  logistic selectivity or knife edge selectivity
#'
#' @param MaxAge maximum age
#' @param ModelType 1=von Bertalanffy growth, 2=Schnute growth
#' @param GrowthParams c(Linf, vbK) von Bertalanffy growth or c(y1, y2, a, b) Schnute growth
#' @param RefnceAges Schnute growth curve reference ages
#' @param ObsRetCatchFreqAtLen observed frequencies in length classes
#' @param MinFreq minimum frequency to include
#' @param lbnd lower bounds of length classes
#' @param midpt mid points of length classes
#' @param ubnd upper bounds of length classes
#'
#' @return plots associated with fitted length-converted catch curve models
#'
#' @examples
#' # Simulate data
#' SampleSize=1000
#' set.seed(123)
#' MaxAge = 30
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = 0.2
#' MaxLen = 1500
#' LenInc = 50
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified is knife-edged at MLL
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 500) # L50, L95 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95 for retention
#' DiscMort = 0
#' # # von Bertalanffy
#' # GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' # Linf = 800
#' # vbK = 0.2
#' # CVSizeAtAge = 0.08
#' # GrowthParams = c(Linf, vbK)
#' # RefnceAges = NA
#' # Schnute
#' GrowthCurveType = 2 # 1 = von Bertalanffy, 2 = Schnute
#' t1 = 1 # growth - Schnute
#' t2 = 10 # growth - Schnute
#' y1 = 400 # growth - Schnute
#' y2 = 1000 # growth - Schnute
#' a = 0.1 # growth - Schnute
#' b = 2.0 # growth - Schnute
#' CVSizeAtAge = 0.08
#' GrowthParams = c(y1, y2, a, b)
#' RefnceAges = c(t1,t2)
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' ObsRetCatchFreqAtLen = as.vector(Res$ObsRetCatchFreqAtLen)
#' MinFreq = 20 # set minimum frequency for larger lengths for analysis
#' # note, this needs to be high enough so that data for ln(n/dt) vs relative age for essentially straight
#' # line - if not, Z will be biased!!!
#' midpt=Res$midpt
#' lbnd=Res$lbnd
#' ubnd=Res$ubnd
#' ModelType = 2 # 1 = von Bertalanffy growth curve (Pauly), 2 = length-converted catch curve - Schnute growth curve)
#' res=GetLenConvCatchCurveResults(ModelType, GrowthParams, RefnceAges, ObsRetCatchFreqAtLen,
#'                                 MinFreq, lbnd, midpt, ubnd)
#' PlotLenConvCatchCurveResults(MaxAge, ModelType, GrowthParams, RefnceAges, ObsRetCatchFreqAtLen, MinFreq,
#'                             lbnd, midpt, ubnd)
#' @export
#'
#'
PlotLenConvCatchCurveResults <- function(MaxAge, ModelType, GrowthParams, RefnceAges, ObsRetCatchFreqAtLen, MinFreq,
                                         lbnd, midpt, ubnd) {

  # .pardefault <- par(no.readonly = TRUE) # store default par settings

    res=GetLenConvCatchCurveResults(ModelType, GrowthParams, RefnceAges, ObsRetCatchFreqAtLen,
                                    MinFreq, lbnd, midpt, ubnd)
    Ages = 1:MaxAge

    if (ModelType == 1) { # Pauly length converted catch curve
      Linf = GrowthParams[1]
      vbK = GrowthParams[2]
      MeanSizeAtAge = Linf * (1 - exp (-vbK * (1:MaxAge))) # von Bertalanffy growth curve
    }

    if (ModelType == 2) { # Length converted catch curve with Schnute growth curve

      MeanSizeAtAge = rep(0,MaxAge)
      t1=RefnceAges[1]
      t2=RefnceAges[2]
      y1=GrowthParams[1]
      y2=GrowthParams[2]
      a=GrowthParams[3]
      b=GrowthParams[4]

      for (t in 1:MaxAge) {
        MeanSizeAtAge[t] = SchnuteGrowthfunction(t, t1, t2, y1, y2, a, b)
      }
    }

    # par(mfrow = c(2,2), mar=c(4,4,0.1,0.1), oma=rep(0.1,4), tck=-0.03, mgp = c(3,1,0))

    # plot data that can be used for the catch curve,
    # relative to full set of length data
    plot(midpt, ObsRetCatchFreqAtLen, ylab=list("Frequency",cex=1.2), xlab=list("Length",cex=1.2),
         cex.lab = 1, frame.plot=F, xlim=c(0,MaxLen), las=1)
    points(midpt, res$ObsCatchFreqAtLen2, pch=16, col="blue")
    legend('topright', col=c("grey","blue"),legend=c("all","analysis"),
           bty='n', cex=0.8,lwd=1.75)

    # plot growth curve, and overlay range of lengths available for catch curve analysis
    plot(Ages,MeanSizeAtAge,"l",frame.plot=F, xlim=c(0,max(Ages)),
         ylim=c(0,MaxLen), cex.lab = 1, col="blue", ylab=list("Length",cex=1.2),xlab=list("Age",cex=1.2), las=1)
    abline(h=lbnd[res$PeakLencl],lty="solid")
    abline(h=ubnd[res$LastLenCl],lty="dotted")
    legend('bottomright', col=c("black","black"),lty=c("solid","dotted"),
           legend=c("start of catch curve","end of catch curve"),bty='n', cex=0.8,lwd=1.75)


    # plot growth curve, and overlay range of lengths available for catch curve analysis
    ymax=1.35*max(res$Obs_ln_n_dt)
    xmax=1.35*max(res$Age_midptlencl)
    plot(res$Age_midptlencl,res$Obs_ln_n_dt,"p",pch=16,frame.plot=F, xlim=c(0,xmax),
         ylim=c(0,ymax), cex.lab = 1, col="black", ylab=list("ln(n/dt)",cex=1.2),xlab=list("Relative age",cex=1.2), las=1)

    Zval = round(res$ParamEst[1,1],2)
    Zest = bquote("Z =" ~ .(Zval) ~ y^-1)
    legend("bottomleft", pch=-1, legend=as.expression(Zest),
           lty="solid",col="black",
           bty='n', cex=0.8,lwd=-1, y.intersp=1, adj=0)
    x = c(res$Age_midptlencl,rev(res$Age_midptlencl)) # using shading for 95% CLs
    y = c(res$Est_ln_n_dtlow, rev(res$Est_ln_n_dtup))
    polygon(x,y,col="pink",border=NA)
    points(res$Age_midptlencl,res$Obs_ln_n_dt,pch=16)
    lines(res$Age_midptlencl, res$Est_ln_n_dt)
    # reset default par options
    # par(.pardefault)

}


#*****************************************
# Catch curve analyses - age-based methods
#*****************************************

#' Determine age at full recruitment into the fishery (assuming knife-edge selection)
#'
#' This function determines the age at full recruitment, based on the age at
#' peak frequency of fish. If RecAssump=0, RecAge is set to the age at
#' peak frequency, but if RecAssump=1, RecAge is set to one year above the age at
#' peak frequency
#'
#' @keywords internal
#'
#' @param RecAssump 0=age at peak frequency, 1=age at peak frequency + 1
#' @param Ages ages for analysis
#' @param ObsAgeFreq observed frequency at age
#'
#' @return age at full recruitment into the fishery (RecAge)
#'
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
#' Res=SimAgeFreqData(SampleSize, MinAge, MaxAge, SelA50, SelA95, NatMort, FMort)
#' ObsAgeFreq = unlist(as.vector(Res$CatchSample))
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


#' Determine age at full recruitment into the fishery (assuming knife-edge selection)
#'
#' This function calculates the oldest age to be included in data for a linear catch curve analysis,
#' according to a specified minimum frequency
#'
#' @keywords internal
#'
#' @param MinFreq minimum frequency of fish for including data for old fish
#' @param RecAge recruitment age
#' @param Ages ages for analysis
#' @param ObsAgeFreq observed age frequency
#'
#' @return oldest age to be included for linear catch curve analysis (LastAgeForLinearCC)
#'
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
#' Res=SimAgeFreqData(SampleSize, MinAge, MaxAge, SelA50, SelA95, NatMort, FMort)
#' ObsAgeFreq = unlist(as.vector(Res$CatchSample))
#' RecAge = CalcRecruitmentAge(RecAssump=0, Ages, ObsAgeFreq)
#' CalcLastAgeForLinearCatchCurve(MinFreq=1, RecAge, Ages, ObsAgeFreq)
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


#' Fit linear catch curve and get results
#'
#' This function fits a linear catch curve to fish age composition data.
#' The age at full recruitment is determined using on the age at
#' peak frequency of fish. If RecAssump=0, RecAge is set to the age at
#' peak frequency, but if RecAssump=1, RecAge is set to one year above the age at
#' peak frequency. The oldest age that should be considered for analysis with
#' the linear catch curve analysis, is calculated given a specified minimum frequency
#' (must be >0). The catch curve is fitted to the natural logarithms of the frequency
#' data using the lm function
#'
#' @param RecAssump 0=age at peak frequency, 1=age at peak frequency + 1
#' @param SpecRecAge specified age at full recruitment. Set to NA if RecAssump = 0 or 1.
#' @param MinFreq minimum frequency of fish for including data for old fish
#' @param Ages ages for analysis
#' @param ObsAgeFreq observed age frequency
#'
#' @return linear catch curve y intercept and associated standard error (yint, yint_se),
#' total mortality and associated standard error (ZMort, ZMort_se), age at full recruitment (RecAge),
#' last age included in analysis (LastAgeForLinearCC), age classes included in analysis (AgeClassesForLinearCC),
#' natural logarithms of age frequencies in analysis (lnObsAgeFreqForLinearCC), natural logarithms of estimated
#' frequencies and associated 95 percent confidence limits (EstlnFreq, EstlnFreq_Zlow, EstlnFreq_Zup), ages
#' for observed frequency data (Ages), estimated frequencies at age, in normal space, with associated 95 percent
#' confidence limits (EstFreq, EstFreq_Zlow, EstFreq_Zup)
#'
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
#' SampleSize = 1000 # required number of fish for age sample
#' Res=SimAgeFreqData(SampleSize, MinAge, MaxAge, SelA50, SelA95, NatMort, FMort)
#' ObsAgeFreq = unlist(as.vector(Res$CatchSample))
#' res=GetLinearCatchCurveResults(RecAssump=0, SpecRecAge=NA, MinFreq=1, Ages, ObsAgeFreq)
#' @export
GetLinearCatchCurveResults <- function(RecAssump, SpecRecAge, MinFreq, Ages, ObsAgeFreq) {

  # get recruitment age, given recruitment assumption
  if (RecAssump==2) {
    RecAge = SpecRecAge
  } else {
    RecAge=CalcRecruitmentAge(RecAssump, Ages, ObsAgeFreq)
  }

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
#' estimator, with associated 95 percent confidence limits. Additional variables are returned for plotting
#'
#' @param RecAssump 0=age at peak frequency, 1=age at peak frequency + 1, 2=specified age at full recruitment
#' @param SpecRecAge specified at at full recruitment, set to NA when RecAssump is set to 0 or 1
#' @param MinAge minimum age
#' @param MaxAge maximum age
#' @param ObsAgeFreq observed age frequency
#'
#' @return total mortality, Z (ZMort), standard error and approximate lower and upper 95 percent confidence
#' limits for Z, (ZMort_se, ZMort_low, ZMort_up), median estimate and lower and upper 95 percent confidence
#' limits for Z from resampling (ZMort_resamp, ZMort_lowresamp, ZMort_upresamp), variables for plotting,
#' including expected frequencies at age with associated confidence limits, calculated from resampling
#' (EstFreq, EstFreq_Zup, EstFreq_Zlow), recruitment age (RecAge),  maximum age in sample (MaxAgeInSample),
#' relative ages used in analysis (CRAges), sample size for relative ages (n) and additional variable
#' calculated in analysis (CR_T)
#'
#' @examples
#' set.seed(123)
#' MinAge = 1
#' MaxAge = 40
#' NatMort <- exp(1.46 - (1.01 * (log(MaxAge)))) # i.e. Hoenig's (1983) eqn for fish
#' FMort = 0.1
#' ZMort = FMort + NatMort
#' SelA50 = 6
#' SelA95 = 8
#' SampleSize = 1000 # required number of fish for age sample
#' Res=SimAgeFreqData(SampleSize, MinAge, MaxAge, SelA50, SelA95, NatMort, FMort)
#' ObsAgeFreq = unlist(as.vector(Res$CatchSample))
#' res=GetChapmanRobsonMortalityResults(RecAssump=0, SpecRecAge=NA, MinAge, MaxAge, ObsAgeFreq)
#' @export
GetChapmanRobsonMortalityResults <- function(RecAssump, SpecRecAge, MinAge, MaxAge, ObsAgeFreq)
{
  Ages = MinAge:MaxAge
  MaxAgeInSample = max(Ages)
  FishAges <- rep(Ages, ObsAgeFreq)

  if (RecAssump==2) {
    RecAge = SpecRecAge
  } else {
    RecAge = CalcRecruitmentAge(RecAssump, Ages, ObsAgeFreq)
  }
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
    FreqFirstAge = 0
    while (FreqFirstAge == 0) { # don't allow any random samples where the frequency for the first recruited age class is zero
      RandCRAges <- sample(CRAges, replace = T)
      RandCRAgeFreq = as.vector(table(factor(RandCRAges, levels=seq(0,MaxAgeInSample-RecAge,1))))
      FreqFirstAge = RandCRAgeFreq[1]
    }
    RandAveAge <- mean(RandCRAges)
    RandCR_T <- n * RandAveAge
    RandZMort[i] <- log(1 + RandAveAge - (n^-1)) - log(RandAveAge) -
      (((n - 1) * (n - 2)) * (n * (RandCR_T + 1) * (n + RandCR_T - 1))^-1)

    # using C&R mortality estimate, calculate an expected frequency by
    # using optimisation to find value of y intercept that provides best
    # fit of expected curve to data
    params = RandCRAgeFreq[1] # set initial value
    tempObjFunc <- function(params) {
      y_intercept = params
      ExpVals = y_intercept * exp(-RandZMort[i] * seq(0, ArrSize - 1,1))
      SumSq = sum((RandCRAgeFreq - ExpVals)^2)
      return(SumSq)
    }
    nlmb = nlminb(params, tempObjFunc)
    RandEstFreq[i,] = nlmb$par * exp(-RandZMort[i] * seq(0, ArrSize - 1,1))
  }

  EstFreq = apply(RandEstFreq,2,median)
  EstFreq_Zlow = apply(RandEstFreq,2,quantile,probs=0.025)
  EstFreq_Zup = apply(RandEstFreq,2,quantile,probs=0.975)
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

#' Return negative log-likelihood for a catch curve with age-based, logistic selectivity
#'
#' This function returns the negative log-likelihood, given age frequency data, a value total mortality,
#' and age-based selectivity parameters, for a catch curve with age-based, logistic selectivity.
#' Function requires an estimate of natural mortality (NatMort), a value for maximum age (MaxAge)
#' and age frequency data (stored in memory in R)
#'
#' @keywords internal
#'
#' @param ln_params model parameters log(c(FMort, SelA50, SelA95)
#'
#' @return negative log-likelihood (NLL)
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
#' with the parameters and data, using nlminb. It provides various statistical outputs including
#' convergence statistics, parameter estimates and associated 95 percent confidence limits and associated
#' variance-covariance matrix, calculated using the MASS package
#'
#' @param ln_params initial values of model parameters log(c(FMort, SelA50, SelA95)
#' @param NatMort natural mortality
#' @param Ages ages in observed data
#' @param ObsAgeFreq observed age frequency data
#'
#' @return negative log-likelihood (nll), nlminb convergence diagnostic (convergence)
#' sample size (SampleSize), parameter estimates with lower and upper 95 percent
#' confidence limits (ParamEst), point estimates for parameters estimated in log space (Estlnparams),
#' variance-covariance matrix for estimated parameters (vcov.Params), standard errors for estimated
#' parameters (lnEstFMort_se, lnEstSelA50_se, lnEstFSelA95_se), selectivity at age (SelAtAge), fishing
#' mortality at age (FAtAge), estimated frequencies at age with associated 95 percent confidence limits
#' (EstFreq, EstFreq_Zlow, EstFreq_Zup), random values of parameters in log space from parametric resampling,
#' using rmultinom function (lnparams.sims), and associated median and lower 2.5 and upper 97.5 percentiles
#' of estimates for frequency at age (EstFreq.sim), selectivity parameters in normal space (SelA50.sim, SelA95.sim),
#' fishing mortality (FMort.sim) and total mortality (EstZMort.sim)
#'
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
#' SampleSize = 1000 # required number of fish for age sample
#' Res=SimAgeFreqData(SampleSize, MinAge, MaxAge, SelA50, SelA95, NatMort, FMort)
#' ObsAgeFreq = unlist(as.vector(Res$CatchSample))
#' Init_FMort = 0.2
#' Init_SelA50 = 5
#' Init_SelA95 = 7
#' ln_params = log(c(Init_FMort, Init_SelA50, Init_SelA95))
#' res=GetLogisticCatchCurveResults(ln_params, NatMort, Ages, ObsAgeFreq)
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
  EstFMort = c(exp(nlmb$par[1]), exp(nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
  EstSelA50 = c(exp(nlmb$par[2]), exp(nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
  EstSelA95 = c(exp(nlmb$par[3]), exp(nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
  SelAtAge = rep(0, length(Ages))
  SelAtAge = 1/(1 + exp(-log(19) * (Ages - EstSelA50[1])/(EstSelA95[1] - EstSelA50[1])))
  FAtAge = SelAtAge * EstFMort[1]
  ZAtAge = NatMort + FAtAge
  ParamEst = t(data.frame(FMort = round(EstFMort, 3), SelA50 = round(EstSelA50,3),
                          SelA95 = round(EstSelA95, 3)))
  colnames(ParamEst) = c("Estimate", "lw_95%CL", "up_95%CL")
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
                 vcov.Params = vcov.Params,
                 lnEstFMort_se = ses[1],
                 lnEstSelA50_se = ses[2],
                 lnEstFSelA95_se = ses[3],
                 SelAtAge = SelAtAge,
                 FAtAge = FAtAge,
                 EstFreq = EstFreq,
                 EstFreq_Zlow = EstFreq_Zlow,
                 EstFreq_Zup = EstFreq_Zup,
                 lnparams.sims = sims,
                 EstFreq.sim = EstFreq.sim,
                 SelA50.sim = SelA50.sim,
                 SelA95.sim = SelA95.sim,
                 FMort.sim = FMort.sim,
                 EstZMort.sim = EstZMort.sim)
  return(results)
}




#' Plot age based catch curve results in normal space
#'
#' This function produces plots of outputs of age-based catch curve analyses in normal space
#'
#' @param RecAssump 0=age at peak frequency, 1=age at peak frequency + 1
#' @param SpecRecAge specified at at full recruitment, set to NA when RecAssump is set to 0 or 1
#' @param MinFreq minimum frequency of fish for including data for old fish
#' @param MinAge minimum age
#' @param MaxAge maximum age
#' @param NatMort natural mortality
#' @param ObsAgeFreq observed age frequency data
#' @param CatchCurveModel 1=C&R, 2=Linear, 3=logist. sel.
#' @param MainLabel plot label
#' @param xaxis_lab y axis label
#' @param yaxis_lab x axis label
#' @param xmax maximum x axis value
#' @param xint x axis interval
#' @param ymax maximum y axis value
#' @param yint y axis interval
#' @param PlotCLs logical (TRUE=plot 95 percent confidence limits for lines)
#'
#' @return plot of expected vs observed proportions at age
#'
#' @examples
#' # Simulate data
#' set.seed(123)
#' MinAge = 1
#' MaxAge = 40
#' Ages = MinAge:MaxAge
#' NatMort <- exp(1.46 - (1.01 * (log(MaxAge)))) # i.e. Hoenig's (1983) eqn for fish
#' FMort = 0.1
#' ZMort = FMort + NatMort
#' SelA50 = 6
#' SelA95 = 8
#' SampleSize = 1000 # required number of fish for age sample
#' Res=SimAgeFreqData(SampleSize, MinAge, MaxAge, SelA50, SelA95, NatMort, FMort)
#' ObsAgeFreq = unlist(as.vector(Res$CatchSample))
#' # Specify catch curve model and required inputs for that model
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
#' Init_FMort = 0.2
#' Init_SelA50 = 5
#' Init_SelA95 = 7
#' ln_params = log(c(Init_FMort, Init_SelA50, SelA95))
#' PlotAgeBasedCatchCurveResults_NormalSpace(RecAssump, SpecRecAge=NA, MinFreq, MinAge, MaxAge, NatMort,
#'                                           ObsAgeFreq, CatchCurveModel, MainLabel=NA,
#'                                           xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA,
#'                                          ymax=NA, yint=NA, PlotCLs=T)
#' @export
PlotAgeBasedCatchCurveResults_NormalSpace <- function(RecAssump, SpecRecAge, MinFreq, MinAge, MaxAge, NatMort,
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
    Res = GetChapmanRobsonMortalityResults(RecAssump, SpecRecAge, MinAge, MaxAge, ObsAgeFreq)
    if (is.na(MainLabel)) MainLabel = "Chapman & Robson"
  }
  # Linear
  if (CatchCurveModel == 2) {
    Res = GetLinearCatchCurveResults(RecAssump, SpecRecAge, MinFreq, Ages, ObsAgeFreq)
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
    xx=length(which(Res$EstFreq_Zup>0))+x # position of last age
    xxx=length(which(Res$EstFreq_Zlow>0))+x # position of last age
    xxxx=length(which(Res$EstFreq>0))+x # position of last age
  }
  # Linear
  if (CatchCurveModel == 2) {
    Z_value = round(Res$ZMort[1],digits=3)
    x=which(Ages==Res$RecAge) # RecAge position
    xx=which(Ages==Res$LastAgeForLinearCC) # position of last age
    xxx=which(Ages==Res$LastAgeForLinearCC) #  # position of last age
    xxxx=which(Ages==Res$LastAgeForLinearCC) #  # position of last age
  } # logistic
  if (CatchCurveModel == 3) {
    Z_value = round(Res$ParamEst[1,1] + NatMort,digits=3)
    x=which(Ages==min(Ages)) # RecAge position
    xx=which(Ages==max(Ages)) # position of last age
    xxx=which(Ages==max(Ages)) #  # position of last age
    xxxx=which(Ages==max(Ages)) #  # position of last age
  }

  j = seq(x,xx,1) # up
  jj = seq(x,xxx,1) # low
  jjj = seq(x,xxxx,1) # est

  # plot catch curve over age frequency data, in normal space
  plot(Ages, ObsAgeFreq, "p", main=MainLabel, cex.main=1.2, pch=16, cex=0.8, xaxt = "n", yaxt = "n",
       xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax)) # observed data (normal space)
  axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
  axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
  axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
  axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
  if (PlotCLs == TRUE) {
    sm1 = spline(Ages[j], Res$EstFreq_Zlow[1:length(jj)], n=100, method="natural")
    sm2 = spline(Ages[j], Res$EstFreq_Zup[1:length(j)], n=100, method="natural")
    x = c(sm1$x, rev(sm2$x)) # using shading for 95% CLs
    y = c(sm1$y, rev(sm2$y))
    polygon(x,y, col="pink",border=NA)
  }
  # lines(Ages[jjj], Res$EstFreq[1:length(jjj)], col="black")
  points(Ages, ObsAgeFreq, pch=16, cex=0.8)
  points(Ages[j], Res$EstFreq[1:length(j)], col="red", pch=1, cex=0.8)

  if (PlotCLs == FALSE) { # if not plotting confidence intervals, can include last age
    xx=which(Ages==max(Ages)) # position of last age
    points(Ages[x:xx], Res$EstFreq[1:length(x:xx)], col="red", pch=1, cex=0.8)
  }

  legend("topright", legend=bquote(paste("Z = ", .(Z_value), " ",y^-1)), y.intersp = 1.5, inset=c(0.13,0),
         lty=1, cex = 1, bty="n",seg.len = 0)
  # legend("topleft", legend=c("Observed","Estimated"), y.intersp = 1.0, inset=c(0.13,0),
  #        lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16,1), col=c("black","red"))
}

#' Plot age based catch curve results in log space
#'
#' This function produces plots of outputs of age-based catch curve analyses in log space
#'
#' @param RecAssump 0=age at peak frequency, 1=age at peak frequency + 1
#' @param SpecRecAge specified at at full recruitment, set to NA when RecAssump is set to 0 or 1
#' @param MinFreq minimum frequency of fish for including data for old fish
#' @param MinAge minimum age
#' @param MaxAge maximum age
#' @param NatMort natural mortality
#' @param ObsAgeFreq observed age frequency data
#' @param CatchCurveModel 1=C&R, 2=Linear, 3=logist. sel.
#' @param MainLabel plot label
#' @param xaxis_lab y axis label
#' @param yaxis_lab x axis label
#' @param xmax maximum x axis value
#' @param xint x axis interval
#' @param ymax maximum y axis value
#' @param yint y axis interval
#' @param PlotCLs logical (TRUE=plot 95 percent confidence limits for lines)
#'
#' @return plot of expected vs observed proportions at age
#'
#' @examples
#' # Simulate data
#' set.seed(123)
#' MinAge = 1
#' MaxAge = 40
#' Ages = MinAge:MaxAge
#' NatMort <- exp(1.46 - (1.01 * (log(MaxAge)))) # i.e. Hoenig's (1983) eqn for fish
#' FMort = 0.1
#' ZMort = FMort + NatMort
#' SelA50 = 6
#' SelA95 = 8
#' SampleSize = 1000 # required number of fish for age sample
#' Res=SimAgeFreqData(SampleSize, MinAge, MaxAge, SelA50, SelA95, NatMort, FMort)
#' ObsAgeFreq = unlist(as.vector(Res$CatchSample))
#' # Specify catch curve model and required inputs for that model
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
#' Init_FMort = 0.2
#' Init_SelA50 = 5
#' Init_SelA95 = 7
#' ln_params = log(c(Init_FMort, Init_SelA50, SelA95))
#' PlotAgeBasedCatchCurveResults_LogSpace(RecAssump, SpecRecAge=NA, MinFreq, MinAge, MaxAge, NatMort,
#'                                        ObsAgeFreq, CatchCurveModel, MainLabel=NA,
#'                                        xaxis_lab=NA, yaxis_lab=NA, ymin=NA, xmax=NA, xint=NA,
#'                                        ymax=NA, yint=NA, PlotCLs=T)
#' @export
PlotAgeBasedCatchCurveResults_LogSpace <- function(RecAssump, SpecRecAge, MinFreq, MinAge, MaxAge, NatMort,
                                                   ObsAgeFreq, CatchCurveModel, MainLabel,
                                                   xaxis_lab, yaxis_lab, ymin, xmax, xint,
                                                   ymax, yint, PlotCLs) {

  if (is.na(xaxis_lab)) xaxis_lab = "Age class"
  if (is.na(yaxis_lab)) yaxis_lab = "ln Frequency"
  ylims = Get_yaxis_scale(log(ObsAgeFreq[which(ObsAgeFreq>0)]))
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
    Res = GetChapmanRobsonMortalityResults(RecAssump, SpecRecAge, MinAge, MaxAge, ObsAgeFreq)
    if (is.na(MainLabel)) MainLabel = "Chapman & Robson"
  }
  # Linear
  if (CatchCurveModel == 2) {
    Res = GetLinearCatchCurveResults(RecAssump, SpecRecAge, MinFreq, Ages, ObsAgeFreq)
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
    Z_value = round(Res$ParamEst[1,1] + NatMort,digits=3)
    # x=which(Ages==min(which(log(Res$EstFreq_Zlow)>0))) # RecAge position
    x=min(which(log(Res$EstFreq_Zlow)>0)) # RecAge position
    xx=length(which(log(Res$EstFreq_Zlow) > -1)) # position of last age
    j = seq(x,xx,1) # up
  }

  plot(Ages, log(ObsAgeFreq), "p", main=MainLabel, pch=16, cex=0.8, cex.main=1.2, xaxt = "n", yaxt = "n",
       xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(ymin,ymax)) # observed data (normal space)
  axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
  axis(2, at = seq(ymin, ymax, yint), line = 0.2, labels = F)
  axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
  axis(2, at = seq(ymin, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
  # Chap-Rob
  if (CatchCurveModel == 1) {
    if (PlotCLs == TRUE) {
      sm1 = spline(Ages[j], log(Res$EstFreq_Zup[k]), n=100, method="natural")
      sm2 = spline(Ages[jj], log(Res$EstFreq_Zlow[kk]), n=100, method="natural")
      x = c(sm1$x, rev(sm2$x)) # using shading for 95% CLs
      y = c(sm1$y, rev(sm2$y))
      polygon(x,y, col="pink",border=NA)
    }
    # lines(Ages[jjj], log(Res$EstFreq[kkk]), col="red")
    points(Ages[jjj], log(Res$EstFreq[kkk]), pch=1, col="red", cex=0.6)
  }
  if (CatchCurveModel == 2) {
    if (PlotCLs == TRUE) {

      sm1 = spline(Ages[j], log(Res$EstFreq_Zup[1:length(j)]), n=100, method="natural")
      sm2 = spline(Ages[j], log(Res$EstFreq_Zlow[1:length(jj)]), n=100, method="natural")
      x = c(sm1$x, rev(sm2$x)) # using shading for 95% CLs
      y = c(sm1$y, rev(sm2$y))
      polygon(x,y, col="pink",border=NA)
    }
    # lines(Ages[jjj], log(Res$EstFreq[1:length(jjj)]), col="red")
    points(Ages[jjj], log(Res$EstFreq[1:length(jjj)]), pch=1, col="red", cex=0.6)
  }
  if (CatchCurveModel == 3) {
    if (PlotCLs == TRUE) {

      sm1 = spline(Ages[j], log(Res$EstFreq_Zup[j]), n=100, method="natural")
      sm2 = spline(Ages[j], log(Res$EstFreq_Zlow[j]), n=100, method="natural")
      x = c(sm1$x, rev(sm2$x)) # using shading for 95% CLs
      y = c(sm1$y, rev(sm2$y))
      polygon(x,y, col="pink",border=NA)
    }
    #lines(Ages[j], log(Res$EstFreq[j]), col="red")
    points(Ages, log(Res$EstFreq), pch=1, col="red", cex=0.6)
  }

  points(Ages, log(ObsAgeFreq), pch=16, cex=0.8)
  legend("topright", legend=bquote(paste("Z = ", .(Z_value), " ",y^-1)), y.intersp = 1.5, inset=c(0.13,0),
         lty=1, cex = 1, bty="n",seg.len = 0)
  # legend("topleft", legend=c("Observed","Estimated"), y.intersp = 1.0, inset=c(0.13,0),
  #        lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16,1), col=c("black","red"))
}


#*********************
# Per recruit analyses
#*********************

#' Get outputs from age-based per recruit analysis for a specified value of fishing mortality
#'
#' This function provides outputs associated with per recruit analysis, and an
#' extended form of this analysis with a Beverton-Holt stock recruitment relationship to account
#' for potential impacts of fishing on recruitment
#'
#' @param MaxModelAge maximum age considered in model
#' @param TimeStep model timestep (e.g. 1 = annual, 1/12 = monthly)
#' @param Linf asymptotic length (von Bertalanffy growth parameter, set to NA if inputting lengths at ages directly from any growth model)
#' @param vbK growth coefficient (von Bertalanffy growth parameter, set to NA if inputting lengths at ages directly from any growth model)
#' @param tzero hypothetical age at zero length (von Bertalanffy growth parameter, set to NA if inputting lengths at ages directly from any growth model)
#' @param EstLenAtAge vector or estimated lengths at ages from any growth model, set to NA if von Bertalanffy growth parameters specified
#' @param lenwt_a weight-length parameter (power relationship)
#' @param ln_lenwt_a weight-length parameter (log-log relationship)
#' @param lenwt_b weight-length parameter (power or log-log relationship)
#' @param WLrel_Type 1=power, 2=log-log relationship (set to NA if inputting weights at ages directly)
#' @param EstWtAtAge vector of weights at ages (set to NA if weight-length growth parameters specified)
#' @param ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' @param ReprodPattern reproductive pattern, 1=gonochoristic, 2=protogynous (female to male sex change) hermaphroditism,
#' 3=protandrous hermaphroditism (male to female sex change)
#' @param InitRatioFem proportion of fish that are females at hatching
#' @param FinalSex_Pmax logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_A50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_A95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param mat_A50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_A95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtAge vector of proportion mature at age (set to NA if using age at maturity parameters)
#' @param sel_A50 logistic parameter for gear selectivity curve
#' @param sel_A95 logistic parameter for gear selectivity curve
#' @param EstSelAtAge vector for gear selectivity at age (set to NA if using age at gear selectivity parameters)
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_A50 logistic parameter for fish retention curve
#' @param ret_A95 logistic parameter for fish retention curve
#' @param EstRetenAtAge vector for fish retention at age (set to NA if using age at fish retention parameters)
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param FMort fishing mortality
#'
#' @return Ages considered in analysis (Ages), yield per recruit for combined sexes (YPR),
#' female, male and combined sex spawning potential ratio (Fem_SPR, Mal_SPR, CombSex_SPR), unfished female
#' survival at age (UnfishFemSurvAtAge, UnfishMalSurvAtAge), unfished female and male biomass mature biomass at age
#' (UnfishFemBiomAtAge, UnfishMalBiomAtAge), fished female survival at age (FishedFemSurvAtAge, FishedMalSurvAtAge),
#' fished female and male biomass mature biomass at age (FishedFemBiomAtAge, FishedMalBiomAtAge), female and male
#' retained catch at age in numbers (FemCatchAtAgeNum, MalCatchAtAgeNum), female and male
#' release catch at age in numbers (FemRelCatchAtAgeNum, MalRelCatchAtAgeNum), female and male retained catch at age in biomass
#' (FemCatchAtAge, MalCatchAtAge), derived Beverton-Holt stock recruitment parameters and associated equilibrium
#' recruitment derived from those parameters (BH_SRRa, BH_SRRb, BH_Equil_Rec), equilibrium recruitment for
#' either Beverton-Holt or Ricker relationship (Equil_Rec), equilibrium catch (Equil_Catch), equilibrium
#' female and male and spawning biomass (Equil_FemSpBiom, Equil_MalSpBiom), equilibrium relative female,
#' male and combined sex spawning biomass (Equilmod_SPR, Equilmod_MalRelBiom, Equilmod_CombSexRelBiom),
#' female and male length at age (FemLenAtAge, MalLenAtAge), female and male weight at age
#' (FemWtAtAge, MalWtAtAge), female, male and combined sex proportion mature at age (FemPropMatAtAge, MalPropMatAtAge, PropFemAtAge),
#' female and male gear selectivity at age (FemGearSelAtAge, MalGearSelAtAge), female and male retention at age probabilities,
#' (FemRetProbAtAge, MalRetProbAtAge), selectivity of female and male fish landings (FemSelLandAtAge, MalSelLandAtAge),
#' female and male selectivity of discards (FemSelDiscAtAge, MalSelDiscAtAge), female and male fishing mortality
#' associated with discarding (FemDiscFAtAge, MalDiscFAtAge), female and male fishing mortality associated with
#' fish landings (FemLandFAtAge, MalLandFAtAge), female and male total fishing mortality at age (FemFAtAge, MalFAtAge),
#' female and male total mortality at age (FemZAtAge, MalZAtAge)
#'
#' @examples
#' # Example 1. Non-hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 20 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be approprriate to use a smaller time step)
#' Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.5, 0.5) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstFemLenAtAge=NA, EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- NA # for log-log relationship
#' lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA, EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' EstSelAtAge <- data.frame(EstFemSelAtAge=NA, EstMalSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
#' EstRetenAtAge <- data.frame(EstFemRetenAtAge=rep(1,MaxModelAge+1), EstMalRetenAtAge=rep(1,MaxModelAge+1)) # fish retention at age (from age 0), inputted as values in data frame
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 0.2 # natural mortality  (year-1)
#' FMort <- 0.4 # estimate of fishing mortality, e.g. from catch curve analysis
#' Res=CalcYPRAndSPRForFMort_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                           lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
#'                           ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
#'                           mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
#'                           EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, FMort)
#' # Example 2: hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 100 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
#' Linf <- c(1000, 1000) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.1, 0.1) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstMalLenAtAge=NA, EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- -11.0 # for log-log relationship
#' lenwt_b <- 3.0 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 2 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA, EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 2 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 1 # Ratio of females to males at age zero
#' FinalSex_Pmax <- 1 # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_A50 <- 35 # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- 60 # Logistic sex change relationship parameters (95% of max probability)
#' mat_A50 <- c(20, 20) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(30, 30) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_A50 <- c(15, 15) # females, males - Logistic age selectivity relationship parameters
#' sel_A95 <- c(25, 25) # females, males - Logistic age selectivity relationship parameters
#' EstSelAtAge <- data.frame(EstFemSelAtAge=NA, EstMalSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
#' EstRetenAtAge <- data.frame(EstFemRetenAtAge=rep(1,MaxModelAge+1), EstMalRetenAtAge=rep(1,MaxModelAge+1)) # fish retention at age (from age 0), inputted as values in data frame
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 0.07 # natural mortality  (year-1)
#' FMort <- 0.07 # estimate of fishing mortality, e.g. from catch curve analysis
#' Res=CalcYPRAndSPRForFMort_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                           lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
#'                           ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
#'                           mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
#'                           EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, FMort)
#' @export
CalcYPRAndSPRForFMort_AB <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                  lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                                  ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                                  mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
                                  EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, FMort) {

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
  if (is.na(EstSelAtAge[1,1])) {
    FemGearSelAtAge <- 1/(1+exp(-log(19) * (Ages - sel_A50[1]) / (sel_A95[1] - sel_A50[1])))
    MalGearSelAtAge <- 1/(1+exp(-log(19) * (Ages - sel_A50[2]) / (sel_A95[2] - sel_A50[2])))
  }
  if (!is.na(EstSelAtAge[1,1])) {
    FemGearSelAtAge <- EstSelAtAge[,1]
    MalGearSelAtAge <- EstSelAtAge[,2]
  }

  # Calculate fish retention at age
  if (is.na(EstRetenAtAge[1,1])) {
    FemRetProbAtAge <- ret_Pmax[1] / (1+exp(-log(19) * (Ages - ret_A50[1]) / (ret_A95[1] - ret_A50[1])))
    MalRetProbAtAge <- ret_Pmax[2] / (1+exp(-log(19) * (Ages - ret_A50[2]) / (ret_A95[2] - ret_A50[2])))
  }
  if (!is.na(EstRetenAtAge[1,1])) {
    FemRetProbAtAge <- EstRetenAtAge[,1]
    MalRetProbAtAge <- EstRetenAtAge[,2]
  }

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
  if (DiscMort >= 0.001) {
    FemFAtAge <- FMort * (FemSelLandAtAge + (DiscMort * FemSelDiscAtAge))
    MalFAtAge <- FMort * (MalSelLandAtAge + (DiscMort * MalSelDiscAtAge))
  } else {
    FemFAtAge <- FMort * FemSelLandAtAge
    MalFAtAge <- FMort * MalSelLandAtAge
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
        (1 - exp(-NatMort))
      UnfishMalSurvAtAge[j] <- (UnfishMalSurvAtAge[j-1] * exp(-NatMort * TimeStep)) /
        (1 - exp(-NatMort))
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
  UnfishFemBiomAtAge <- UnfishFemSurvAtAge * FemPropMatAtAge * (((FemWtAtAge * 1000) ^ ReprodScale) / 1000)
  UnfishMalBiomAtAge <- UnfishMalSurvAtAge * MalPropMatAtAge * (((MalWtAtAge * 1000) ^ ReprodScale) / 1000)

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
        (1 - exp(-(FMort+NatMort)))
      FishedMalSurvAtAge[j] <- FishedMalSurvAtAge[j-1] * exp(-MalZAtAge[j-1] * TimeStep) /
        (1 - exp(-(FMort+NatMort)))
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
  FemPropMatAtAge * FemWtAtAge
  FishedFemBiomAtAge <- FishedFemSurvAtAge * FemPropMatAtAge * (((FemWtAtAge * 1000) ^ ReprodScale) / 1000)
  FishedMalBiomAtAge <- FishedMalSurvAtAge * MalPropMatAtAge * (((MalWtAtAge * 1000) ^ ReprodScale) / 1000)

  # calculate female and male catch at age (in numbers) - Baranov catch equation
  FemCatchAtAgeNum <- FishedFemSurvAtAge * (FemLandFAtAge/FemZAtAge) *
    (1 - exp(-(FemZAtAge * TimeStep)))
  MalCatchAtAgeNum <- FishedMalSurvAtAge * (MalLandFAtAge/MalZAtAge) *
    (1 - exp(-(MalZAtAge * TimeStep)))

  # calculate female and male discarded catch at age (in numbers) - Baranov catch equation
  # (including fish that either survive or die after release)
  FemRelCatchAtAgeNum <- FishedFemSurvAtAge * ((FMort * FemSelDiscAtAge)/FemZAtAge) *
    (1 - exp(-(FemZAtAge * TimeStep)))
  MalRelCatchAtAgeNum <- FishedMalSurvAtAge * ((FMort * MalSelDiscAtAge)/MalZAtAge) *
    (1 - exp(-(MalZAtAge * TimeStep)))

  # calculate female and male catch at age (in biomass)
  FemCatchAtAge <- FemCatchAtAgeNum * FemWtAtAge
  MalCatchAtAge <- MalCatchAtAgeNum * MalWtAtAge

  # calculate yield per recruit in kg
  YPR <- max(0,sum(FemCatchAtAge) + sum(MalCatchAtAge))

  # calculate female and male spawning biomass per recruit in kg
  FishFemSpawnBiom <- sum(FishedFemBiomAtAge)
  FishMalSpawnBiom <- sum(FishedMalBiomAtAge)
  FishCombSexSpawnBiom <- FishFemSpawnBiom + FishMalSpawnBiom

  # calculate spawning potential ratio (SPR)
  Fem_SPR <- max(0,FishFemSpawnBiom / UnfishFemSpawnBiom)
  Mal_SPR <- max(0,FishMalSpawnBiom / UnfishMalSpawnBiom)
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
  Equil_Catch <- max(0,Equil_Rec * YPR)

  # calculate equilibrium female spawning biomass
  Equil_FemSpBiom <- Equil_Rec * FishFemSpawnBiom
  Equil_MalSpBiom <- Equil_Rec * FishMalSpawnBiom

  # calculate equilibrium model SPR
  Equilmod_FemRelBiom <- max(0,Equil_FemSpBiom / UnfishFemSpawnBiom)
  Equilmod_MalRelBiom <- max(0,Equil_MalSpBiom / UnfishMalSpawnBiom)
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
                 FemRelCatchAtAgeNum = FemRelCatchAtAgeNum,
                 MalRelCatchAtAgeNum = MalRelCatchAtAgeNum,
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


#' Calculates mean size at age for each sex for length-based per recruit analysis
#'
#' This function calculates mean size at age for each sex applying the von Bertalanffy growth function or Schnute function,
#' for length-based per recruit analysis,
#'
#' @keywords internal
#'
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
#' @param MaxModelAge maximum age considered by model
#' @param Ages ages considered by model
#' @param GrowthParams growth parameters of either von Bertalanffy or Schnute model
#' @param RefnceAges reference ages for Schnute model, set to NA if using von Bertalanffy model
#'
#' @return MeanSizeAtAge
CalcMeanSizeAtAgeForPerRecruitAnalysis <- function(GrowthCurveType, MaxModelAge, Ages, GrowthParams, RefnceAges) {

  MeanSizeAtAge <- data.frame(matrix(nrow = 2, ncol = length(Ages)))
  colnames(MeanSizeAtAge) <- round(Ages,2)
  MeanSizeAtAge <- as.matrix(MeanSizeAtAge)

  # calculate length at age - von Bertalanffy growth curve
  if (GrowthCurveType == 1) {
    FemLenAtAge <- GrowthParams[1,1] * (1 - exp(-GrowthParams[1,2] * (Ages - GrowthParams[1,3])))
    MalLenAtAge <- GrowthParams[2,1] * (1 - exp(-GrowthParams[2,2] * (Ages - GrowthParams[2,3])))
  }
  # Schnute growth curve
  if (GrowthCurveType == 2) {
    t1=RefnceAges[1]
    t2=RefnceAges[2]
    for (i in 1:2) {
      y1=GrowthParams[i,1]
      y2=GrowthParams[i,2]
      a=GrowthParams[i,3]
      b=GrowthParams[i,4]
      for (t in 1:MaxModelAge) {
        if (i==1) {
          FemLenAtAge[t] = SchnuteGrowthfunction(t, t1, t2, y1, y2, a, b)
        }
        if (i == 2) {
          MalLenAtAge[t] = SchnuteGrowthfunction(t, t1, t2, y1, y2, a, b)
        }
      }
    }
  }

  # set negative lengths to zero
  FemLenAtAge[which(FemLenAtAge<0)] = 0
  MalLenAtAge[which(MalLenAtAge<0)] = 0
  MeanSizeAtAge[1,] <- FemLenAtAge
  MeanSizeAtAge[2,] <- MalLenAtAge

  return(MeanSizeAtAge)

}


#' Calculates weight at length for each sex for length-based per recruit analysis
#'
#' This function calculates mean size at age for each sex, based on power curve,
#' log-log relationship, or allows for these to be used as vectors inputted by user,
#' for length-based per recruit analysis
#'
#' @keywords internal
#'
#' @param lenwt_a weight-length parameter
#' @param ln_lenwt_a weight-length parameter
#' @param lenwt_b weight-length parameter
#' @param WLrel_Type 1=power, 2=log-log
#' @param EstWtAtLen user-specified weights at lengths
#' @param midpt length class mid-points
#'
#' @return FemWtAtLen
#' @return MalWtAtLen
CalcWeightAtLengthForPerRecruitAnalysis <- function(lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen, midpt) {


  # calculate weight at length (i.e. for corresponding relative age)
  if (is.na(EstWtAtLen[1,1])) {
    if (WLrel_Type == 1) { # power relationship
      FemWtAtLen <- (lenwt_a * midpt ^ lenwt_b) / 1000 # weight at length, kg
      MalWtAtLen <- (lenwt_a * midpt ^ lenwt_b) / 1000 # weight at length, kg
    }
    if (WLrel_Type == 2) { # log-log relationship
      FemWtAtLen <- exp(ln_lenwt_a + lenwt_b * log(midpt)) / 1000 # weight at length, kg
      MalWtAtLen <- exp(ln_lenwt_a + lenwt_b * log(midpt)) / 1000 # weight at length, kg
    }
  }
  if (!is.na(EstWtAtLen[1,1])) {
    FemWtAtLen <- EstWtAtLen[,1]
    MalWtAtLen <- EstWtAtLen[,2]
  }

  results = list(FemWtAtLen=FemWtAtLen,
                 MalWtAtLen=MalWtAtLen)

}

#' Calculates inputs required to calculate length transition matrices, for length-based per recruit analysis
#'
#' This function calculates inputs required to calculate length transition matrices, for length-based per
#' recruit analysis, including the mean size to which a fish grows, for a  model timestep, given its initial length and specified
#' growth curve (MeanEndingLength), and the amount by which it grows (TimestepGrowthSizeInc)
#'
#' @keywords internal
#'
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
#' @param Ages ages considered by model
#' @param GrowthParams growth parameters of either von Bertalanffy or Schnute model
#' @param RefnceAges reference ages for Schnute model, set to NA if using von Bertalanffy model
#' @param midpt length class mid-points
#'
#' @return MeanEndingLength
#' @return TimestepGrowthSizeInc
GetLTMInputsForPerRecruitAnalysis <- function(GrowthCurveType, TimeStep, MaxModelAge, GrowthParams, RefnceAges, midpt) {

  nLenCl = length(midpt)
  MeanEndingLength <- data.frame(matrix(nrow = 2, ncol = nLenCl))
  colnames(MeanEndingLength) <- midpt
  MeanEndingLength <- as.matrix(MeanEndingLength)
  TimestepGrowthSizeInc <- MeanEndingLength

  # get key inputs for length transition matrices
  if (GrowthCurveType == 1) { # von Bertalanffy
    for (i in 1:2) {
      MeanEndingLength[i,] = midpt + (GrowthParams[i,1] - midpt) * (1 - exp(-GrowthParams[i,2]*TimeStep))
      TimestepGrowthSizeInc[i,] = MeanEndingLength[i,] - midpt # amount of annual growth with respect to initial length
    }
  }
  if (GrowthCurveType == 2) { # Schnute
    for (i in 1:2) {
      y1=GrowthParams[i,1]
      y2=GrowthParams[i,2]
      a=GrowthParams[i,3]
      b=GrowthParams[i,4]
      t1=RefnceAges[1]
      t2=RefnceAges[2]

      for (t in 1:MaxModelAge) {
        MeanSizeAtAge[i,t] =  SchnuteGrowthfunction(t, t1, t2, y1, y2, a, b)
      }
      GrowthParamsForSex = c(y1, y2, a, b)
      RefnceAgesForSex = c(t1, t2)
      MeanEndingLength[i,] = CalcLengthAfterGrowthForTimetep(GrowthCurveType=2, TimeStep, GrowthParamsForSex, RefnceAgesForSex, midpt)
      TimestepGrowthSizeInc[i,] = MeanEndingLength[i,] - midpt
    }
  }

  results = list(MeanEndingLength = MeanEndingLength,
                 TimestepGrowthSizeInc = TimestepGrowthSizeInc)

  return(results)

}

#' Calculates expected lengths and associated statistics from length-based per recruit analysis
#'
#' This function outputs random fish lengths, calculated from length-based per recruit analysis
#' with specified fishing mortality (PerRecFishLen), expected mean length (PerRecLenFreq),
#' 95 percent confidence limits for mean length (MeanCatchLen.95ci), and 60 and 95 percent
#' prediction intervals for mean length (MeanCatchLen.60pi, MeanCatchLen.95pi). Values are calculated
#' for females, males and combined sexes.
#'
#' @keywords internal
#'
#' @param midpt length class mid-points
#' @param FemCatchNum female per recruit catch numbers at length
#' @param MalCatchNum male per recruit catch numbers at length
#' @param CombSexCatchNum
#'
#' @return PerRecRandFishLen
#' @return PerRecMeanCatchLen
#' @return PerRecLenFreq
#' @return MeanCatchLen.95ci
#' @return MeanCatchLen.60pi
#' @return MeanCatchLen.95pi
#' @return PerRecRandFishWt
#' @return PerRecMeanCatchWt
#' @return MeanCatchWt.95ci
#' @return MeanCatchWt.60pi
#' @return MeanCatchWt.95pi
#'
CalcPerRecruitFishLenAndFishWtStats<- function(midpt, FemCatchNum, MalCatchNum, CombSexCatchNum,
                                               FemWtAtLen, MalWtAtLen) {

  # get approx. 95 and 60% confidence intervals for mean catch length
  FemProbs = FemCatchNum/sum(FemCatchNum)
  FemLenFreq = as.vector(rmultinom(1,10000,FemProbs))
  FemFishLen = rep(midpt, FemLenFreq)
  FemFishWt = rep(FemWtAtLen, FemLenFreq)
  MalProbs = MalCatchNum/sum(MalCatchNum)
  MalLenFreq = as.vector(rmultinom(1,10000,MalProbs))
  MalFishLen = rep(midpt, MalLenFreq)
  MalFishWt = rep(MalWtAtLen, MalLenFreq)
  CombSexProbs = CombSexCatchNum/sum(CombSexCatchNum)
  CombSexLenFreq = as.vector(rmultinom(1,10000,CombSexProbs))
  CombSexFishLen = rep(midpt, CombSexLenFreq)
  FishWtsBothSexes = c(FemFishWt, MalFishWt)
  CombSexFishWt = sample(FishWtsBothSexes, 10000, replace = FALSE, prob=NULL)

  # random fish lengths
  PerRecRandFishLen = data.frame(FemFishLen=FemFishLen,
                                 MalFishLen=MalFishLen,
                                 CombSexFishLen=CombSexFishLen)

  # random fish weights
  PerRecRandFishWt = data.frame(FemFishWt=FemFishWt,
                                MalFishWt=MalFishWt,
                                CombSexFishWt=CombSexFishWt)

  # random fish length frequency distribution
  PerRecLenFreq = data.frame(FemLenFreq=FemLenFreq,
                             MalLenFreq=MalLenFreq,
                             CombSexLenFreq=CombSexLenFreq)

  # calculate mean catch length
  # FemMeanCatchLen = sum(FemCatchNum * midpt) / sum(FemCatchNum)
  # MalMeanCatchLen = sum(MalCatchNum * midpt) / sum(MalCatchNum)
  # CombSexMeanCatchLen = sum((CombSexCatchNum * midpt)) / sum(CombSexCatchNum)
  FemMeanCatchLen = mean(FemFishLen)
  MalMeanCatchLen = mean(MalFishLen)
  CombSexMeanCatchLen = mean(CombSexFishLen)
  PerRecMeanCatchLen = data.frame(FemMeanCatchLen=FemMeanCatchLen,
                                  MalMeanCatchLen=MalMeanCatchLen,
                                  CombSexMeanCatchLen=CombSexMeanCatchLen)

  # calculate mean catch weight
  # FemMeanCatchWt = sum(FemCatchNum * FemWtAtLen) / sum(FemCatchNum)
  # MalMeanCatchWt = sum(MalCatchNum * MalWtAtLen) / sum(MalCatchNum)
  FemMeanCatchWt = mean(FemFishWt)
  MalMeanCatchWt = mean(MalFishWt)
  CombSexMeanCatchWt = mean(CombSexFishWt)
  PerRecMeanCatchWt = data.frame(FemMeanCatchWt=FemMeanCatchWt,
                                 MalMeanCatchWt=MalMeanCatchWt,
                                 CombSexMeanCatchWt=CombSexMeanCatchWt)

  # 95 percent confidence limits for mean length
  FemMeanCatchLen.95ci = FemMeanCatchLen + c(-1.96,1.96) * (sd(FemFishLen)/sqrt(10000))
  MalMeanCatchLen.95ci = MalMeanCatchLen + c(-1.96,1.96) * (sd(MalFishLen)/sqrt(10000))
  CombSexMeanCatchLen.95ci = CombSexMeanCatchLen + c(-1.96,1.96) * (sd(CombSexFishLen)/sqrt(10000))
  MeanCatchLen.95ci = data.frame(FemMeanCatchLen.95ci=FemMeanCatchLen.95ci,
                                 MalMeanCatchLen.95ci=MalMeanCatchLen.95ci,
                                 CombSexMeanCatchLen.95ci=CombSexMeanCatchLen.95ci)

  # 95 percent confidence limits for mean weight
  FemMeanCatchWt.95ci = FemMeanCatchWt + c(-1.96,1.96) * (sd(FemFishWt)/sqrt(10000))
  MalMeanCatchWt.95ci = MalMeanCatchWt + c(-1.96,1.96) * (sd(MalFishWt)/sqrt(10000))
  CombSexMeanCatchWt.95ci = CombSexMeanCatchWt + c(-1.96,1.96) * (sd(CombSexFishWt)/sqrt(10000))
  MeanCatchWt.95ci = data.frame(FemMeanCatchWt.95ci=FemMeanCatchWt.95ci,
                                MalMeanCatchWt.95ci=MalMeanCatchWt.95ci,
                                CombSexMeanCatchWt.95ci=CombSexMeanCatchWt.95ci)

  # calculate standard error of the prediction
  FemMSE = sqrt(sum((FemFishLen-FemMeanCatchLen)^2) / (10000 - 2))
  MalMSE = sqrt(sum((MalFishLen-MalMeanCatchLen)^2) / (10000 - 2))
  CombSexMSE = sqrt(sum((CombSexFishLen-CombSexMeanCatchLen)^2) / (10000 - 2))

  # 60 and 95 percent prediction intervals for mean catch length
  FemMeanCatchLen.95pi = FemMeanCatchLen + c(-1.96,1.96) * FemMSE * (sqrt(1 + (1/10000)))
  MalMeanCatchLen.95pi = MalMeanCatchLen + c(-1.96,1.96) * MalMSE * (sqrt(1 + (1/10000)))
  CombSexMeanCatchLen.95pi = CombSexMeanCatchLen + c(-1.96,1.96) * CombSexMSE * (sqrt(1 + (1/10000)))
  MeanCatchLen.95pi = data.frame(FemMeanCatchLen.95pi=FemMeanCatchLen.95pi,
                                 MalMeanCatchLen.95pi=MalMeanCatchLen.95pi,
                                 CombSexMeanCatchLen.95pi=CombSexMeanCatchLen.95pi)

  FemMeanCatchLen.60pi = FemMeanCatchLen + c(-0.84,0.84) * FemMSE * (sqrt(1 + (1/10000)))
  MalMeanCatchLen.60pi = MalMeanCatchLen + c(-0.84,0.84) * MalMSE * (sqrt(1 + (1/10000)))
  CombSexMeanCatchLen.60pi = CombSexMeanCatchLen + c(-0.84,0.84) * CombSexMSE * (sqrt(1 + (1/10000)))
  MeanCatchLen.60pi = data.frame(FemMeanCatchLen.60pi=FemMeanCatchLen.60pi,
                                 MalMeanCatchLen.60pi=MalMeanCatchLen.60pi,
                                 CombSexMeanCatchLen.60pi=CombSexMeanCatchLen.60pi)

  # calculate standard error of the prediction, for weight
  FemMSEwt = sqrt(sum((FemFishWt - FemMeanCatchWt)^2) / (10000 - 2))
  MalMSEwt = sqrt(sum((MalFishWt - MalMeanCatchWt)^2) / (10000 - 2))
  CombSexMSEwt = sqrt(sum((CombSexFishWt - CombSexMeanCatchWt)^2) / (10000 - 2))

  # 60 and 95 percent prediction intervals for mean catch weight
  FemMeanCatchWt.95pi = FemMeanCatchWt + c(-1.96,1.96) * FemMSEwt * (sqrt(1 + (1/10000)))
  MalMeanCatchWt.95pi = MalMeanCatchWt + c(-1.96,1.96) * MalMSEwt * (sqrt(1 + (1/10000)))
  CombSexMeanCatchWt.95pi = CombSexMeanCatchWt + c(-1.96,1.96) * CombSexMSEwt * (sqrt(1 + (1/10000)))
  MeanCatchWt.95pi = data.frame(FemMeanCatchWt.95pi=FemMeanCatchWt.95pi,
                                MalMeanCatchWt.95pi=MalMeanCatchWt.95pi,
                                CombSexMeanCatchWt.95pi=CombSexMeanCatchWt.95pi)

  FemMeanCatchWt.60pi = FemMeanCatchWt + c(-0.84,0.84) * FemMSEwt * (sqrt(1 + (1/10000)))
  MalMeanCatchWt.60pi = MalMeanCatchWt + c(-0.84,0.84) * MalMSEwt * (sqrt(1 + (1/10000)))
  CombSexMeanCatchWt.60pi = CombSexMeanCatchWt + c(-0.84,0.84) * CombSexMSEwt * (sqrt(1 + (1/10000)))
  MeanCatchWt.60pi = data.frame(FemMeanCatchWt.60pi=FemMeanCatchWt.60pi,
                                MalMeanCatchWt.60pi=MalMeanCatchWt.60pi,
                                CombSexMeanCatchWt.60pi=CombSexMeanCatchWt.60pi)

  # save random fish lengths and ages
  RandFishLen = data.frame(FemFishLen=FemFishLen,
                           MalFishLen=MalFishLen,
                           CombSexFishLen=CombSexFishLen)

  RandFishWt = data.frame(FemFishWt=FemFishWt,
                          MalFishWt=MalFishWt,
                          CombSexFishWt=CombSexFishWt)

  Results = list(PerRecRandFishLen=PerRecRandFishLen,
                 PerRecLenFreq=PerRecLenFreq,
                 PerRecMeanCatchLen=PerRecMeanCatchLen,
                 MeanCatchLen.95ci=MeanCatchLen.95ci,
                 MeanCatchLen.60pi=MeanCatchLen.60pi,
                 MeanCatchLen.95pi=MeanCatchLen.95pi,
                 PerRecRandFishWt=PerRecRandFishWt,
                 PerRecMeanCatchWt=PerRecMeanCatchWt,
                 MeanCatchWt.95ci=MeanCatchWt.95ci,
                 MeanCatchWt.60pi=MeanCatchWt.60pi,
                 MeanCatchWt.95pi=MeanCatchWt.95pi,
                 RandFishLen = RandFishLen,
                 RandFishWt = RandFishWt)

  return(Results)

}

#' Get outputs from length-based per recruit analysis for a specified value of fishing mortality
#'
#' This function provides outputs associated with per recruit analysis, and an
#' extended form of this analysis with a stock recruitment relationship to account
#' for potential impacts of fishing on recruitment
#'
#' @param MaxModelAge maximum age considered by model
#' @param TimeStep model time step (in y)
#' @param lbnd lower bounds of length classes
#' @param ubnd upper bounds of length classes
#' @param midpt mid points bounds of length classes
#' @param nLenCl number of length classes
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
#' @param GrowthParams growth parameters of either von Bertalanffy or Schnute model
#' @param RefnceAges reference ages for Schnute model, set to NA if using von Bertalanffy model
#' @param CVSizeAtAge coefficient of variation for size at age
#' @param lenwt_a weight-length parameter
#' @param ln_lenwt_a weight-length parameter
#' @param lenwt_b weight-length parameter
#' @param WLrel_Type 1=power, 2=log-log
#' @param EstWtAtLen user-specified weights at lengths
#' @param ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' @param ReprodPattern reproductive pattern, 1=gonochoristic, 2=protogynous (female to male sex change) hermaphroditism,
#' 3=protandrous hermaphroditism (male to female sex change)
#' @param InitRatioFem proportion of fish that are females at hatching
#' @param FinalSex_Pmax logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param mat_L50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_L95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtLen vector of proportion mature at length (set to NA if using age at maturity parameters)
#' @param sel_L50 logistic parameter for gear selectivity curve
#' @param sel_L95 logistic parameter for gear selectivity curve
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_L50 logistic parameter for fish retention curve
#' @param ret_L95 logistic parameter for fish retention curve
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param FMort fishing mortality
#'
#' @return length distribution of recruits (RecLenDist), unfished female and male numbers per recruit at age and length
#' (Unfish_FemNPerRec, Unfish_MalNPerRec), fished female and male numbers per recruit at age and length (Fish_FemNPerRec,
#' Fish_MalNPerRec), unfished female and male biomass per recruit at age and length (Unfish_FemBiomPerRecAtAge,
#' Unfish_MalBiomPerRecAtAge), fished female and male biomass per recruit at age and length (Fish_FemBiomPerRecAtAge,
#' Fish_MalBiomPerRecAtAge), unfished female and male biomass at age (Unfish_FemBiomAtAge, Unfish_MalBiomAtAge),
#' fished female and male biomass at age (Fish_FemBiomAtAge, Fish_MalBiomAtAge), female and male catch biomass (FemCatchBiom, MalCatchBiom),
#' yield per recruit for combined sexes (YPR), female, male and combined sex spawning potential ratio (Fem_SPR, Mal_SPR, CombSex_SPR),
#' derived Beverton-Holt stock recruitment parameters and associated equilibrium recruitment derived from those
#' parameters (BH_SRRa, BH_SRRb, BH_Equil_Rec), equilibrium recruitment for either Beverton-Holt or Ricker relationship (Equil_Rec),
#' equilibrium catch (Equil_Catch), female and male and spawning biomass (Equil_FemSpBiom, Equil_MalSpBiom), equilibrium relative
#' female, male and combined sex spawning biomass (Equilmod_FemRelBiom, Equilmod_MalRelBiom, Equilmod_CombSexRelBiom), mean size at
#' age of females and males (MeanSizeAtAge), female and male weight at length (FemWtAtLen, MalWtAtLen), female and male proportion
#' mature at length (FemPropMatAtLen, MalPropMatAtLen), female and male gear selectivity at length (FemGearSelAtLen, MalGearSelAtLen),
#' female and male retention at length (FemRetProbAtLen, MalRetProbAtLen), female and male selectivity of landings at length
#' (FemSelLandAtLen, MalSelLandAtLen), female and male selectivity of discards at length (FemSelDiscAtLen, MalSelDiscAtLen),
#' female and male discard mortality at length (FemDiscFAtLen, MalDiscFAtLen), female and male mortality associated with landings
#' at length (FemLandFAtLen, MalLandFAtLen), female and male fishing mortality at length (FemFAtLen, MalFAtLen),
#' female and male total mortality at length (FemZAtLen, MalZAtLen)
#'
#' @examples
#' MaxModelAge <- 20 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
#' MaxLen = 800
#' LenInc = 20
#' lbnd = seq(0,MaxLen - LenInc, LenInc)
#' ubnd = lbnd + LenInc
#' midpt = lbnd + (LenInc/2)
#' nLenCl = length(midpt)
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.2, 0.2) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' GrowthParams = data.frame(Linf=Linf, vbK=vbK, tzero=tzero)
#' RefnceAges = NA
#' # GrowthParams c(Linf, vbK, tzero) single sex von Bertalanffy, or data.frame(Linf=Linf, vbK=vbK, tzero=tzero),
#' #' both sexes von Bertalanffy, or c(y1, y2, a, b) single sex Schnute, or data.frame(y1=y1, y2=y2, a=a, b=b), both sexes Schnute
#' CVSizeAtAge = c(0.05,0.05)
#' lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- NA # for log-log relationship
#' lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' EstWtAtLen <- data.frame(EstFemWtAtLen=NA,
#'                          EstMalWtAtLen=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_L50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_L95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' mat_L50 <- c(250, 250) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_L95 <- c(300, 300) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
#'                           EstMalMatAtLen=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_L50 <- c(250, 250) # females, males - Logistic age selectivity relationship parameters
#' sel_L95 <- c(300, 300) # females, males - Logistic age selectivity relationship parameters
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_L50 <- c(250, 250) # females, males - Logistic age fish retention at age parameters
#' ret_L95 <- c(300, 300) # females, males - Logistic age fish retention at age parameters
#' DiscMort <- 0.25 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 4.22 / MaxModelAge # natural mortality  (year-1)
#' FMort = 0.2
#' Res=CalcYPRAndSPRForFMort_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
#'                              RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
#'                              EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
#'                              FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, ret_Pmax,
#'                              ret_L50, ret_L95, DiscMort, Steepness, SRrel_Type, NatMort, FMort)
#' @export
CalcYPRAndSPRForFMort_LB<- function(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                    RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
                                    EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
                                    FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, ret_Pmax,
                                    ret_L50, ret_L95, DiscMort, Steepness, SRrel_Type, NatMort, FMort) {

  if (FMort < 0.0001) FMort = 0.0001 # set negligible fishing mortality for unfished stock,
  # to allow calculations associated with expected catch size distribtion, i.e. from a survey

  # number of model time steps
  nTimeSteps <- 1 + (MaxModelAge / TimeStep)
  Ages <- seq(TimeStep,MaxModelAge, TimeStep)

  # set up storage
  RecLenDist <- data.frame(matrix(nrow = 2, ncol = nLenCl))
  colnames(RecLenDist) <- round(midpt,2)
  RecLenDist <- as.matrix(RecLenDist)

  # calculate mean size at age
  MeanSizeAtAge = CalcMeanSizeAtAgeForPerRecruitAnalysis(GrowthCurveType, MaxModelAge, Ages, GrowthParams, RefnceAges)

  # calculate weight at length
  res=CalcWeightAtLengthForPerRecruitAnalysis(lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen, midpt)
  FemWtAtLen = res$FemWtAtLen
  MalWtAtLen = res$MalWtAtLen

  # for sex-changing species, calculate proportion female at age
  if (ReprodPattern == 1) { # gonochorism
    PropFemAtLen = NA
  }
  if (ReprodPattern == 2) { # protogynous hermaphroditism (female to male sex change)
    PropFemAtLen = 1 - (FinalSex_Pmax / (1 + exp(-log(19) * (midpt - FinalSex_L50) / (FinalSex_L95 - FinalSex_L50))))
  }
  if (ReprodPattern == 3) { # protandrous hermaphroditism (male to female sex change)
    PropFemAtLen = FinalSex_Pmax / (1 + exp(-log(19) * (midpt - FinalSex_L50) / (FinalSex_L95 - FinalSex_L50)))
  }

  # calculate proportion mature at length
  if (is.na(EstMatAtLen[1,1])) {
    FemPropMatAtLen <- 1 / (1 + exp(-log(19) * (midpt - mat_L50[1]) / (mat_L95[1] - mat_L50[1])))
    MalPropMatAtLen <- 1 / (1 + exp(-log(19) * (midpt - mat_L50[2]) / (mat_L95[2] - mat_L50[2])))
  }
  if (!is.na(EstMatAtLen[1,1])) {
    FemPropMatAtLen <- EstMatAtLen[,1]
    MalPropMatAtLen <- EstMatAtLen[,2]
  }


  # Calculate gear selectivity at age
  FemGearSelAtLen <- 1/(1+exp(-log(19) * (midpt - sel_L50[1]) / (sel_L95[1] - sel_L50[1])))
  MalGearSelAtLen <- 1/(1+exp(-log(19) * (midpt - sel_L50[2]) / (sel_L95[2] - sel_L50[2])))

  # Calculate fish retention at age
  FemRetProbAtLen <- ret_Pmax[1] / (1+exp(-log(19) * (midpt - ret_L50[1]) / (ret_L95[1] - ret_L50[1])))
  MalRetProbAtLen <- ret_Pmax[2] / (1+exp(-log(19) * (midpt - ret_L50[2]) / (ret_L95[2] - ret_L50[2])))

  # Calculate selectivity of landings at age
  FemSelLandAtLen <- FemGearSelAtLen * FemRetProbAtLen
  MalSelLandAtLen <- MalGearSelAtLen * MalRetProbAtLen

  # Calculate selectivity of discards at age
  FemSelDiscAtLen <- FemGearSelAtLen * (1 - FemRetProbAtLen)
  MalSelDiscAtLen <- MalGearSelAtLen * (1 - MalRetProbAtLen)

  # calculate female and male mortality at age associated with discarding of undersize fish
  FemDiscFAtLen <- FMort * DiscMort * FemSelDiscAtLen
  MalDiscFAtLen <- FMort * DiscMort * MalSelDiscAtLen

  # calculate female and male mortality at age associated with landings
  FemLandFAtLen <- FMort * FemSelLandAtLen
  MalLandFAtLen <- FMort * MalSelLandAtLen

  # calculate (total) female and male fishing mortality at age
  if (DiscMort >= 0.001) {
    FemFAtLen <- FMort * (FemSelLandAtLen + (DiscMort * FemSelDiscAtLen))
    MalFAtLen <- FMort * (MalSelLandAtLen + (DiscMort * MalSelDiscAtLen))
  } else {
    FemFAtLen <- FMort * FemSelLandAtLen
    MalFAtLen <- FMort * MalSelLandAtLen
  }

  # calculate female and male total mortality at age
  FemZAtLen <- FemFAtLen + NatMort
  MalZAtLen <- MalFAtLen + NatMort

  # calculate initial size distribution of female and male recruits
  RecLenDist = CalcSizeDistOfRecruits(MeanSizeAtAge, CVSizeAtAge, lbnd, ubnd, midpt, nLenCl)

  # get required inputs to calculate length transition matrices
  results = GetLTMInputsForPerRecruitAnalysis(GrowthCurveType, TimeStep, MaxModelAge, GrowthParams, RefnceAges, midpt)
  MeanEndingLength = results$MeanEndingLength
  TimestepGrowthSizeInc = results$TimestepGrowthSizeInc

  # calculate female and male length transition matrices
  TimestepGrowthSizeInc_Fem = TimestepGrowthSizeInc[1,]
  CV_Fem = CVSizeAtAge[1]
  LTM_Fem = CalcLTM_cpp(TimestepGrowthSizeInc_Fem, CV_Fem, lbnd, midpt, ubnd, nLenCl)
  CV_Mal = CVSizeAtAge[2]
  TimestepGrowthSizeInc_Mal = TimestepGrowthSizeInc[2,]

  LTM_Mal = CalcLTM_cpp(TimestepGrowthSizeInc_Mal, CV_Mal, lbnd, midpt, ubnd, nLenCl)

  # update survival and growth
  PopnRes = UpdateGrowthAndSurvival_cpp(ReprodPattern, TimeStep, nTimeSteps, nLenCl, InitRatioFem, RecLenDist, NatMort, FemZAtLen,
                                        MalZAtLen, PropFemAtLen, LTM_Fem, LTM_Mal, FemWtAtLen, MalWtAtLen, ReprodScale,
                                        FemPropMatAtLen, MalPropMatAtLen)

  Unfish_FemNPerRecAtAge=PopnRes$Unfish_FemNPerRecAtAge
  Fish_FemNPerRecAtAge=PopnRes$Fish_FemNPerRecAtAge
  Unfish_MalNPerRecAtAge=PopnRes$Unfish_MalNPerRecAtAge
  Fish_MalNPerRecAtAge=PopnRes$Fish_MalNPerRecAtAge
  Unfish_FemNPerRec=PopnRes$Unfish_FemNPerRec
  Fish_FemNPerRec=PopnRes$Fish_FemNPerRec
  Unfish_MalNPerRec=PopnRes$Unfish_MalNPerRec
  Fish_MalNPerRec=PopnRes$Fish_MalNPerRec
  Unfish_FemBiomPerRecAtAge = PopnRes$Unfish_FemBiomPerRecAtAge
  Fish_FemBiomPerRecAtAge = PopnRes$Fish_FemBiomPerRecAtAge
  Unfish_MalBiomPerRecAtAge = PopnRes$Unfish_MalBiomPerRecAtAge
  Fish_MalBiomPerRecAtAge = PopnRes$Fish_MalBiomPerRecAtAge
  Unfish_FemBiomAtAge = PopnRes$Unfish_FemBiomAtAge
  Fish_FemBiomAtAge = PopnRes$Fish_FemBiomAtAge
  Unfish_MalBiomAtAge = PopnRes$Unfish_MalBiomAtAge
  Fish_MalBiomAtAge = PopnRes$Fish_MalBiomAtAge

  # Unfished spawning biomass in kg
  UnfishFemSpawnBiom <- sum(Unfish_FemBiomAtAge)
  UnfishMalSpawnBiom <- sum(Unfish_MalBiomAtAge)
  UnfishCombSexSpawnBiom <- UnfishFemSpawnBiom + UnfishMalSpawnBiom

  # Fished spawning biomass in kg
  FishFemSpawnBiom <- sum(Fish_FemBiomAtAge)
  FishMalSpawnBiom <- sum(Fish_MalBiomAtAge)
  FishCombSexSpawnBiom <- FishFemSpawnBiom + FishMalSpawnBiom

  # calculate catch yield (biomass)
  FemCatchBiom = (FemLandFAtLen/FemZAtLen) * (1-exp(-FemZAtLen)) * Fish_FemNPerRec * FemWtAtLen
  MalCatchBiom = (MalLandFAtLen/MalZAtLen) * (1-exp(-MalZAtLen)) * Fish_MalNPerRec * MalWtAtLen
  CombSexCatchBiom = FemCatchBiom + MalCatchBiom

  # calculate catch yield (numbers)
  FemCatchNum = (FemLandFAtLen/FemZAtLen) * (1-exp(-FemZAtLen)) * Fish_FemNPerRec
  MalCatchNum = (MalLandFAtLen/MalZAtLen) * (1-exp(-MalZAtLen)) * Fish_MalNPerRec
  CombSexCatchNum = FemCatchNum + MalCatchNum

  # Random data and stats for expected lengths and weights in catches
  CatchLenWtRes = CalcPerRecruitFishLenAndFishWtStats(midpt, FemCatchNum, MalCatchNum, CombSexCatchNum,
                                                 FemWtAtLen, MalWtAtLen)

  PerRecRandFishLen = CatchLenWtRes$PerRecRandFishLen
  PerRecLenFreq = CatchLenWtRes$PerRecLenFreq
  PerRecMeanCatchLen = CatchLenWtRes$PerRecMeanCatchLen
  MeanCatchLen.95ci = CatchLenWtRes$MeanCatchLen.95ci
  MeanCatchLen.60pi = CatchLenWtRes$MeanCatchLen.60pi
  MeanCatchLen.95pi = CatchLenWtRes$MeanCatchLen.95pi
  PerRecRandFishWt = CatchLenWtRes$PerRecRandFishWt
  PerRecMeanCatchWt = CatchLenWtRes$PerRecMeanCatchWt
  MeanCatchWt.95ci = CatchLenWtRes$MeanCatchWt.95ci
  MeanCatchWt.60pi = CatchLenWtRes$MeanCatchWt.60pi
  MeanCatchWt.95pi = CatchLenWtRes$MeanCatchWt.95pi
  RandFishLen = CatchLenWtRes$RandFishLen
  RandFishWt = CatchLenWtRes$RandFishWt

  # calculate yield per recruit in kg
  YPR <- max(0,sum(FemCatchBiom) + sum(MalCatchBiom))

  # calculate spawning potential ratio (SPR)
  Fem_SPR <- max(0,FishFemSpawnBiom / UnfishFemSpawnBiom)
  Mal_SPR <- max(0,FishMalSpawnBiom / UnfishMalSpawnBiom)
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
  #check = (Unfish_FemBiomPerRecSpawnSeas-SR_alpha)/(SR_beta*Unfish_FemBiomPerRecSpawnSeas)

  # calculate equilibrium catch
  Equil_Catch <- max(0,Equil_Rec * YPR)

  # calculate equilibrium female spawning biomass
  Equil_FemSpBiom <- Equil_Rec * FishFemSpawnBiom
  Equil_MalSpBiom <- Equil_Rec * FishMalSpawnBiom

  # calculate equilibrium model SPR
  Equilmod_FemRelBiom <- max(0,Equil_FemSpBiom / UnfishFemSpawnBiom)
  Equilmod_MalRelBiom <- max(0,Equil_MalSpBiom / UnfishMalSpawnBiom)
  Equilmod_CombSexRelBiom <- max(0,(Equil_FemSpBiom + Equil_MalSpBiom) / (UnfishFemSpawnBiom + UnfishMalSpawnBiom))

  Results = list(RecLenDist=RecLenDist,
                 PropFemAtLen=PropFemAtLen,
                 Unfish_FemNPerRec=Unfish_FemNPerRec,
                 Unfish_MalNPerRec=Unfish_MalNPerRec,
                 Fish_FemNPerRec=Fish_FemNPerRec,
                 Fish_MalNPerRec=Fish_MalNPerRec,
                 Unfish_FemNPerRecAtAge=Unfish_FemNPerRecAtAge,
                 Unfish_MalNPerRecAtAge=Unfish_MalNPerRecAtAge,
                 Fish_FemNPerRecAtAge=Fish_FemNPerRecAtAge,
                 Fish_MalNPerRecAtAge=Fish_MalNPerRecAtAge,
                 Unfish_FemBiomPerRecAtAge=Unfish_FemBiomPerRecAtAge,
                 Unfish_MalBiomPerRecAtAge=Unfish_MalBiomPerRecAtAge,
                 Fish_FemBiomPerRecAtAge=Fish_FemBiomPerRecAtAge,
                 Fish_MalBiomPerRecAtAge=Fish_MalBiomPerRecAtAge,
                 Unfish_FemBiomAtAge=Unfish_FemBiomAtAge,
                 Unfish_MalBiomAtAge=Unfish_MalBiomAtAge,
                 Fish_FemBiomAtAge=Fish_FemBiomAtAge,
                 Fish_MalBiomAtAge=Fish_MalBiomAtAge,
                 FemCatchBiom=FemCatchBiom,
                 MalCatchBiom=MalCatchBiom,
                 FemCatchBiom=FemCatchBiom,
                 FemCatchNum=FemCatchNum,
                 MalCatchNum=MalCatchNum,
                 CombSexCatchNum=CombSexCatchNum,
                 PerRecRandFishLen=PerRecRandFishLen,
                 PerRecLenFreq=PerRecLenFreq,
                 PerRecMeanCatchLen=PerRecMeanCatchLen,
                 MeanCatchLen.95ci=MeanCatchLen.95ci,
                 MeanCatchLen.60pi=MeanCatchLen.60pi,
                 MeanCatchLen.95pi=MeanCatchLen.95pi,
                 RandFishLen = RandFishLen,
                 PerRecRandFishWt=PerRecRandFishWt,
                 PerRecMeanCatchWt=PerRecMeanCatchWt,
                 MeanCatchWt.95ci=MeanCatchWt.95ci,
                 MeanCatchWt.60pi=MeanCatchWt.60pi,
                 MeanCatchWt.95pi=MeanCatchWt.95pi,
                 RandFishWt = RandFishWt,
                 YPR = YPR,
                 Fem_SPR = Fem_SPR,
                 Mal_SPR = Mal_SPR,
                 CombSex_SPR = CombSex_SPR,
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
                 MeanSizeAtAge = MeanSizeAtAge,
                 FemWtAtLen = FemWtAtLen,
                 MalWtAtLen = MalWtAtLen,
                 FemPropMatAtLen = FemPropMatAtLen,
                 MalPropMatAtLen = MalPropMatAtLen,
                 FemGearSelAtLen = FemGearSelAtLen,
                 MalGearSelAtLen = MalGearSelAtLen,
                 FemRetProbAtLen = FemRetProbAtLen,
                 MalRetProbAtLen = MalRetProbAtLen,
                 FemSelLandAtLen = FemSelLandAtLen,
                 MalSelLandAtLen = MalSelLandAtLen,
                 FemSelDiscAtLen = FemSelDiscAtLen,
                 MalSelDiscAtLen = MalSelDiscAtLen,
                 FemDiscFAtLen = FemDiscFAtLen,
                 MalDiscFAtLen = MalDiscFAtLen,
                 FemLandFAtLen = FemLandFAtLen,
                 MalLandFAtLen = MalLandFAtLen,
                 FemFAtLen=FemFAtLen,
                 MalFAtLen=MalFAtLen,
                 FemZAtLen=FemZAtLen,
                 MalZAtLen=MalZAtLen)

  return(Results)
}



#' Get outputs from age-based per recruit analysis across a range of fishing mortality values
#'
#' This function provides outputs associated with per recruit analysis, and an
#' extended form of this analysis with a Beverton-Holt stock recruitment relationship to account
#' for potential impacts of fishing on recruitment. Outputs are provided for a range of
#' fishing mortality values, including the current, estimated value.
#'
#' @param MaxModelAge maximum age considered in model
#' @param TimeStep model timestep (e.g. 1 = annual, 1/12 = monthly)
#' @param Linf asymptotic length (von Bertalanffy growth parameter, set to NA if inputting lengths at ages directly from any growth model)
#' @param vbK growth coefficient (von Bertalanffy growth parameter, set to NA if inputting lengths at ages directly from any growth model)
#' @param tzero hypothetical age at zero length (von Bertalanffy growth parameter, set to NA if inputting lengths at ages directly from any growth model)
#' @param EstLenAtAge vector or estimated lengths at ages from any growth model, set to NA if von Bertalanffy growth parameters specified
#' @param lenwt_a weight-length parameter (power relationship)
#' @param ln_lenwt_a weight-length parameter (log-log relationship)
#' @param lenwt_b weight-length parameter (power or log-log relationship)
#' @param WLrel_Type 1=power, 2=log-log relationship (set to NA if inputting weights at ages directly)
#' @param EstWtAtAge vector of weights at ages (set to NA if weight-length growth parameters specified)
#' @param ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' @param ReprodPattern reproductive pattern, 1=gonochoristic, 2=protogynous (female to male sex change) hermaphroditism,
#' 3=protandrous hermaphroditism (male to female sex change)
#' @param InitRatioFem proportion of fish that are females at hatching
#' @param FinalSex_Pmax logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_A50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_A95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param mat_A50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_A95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtAge vector of proportion mature at age (set to NA if using age at maturity parameters)
#' @param sel_A50 logistic parameter for gear selectivity curve
#' @param sel_A95 logistic parameter for gear selectivity curve
#' @param EstSelAtAge vector for gear selectivity at age (set to NA if using age at selectivity parameters)
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_A50 logistic parameter for fish retention curve
#' @param ret_A95 logistic parameter for fish retention curve
#' @param EstRetenAtAge vector for fish retention at age (set to NA if using age at retention parameters)
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param Current_F estimate of current fishing mortality
#'
#' @return Age at each timestep (Ages), female and male length at age (FemLenAtAge, MalLenAtAge),
#' female and male weight at age (FemWtAtAge, FemWtAtAge), proportion female at age (PropFemAtAge), female and male
#' proportion mature at age (FemPropMatAtAge, MalPropMatAtAge), female and male selectivity at age (FemSelAtAge, MalSelAtAge),
#' yield-per-recruit (YPR), unfished female, male and combined sex spawning potential ratio (Fem_SPR, Mal_SPR,
#' CombSex_SPR), unfished female and male per recruit survival at age (UnfishFemSurvAtAge, UnfishMalSurvAtAge),
#' unfished female and male mature biomass at age (UnfishFemBiomAtAge, UnfishMalBiomAtAge), fished female and
#' male spawning potential ratio (Fem_SPR, Mal_SPR), fished female and male per recruit survival at age
#' (FishedFemSurvAtAge, FishedMalSurvAtAge), fished female and male mature biomass at age (FishedFemBiomAtAge,
#' FishedMalBiomAtAge), female and male per recruit retained catches at age in numbers (FemCatchAtAgeNum, MalCatchAtAgeNum),
#' female and male per recruit released catches at age in numbers, (FemRelCatchAtAgeNum, MalRelCatchAtAgeNum),
#' female and male per recruit retained catches at age in biomass (FemCatchAtAge, MalCatchAtAge), and for the extended model, equilibrium recruitment (Equil_Rec),
#' equilibrium catch (Equil_Catch), equilibrium female and male spawning biomass (Equil_FemSpBiom, Equil_MalSpBiom) and
#' relative female, male and combined sex spawning biomass (Equilmod_FemRelBiom, Equilmod_MalRelBiom, Equilmod_CombSexRelBiom),
#' Beverton-Holt stock-recruitment parameters (BH_SRRa, BH_SRRb), equilibrium recruitment from Beverton-Holt relationship (BH_Equil_Rec),
#' female and male gear selectivity at age (FemGearSelAtAge, MalGearSelAtAge), female and male retention at age probabilities,
#' (FemRetProbAtAge, MalRetProbAtAge), selectivity of female and male fish landings (FemSelLandAtAge, MalSelLandAtAge),
#' female and male selectivity of discards (FemSelDiscAtAge, MalSelDiscAtAge), female and male fishing mortality
#' associated with discarding (FemDiscFAtAge, MalDiscFAtAge), female and male fishing mortality associated with
#' fish landings (FemLandFAtAge, MalLandFAtAge), female and male total fishing mortality at age (FemFAtAge, MalFAtAge),
#' female and male total mortality at age (FemZAtAge, MalZAtAge)range of fishing mortality values applied for which per quantities are calculated
#' (FishMort = seq(0,2,0.01)), equilibrium recruitment vs FMort (Equil_Rec), equilibrium catch vs FMort (Equil_Catch),
#' equilibrium female spawning biomass vs FMort (Equil_FemSpBiom), equilibrium spawning potential ratio vs FMort
#' (Equilmod_SPR), maximum yield per recruit (maxypr), maximum equilibrium catch (maxeqCatch), fishing mortality
#' associated with maxypr (Fmax), fishing mortality associated with maxeqCatch (F_MSY or FmaxeqCatch), biomass target at 1.2maxeqCatch,
#' (BMSY_Targ), biomass threshold at maxeqCatch (BMSY_Thresh), biomass threshold at 0.5maxeqCatch (BMSY_Lim),
#' YPR vs FishMort (YPRResults), Equil_Catch vs FishMort (EquilCatchResults), Equilmod_FemRelBiom vs FishMort (Fem_SPRResults),
#' Equilmod_MalRelBiom vs FishMort (Mal_SPRResults), Equilmod_FemRelBiom vs FishMort (Equilmod_FemRelBiomResults),
#' Equilmod_MalRelBiom vs FishMort (Equilmod_MalRelBiomResults), Equilmod_CombSexRelBiom vs FishMort (CombSex_SPRResults),
#' Equil_Rec vs FishMort (EquilRecResults)
#'
#' @examples
#' # Example 1. Non-hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 20 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
#' Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.5, 0.5) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstFemLenAtAge=NA, EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- NA # for log-log relationship
#' lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA, EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' EstSelAtAge <- data.frame(EstFemSelAtAge=NA, EstMalSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
#' EstRetenAtAge <- data.frame(EstFemRetenAtAge=rep(1,MaxModelAge+1), EstMalRetenAtAge=rep(1,MaxModelAge+1)) # fish retention at age (from age 0), inputted as values in data frame
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 0.2 # natural mortality  (year-1)
#' Current_F <- 0.07 # estimate of fishing mortality, e.g. from catch curve analysis
#' Res = GetPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                            lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
#'                            ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
#'                            mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
#'                            EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, Current_F)
#' # Example 2: hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 100 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
#' Linf <- c(1000, 1000) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.1, 0.1) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstMalLenAtAge=NA, EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- -11.0 # for log-log relationship
#' lenwt_b <- 3.0 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 2 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA, EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 2 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 1 # Ratio of females to males at age zero
#' FinalSex_Pmax <- 1 # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_A50 <- 35 # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- 60 # Logistic sex change relationship parameters (95% of max probability)
#' mat_A50 <- c(20, 20) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(30, 30) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_A50 <- c(15, 15) # females, males - Logistic age selectivity relationship parameters
#' sel_A95 <- c(25, 25) # females, males - Logistic age selectivity relationship parameters
#' EstSelAtAge <- data.frame(EstFemSelAtAge=NA, EstMalSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
#' EstRetenAtAge <- data.frame(EstFemRetenAtAge=rep(1,MaxModelAge+1), EstMalRetenAtAge=rep(1,MaxModelAge+1)) # fish retention at age (from age 0), inputted as values in data frame
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 0.07 # natural mortality  (year-1)
#' Current_F <- 0.07 # estimate of fishing mortality, e.g. from catch curve analysis
#' Res = GetPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                            lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
#'                            ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
#'                            mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
#'                            EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, Current_F)
#' @export
GetPerRecruitResults_AB <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                 lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                                 ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                                 mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge,
                                 ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness,
                                 SRrel_Type, NatMort, Current_F) {

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
    Res = CalcYPRAndSPRForFMort_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                                ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                                mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge,
                                ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness,
                                SRrel_Type, NatMort, FMort)
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
  Res2 = CalcYPRAndSPRForFMort_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                               lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                               ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                               mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
                               EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, FMort)



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
                 FemRelCatchAtAgeNum = Res2$FemRelCatchAtAgeNum,
                 MalRelCatchAtAgeNum = Res2$MalRelCatchAtAgeNum,
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

#' Get outputs from length-based per recruit analysis across a range of fishing mortality values
#'
#' This function provides outputs associated with per recruit analysis, and an
#' extended form of this analysis with a Beverton-Holt stock recruitment relationship to account
#' for potential impacts of fishing on recruitment. Outputs are provided for a range of
#' fishing mortality values, including the current, estimated value.
#'
#' @param MaxModelAge maximum age considered by model
#' @param TimeStep model time step (in y)
#' @param lbnd lower bounds of length classes
#' @param ubnd upper bounds of length classes
#' @param midpt mid points bounds of length classes
#' @param nLenCl number of length classes
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
#' @param GrowthParams growth parameters of either von Bertalanffy or Schnute model
#' @param RefnceAges reference ages for Schnute model, set to NA if using von Bertalanffy model
#' @param CVSizeAtAge coefficient of variation for size at age
#' @param lenwt_a weight-length parameter
#' @param ln_lenwt_a weight-length parameter
#' @param lenwt_b weight-length parameter
#' @param WLrel_Type 1=power, 2=log-log
#' @param EstWtAtLen user-specified weights at lengths
#' @param ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)'
#' @param ReprodPattern reproductive pattern, 1=gonochoristic, 2=protogynous (female to male sex change) hermaphroditism,
#' 3=protandrous hermaphroditism (male to female sex change)
#' @param InitRatioFem proportion of fish that are females at hatching
#' @param FinalSex_Pmax logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param mat_L50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_L95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtLen vector of proportion mature at length (set to NA if using age at maturity parameters)
#' @param sel_L50 logistic parameter for gear selectivity curve
#' @param sel_L95 logistic parameter for gear selectivity curve
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_L50 logistic parameter for fish retention curve
#' @param ret_L95 logistic parameter for fish retention curve
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param Current_F estimated current fishing mortality
#'
#' @return yield per recruit for combined sexes (YPR), female, male and combined sex spawning potential ratio
#' (Fem_SPR, Mal_SPR, CombSex_SPR), female and male length at age (FemLenAtAge, MalLenAtAge),
#' female and male weight at length (FemWtAtLen, MalWtAtLen), female and male proportion mature at length
#' (FemPropMatAtLen, MalPropMatAtLen), unfished female and male numbers per recruit at age and length
#' (Unfish_FemNPerRec, Unfish_MalNPerRec), fished female and male numbers per recruit at age and length (Fish_FemNPerRec,
#' Fish_MalNPerRec), unfished female and male biomass per recruit at age and length (Unfish_FemBiomPerRecAtAge,
#' Unfish_MalBiomPerRecAtAge), fished female and male biomass per recruit at age and length (Fish_FemBiomPerRecAtAge,
#' Fish_MalBiomPerRecAtAge), unfished female and male biomass at age (Unfish_FemBiomAtAge, Unfish_MalBiomAtAge),
#' fished female and male biomass at age (Fish_FemBiomAtAge, Fish_MalBiomAtAge), female and male catch biomass
#' (FemCatchBiom, MalCatchBiom), equilibrium recruitment for either Beverton-Holt or Ricker relationship (Equil_Rec),
#' equilibrium catch (Equil_Catch), female and male and spawning biomass (Equil_FemSpBiom, Equil_MalSpBiom), equilibrium relative
#' female, male and combined sex spawning biomass (Equilmod_FemRelBiom, Equilmod_MalRelBiom, Equilmod_CombSexRelBiom),
#' parameters and equilbrium recruitment from Beverton-Holt stock-recruitment relationship (BH_SRRa, BH_SRRb, BH_Equil_Rec),
#' female and male gear selectivity at length (FemGearSelAtLen, MalGearSelAtLen), female and male retention at length
#' (FemRetProbAtLen, MalRetProbAtLen), female and male selectivity of landings at length (FemSelLandAtLen, MalSelLandAtLen),
#' female and male selectivity of discards at length (FemSelDiscAtLen, MalSelDiscAtLen),
#' female and male discard mortality at length (FemDiscFAtLen, MalDiscFAtLen), female and male mortality associated with landings
#' at length (FemLandFAtLen, MalLandFAtLen), female and male fishing mortality at length (FemFAtLen, MalFAtLen),
#' female and male total mortality at length (FemZAtLen, MalZAtLen), fishing mortality values evaluated (FishMort = seq(0,2,0.01)),
#' maximum yield per recruit (maxypr), fishing mortality associated with maxeqCatch (F_MSY or FmaxeqCatch),
#' biomass target at 1.2maxeqCatch, (BMSY_Targ), biomass threshold at maxeqCatch (BMSY_Thresh), biomass threshold at 0.5maxeqCatch (BMSY_Lim),
#' YPR vs FishMort (YPRResults), Equil_Catch vs FishMort (EquilCatchResults), Equilmod_FemRelBiom vs FishMort (Fem_SPRResults),
#' Equilmod_MalRelBiom vs FishMort (Mal_SPRResults), Equilmod_FemRelBiom vs FishMort (Equilmod_FemRelBiomResults),
#' Equilmod_MalRelBiom vs FishMort (Equilmod_MalRelBiomResults), Equilmod_CombSexRelBiom vs FishMort (Equilmod_CombSexRelBiomResults),
#' Equil_Rec vs FishMort (EquilRecResults)
#'
#' @examples
#' MaxModelAge <- 20 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
#' MaxLen = 800
#' LenInc = 20
#' lbnd = seq(0,MaxLen - LenInc, LenInc)
#' ubnd = lbnd + LenInc
#' midpt = lbnd + (LenInc/2)
#' nLenCl = length(midpt)
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.2, 0.2) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' GrowthParams = data.frame(Linf=Linf, vbK=vbK, tzero=tzero)
#' RefnceAges = NA
#' # GrowthParams c(Linf, vbK, tzero) single sex von Bertalanffy, or data.frame(Linf=Linf, vbK=vbK, tzero=tzero),
#' #' both sexes von Bertalanffy, or c(y1, y2, a, b) single sex Schnute, or data.frame(y1=y1, y2=y2, a=a, b=b), both sexes Schnute
#' CVSizeAtAge = c(0.05,0.05)
#' lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- NA # for log-log relationship
#' lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' EstWtAtLen <- data.frame(EstFemWtAtLen=NA,
#'                          EstMalWtAtLen=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_L50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_L95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' mat_L50 <- c(250, 250) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_L95 <- c(300, 300) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
#'                           EstMalMatAtLen=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_L50 <- c(250, 250) # females, males - Logistic age selectivity relationship parameters
#' sel_L95 <- c(300, 300) # females, males - Logistic age selectivity relationship parameters
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_L50 <- c(250, 250) # females, males - Logistic age fish retention at age parameters
#' ret_L95 <- c(300, 300) # females, males - Logistic age fish retention at age parameters
#' DiscMort <- 0.25 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 4.22 / MaxModelAge # natural mortality  (year-1)
#' Current_F = 0.2
#' Res=GetPerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
#'                             RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
#'                             EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
#'                             FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, ret_Pmax,
#'                             ret_L50, ret_L95, DiscMort, Steepness, SRrel_Type, NatMort, Current_F)
#' @export
GetPerRecruitResults_LB <- function(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                    RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
                                    EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
                                    FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, ret_Pmax,
                                    ret_L50, ret_L95, DiscMort, Steepness, SRrel_Type, NatMort, Current_F) {

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

  FemMeanCatchLen <- rep(0,nFVals) # mean length of catch
  MalMeanCatchLen <- rep(0,nFVals)
  CombSexMeanCatchLen <- rep(0,nFVals)
  FemMeanCatchLen.lw95ci <- rep(0,nFVals) # 95 percent confidence limits
  FemMeanCatchLen.up95ci <- rep(0,nFVals)
  MalMeanCatchLen.lw95ci <- rep(0,nFVals)
  MalMeanCatchLen.up95ci <- rep(0,nFVals)
  CombSexMeanCatchLen.lw95ci <- rep(0,nFVals)
  CombSexMeanCatchLen.up95ci <- rep(0,nFVals)
  FemMeanCatchLen.lw60pi <- rep(0,nFVals) # 60 percent confidence prediction limits
  FemMeanCatchLen.up60pi <- rep(0,nFVals)
  MalMeanCatchLen.lw60pi <- rep(0,nFVals)
  MalMeanCatchLen.up60pi <- rep(0,nFVals)
  CombSexMeanCatchLen.lw60pi <- rep(0,nFVals)
  CombSexMeanCatchLen.up60pi <- rep(0,nFVals)
  FemMeanCatchLen.lw95pi <- rep(0,nFVals) # 95 percent confidence prediction limits
  FemMeanCatchLen.up95pi <- rep(0,nFVals)
  MalMeanCatchLen.lw95pi <- rep(0,nFVals)
  MalMeanCatchLen.up95pi <- rep(0,nFVals)
  CombSexMeanCatchLen.lw95pi <- rep(0,nFVals)
  CombSexMeanCatchLen.up95pi <- rep(0,nFVals)

  FemMeanCatchWt <- rep(0,nFVals) # mean weight of catch
  MalMeanCatchWt <- rep(0,nFVals)
  CombSexMeanCatchWt <- rep(0,nFVals)
  FemMeanCatchWt.lw95ci <- rep(0,nFVals) # 95 percent confidence limits
  FemMeanCatchWt.up95ci <- rep(0,nFVals)
  MalMeanCatchWt.lw95ci <- rep(0,nFVals)
  MalMeanCatchWt.up95ci <- rep(0,nFVals)
  CombSexMeanCatchWt.lw95ci <- rep(0,nFVals)
  CombSexMeanCatchWt.up95ci <- rep(0,nFVals)
  FemMeanCatchWt.lw60pi <- rep(0,nFVals) # 60 percent confidence prediction limits
  FemMeanCatchWt.up60pi <- rep(0,nFVals)
  MalMeanCatchWt.lw60pi <- rep(0,nFVals)
  MalMeanCatchWt.up60pi <- rep(0,nFVals)
  CombSexMeanCatchWt.lw60pi <- rep(0,nFVals)
  CombSexMeanCatchWt.up60pi <- rep(0,nFVals)
  FemMeanCatchWt.lw95pi <- rep(0,nFVals) # 95 percent confidence prediction limits
  FemMeanCatchWt.up95pi <- rep(0,nFVals)
  MalMeanCatchWt.lw95pi <- rep(0,nFVals)
  MalMeanCatchWt.up95pi <- rep(0,nFVals)
  CombSexMeanCatchWt.lw95pi <- rep(0,nFVals)
  CombSexMeanCatchWt.up95pi <- rep(0,nFVals)


  for (k in 1:nFVals) {
    FMort = FishMort[k]
    # cat("k",k,"FMort",FMort,'\n')
    Res = CalcYPRAndSPRForFMort_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                   RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
                                   EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
                                   FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, ret_Pmax,
                                   ret_L50, ret_L95, DiscMort, Steepness, SRrel_Type, NatMort, FMort)
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

    FemMeanCatchLen[k] <- Res$PerRecMeanCatchLen$FemMeanCatchLen
    MalMeanCatchLen[k] <- Res$PerRecMeanCatchLen$MalMeanCatchLen
    CombSexMeanCatchLen[k] <- Res$PerRecMeanCatchLen$CombSexMeanCatchLen
    FemMeanCatchLen.lw95ci[k] <- Res$MeanCatchLen.95ci$FemMeanCatchLen.95ci[1]   # 95 percent confidence limits
    FemMeanCatchLen.up95ci[k] <- Res$MeanCatchLen.95ci$FemMeanCatchLen.95ci[2]
    MalMeanCatchLen.lw95ci[k] <- Res$MeanCatchLen.95ci$MalMeanCatchLen.95ci[1]
    MalMeanCatchLen.up95ci[k] <- Res$MeanCatchLen.95ci$MalMeanCatchLen.95ci[2]
    CombSexMeanCatchLen.lw95ci[k] <- Res$MeanCatchLen.95ci$CombSexMeanCatchLen.95ci[1]
    CombSexMeanCatchLen.up95ci[k] <- Res$MeanCatchLen.95ci$CombSexMeanCatchLen.95ci[2]
    FemMeanCatchLen.lw60pi[k] <- Res$MeanCatchLen.60pi$FemMeanCatchLen.60pi[1] # 60 percent confidence prediction limits
    FemMeanCatchLen.up60pi[k] <- Res$MeanCatchLen.60pi$FemMeanCatchLen.60pi[2]
    MalMeanCatchLen.lw60pi[k] <- Res$MeanCatchLen.60pi$MalMeanCatchLen.60pi[1]
    MalMeanCatchLen.up60pi[k] <- Res$MeanCatchLen.60pi$MalMeanCatchLen.60pi[2]
    CombSexMeanCatchLen.lw60pi[k] <- Res$MeanCatchLen.60pi$CombSexMeanCatchLen.60pi[1]
    CombSexMeanCatchLen.up60pi[k] <- Res$MeanCatchLen.60pi$CombSexMeanCatchLen.60pi[2]
    FemMeanCatchLen.lw95pi[k] <- Res$MeanCatchLen.95pi$FemMeanCatchLen.95pi[1] # 95 percent confidence prediction limits
    FemMeanCatchLen.up95pi[k] <- Res$MeanCatchLen.95pi$FemMeanCatchLen.95pi[2]
    MalMeanCatchLen.lw95pi[k] <- Res$MeanCatchLen.95pi$MalMeanCatchLen.95pi[1]
    MalMeanCatchLen.up95pi[k] <- Res$MeanCatchLen.95pi$MalMeanCatchLen.95pi[2]
    CombSexMeanCatchLen.lw95pi[k] <- Res$MeanCatchLen.95pi$CombSexMeanCatchLen.95pi[1]
    CombSexMeanCatchLen.up95pi[k] <- Res$MeanCatchLen.95pi$CombSexMeanCatchLen.95pi[2]

    FemMeanCatchWt[k] <- Res$PerRecMeanCatchWt$FemMeanCatchWt # mean weight of catch
    MalMeanCatchWt[k] <- Res$PerRecMeanCatchWt$MalMeanCatchWt
    CombSexMeanCatchWt[k] <- Res$PerRecMeanCatchWt$CombSexMeanCatchWt
    FemMeanCatchWt.lw95ci[k] <- Res$MeanCatchWt.95ci$FemMeanCatchWt.95ci[1] # 95 percent confidence limits
    FemMeanCatchWt.up95ci[k] <- Res$MeanCatchWt.95ci$FemMeanCatchWt.95ci[2]
    MalMeanCatchWt.lw95ci[k] <- Res$MeanCatchWt.95ci$MalMeanCatchWt.95ci[1]
    MalMeanCatchWt.up95ci[k] <- Res$MeanCatchWt.95ci$MalMeanCatchWt.95ci[2]
    CombSexMeanCatchWt.lw95ci[k] <- Res$MeanCatchWt.95ci$CombSexMeanCatchWt.95ci[1]
    CombSexMeanCatchWt.up95ci[k] <- Res$MeanCatchWt.95ci$CombSexMeanCatchWt.95ci[2]
    FemMeanCatchWt.lw60pi[k] <- Res$MeanCatchWt.60pi$FemMeanCatchWt.60pi[1] # 60 percent confidence prediction limits
    FemMeanCatchWt.up60pi[k] <- Res$MeanCatchWt.60pi$FemMeanCatchWt.60pi[2]
    MalMeanCatchWt.lw60pi[k] <- Res$MeanCatchWt.60pi$MalMeanCatchWt.60pi[1]
    MalMeanCatchWt.up60pi[k] <- Res$MeanCatchWt.60pi$MalMeanCatchWt.60pi[2]
    CombSexMeanCatchWt.lw60pi[k] <- Res$MeanCatchWt.60pi$CombSexMeanCatchWt.60pi[1]
    CombSexMeanCatchWt.up60pi[k] <- Res$MeanCatchWt.60pi$CombSexMeanCatchWt.60pi[2]
    FemMeanCatchWt.lw95pi[k] <- Res$MeanCatchWt.95pi$FemMeanCatchWt.95pi[1] # 95 percent confidence prediction limits
    FemMeanCatchWt.up95pi[k] <- Res$MeanCatchWt.95pi$FemMeanCatchWt.95pi[2]
    MalMeanCatchWt.lw95pi[k] <- Res$MeanCatchWt.95pi$MalMeanCatchWt.95pi[1]
    MalMeanCatchWt.up95pi[k] <- Res$MeanCatchWt.95pi$FemMeanCatchWt.95pi[2]
    CombSexMeanCatchWt.lw95pi[k] <- Res$MeanCatchWt.95pi$CombSexMeanCatchWt.95pi[1]
    CombSexMeanCatchWt.up95pi[k] <- Res$MeanCatchWt.95pi$CombSexMeanCatchWt.95pi[2]
  }

  # Length results
  MeanCatchLenResults = data.frame(FishMort=FishMort,
                                   FemMeanCatchLen=FemMeanCatchLen,
                                   MalMeanCatchLen=MalMeanCatchLen,
                                   CombSexMeanCatchLen=CombSexMeanCatchLen,
                                   FemMeanCatchLen.lw95ci=FemMeanCatchLen.lw95ci,
                                   MalMeanCatchLen.lw95ci=MalMeanCatchLen.lw95ci,
                                   CombSexMeanCatchLen.lw95ci=CombSexMeanCatchLen.lw95ci,
                                   FemMeanCatchLen.up95ci=FemMeanCatchLen.up95ci,
                                   MalMeanCatchLen.up95ci=MalMeanCatchLen.up95ci,
                                   CombSexMeanCatchLen.up95ci=CombSexMeanCatchLen.up95ci,
                                   FemMeanCatchLen.lw60pi=FemMeanCatchLen.lw60pi,
                                   MalMeanCatchLen.lw60pi=MalMeanCatchLen.lw60pi,
                                   CombSexMeanCatchLen.lw60pi=CombSexMeanCatchLen.lw60pi,
                                   FemMeanCatchLen.up60pi=FemMeanCatchLen.up60pi,
                                   MalMeanCatchLen.up60pi=MalMeanCatchLen.up60pi,
                                   CombSexMeanCatchLen.up60pi=CombSexMeanCatchLen.up60pi,
                                   FemMeanCatchLen.lw95pi=FemMeanCatchLen.lw95pi,
                                   MalMeanCatchLen.lw95pi=MalMeanCatchLen.lw95pi,
                                   CombSexMeanCatchLen.lw95pi=CombSexMeanCatchLen.lw95pi,
                                   FemMeanCatchLen.up95pi=FemMeanCatchLen.up95pi,
                                   MalMeanCatchLen.up95pi=MalMeanCatchLen.up95pi,
                                   CombSexMeanCatchLen.up95pi=CombSexMeanCatchLen.up95pi)


  # weight results
  MeanCatchWtResults = data.frame(FishMort=FishMort,
                                  FemMeanCatchWt=FemMeanCatchWt,
                                  MalMeanCatchWt=MalMeanCatchWt,
                                  CombSexMeanCatchWt=CombSexMeanCatchWt,
                                  FemMeanCatchWt.lw95ci=FemMeanCatchWt.lw95ci,
                                  MalMeanCatchWt.lw95ci=MalMeanCatchWt.lw95ci,
                                  CombSexMeanCatchWt.lw95ci=CombSexMeanCatchWt.lw95ci,
                                  FemMeanCatchWt.up95ci=FemMeanCatchWt.up95ci,
                                  MalMeanCatchWt.up95ci=MalMeanCatchWt.up95ci,
                                  CombSexMeanCatchWt.up95ci=CombSexMeanCatchWt.up95ci,
                                  FemMeanCatchWt.lw60pi=FemMeanCatchWt.lw60pi,
                                  MalMeanCatchWt.lw60pi=MalMeanCatchWt.lw60pi,
                                  CombSexMeanCatchWt.lw60pi=CombSexMeanCatchWt.lw60pi,
                                  FemMeanCatchWt.up60pi=FemMeanCatchWt.up60pi,
                                  MalMeanCatchWt.up60pi=MalMeanCatchWt.up60pi,
                                  CombSexMeanCatchWt.up60pi=CombSexMeanCatchWt.up60pi,
                                  FemMeanCatchWt.lw95pi=FemMeanCatchWt.lw95pi,
                                  MalMeanCatchWt.lw95pi=MalMeanCatchWt.lw95pi,
                                  CombSexMeanCatchWt.lw95pi=CombSexMeanCatchWt.lw95pi,
                                  FemMeanCatchWt.up95pi=FemMeanCatchWt.up95pi,
                                  MalMeanCatchWt.up95pi=MalMeanCatchWt.up95pi,
                                  CombSexMeanCatchWt.up95pi=CombSexMeanCatchWt.up95pi)

  maxypr <- max(YPRResults) # maximum yield per recruit
  maxeqCatch <- max(EquilCatchResults) # maximum equilibrium catch
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
  # cat("2: FMort",FMort,'\n')
  Res2 = CalcYPRAndSPRForFMort_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                  RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
                                  EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
                                  FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, ret_Pmax,
                                  ret_L50, ret_L95, DiscMort, Steepness, SRrel_Type, NatMort, FMort)

  Results = list(PropFemAtLen = Res2$PropFemAtLen,
                 YPR = Res2$YPR,
                 Fem_SPR = Res2$Fem_SPR,
                 Mal_SPR = Res2$Mal_SPR,
                 CombSex_SPR = Res2$CombSex_SPR,
                 FemLenAtAge = as.vector(Res2$MeanSizeAtAge[1,]),
                 MalLenAtAge = as.vector(Res2$MeanSizeAtAge[2,]),
                 FemWtAtLen = Res2$FemWtAtLen,
                 MalWtAtLen = Res2$MalWtAtLen,
                 FemPropMatAtLen = Res2$FemPropMatAtLen,
                 MalPropMatAtLen = Res2$MalPropMatAtLen,
                 Unfish_FemNPerRec = Res2$Unfish_FemNPerRec,
                 Unfish_MalNPerRec = Res2$Unfish_MalNPerRec,
                 Fish_FemNPerRec = Res2$Fish_FemNPerRec,
                 Fish_MalNPerRec = Res2$Fish_MalNPerRec,
                 Unfish_FemBiomPerRecAtAge = Res2$Unfish_FemBiomPerRecAtAge,
                 Unfish_MalBiomPerRecAtAge = Res2$Unfish_MalBiomPerRecAtAge,
                 Fish_FemBiomPerRecAtAge = Res2$Fish_FemBiomPerRecAtAge,
                 Fish_MalBiomPerRecAtAge = Res2$Fish_MalBiomPerRecAtAge,
                 Unfish_FemBiomAtAge=Res2$Unfish_FemBiomAtAge,
                 Unfish_MalBiomAtAge=Res2$Unfish_MalBiomAtAge,
                 Fish_FemBiomAtAge=Res2$Fish_FemBiomAtAge,
                 Fish_MalBiomAtAge=Res2$Fish_MalBiomAtAge,
                 FemCatchBiom = Res2$FemCatchBiom,
                 MalCatchBiom = Res2$MalCatchBiom,
                 PerRecRandFishLen = Res2$PerRecRandFishLen,
                 PerRecLenFreq = Res2$PerRecLenFreq,
                 PerRecMeanCatchLen = Res2$PerRecMeanCatchLen,
                 MeanCatchLen.95ci = Res2$MeanCatchLen.95ci,
                 MeanCatchLen.60pi = Res2$MeanCatchLen.60pi,
                 MeanCatchLen.95pi = Res2$MeanCatchLen.95pi,
                 PerRecRandFishWt = Res2$PerRecRandFishWt,
                 PerRecMeanCatchWt = Res2$PerRecMeanCatchWt,
                 MeanCatchWt.95ci = Res2$MeanCatchWt.95ci,
                 MeanCatchWt.60pi = Res2$MeanCatchWt.60pi,
                 MeanCatchWt.95pi = Res2$MeanCatchWt.95pi,
                 MeanCatchLenResults = MeanCatchLenResults,
                 MeanCatchWtResults=MeanCatchWtResults,
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
                 FemGearSelAtLen = Res2$FemGearSelAtLen,
                 MalGearSelAtLen = Res2$MalGearSelAtLen,
                 FemRetProbAtLen = Res2$FemRetProbAtLen,
                 MalRetProbAtLen = Res2$MalRetProbAtLen,
                 FemSelLandAtLen = Res2$FemSelLandAtLen,
                 MalSelLandAtLen = Res2$MalSelLandAtLen,
                 FemSelDiscAtLen = Res2$FemSelDiscAtLen,
                 MalSelDiscAtLen = Res2$MalSelDiscAtLen,
                 FemDiscFAtLen = Res2$FemDiscFAtLen,
                 MalDiscFAtLen = Res2$MalDiscFAtLen,
                 FemLandFAtLen = Res2$FemLandFAtLen,
                 MalLandFAtLen = Res2$MalLandFAtLen,
                 FemFAtLen = Res2$FemFAtLen,
                 MalFAtLen = Res2$MalFAtLen,
                 FemZAtLen = Res2$FemZAtLen,
                 MalZAtLen = Res2$MalZAtLen,
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


#' Get plots associated with age-based per recruit analysis and extended analysis
#' with a stock-recruitment relationship
#'
#' This function provides a range of plots associated with age-based per recruit analysis, and an
#' extended form of this analysis with a Beverton-Holt stock recruitment relationship to account
#' for potential impacts of fishing on recruitment. Outputs include information on biology
#' inputs and various per recruit outputs
#'
#' @param MaxModelAge maximum age considered in model
#' @param TimeStep model timestep (e.g. 1 = annual, 1/12 = monthly)
#' @param Linf asymptotic length (von Bertalanffy growth parameter, set to NA if inputting lengths at ages directly from any growth model)
#' @param vbK growth coefficient (von Bertalanffy growth parameter, set to NA if inputting lengths at ages directly from any growth model)
#' @param tzero hypothetical age at zero length (von Bertalanffy growth parameter, set to NA if inputting lengths at ages directly from any growth model)
#' @param EstLenAtAge vector or estimated lengths at ages from any growth model, set to NA if von Bertalanffy growth parameters specified
#' @param lenwt_a weight-length parameter (power relationship)
#' @param ln_lenwt_a weight-length parameter (log-log relationship)
#' @param lenwt_b weight-length parameter (power or log-log relationship)
#' @param WLrel_Type 1=power, 2=log-log relationship (set to NA if inputting weights at ages directly)
#' @param EstWtAtAge vector of weights at ages (set to NA if weight-length growth parameters specified)
#' @param ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' @param ReprodPattern reproductive pattern, 1=gonochoristic, 2=protogynous (female to male sex change) hermaphroditism,
#' 3=protandrous hermaphroditism (male to female sex change)
#' @param InitRatioFem proportion of fish that are females at hatching
#' @param FinalSex_Pmax logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_A50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_A95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param mat_A50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_A95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtAge vector of proportion mature at age (set to NA if using age at maturity parameters)
#' @param sel_A50 logistic parameter for gear selectivity curve
#' @param sel_A95 logistic parameter for gear selectivity curve
#' @param EstSelAtAge vector of gear selectivity at age (set to NA if using age at gear selectivity parameters)
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_A50 logistic parameter for fish retention curve
#' @param ret_A95 logistic parameter for fish retention curve
#' @param EstRetenAtAge vector of fish retention at age (set to NA if using age at fish retention parameters)
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param RefPointPlotOpt plotting option for reference points, 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' @param Current_F estimate of current fishing mortality
#'
#' @return
#' plots including sex specific growth curves, weight at age curves, selectivity vs maturity at age curves,
#' total and natural mortality at age curves, fished vs unfished survival at age, fished vs unfished
#' female mature biomass at age, catch at age, at current estimated fishing mortality, yield per recruit
#' and equilibrium catch vs fishing mortality values, i.e. FishMort = seq(0,3,0.01), spawning potential
#' ratio and relative  spawning biomass vs FishMort, equilibrium recruitment vs FishMort
#'
#' @examples
#' # Example 1. Non-hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 20 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
#' Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.5, 0.5) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstMalLenAtAge=NA, EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- NA # for log-log relationship
#' lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA, EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' EstSelAtAge <- data.frame(EstFemSelAtAge=NA, EstMalSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
#' EstRetenAtAge <- data.frame(EstFemRetenAtAge=rep(1,MaxModelAge+1), EstMalRetenAtAge=rep(1,MaxModelAge+1)) # fish retention at age (from age 0), inputted as values in data frame
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 0.2 # natural mortality  (year-1)
#' Current_F <- 0.07 # estimate of fishing mortality, e.g. from catch curve analysis
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' PlotPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                       lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
#'                       ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
#'                       mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
#'                       EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F)
#' # # Example 2: hermaphroditic species
#' # InitRecruit <- 1 # Initial recruitment
#' # MaxModelAge <- 100 # maximum age considered by model, years
#' # TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be approprriate to use a smaller time step)
#' #  Linf <- c(1000, 1000) # mm - von Bertalanffy growth model parameters - Females, males
#' # vbK <- c(0.1, 0.1) # year-1 - von Bertalanffy growth model parameters - Females, males
#' # tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' # EstLenAtAge <- data.frame(EstMalLenAtAge=NA,
#' #                           EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' # lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' # ln_lenwt_a <- -11.0 # for log-log relationship
#' # lenwt_b <- 3.0 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' # WLrel_Type <- 2 # 1=power, 2=log-log relationship
#' # EstWtAtAge <- data.frame(EstFemWtAtAge=NA,
#' #                          EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' # ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' # ReprodPattern <- 2 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' # InitRatioFem <- 1 # Ratio of females to males at age zero
#' # FinalSex_Pmax <- 1 # Logistic sex change relationship parameters (max probability of final sex)
#' # FinalSex_A50 <- 35 # Logistic sex change relationship parameters (inflection point)
#' # FinalSex_A95 <- 60 # Logistic sex change relationship parameters (95% of max probability)
#' # mat_A50 <- c(20, 20) # females, males - Logistic length (mm) at maturity relationship parameters
#' # mat_A95 <- c(30, 30) # females, males - Logistic length (mm) at maturity relationship parameters
#' # EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' # sel_A50 <- c(15, 15) # females, males - Logistic age selectivity relationship parameters
#' # sel_A95 <- c(25, 25) # females, males - Logistic age selectivity relationship parameters
#' # EstSelAtAge <- data.frame(EstFemSelAtAge=NA, EstMalSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' # ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' # ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
#' # ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
#' # EstRetenAtAge <- data.frame(EstFemRetenAtAge=rep(1,MaxModelAge+1), EstMalRetenAtAge=rep(1,MaxModelAge+1)) # fish retention at age (from age 0), inputted as values in data frame
#' # DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' # Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' # SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' # NatMort = 0.07 # natural mortality  (year-1)
#' # Current_F <- 0.07 # estimate of fishing mortality, e.g. from catch curve analysis
#' # RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' # PlotPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#' #                       lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
#' #                       ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
#' #                       mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
#' #                       EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F)
#' @export
PlotPerRecruitResults_AB <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                  lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                                  ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                                  mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
                                  EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F) {

  .pardefault <- par(no.readonly = TRUE) # store current par settings

  Res = GetPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                             ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                             mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge,
                             ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness,
                             SRrel_Type, NatMort, Current_F)


  #Plot 1:
  par(mfrow = c(2,2), mar=c(3.5,4,2,2),
      oma=c(1,1,1,1),tck=-0.03,mgp = c(3, 0.5, 0))

  # plot growth curve
  y1=max(Res$FemLenAtAge)
  y2=max(Res$MalLenAtAge)
  if (y1 > y2) {
    ylims = Get_yaxis_scale(c(0,Res$FemLenAtAge))
  } else {
    ylims = Get_yaxis_scale(c(0,Res$MalLenAtAge))
  }
  ymax = ylims$ymax; yint = ylims$yint
  xlims = Get_xaxis_scale(Res$Ages)
  xmax = xlims$xmax; xint = xlims$xint
  plot(Res$Ages,Res$FemLenAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),
       col="red",yaxt="n",xaxt="n",ylab="",xlab="")
  lines(Res$Ages, Res$MalLenAtAge,col="blue")
  axis(1,at=seq(0,xmax,xint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,xmax,xint), labels=seq(0,xmax,xint),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
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
  xlims = Get_xaxis_scale(Res$Ages)
  xmax = xlims$xmax; xint = xlims$xint
  plot(Res$Ages,Res$FemWtAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),col="red",yaxt="n",xaxt="n",
       ylab="",xlab="")
  lines(Res$Ages,Res$MalWtAtAge,col="blue")
  axis(1,at=seq(0,xmax,xint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,xmax,xint), labels=seq(0,xmax,xint),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Total weight (kg"),plain(")"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('bottomright', col=c("red","blue"),legend=c("Female","Male"),bty='n', cex=0.8,lwd=1.75)

  # plot female maturity and selectivity at age
  xlims = Get_xaxis_scale(Res$Ages)
  xmax = xlims$xmax; xint = xlims$xint
  plot(Res$Ages,Res$FemPropMatAtAge,"l", pch=16,frame.plot=F,ylim=c(0,1),xlim=c(0,xmax),col="red",yaxt="n",xaxt="n",
       ylab="",xlab="",cex=0.8)
  if (DiscMort == 0) {
    lines(Res$Ages, Res$FemSelLandAtAge, "l", col="red",lty="dotted", cex=0.8)
    legend('bottomright', col=c("red","red"),legend=c("Fem. mature","Fem. reten."),
           lty=c("solid","dotted"),bty='n', cex=0.8,lwd=1.75)
  } else {
    lines(Res$Ages,Res$FemSelLandAtAge, "l", col="red",lty="dotted",cex=0.8)
    lines(Res$Ages,Res$FemSelDiscAtAge, "l", col="brown",lty="dotted",cex=0.8)
    legend('bottomright', col=c("red","red","brown"),legend=c("Fem. mature","Fem. land.","Fem. disc."),
           lty=c("solid","dotted","dotted"),bty='n', cex=0.8,lwd=1.75)
  }
  axis(1,at=seq(0,xmax,xint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,1,0.5), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,xmax,xint), labels=seq(0,xmax,xint),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,1,0.5), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Proportion"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  if (ReprodPattern > 1) {
    lines(Res$Ages,Res$PropFemAtAge, "l", col="black",lty="solid",cex=0.8)
    legend('topright', col="black",legend="Prop. Fem.",
           lty="solid",bty='n', cex=0.8,lwd=1.75)
  }

  # plot male maturity and selectivity at age
  xlims = Get_xaxis_scale(Res$Ages)
  xmax = xlims$xmax; xint = xlims$xint
  plot(Res$Ages, Res$MalPropMatAtAge,"l", pch=16, frame.plot=F,ylim=c(0,1),xlim=c(0,xmax),col="blue",yaxt="n",xaxt="n",
       ylab="",xlab="", cex=0.8)
  if (DiscMort == 0) {
    lines(Res$Ages, Res$FemSelLandAtAge, "l", col="blue",lty="dotted", cex=0.8)
    legend('bottomright', col=c("blue","blue"),legend=c("Male mature","Male reten."),
           lty=c("solid","dotted"),bty='n', cex=0.8,lwd=1.75)
  } else {
    lines(Res$Ages, Res$FemSelLandAtAge, "l", col="blue",lty="dotted", cex=0.8)
    lines(Res$Ages, Res$FemSelDiscAtAge, "l", col="purple",lty="dotted", cex=0.8)
    legend('bottomright', col=c("blue","blue","purple"),legend=c("Male mature","Male land.", "Male disc."),
           lty=c("solid","dotted","dotted"),bty='n', cex=0.8,lwd=1.75)
  }
  axis(1,at=seq(0,xmax,xint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,1,0.5), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,xmax,xint), labels=seq(0,xmax,xint),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
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
  xlims = Get_xaxis_scale(Res$Ages)
  xmax = xlims$xmax; xint = xlims$xint
  plot(Res$Ages, Res$FemFAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),col="red",yaxt="n",xaxt="n",
       ylab="",xlab="")
  # lines(Res$Ages,Res$FemZAtAge,lty="dotted",col="red")
  lines(Res$Ages,Res$FemDiscFAtAge,lty="dotted",col="brown")
  lines(Res$Ages,Res$FemLandFAtAge,lty="dotted",col="purple")
  lines(Res$Ages,rep(NatMort,length(Res$Ages)),lty="dashed")
  axis(1,at=seq(0,xmax,xint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,xmax,xint), labels=seq(0,xmax,xint),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Mortality") ~ (year^{-1}))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col=c("red","brown","purple","black"),lty=c("solid","dotted","dotted","dashed"),
         legend=c("Fem. F","Fem. DiscF","Fem. LandF","Fem. M"),bty='n', cex=1.0,lwd=1.75)

  # plot male mortality at age
  ylims = Get_yaxis_scale(Res$MalFAtAge)
  ymax = ylims$ymax; yint = ylims$yint
  xlims = Get_xaxis_scale(Res$Ages)
  xmax = xlims$xmax; xint = xlims$xint
  plot(Res$Ages, Res$MalFAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),col="blue",yaxt="n",xaxt="n",
       ylab="",xlab="")
  # lines(Res$Ages, Res$MalZAtAge,lty="dotted",col="blue")
  lines(Res$Ages,Res$MalDiscFAtAge,lty="dotted",col="brown")
  lines(Res$Ages,Res$MalLandFAtAge,lty="dotted",col="purple")
  lines(Res$Ages,rep(NatMort,length(Res$Ages)),lty="dashed")
  axis(1,at=seq(0,xmax,xint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,xmax,xint), labels=seq(0,xmax,xint),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Mortality") ~ (year^{-1}))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col=c("blue","brown","purple","black"),lty=c("solid","dotted","dotted","dashed"),
         legend=c("Mal. F","Mal. DiscF","Mal.LandF","Mal. M"),bty='n', cex=1.0,lwd=1.75)

  # plot fished and unfished female survival in terms of numbers given specified current fully-selected fishing mortality
  ylims = Get_yaxis_scale(Res$UnfishFemSurvAtAge)
  ymax = ylims$ymax; yint = ylims$yint
  if (ymax > 1) {
    ymax=1
    yint=0.2
  }
  xlims = Get_xaxis_scale(Res$Ages)
  xmax = xlims$xmax; xint = xlims$xint
  plot(Res$Ages, Res$UnfishFemSurvAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),col="red",yaxt="n",xaxt="n",
       ylab="",xlab="")
  lines(Res$Ages, Res$FishedFemSurvAtAge,col="red",lty="dotted")
  axis(1,at=seq(0,xmax,xint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,xmax,xint), labels=seq(0,xmax,xint),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Rel. survival"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col="red",lty=c("solid","dotted"),
         legend=c("Fem. unfish","Fem. fish"),bty='n', cex=1.0,lwd=1.75)

  # plot fished and unfished male survival in terms of numbers given specified current fully-selected fishing mortality
  ylims = Get_yaxis_scale(Res$UnfishMalSurvAtAge)
  ymax = ylims$ymax; yint = ylims$yint
  if (ymax > 1) {
    ymax=1
    yint=0.2
  }
  xlims = Get_xaxis_scale(Res$Ages)
  xmax = xlims$xmax; xint = xlims$xint
  plot(Res$Ages, Res$UnfishMalSurvAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),col="blue",yaxt="n",xaxt="n",
       ylab="",xlab="")
  lines(Res$Ages, Res$FishedMalSurvAtAge,col="blue",lty="dotted")
  axis(1,at=seq(0,xmax,xint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,xmax,xint), labels=seq(0,xmax,xint),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Rel. survival"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col="blue",lty=c("solid","dotted"),
         legend=c("Mal. unfish","Mal. fish"),bty='n', cex=1.0,lwd=1.75)

  # plot fished and unfished mature female biomass at age given specified current fully-selected fishing mortality
  ylims = Get_yaxis_scale(Res$UnfishFemBiomAtAge)
  ymax = ylims$ymax; yint = ylims$yint
  if (ymax == 0) {
    ymax = round(1.4 * max(Res$UnfishFemBiomAtAge),1)
    yint = ymax/4
  }
  xlims = Get_xaxis_scale(Res$Ages)
  xmax = xlims$xmax; xint = xlims$xint
  plot(Res$Ages, Res$UnfishFemBiomAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),
       col="red",yaxt="n",xaxt="n",ylab="",xlab="")
  lines(Res$Ages, Res$FishedFemBiomAtAge,col="red",lty="dotted")
  axis(1,at=seq(0,xmax,xint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,xmax,xint), labels=seq(0,xmax,xint),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Mat. biom. at age (kg"),plain(")"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col=c("red","red"),lty=c("solid","dotted"),legend=c("Fem. unfish","Fem unfish"),bty='n', cex=1.0,lwd=1.75)

  # plot fished and unfished mature male biomass at age given specified current fully-selected fishing mortality
  ylims = Get_yaxis_scale(Res$UnfishMalBiomAtAge)
  ymax = ylims$ymax; yint = ylims$yint
  xlims = Get_xaxis_scale(Res$Ages)
  xmax = xlims$xmax; xint = xlims$xint
  if (ymax == 0) {
    ymax = round(1.4 * max(Res$UnfishFemBiomAtAge),1)
    yint = ymax/4
  }
  plot(Res$Ages, Res$UnfishMalBiomAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),
       col="blue",yaxt="n",xaxt="n",ylab="",xlab="")
  lines(Res$Ages, Res$FishedMalBiomAtAge,col="blue",lty="dotted")
  axis(1,at=seq(0,xmax,xint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,xmax,xint), labels=seq(0,xmax,xint),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Mat. biom. at age (kg"),plain(")"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col=c("blue","blue"),lty=c("solid","dotted"),legend=c("Mal. unfish","Mal. unfish"),bty='n', cex=1.0,lwd=1.75)
  #Plot 3:

  # plot female and male catch at age, given specified current fully-selected fishing mortality
  par(mfrow = c(3,2), mar=c(3.5,4,2,2),
      oma=c(1,1,1,1),tck=-0.03,mgp = c(3, 0.5, 0))

  FemCatchNumAtAgeProp <- Res$FemCatchAtAgeNum / sum(Res$FemCatchAtAgeNum)
  MalCatchNumAtAgeProp <- Res$MalCatchAtAgeNum / sum(Res$MalCatchAtAgeNum)
  ylims = Get_yaxis_scale(FemCatchNumAtAgeProp)
  ymax = ylims$ymax; yint = ylims$yint
  xlims = Get_xaxis_scale(Res$Ages)
  xmax = xlims$xmax; xint = xlims$xint
  plot(Res$Ages, FemCatchNumAtAgeProp,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),
       col="red",yaxt="n",xaxt="n",ylab="",xlab="")
  lines(Res$Ages, MalCatchNumAtAgeProp,col="blue","l")
  axis(1,at=seq(0,xmax,xint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,xmax,xint), labels=seq(0,xmax,xint),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
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

  # reset default par options
  par(.pardefault)

}


#' Get plots associated with length-based per recruit analysis and extended analysis
#' with a stock-recruitment relationship
#'
#' This function provides a range of plots associated with length-based per recruit analysis, and an
#' extended form of this analysis with a Beverton-Holt stock recruitment relationship to account
#' for potential impacts of fishing on recruitment. Outputs include information on biology
#' inputs and various per recruit outputs
#'
#' @param MaxModelAge maximum age considered by model
#' @param TimeStep model time step (in y)
#' @param lbnd lower bounds of length classes
#' @param ubnd upper bounds of length classes
#' @param midpt mid points bounds of length classes
#' @param nLenCl number of length classes
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
#' @param GrowthParams growth parameters of either von Bertalanffy or Schnute model
#' @param RefnceAges reference ages for Schnute model, set to NA if using von Bertalanffy model
#' @param CVSizeAtAge coefficient of variation for size at age
#' @param lenwt_a weight-length parameter
#' @param ln_lenwt_a weight-length parameter
#' @param lenwt_b weight-length parameter
#' @param WLrel_Type 1=power, 2=log-log
#' @param EstWtAtLen user-specified weights at lengths
#' @param ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' @param ReprodPattern reproductive pattern, 1=gonochoristic, 2=protogynous (female to male sex change) hermaphroditism,
#' 3=protandrous hermaphroditism (male to female sex change)
#' @param InitRatioFem proportion of fish that are females at hatching
#' @param FinalSex_Pmax logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param mat_L50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_L95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtLen vector of proportion mature at length (set to NA if using age at maturity parameters)
#' @param sel_L50 logistic parameter for gear selectivity curve
#' @param sel_L95 logistic parameter for gear selectivity curve
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_L50 logistic parameter for fish retention curve
#' @param ret_L95 logistic parameter for fish retention curve
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param RefPointPlotOpt plotting option for reference points, 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' @param Current_F estimated current fishing mortality
#'
#' @examples
#' MaxModelAge <- 20 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
#' MaxLen = 800
#' LenInc = 20
#' lbnd = seq(0,MaxLen - LenInc, LenInc)
#' ubnd = lbnd + LenInc
#' midpt = lbnd + (LenInc/2)
#' nLenCl = length(midpt)
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.2, 0.2) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' GrowthParams = data.frame(Linf=Linf, vbK=vbK, tzero=tzero)
#' RefnceAges = NA
#' # GrowthParams c(Linf, vbK, tzero) single sex von Bertalanffy, or data.frame(Linf=Linf, vbK=vbK, tzero=tzero),
#' #' both sexes von Bertalanffy, or c(y1, y2, a, b) single sex Schnute, or data.frame(y1=y1, y2=y2, a=a, b=b), both sexes Schnute
#' CVSizeAtAge = c(0.05,0.05)
#' lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- NA # for log-log relationship
#' lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' EstWtAtLen <- data.frame(EstFemWtAtLen=NA,
#'                          EstMalWtAtLen=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_L50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_L95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' mat_L50 <- c(250, 250) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_L95 <- c(300, 300) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
#'                           EstMalMatAtLen=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_L50 <- c(250, 250) # females, males - Logistic age selectivity relationship parameters
#' sel_L95 <- c(300, 300) # females, males - Logistic age selectivity relationship parameters
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_L50 <- c(250, 250) # females, males - Logistic age fish retention at age parameters
#' ret_L95 <- c(300, 300) # females, males - Logistic age fish retention at age parameters
#' DiscMort <- 0.25 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 4.22 / MaxModelAge # natural mortality  (year-1)
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' Current_F = 0.2
#' PlotPerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
#'                                      RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
#'                                      EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
#'                                      FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, ret_Pmax,
#'                                      ret_L50, ret_L95, DiscMort, Steepness, SRrel_Type, NatMort, Current_F)
#' @export
PlotPerRecruitResults_LB <- function(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                     RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
                                     EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
                                     FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, ret_Pmax,
                                     ret_L50, ret_L95, DiscMort, Steepness, SRrel_Type, NatMort, Current_F) {

  .pardefault <- par(no.readonly = TRUE) # store current par settings

  Res = GetPerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
                                EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
                                FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, ret_Pmax,
                                ret_L50, ret_L95, DiscMort, Steepness, SRrel_Type, NatMort, Current_F)

  nTimeSteps <- 1 + (MaxModelAge / TimeStep)
  Ages <- seq(TimeStep,MaxModelAge,TimeStep)


  #Plot 1:
  par(mfrow = c(2,2), mar=c(3.5,4,2,2),
      oma=c(1,1,1,1),tck=-0.03,mgp = c(3, 0.5, 0))

  # plot growth curve
  names(Res)
  y1=max(Res$FemLenAtAge)
  y2=max(Res$MalLenAtAge)
  if (y1 > y2) {
    ylims = Get_yaxis_scale(c(0,Res$FemLenAtAge))
  } else {
    ylims = Get_yaxis_scale(c(0,Res$MalLenAtAge))
  }
  ymax = ylims$ymax; yint = ylims$yint
  xlims = Get_xaxis_scale(Ages)
  xmax = xlims$xmax; xint = xlims$xint
  plot(Ages,Res$FemLenAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),
       col="red",yaxt="n",xaxt="n",ylab="",xlab="")
  lines(Ages, Res$MalLenAtAge,col="blue")
  axis(1,at=seq(0,xmax,xint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,xmax,xint), labels=seq(0,xmax,xint),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Total length (mm"),plain(")"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('bottomright', col=c("red","blue"),legend=c("Female","Male"),bty='n', cex=0.8,lwd=1.75)

  # plot weight at length
  y1=max(Res$FemWtAtLen)
  y2=max(Res$MalWtAtLen)
  if (y1 > y2) {
    ylims = Get_yaxis_scale(c(0,Res$FemWtAtLen))
  } else {
    ylims = Get_yaxis_scale(c(0,Res$MalWtAtLen))
  }
  ymax = ylims$ymax; yint = ylims$yint
  xlims=Get_xaxis_scale(0:MaxLen)
  xmax = xlims$xmax; xint = xlims$xint

  plot(midpt,Res$FemWtAtLen,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),col="red",yaxt="n",xaxt="n",
       ylab="",xlab="")
  lines(midpt,Res$MalWtAtLen,col="blue")
  axis(1,at=seq(0,xmax,xint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,xmax,xint), labels=seq(0,xmax,xint),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Total weight (kg"),plain(")"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Length (mm"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topleft', col=c("red","blue"),legend=c("Female","Male"),bty='n', cex=0.8,lwd=1.75)

  # plot female maturity and selectivity at length
  xlims=Get_xaxis_scale(0:MaxLen)
  xmax = xlims$xmax; xint = xlims$xint
  plot(midpt,Res$FemPropMatAtLen,"l", pch=16,frame.plot=F,ylim=c(0,1),xlim=c(0,xmax),col="red",yaxt="n",xaxt="n",
       ylab="",xlab="",cex=0.8)
  if (DiscMort == 0) {
    lines(midpt, Res$FemSelLandAtLen, "l", col="red",lty="dotted", cex=0.8)
    legend('topleft', col=c("red","red"),legend=c("Fem. mature","Fem. reten."),
           lty=c("solid","dotted"),bty='n', cex=0.8,lwd=1.75)
  } else {
    lines(midpt,Res$FemSelLandAtLen, "l", col="red",lty="dotted",cex=0.8)
    lines(midpt,Res$FemSelDiscAtLen, "l", col="brown",lty="dotted",cex=0.8)
    legend('topleft', col=c("red","red","brown"),legend=c("Fem. mature","Fem. land.","Fem. disc."),
           lty=c("solid","dotted","dotted"),bty='n', cex=0.8,lwd=1.75)
  }
  axis(1,at=seq(0,xmax,xint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,1,0.5), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,xmax,xint), labels=seq(0,xmax,xint),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,1,0.5), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Proportion"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Length (mm"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  if (ReprodPattern > 1) {
    lines(midpt,Res$PropFemAtLen, "l", col="black",lty="solid",cex=0.8)
    legend('topright', col="black",legend="Prop. Fem.",
           lty="solid",bty='n', cex=0.8,lwd=1.75)
  }

  # plot male maturity and selectivity at length
  xlims=Get_xaxis_scale(0:MaxLen)
  xmax = xlims$xmax; xint = xlims$xint
  plot(midpt, Res$MalPropMatAtLen,"l", pch=16, frame.plot=F,ylim=c(0,1),xlim=c(0,xmax),col="blue",yaxt="n",xaxt="n",
       ylab="",xlab="", cex=0.8)
  if (DiscMort == 0) {
    lines(midpt, Res$FemSelLandAtLen, "l", col="blue",lty="dotted", cex=0.8)
    legend('topleft', col=c("blue","blue"),legend=c("Male mature","Male reten."),
           lty=c("solid","dotted"),bty='n', cex=0.8,lwd=1.75)
  } else {
    lines(midpt, Res$FemSelLandAtLen, "l", col="blue",lty="dotted", cex=0.8)
    lines(midpt, Res$FemSelDiscAtLen, "l", col="purple",lty="dotted", cex=0.8)
    legend('topleft', col=c("blue","blue","purple"),legend=c("Male mature","Male land.", "Male disc."),
           lty=c("solid","dotted","dotted"),bty='n', cex=0.8,lwd=1.75)
  }
  axis(1,at=seq(0,xmax,xint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,1,0.5), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,xmax,xint), labels=seq(0,xmax,xint),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,1,0.5), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Proportion"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  if (ReprodPattern > 1) {
    lines(midpt,1-Res$PropFemAtLen, "l", col="black",lty="solid",cex=0.8)
    legend('topright', col="black",legend="Prop. Male",
           lty="solid",bty='n', cex=0.8,lwd=1.75)
  }

  #plot 2:
  par(mfrow = c(3,2), mar=c(3.5,4,2,2),
      oma=c(1,1,1,1),tck=-0.03,mgp = c(3, 0.5, 0))

  # plot female mortality at length
  ylims = Get_yaxis_scale(Res$FemFAtLen)
  ymax = ylims$ymax; yint = ylims$yint
  xlims=Get_xaxis_scale(0:MaxLen)
  xmax = xlims$xmax; xint = xlims$xint
  plot(midpt, Res$FemFAtLen,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),col="red",yaxt="n",xaxt="n",
       ylab="",xlab="")
  # lines(Ages,Res$FemZAtAge,lty="dotted",col="red")
  lines(midpt,Res$FemDiscFAtLen,lty="dotted",col="brown")
  lines(midpt,Res$FemLandFAtLen,lty="dotted",col="purple")
  lines(midpt,rep(NatMort,length(midpt)),lty="dashed")
  axis(1,at=seq(0,xmax,xint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,xmax,xint), labels=seq(0,xmax,xint),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Mortality") ~ (year^{-1}))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Length (mm"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col=c("red","brown","purple","black"),lty=c("solid","dotted","dotted","dashed"),
         legend=c("Fem. F","Fem. DiscF","Fem. LandF","Fem. M"),bty='n', cex=1.0,lwd=1.75)

  # plot male mortality at age
  ylims = Get_yaxis_scale(Res$MalFAtLen)
  ymax = ylims$ymax; yint = ylims$yint
  xlims=Get_xaxis_scale(0:MaxLen)
  xmax = xlims$xmax; xint = xlims$xint
  plot(midpt, Res$MalFAtLen,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),col="blue",yaxt="n",xaxt="n",
       ylab="",xlab="")
  # lines(Ages, Res$MalZAtAge,lty="dotted",col="blue")
  lines(midpt,Res$MalDiscFAtLen,lty="dotted",col="brown")
  lines(midpt,Res$MalLandFAtLen,lty="dotted",col="purple")
  lines(midpt,rep(NatMort,length(midpt)),lty="dashed")
  axis(1,at=seq(0,xmax,xint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,xmax,xint), labels=seq(0,xmax,xint),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Mortality") ~ (year^{-1}))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Length (mm"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col=c("blue","brown","purple","black"),lty=c("solid","dotted","dotted","dashed"),
         legend=c("Mal. F","Mal. DiscF","Mal.LandF","Mal. M"),bty='n', cex=1.0,lwd=1.75)

  # plot fished and unfished female survival in terms of numbers given specified current fully-selected fishing mortality
  ylims = Get_yaxis_scale(Res$Unfish_FemNPerRec)
  ymax = ylims$ymax; yint = ylims$yint
  if (ymax <= 0) {
    ymax = round(max(1.2*Res$Unfish_FemNPerRec),1)
    yint = ymax
  }
  xlims=Get_xaxis_scale(0:MaxLen)
  xmax = xlims$xmax; xint = xlims$xint
  plot(midpt, Res$Unfish_FemNPerRec,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),col="red",yaxt="n",xaxt="n",
       ylab="",xlab="")
  lines(midpt, Res$Fish_FemNPerRec,col="red",lty="dotted")
  axis(1,at=seq(0,xmax,xint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,xmax,xint), labels=seq(0,xmax,xint),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Rel. survival"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Length (mm"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col="red",lty=c("solid","dotted"),
         legend=c("Fem. unfish","Fem. fish"),bty='n', cex=1.0,lwd=1.75)

  # plot fished and unfished male survival in terms of numbers given specified current fully-selected fishing mortality
  ylims = Get_yaxis_scale(Res$Unfish_MalNPerRec)
  ymax = ylims$ymax; yint = ylims$yint
  if (ymax <= 0) {
    ymax = round(max(1.2*Res$Unfish_FemNPerRec),1)
    yint = ymax
  }
  xlims=Get_xaxis_scale(0:MaxLen)
  xmax = xlims$xmax; xint = xlims$xint
  plot(midpt, Res$Unfish_MalNPerRec,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),col="blue",yaxt="n",xaxt="n",
       ylab="",xlab="")
  lines(midpt, Res$Fish_MalNPerRec,col="blue",lty="dotted")
  axis(1,at=seq(0,xmax,xint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,xmax,xint), labels=seq(0,xmax,xint),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Rel. survival"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Length (mm"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col="blue",lty=c("solid","dotted"),
         legend=c("Mal. unfish","Mal. fish"),bty='n', cex=1.0,lwd=1.75)


  # plot fished and unfished mature female biomass at age given specified current fully-selected fishing mortality
  ylims = Get_yaxis_scale(Res$Unfish_FemBiomAtAge)
  ymax = ylims$ymax; yint = ylims$yint
  if (max(Res$Unfish_FemBiomAtAge) < 0.2 * ymax) {
    ymax = ymax/5
    yint = yint/5
  }
  xlims=Get_xaxis_scale(c(0,Ages))
  xmax = xlims$xmax; xint = xlims$xint
  plot(c(0,Ages), Res$Unfish_FemBiomAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),
       col="red",yaxt="n",xaxt="n",ylab="",xlab="")
  lines(c(0,Ages), Res$Fish_FemBiomAtAge,col="red",lty="dotted")
  axis(1,at=seq(0,xmax,xint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,xmax,xint), labels=seq(0,xmax,xint),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Biom. at age (kg"),plain(")"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col=c("red","red"),lty=c("solid","dotted"),legend=c("Fem. unfish","Fem fish"),bty='n', cex=1.0,lwd=1.75)


  # plot fished and unfished mature male biomass at age given specified current fully-selected fishing mortality
  ylims = Get_yaxis_scale(Res$Unfish_MalBiomAtAge)
  ymax = ylims$ymax; yint = ylims$yint
  if (max(Res$Unfish_MalBiomAtAge) < 0.2 * ymax) {
    ymax = ymax/5
    yint = yint/5
  }
  xlims=Get_xaxis_scale(c(0,Ages))
  xmax = xlims$xmax; xint = xlims$xint
  plot(c(0,Ages), Res$Unfish_MalBiomAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),
       col="blue",yaxt="n",xaxt="n",ylab="",xlab="")
  lines(c(0,Ages), Res$Fish_MalBiomAtAge,col="blue",lty="dotted")
  axis(1,at=seq(0,xmax,xint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8, lwd=1.75,lab=F) # y axis
  axis(1,at=seq(0,xmax,xint), labels=seq(0,xmax,xint),cex.axis=0.8,line=0,las=1,lwd=1.5,tick=F) #add y labels
  axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
  mtext(expression(paste(plain("Biom. at age (kg"),plain(")"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
  mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  legend('topright', col=c("blue","blue"),lty=c("solid","dotted"),legend=c("Mal unfish","Mal. fish"),bty='n', cex=1.0,lwd=1.75)

  #Plot 3:

  # plot female and male catch at age, given specified current fully-selected fishing mortality
  par(mfrow = c(2,2), mar=c(3.5,4,2,2),
      oma=c(1,1,1,1),tck=-0.03,mgp = c(3, 0.5, 0))


  # plot yield per recruit (per recruit analysis) and equilibrium catch (equilibrium age-structured model)
  # given specified current fully-selected fishing mortality
  ylims = Get_yaxis_scale(Res$YPRResults)
  ymax = ylims$ymax; yint = ylims$yint
  if (max(Res$YPRResults) < 0.2 * ymax) {
    ymax = ymax/5
    yint = yint/5
  }
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

  # plot equilibrium recruitment vs F
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

  # reset default par options
  par(.pardefault)

}

#' Deterministic plot of equilibrium yield vs F from extended age-based per recruit analysis
#' with a stock-recruitment relationship
#'
#' This function provides a deterministic plot of equilibrium yield vs F from extended age-based per recruit analysis
#' with a stock-recruitment relationship
#'
#' @param MaxModelAge maximum age considered in model
#' @param TimeStep model timestep (e.g. 1 = annual, 1/12 = monthly)
#' @param Linf asymptotic length (von Bertalanffy growth parameter, set to NA if inputting lengths at ages directly from any growth model)
#' @param vbK growth coefficient (von Bertalanffy growth parameter, set to NA if inputting lengths at ages directly from any growth model)
#' @param tzero hypothetical age at zero length (von Bertalanffy growth parameter, set to NA if inputting lengths at ages directly from any growth model)
#' @param EstLenAtAge vector or estimated lengths at ages from any growth model, set to NA if von Bertalanffy growth parameters specified
#' @param lenwt_a weight-length parameter (power relationship)
#' @param ln_lenwt_a weight-length parameter (log-log relationship)
#' @param lenwt_b weight-length parameter (power or log-log relationship)
#' @param WLrel_Type 1=power, 2=log-log relationship (set to NA if inputting weights at ages directly)
#' @param EstWtAtAge vector of weights at ages (set to NA if weight-length growth parameters specified)
#' @param ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' @param ReprodPattern reproductive pattern, 1=gonochoristic, 2=protogynous (female to male sex change) hermaphroditism,
#' 3=protandrous hermaphroditism (male to female sex change)
#' @param InitRatioFem proportion of fish that are females at hatching
#' @param FinalSex_Pmax logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_A50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_A95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param mat_A50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_A95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtAge vector of proportion mature at age (set to NA if using age at maturity parameters)
#' @param sel_A50 logistic parameter for gear selectivity curve
#' @param sel_A95 logistic parameter for gear selectivity curve
#' @param EstSelAtAge vector of gear selectivity at age (set to NA if using age at gear selectivity parameters)
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_A50 logistic parameter for fish retention curve
#' @param ret_A95 logistic parameter for fish retention curve
#' @param EstRetenAtAge vector of fish retention at age (set to NA if using age at fish retention parameters)
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param RefPointPlotOpt plotting option for reference points, 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' @param Current_F estimate of current fishing mortality
#' @param xaxis_lab y axis label
#' @param yaxis_lab x axis label
#' @param xmax maximum x axis value
#' @param xint x axis interval
#' @param ymax maximum y axis value
#' @param yint y axis interval
#'
#' @return
#' deterministic plot of yield vs F from per recruit analysis and extended analysis
#' with a stock-recruitment relationship.
#' @examples
#' # Example 1. Non-hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 20 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
#' Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.5, 0.5) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstMalLenAtAge=NA, EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- NA # for log-log relationship
#' lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA, EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' EstSelAtAge <- data.frame(EstFemSelAtAge=NA, EstMalSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
#' EstRetenAtAge <- data.frame(EstFemRetenAtAge=rep(1,MaxModelAge+1), EstMalRetenAtAge=rep(1,MaxModelAge+1)) # fish retention at age (from age 0), inputted as values in data frame
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 0.2 # natural mortality  (year-1)
#' Current_F <- 0.07 # estimate of fishing mortality, e.g. from catch curve analysis
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' PlotPerRecruit_EqYield_no_err_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                      lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
#'                      ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
#'                      mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
#'                      EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F,
#'                      MainLabel=NA, xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA, ymax=NA, yint=NA)
#' # # Example 2: hermaphroditic species
#' # InitRecruit <- 1 # Initial recruitment
#' # MaxModelAge <- 100 # maximum age considered by model, years
#' # TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be approprriate to use a smaller time step)
#' # Linf <- c(1000, 1000) # mm - von Bertalanffy growth model parameters - Females, males
#' # vbK <- c(0.1, 0.1) # year-1 - von Bertalanffy growth model parameters - Females, males
#' # tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' # EstLenAtAge <- data.frame(EstMalLenAtAge=NA, EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' # lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' # ln_lenwt_a <- -11.0 # for log-log relationship
#' # lenwt_b <- 3.0 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' # WLrel_Type <- 2 # 1=power, 2=log-log relationship
#' # EstWtAtAge <- data.frame(EstFemWtAtAge=NA, EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' # ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' # ReprodPattern <- 2 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' # InitRatioFem <- 1 # Ratio of females to males at age zero
#' # FinalSex_Pmax <- 1 # Logistic sex change relationship parameters (max probability of final sex)
#' # FinalSex_A50 <- 35 # Logistic sex change relationship parameters (inflection point)
#' # FinalSex_A95 <- 60 # Logistic sex change relationship parameters (95% of max probability)
#' # mat_A50 <- c(20, 20) # females, males - Logistic length (mm) at maturity relationship parameters
#' # mat_A95 <- c(30, 30) # females, males - Logistic length (mm) at maturity relationship parameters
#' # EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' # sel_A50 <- c(15, 15) # females, males - Logistic age selectivity relationship parameters
#' # sel_A95 <- c(25, 25) # females, males - Logistic age selectivity relationship parameters
#' # EstSelAtAge <- data.frame(EstFemSelAtAge=NA, EstMalSelAtAge=NA) # fish retention at age (from age 0), inputted as values in data frame
#' # ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' # ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
#' # ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
#' # EstRetenAtAge <- data.frame(EstFemRetenAtAge=rep(1,MaxModelAge+1), EstMalRetenAtAge=rep(1,MaxModelAge+1)) # fish retention at age (from age 0), inputted as values in data frame
#' # DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' # Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' # SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' # NatMort = 0.07 # natural mortality  (year-1)
#' # Current_F <- 0.07 # estimate of fishing mortality, e.g. from catch curve analysis
#' # RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' # PlotPerRecruit_EqYield_no_err_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#' #                       lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
#' #                       ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
#' #                       mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
#' #                       EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F,
#' #                       MainLabel=NA, xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA, ymax=NA, yint=NA)
#' @export
PlotPerRecruit_EqYield_no_err_AB <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                                             ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                                             mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
                                             EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F,
                                             MainLabel, xaxis_lab, yaxis_lab, xmax, xint, ymax, yint) {

  Res = GetPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                                ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                                mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge,
                                ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness,
                                SRrel_Type, NatMort, Current_F)

  if (is.na(yaxis_lab)) yaxis_lab = "Catch (kg)"
  if (is.na(xaxis_lab)) xaxis_lab = expression(paste(italic("F") ~ (year^{-1})))
  if (is.na(xmax)) xmax = max(Res$FishMort)
  if (is.na(xint)) xint = 0.5
  if (is.na(ymax)) {
    ylims = Get_yaxis_scale(Res$YPRResults)
    ymax = ylims$ymax; yint = ylims$yint
    if (max(Res$YPRResults) < 0.2 * ymax) {
      ymax = ymax/5
      yint = yint/5
    }
  }

  # F vs YPR and equil catch
  plot(Res$FishMort, Res$EquilCatchResults, "l", frame.plot=F, ylim=c(0, ymax), xlim=c(0, xmax),
       col="black", yaxt="n", xaxt="n", ylab="", xlab="")
  lines(Res$FishMort, Res$YPRResults, lty="dotted")
  points(Current_F, Res$Equil_Catch, cex=1.2, col="black", pch=16)
  points(Current_F, Res$YPR, cex=1.2, col="black", pch=1)
  axis(1, at=seq(0, xmax, xint), cex.axis=1, lwd=1.75, lab=F, line=-0.3)
  axis(2, at=seq(0, ymax, yint), cex.axis=1, lwd=1.75, lab=F, line=-0.3)
  axis(1, at=seq(0, xmax, xint), labels=seq(0, xmax, xint),
       cex.axis=1, line=-0.5, las=1, lwd=1.5, tick=F)
  axis(2, at=seq(0, ymax, yint), cex.axis=1, line=-0.5, las=1, lwd=1.5, tick=F)
  mtext(yaxis_lab, las=3, side=2, line=3, cex=1.2, lwd=1.75)
  mtext(xaxis_lab, las=1, side=1, line=3, cex=1.2, lwd=1.75)
  legend("topright", col=c("black", "black"), lty=c("solid", "dotted"), legend=c("Eq. catch","YPR"),
         bty="n", cex=0.8, lwd=1.75, inset = 0.05)
  legend("topleft", col="black", pch = c(16,1), lty=0, legend=c("Est_EqYield","Est_YPR"),
         bty="n", cex=0.8, inset = 0.05)
}


#' Deterministic plot of equilibrium yield vs F from extended age-based per recruit analysis
#' with a stock-recruitment relationship
#'
#' This function provides a deterministic plot of equilibrium yield vs F from extended age-based per recruit analysis
#' with a stock-recruitment relationship
#'
#' Get outputs from length-based per recruit analysis across a range of fishing mortality values
#'
#' This function provides outputs associated with per recruit analysis, and an
#' extended form of this analysis with a Beverton-Holt stock recruitment relationship to account
#' for potential impacts of fishing on recruitment. Outputs are provided for a range of
#' fishing mortality values, including the current, estimated value.
#'
#' @param MaxModelAge maximum age considered by model
#' @param TimeStep model time step (in y)
#' @param lbnd lower bounds of length classes
#' @param ubnd upper bounds of length classes
#' @param midpt mid points bounds of length classes
#' @param nLenCl number of length classes
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
#' @param GrowthParams growth parameters of either von Bertalanffy or Schnute model
#' @param RefnceAges reference ages for Schnute model, set to NA if using von Bertalanffy model
#' @param CVSizeAtAge coefficient of variation for size at age
#' @param lenwt_a weight-length parameter
#' @param ln_lenwt_a weight-length parameter
#' @param lenwt_b weight-length parameter
#' @param WLrel_Type 1=power, 2=log-log
#' @param EstWtAtLen user-specified weights at lengths
#' @param ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' @param ReprodPattern reproductive pattern, 1=gonochoristic, 2=protogynous (female to male sex change) hermaphroditism,
#' 3=protandrous hermaphroditism (male to female sex change)
#' @param InitRatioFem proportion of fish that are females at hatching
#' @param FinalSex_Pmax logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param mat_L50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_L95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtLen vector of proportion mature at length (set to NA if using age at maturity parameters)
#' @param sel_L50 logistic parameter for gear selectivity curve
#' @param sel_L95 logistic parameter for gear selectivity curve
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_L50 logistic parameter for fish retention curve
#' @param ret_L95 logistic parameter for fish retention curve
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param RefPointPlotOpt plotting option for reference points, 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' @param Current_F estimated current fishing mortality
#' @param xaxis_lab y axis label
#' @param yaxis_lab x axis label
#' @param xmax maximum x axis value
#' @param xint x axis interval
#' @param ymax maximum y axis value
#' @param yint y axis interval
#'
#' @examples
#' MaxModelAge <- 20 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
#' MaxLen = 800
#' LenInc = 20
#' lbnd = seq(0,MaxLen - LenInc, LenInc)
#' ubnd = lbnd + LenInc
#' midpt = lbnd + (LenInc/2)
#' nLenCl = length(midpt)
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.2, 0.2) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' GrowthParams = data.frame(Linf=Linf, vbK=vbK, tzero=tzero)
#' RefnceAges = NA
#' # GrowthParams c(Linf, vbK, tzero) single sex von Bertalanffy, or data.frame(Linf=Linf, vbK=vbK, tzero=tzero),
#' #' both sexes von Bertalanffy, or c(y1, y2, a, b) single sex Schnute, or data.frame(y1=y1, y2=y2, a=a, b=b), both sexes Schnute
#' CVSizeAtAge = c(0.05,0.05)
#' lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- NA # for log-log relationship
#' lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' EstWtAtLen <- data.frame(EstFemWtAtLen=NA,
#'                          EstMalWtAtLen=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_L50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_L95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' mat_L50 <- c(250, 250) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_L95 <- c(300, 300) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
#'                           EstMalMatAtLen=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_L50 <- c(250, 250) # females, males - Logistic age selectivity relationship parameters
#' sel_L95 <- c(300, 300) # females, males - Logistic age selectivity relationship parameters
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_L50 <- c(250, 250) # females, males - Logistic age fish retention at age parameters
#' ret_L95 <- c(300, 300) # females, males - Logistic age fish retention at age parameters
#' DiscMort <- 0.25 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 4.22 / MaxModelAge # natural mortality  (year-1)
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' Current_F = 0.2
#' PlotPerRecruit_EqYield_no_err_LB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                                              lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
#'                                              ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
#'                                              mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
#'                                              EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F,
#'                                              MainLabel=NA, xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA, ymax=NA, yint=NA)
#' @export
PlotPerRecruit_EqYield_no_err_LB <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                                             ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                                             mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
                                             EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F,
                                             MainLabel, xaxis_lab, yaxis_lab, xmax, xint, ymax, yint) {

Res=GetPerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                            RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
                            EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
                            FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, ret_Pmax,
                            ret_L50, ret_L95, DiscMort, Steepness, SRrel_Type, NatMort, Current_F)

if (is.na(yaxis_lab)) yaxis_lab = "Catch (kg)"
if (is.na(xaxis_lab)) xaxis_lab = expression(paste(italic("F") ~ (year^{-1})))
if (is.na(xmax)) xmax = max(Res$FishMort)
if (is.na(xint)) xint = 0.5
if (is.na(ymax)) {
  ylims = Get_yaxis_scale(Res$Unfish_FemBiomAtAge)
  ymax = ylims$ymax; yint = ylims$yint
  if (max(Res$Unfish_FemBiomAtAge) < 0.2 * ymax) {
    ymax = ymax/5
    yint = yint/5
  }
}

# F vs YPR and equil catch
plot(Res$FishMort, Res$EquilCatchResults, "l", frame.plot=F, ylim=c(0, ymax), xlim=c(0, xmax),
     col="black", yaxt="n", xaxt="n", ylab="", xlab="")
lines(Res$FishMort, Res$YPRResults, lty="dotted")
points(Current_F, Res$Equil_Catch, cex=1.2, col="black", pch=16)
points(Current_F, Res$YPR, cex=1.2, col="black", pch=1)
axis(1, at=seq(0, xmax, xint), cex.axis=1, lwd=1.75, lab=F, line = -0.3)
axis(2, at=seq(0, ymax, yint), cex.axis=1, lwd=1.75, lab=F, line = -0.3)
axis(1, at=seq(0, xmax, xint), labels=seq(0, xmax, xint),
     cex.axis=1, line=-0.5, las=1, lwd=1.5, tick=F)
axis(2, at=seq(0, ymax, yint), cex.axis=1, line=-0.5, las=1, lwd=1.5, tick=F)
mtext(yaxis_lab, las=3, side=2, line=3, cex=1.2, lwd=1.75)
mtext(yaxis_lab, las=1, side=1, line=3, cex=1.2, lwd=1.75)
legend("topright", col=c("black", "black"), lty=c("solid", "dotted"), legend=c("Eq. catch","YPR"),
       bty="n", cex=0.8, lwd=1.75, inset = 0.05)
legend("topleft", col="black", pch = c(16,1), lty=0, legend=c("Est_EqYield","Est_YPR"),
       bty="n", cex=0.8, inset = 0.05)
}

#' Deterministic plot of SPR and relative equilibrium biomass vs F from age-based per recruit analysis and extended analysis
#' with a stock-recruitment relationship
#'
#' This function provides a deterministic plot of SPR and relative equilibrium biomass  vs F from age-based per recruit analysis and extended analysis
#' with a stock-recruitment relationship
#'
#' @param MaxModelAge maximum age considered in model
#' @param TimeStep model timestep (e.g. 1 = annual, 1/12 = monthly)
#' @param Linf asymptotic length (von Bertalanffy growth parameter, set to NA if inputting lengths at ages directly from any growth model)
#' @param vbK growth coefficient (von Bertalanffy growth parameter, set to NA if inputting lengths at ages directly from any growth model)
#' @param tzero hypothetical age at zero length (von Bertalanffy growth parameter, set to NA if inputting lengths at ages directly from any growth model)
#' @param EstLenAtAge vector or estimated lengths at ages from any growth model, set to NA if von Bertalanffy growth parameters specified
#' @param lenwt_a weight-length parameter (power relationship)
#' @param ln_lenwt_a weight-length parameter (log-log relationship)
#' @param lenwt_b weight-length parameter (power or log-log relationship)
#' @param WLrel_Type 1=power, 2=log-log relationship (set to NA if inputting weights at ages directly)
#' @param EstWtAtAge vector of weights at ages (set to NA if weight-length growth parameters specified)
#' @param ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' @param ReprodPattern reproductive pattern, 1=gonochoristic, 2=protogynous (female to male sex change) hermaphroditism,
#' 3=protandrous hermaphroditism (male to female sex change)
#' @param InitRatioFem proportion of fish that are females at hatching
#' @param FinalSex_Pmax logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_A50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_A95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param mat_A50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_A95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtAge vector of proportion mature at age (set to NA if using age at maturity parameters)
#' @param sel_A50 logistic parameter for gear selectivity curve
#' @param sel_A95 logistic parameter for gear selectivity curve
#' @param EstSelAtAge vector of gear selectivity at age (set to NA if using age at gear selectivity parameters)
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_A50 logistic parameter for fish retention curve
#' @param ret_A95 logistic parameter for fish retention curve
#' @param EstRetenAtAge vector of fish retention at age (set to NA if using age at fish retention parameters)
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param RefPointPlotOpt plotting option for reference points, 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' @param Current_F estimate of current fishing mortality
#'
#' @return
#' deterministic plot of SPR and relative equilibrium biomass  vs F from per recruit analysis and extended analysis
#' with a stock-recruitment relationship
#'
#' @examples
#' # Example 1. Non-hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 20 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
#' Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.5, 0.5) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstMalLenAtAge=NA, EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- NA # for log-log relationship
#' lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA, EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' EstSelAtAge <- data.frame(EstFemSelAtAge=NA, EstMalSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
#' EstRetenAtAge <- data.frame(EstFemRetenAtAge=rep(1,MaxModelAge+1), EstMalRetenAtAge=rep(1,MaxModelAge+1)) # fish retention at age (from age 0), inputted as values in data frame
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 0.2 # natural mortality  (year-1)
#' Current_F <- 0.2 # estimate of fishing mortality, e.g. from catch curve analysis
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' PlotPerRecruit_Biom_no_err_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                       lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
#'                       ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
#'                       mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
#'                       EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F)
#' # # Example 2: hermaphroditic species
#' # InitRecruit <- 1 # Initial recruitment
#' # MaxModelAge <- 100 # maximum age considered by model, years
#' # TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
#' # Linf <- c(1000, 1000) # mm - von Bertalanffy growth model parameters - Females, males
#' # vbK <- c(0.1, 0.1) # year-1 - von Bertalanffy growth model parameters - Females, males
#' # tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' # EstLenAtAge <- data.frame(EstMalLenAtAge=NA, EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' # lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' # ln_lenwt_a <- -11.0 # for log-log relationship
#' # lenwt_b <- 3.0 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' # WLrel_Type <- 2 # 1=power, 2=log-log relationship
#' # EstWtAtAge <- data.frame(EstFemWtAtAge=NA, EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' # ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' # ReprodPattern <- 2 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' # InitRatioFem <- 1 # Ratio of females to males at age zero
#' # FinalSex_Pmax <- 1 # Logistic sex change relationship parameters (max probability of final sex)
#' # FinalSex_A50 <- 35 # Logistic sex change relationship parameters (inflection point)
#' # FinalSex_A95 <- 60 # Logistic sex change relationship parameters (95% of max probability)
#' # mat_A50 <- c(20, 20) # females, males - Logistic length (mm) at maturity relationship parameters
#' # mat_A95 <- c(30, 30) # females, males - Logistic length (mm) at maturity relationship parameters
#' # EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' # sel_A50 <- c(15, 15) # females, males - Logistic age selectivity relationship parameters
#' # sel_A95 <- c(25, 25) # females, males - Logistic age selectivity relationship parameters
#' # EstSelAtAge <- data.frame(EstFemSelAtAge=NA, EstMalSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' # ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' # ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
#' # ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
#' # EstRetenAtAge <- data.frame(EstFemRetenAtAge=rep(1,MaxModelAge+1), EstMalRetenAtAge=rep(1,MaxModelAge+1)) # fish retention at age (from age 0), inputted as values in data frame
#' # DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' # Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' # SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' # NatMort = 0.07 # natural mortality  (year-1)
#' # Current_F <- 0.07 # estimate of fishing mortality, e.g. from catch curve analysis
#' # RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' # PlotPerRecruit_Biom_no_err_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#' #                       lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
#' #                       ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
#' #                       mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
#' #                       EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F)
#' @export
PlotPerRecruit_Biom_no_err_AB <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                       lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                                       ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                                       mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
                                       EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F) {

  Res = GetPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                             ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                             mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge,
                             ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness,
                             SRrel_Type, NatMort, Current_F)

  # F vs SPR and Brel
  xmax = max(Res$FishMort)
  if (NatMort < 0.15) xmax = 1.0
  plot(Res$FishMort, Res$Equilmod_FemRelBiomResults, "l", frame.plot=F, ylim=c(0, 1), xlim=c(0, xmax),
       col="black", yaxt="n", xaxt="n", ylab="", xlab="")
  points(Current_F, Res$Equilmod_FemRelBiom, cex=1.2, col="black", pch=16)
  axis(1, at=seq(0, xmax, 0.5), cex.axis=1, lwd=1.75, lab=F)
  axis(2, at=seq(0, 1, 0.2), cex.axis=1, lwd=1.75, lab=F)
  axis(1, at=seq(0, xmax, 0.5), labels = seq(0, xmax, 0.5),
       cex.axis=1, line=0.5, las=1, lwd=1.5, tick=F)
  axis(2, at=seq(0, 1, 0.2), cex.axis=1, line=0.5, las=1, lwd=1.5, tick=F)
  mtext(expression(paste(plain("Relative spawning biomass"))), las=3, side=2, line=3, cex=1, lwd=1.75)
  mtext(expression(paste(italic("F") ~ (year^{-1}))), las=1, side=1, line=3, cex=1, lwd=1.75)
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


#' Deterministic plot of SPR and relative equilibrium biomass vs F from length-based per recruit analysis and extended analysis
#' with a stock-recruitment relationship
#'
#' This function provides a deterministic plot of SPR and relative equilibrium biomass  vs F from length-based per recruit analysis and extended analysis
#' with a stock-recruitment relationship
#'
#' @param MaxModelAge maximum age considered by model
#' @param TimeStep model time step (in y)
#' @param lbnd lower bounds of length classes
#' @param ubnd upper bounds of length classes
#' @param midpt mid points bounds of length classes
#' @param nLenCl number of length classes
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
#' @param GrowthParams growth parameters of either von Bertalanffy or Schnute model
#' @param RefnceAges reference ages for Schnute model, set to NA if using von Bertalanffy model
#' @param CVSizeAtAge coefficient of variation for size at age
#' @param lenwt_a weight-length parameter
#' @param ln_lenwt_a weight-length parameter
#' @param lenwt_b weight-length parameter
#' @param WLrel_Type 1=power, 2=log-log
#' @param EstWtAtLen user-specified weights at lengths
#' @param ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' @param ReprodPattern reproductive pattern, 1=gonochoristic, 2=protogynous (female to male sex change) hermaphroditism,
#' 3=protandrous hermaphroditism (male to female sex change)
#' @param InitRatioFem proportion of fish that are females at hatching
#' @param FinalSex_Pmax logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param mat_L50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_L95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtLen vector of proportion mature at length (set to NA if using age at maturity parameters)
#' @param sel_L50 logistic parameter for gear selectivity curve
#' @param sel_L95 logistic parameter for gear selectivity curve
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_L50 logistic parameter for fish retention curve
#' @param ret_L95 logistic parameter for fish retention curve
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param RefPointPlotOpt plotting option for reference points, 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' @param Current_F estimated current fishing mortality
#'
#' @examples
#' MaxModelAge <- 20 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
#' MaxLen = 800
#' LenInc = 20
#' lbnd = seq(0,MaxLen - LenInc, LenInc)
#' ubnd = lbnd + LenInc
#' midpt = lbnd + (LenInc/2)
#' nLenCl = length(midpt)
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.2, 0.2) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' GrowthParams = data.frame(Linf=Linf, vbK=vbK, tzero=tzero)
#' RefnceAges = NA
#' # GrowthParams c(Linf, vbK, tzero) single sex von Bertalanffy, or data.frame(Linf=Linf, vbK=vbK, tzero=tzero),
#' #' both sexes von Bertalanffy, or c(y1, y2, a, b) single sex Schnute, or data.frame(y1=y1, y2=y2, a=a, b=b), both sexes Schnute
#' CVSizeAtAge = c(0.05, 0.05)
#' lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- NA # for log-log relationship
#' lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' EstWtAtLen <- data.frame(EstFemWtAtLen=NA,
#'                          EstMalWtAtLen=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_L50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_L95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' mat_L50 <- c(250, 250) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_L95 <- c(300, 300) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
#'                           EstMalMatAtLen=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_L50 <- c(250, 250) # females, males - Logistic age selectivity relationship parameters
#' sel_L95 <- c(300, 300) # females, males - Logistic age selectivity relationship parameters
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_L50 <- c(250, 250) # females, males - Logistic age fish retention at age parameters
#' ret_L95 <- c(300, 300) # females, males - Logistic age fish retention at age parameters
#' DiscMort <- 0.25 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 4.22 / MaxModelAge # natural mortality  (year-1)
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' Current_F = 0.2
#' PlotPerRecruit_Biom_no_err_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
#'                                           RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
#'                                           EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
#'                                           FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, ret_Pmax,
#'                                           ret_L50, ret_L95, DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F)
#' @export
PlotPerRecruit_Biom_no_err_LB <- function(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                          RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
                                          EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
                                          FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, ret_Pmax,
                                          ret_L50, ret_L95, DiscMort, Steepness, SRrel_Type, NatMort, RefPointPlotOpt, Current_F) {

  Res = GetPerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
                                EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
                                FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, ret_Pmax,
                                ret_L50, ret_L95, DiscMort, Steepness, SRrel_Type, NatMort, Current_F)

  # F vs SPR and Brel
  xmax = max(Res$FishMort)
  if (NatMort < 0.15) xmax = 1.0
  plot(Res$FishMort, Res$Equilmod_FemRelBiomResults, "l", frame.plot=F, ylim=c(0, 1), xlim=c(0, xmax),
       col="black", yaxt="n", xaxt="n", ylab="", xlab="")
  points(Current_F, Res$Equilmod_FemRelBiom, cex=1.2, col="black", pch=16)
  axis(1, at=seq(0, xmax, 0.5), cex.axis=1, lwd=1.75, lab=F)
  axis(2, at=seq(0, 1, 0.2), cex.axis=1, lwd=1.75, lab=F)
  axis(1, at=seq(0, xmax, 0.5), labels = seq(0, xmax, 0.5),
       cex.axis=1, line=0.5, las=1, lwd=1.5, tick=F)
  axis(2, at=seq(0, 1, 0.2), cex.axis=1, line=0.5, las=1, lwd=1.5, tick=F)
  mtext(expression(paste(plain("Relative spawning biomass"))), las=3, side=2, line=3, cex=1, lwd=1.75)
  mtext(expression(paste(italic("F") ~ (year^{-1}))), las=1, side=1, line=3, cex=1, lwd=1.75)
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


#' Get outputs from per age-based recruit analysis across a range of fishing mortality values, allowing for error in
#' M, h and F
#'
#' This function provides outputs associated with age-based per recruit analysis, and an
#' extended form of this analysis with a Beverton-Holt stock recruitment relationship to account
#' for potential impacts of fishing on recruitment. Outputs are provided for a range of
#' fishing mortality values, including the current, estimated value, with specified error in M, h and F.
#'
#' @param MaxModelAge maximum age considered in model
#' @param TimeStep model timestep (e.g. 1 = annual, 1/12 = monthly)
#' @param Linf asymptotic length (von Bertalanffy growth parameter, set to NA if inputting lengths at ages directly from any growth model)
#' @param vbK growth coefficient (von Bertalanffy growth parameter, set to NA if inputting lengths at ages directly from any growth model)
#' @param tzero hypothetical age at zero length (von Bertalanffy growth parameter, set to NA if inputting lengths at ages directly from any growth model)
#' @param EstLenAtAge vector or estimated lengths at ages from any growth model, set to NA if von Bertalanffy growth parameters specified
#' @param lenwt_a weight-length parameter (power relationship)
#' @param ln_lenwt_a weight-length parameter (log-log relationship)
#' @param lenwt_b weight-length parameter (power or log-log relationship)
#' @param WLrel_Type 1=power, 2=log-log relationship (set to NA if inputting weights at ages directly)
#' @param EstWtAtAge vector of weights at ages (set to NA if weight-length growth parameters specified)
#' @param ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' @param ReprodPattern reproductive pattern, 1=gonochoristic, 2=protogynous (female to male sex change) hermaphroditism,
#' 3=protandrous hermaphroditism (male to female sex change)
#' @param InitRatioFem proportion of fish that are females at hatching
#' @param FinalSex_Pmax logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_A50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_A95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param mat_A50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_A95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtAge vector of proportion mature at age (set to NA if using age at maturity parameters)
#' @param sel_A50 logistic parameter for gear selectivity curve
#' @param sel_A95 logistic parameter for gear selectivity curve
#' @param EstSelAtAge vector of gear selectivity at age (set to NA if using age at gear selectivity parameters)
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_A50 logistic parameter for fish retention curve
#' @param ret_A95 logistic parameter for fish retention curve
#' @param EstRetenAtAge vector of fish retention at age (set to NA if using age at fish retention parameters)
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param Steepness_sd specified error for steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param NatMort_sd specified error for natural mortality
#' @param RefPointPlotOpt plotting option for reference points, 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' @param Current_F estimate of current fishing mortality
#' @param Current_F_sd error for estimate of current fishing mortality
#' @param nReps number of random parameter sets from parametric resampling to generate per recruit outputs with error
#'
#' @return range of fishing mortality values for calculating per recruit quantities (PerRec_FValues)
#' estimated female SPR and equilibrium relative biomass values for trials with random M, h and F,
#' (Fem_SPR_Vals, Equil_RelFemSpBiom_Vals), estimated relative biomass at maximum sustainable yield
#' by trial (BMSY_Vals), estimated female spawning potential ratio values by F for each trial (Sim_FemSPR)
#' estimated female relative biomass values by F for each trial (Sim_Equil_RelFemSpBiom), median,
#' lower 95 and upper 95 percentile values for SPR for current F (EstFemSPR, EstLow95FemSPR, EstUp95FemSPR), median,
#' lower 95 and upper 95 percentile values for equilibrium relative female spawning biomass, median
#' BMSY ratio (EstBMSYratio), summary matrix containing key outputs from analysis (ResSummary_with_err)
#'
#' @examples
#' # Example. Non-hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 20 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be approprriate to use a smaller time step)
#' Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.5, 0.5) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstFemLenAtAge=NA, EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- NA # for log-log relationship
#' lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA, EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' EstSelAtAge <- data.frame(EstFemSelAtAge=NA, EstMalSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
#' EstRetenAtAge <- data.frame(EstFemRetenAtAge=rep(1,MaxModelAge+1), EstMalRetenAtAge=rep(1,MaxModelAge+1)) # fish retention at age (from age 0), inputted as values in data frame
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
#' GetPerRecruitResults_AB_with_err(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                               lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern,
#'                               InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95, mat_A50, mat_A95,
#'                               EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95, EstRetenAtAge,
#'                               DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd, Current_F,
#'                               Current_F_sd, nReps)
#' @export
GetPerRecruitResults_AB_with_err <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                          lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern,
                                          InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95, mat_A50, mat_A95,
                                          EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
                                          EstRetenAtAge, DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd,
                                          Current_F, Current_F_sd, nReps) {



  FValues = rnorm(nReps, Current_F, Current_F_sd)
  hValues = rnorm(nReps, Steepness, Steepness_sd)
  MValues = rnorm(nReps, NatMort, NatMort_sd)


  for (i in 1:nReps) {
    FMort = FValues[i]
    Steepness = hValues[i]
    NatMort = MValues[i]
    PREst = GetPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                 lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                                 ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                                 mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge,
                                 ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness,
                                 SRrel_Type, NatMort, Current_F)

    if (i == 1) {
      Fem_SPR_Vals = rep(0, nReps)
      Equil_RelFemSpBiom_Vals = rep(0, nReps)
      BMSY_Vals = rep(0,nReps)
      Sim_FemSPR = data.frame(matrix(nrow=nReps, ncol=length(PREst$FishMort)))
      colnames(Sim_FemSPR) = PREst$FishMort
      Sim_FemSPR = as.matrix(Sim_FemSPR)
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
  colnames(ResSummary_with_err)=c("Fem_SPR", "Fem_LowSPR", "Fem_UppSPR", "EquilSB", "LowEquilSB", "UppEquilSB", "BMSYratio")

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


#' Get outputs from length-based per recruit analysis across a range of fishing mortality values, allowing for error in
#' M, h and F
#'
#' This function provides outputs associated with length-based per recruit analysis, and an
#' extended form of this analysis with a Beverton-Holt stock recruitment relationship to account
#' for potential impacts of fishing on recruitment. Outputs are provided for a range of
#' fishing mortality values, including the current, estimated value, with specified error in M, h and F.
#'
#' @param MaxModelAge maximum age considered by model
#' @param TimeStep model time step (in y)
#' @param lbnd lower bounds of length classes
#' @param ubnd upper bounds of length classes
#' @param midpt mid points bounds of length classes
#' @param nLenCl number of length classes
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
#' @param GrowthParams growth parameters of either von Bertalanffy or Schnute model
#' @param RefnceAges reference ages for Schnute model, set to NA if using von Bertalanffy model
#' @param CVSizeAtAge coefficient of variation for size at age
#' @param lenwt_a weight-length parameter
#' @param ln_lenwt_a weight-length parameter
#' @param lenwt_b weight-length parameter
#' @param WLrel_Type 1=power, 2=log-log
#' @param EstWtAtLen user-specified weights at lengths
#' @param ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' @param ReprodPattern reproductive pattern, 1=gonochoristic, 2=protogynous (female to male sex change) hermaphroditism,
#' 3=protandrous hermaphroditism (male to female sex change)
#' @param InitRatioFem proportion of fish that are females at hatching
#' @param FinalSex_Pmax logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param mat_L50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_L95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtLen vector of proportion mature at length (set to NA if using age at maturity parameters)
#' @param sel_L50 logistic parameter for gear selectivity curve
#' @param sel_L95 logistic parameter for gear selectivity curve
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_L50 logistic parameter for fish retention curve
#' @param ret_L95 logistic parameter for fish retention curve
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param Steepness_sd standard deviation for steepness parameter
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param NatMort_sd standard deviation for natural mortality
#' @param Current_F estimated current fishing mortality
#' @param Current_F_sd standard deviation for estimated current fishing mortality
#' @param nReps number of resamping trials
#'
#' @return fishing mortality values for analysis (PerRec_FValues), female SPR estimates vs F (Fem_SPR_Vals)
#' estimates of equilbrium female spawning biomass vs F (Equil_RelFemSpBiom_Vals), estimated BMSY values (BMSY_Vals)
#' resampled female SPR values vs F (Sim_FemSPR=Sim_FemSPR), resampled equilbrium female spawning biomass values vs F (Sim_Equil_RelFemSpBiom),
#' median female SPR values vs F (EstFemSPR), low 95 percent female SPR values vs F (EstLow95FemSPR), upper 95 percent
#' female SPR values vs F (EstUp95FemSPR), median equilibrium relative female spawning biomass vs F (EstEquilRelFemSpBiom)
#' low 95 percent female relative female spawning biomass vs F (Low95EquilRelFemSpBiom), upper 95 percent relative female
#' spawning biomass vs F (Upp95EquilRelFemSpBiom), BMSY ratio estimates (EstBMSYratio), summary data frame with estimates
#' and associated 95 percent confidence limits.
#'
#' @examples
#' MaxModelAge <- 20 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
#' MaxLen = 800
#' LenInc = 20
#' lbnd = seq(0,MaxLen - LenInc, LenInc)
#' ubnd = lbnd + LenInc
#' midpt = lbnd + (LenInc/2)
#' nLenCl = length(midpt)
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.2, 0.2) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' GrowthParams = data.frame(Linf=Linf, vbK=vbK, tzero=tzero)
#' RefnceAges = NA
#' # GrowthParams c(Linf, vbK, tzero) single sex von Bertalanffy, or data.frame(Linf=Linf, vbK=vbK, tzero=tzero),
#' #' both sexes von Bertalanffy, or c(y1, y2, a, b) single sex Schnute, or data.frame(y1=y1, y2=y2, a=a, b=b), both sexes Schnute
#' CVSizeAtAge = c(0.05,0.05)
#' lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- NA # for log-log relationship
#' lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' EstWtAtLen <- data.frame(EstFemWtAtLen=NA,
#'                          EstMalWtAtLen=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_L50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_L95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' mat_L50 <- c(250, 250) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_L95 <- c(300, 300) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
#'                           EstMalMatAtLen=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_L50 <- c(250, 250) # females, males - Logistic age selectivity relationship parameters
#' sel_L95 <- c(300, 300) # females, males - Logistic age selectivity relationship parameters
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_L50 <- c(250, 250) # females, males - Logistic age fish retention at age parameters
#' ret_L95 <- c(300, 300) # females, males - Logistic age fish retention at age parameters
#' DiscMort <- 0.25 # discard mortality (e.g. 50% released fish die = 0.5)
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' Steepness_sd <- 0.025
#' NatMort <- 0.2 # natural mortality  (year-1)
#' NatMort_sd <- 0.025
#' Current_F <- 0.1 # estimate of fishing mortality, e.g. from catch curve analysis
#' Current_F_sd <- 0.005
#' nReps = 10 # number of resampling trials. Set to low number to test, then much higher for final analysis.
#' GetPerRecruitResults_LB_with_err(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
#'                                  RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
#'                                  EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
#'                                  FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, ret_Pmax,
#'                                  ret_L50, ret_L95, DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd,
#'                                  Current_F, Current_F_sd, nReps)
#' @export
GetPerRecruitResults_LB_with_err <- function(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                             RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
                                             EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
                                             FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, ret_Pmax,
                                             ret_L50, ret_L95, DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd,
                                             Current_F, Current_F_sd, nReps) {



  FValues = rnorm(nReps, Current_F, Current_F_sd)
  hValues = rnorm(nReps, Steepness, Steepness_sd)
  MValues = rnorm(nReps, NatMort, NatMort_sd)

  for (i in 1:nReps) {
    FMort = FValues[i]
    Steepness = hValues[i]
    NatMort = MValues[i]
    PREst = GetPerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                    RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
                                    EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
                                    FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, ret_Pmax,
                                    ret_L50, ret_L95, DiscMort, Steepness, SRrel_Type, NatMort, Current_F)

    if (i == 1) {
      Fem_SPR_Vals = rep(0, nReps)
      Mal_SPR_Vals = rep(0, nReps)
      CombSex_SPR_Vals = rep(0, nReps)

      Equil_RelFemSpBiom_Vals = rep(0, nReps)
      Equil_RelMalSpBiom_Vals = rep(0, nReps)
      Equil_RelCombSexSpBiom_Vals = rep(0, nReps)

      BMSY_Vals = rep(0,nReps)
      Sim_FemSPR = data.frame(matrix(nrow=nReps, ncol=length(PREst$FishMort)))
      colnames(Sim_FemSPR) = PREst$FishMort
      Sim_Equil_RelFemSpBiom = as.matrix(Sim_FemSPR)
      Sim_Equil_RelMalSpBiom = Sim_Equil_RelFemSpBiom
      Sim_Equil_RelCombSexSpBiom = Sim_Equil_RelFemSpBiom
    }

    Fem_SPR_Vals[i] = PREst$Fem_SPR
    Mal_SPR_Vals[i] = PREst$Mal_SPR
    CombSex_SPR_Vals[i] = PREst$CombSex_SPR

    Equil_RelFemSpBiom_Vals[i] = PREst$Equilmod_FemRelBiom
    Equil_RelMalSpBiom_Vals[i] = PREst$Equilmod_MalRelBiom
    Equil_RelCombSexSpBiom_Vals[i] = PREst$Equilmod_CombSexRelBiom

    BMSY_Vals[i] = PREst$BMSY_Thresh
    Sim_FemSPR[i,] = PREst$Fem_SPRResults
    Sim_Equil_RelFemSpBiom[i,] = PREst$Equilmod_FemRelBiomResults
    Sim_Equil_RelMalSpBiom[i,] = PREst$Equilmod_MalRelBiomResults
    Sim_Equil_RelCombSexSpBiom[i,] = PREst$Equilmod_CombSexRelBiomResults

    cat("i",i,'\n')
  }

  # Save key outputs (median and 95% confidence levels)
  EstMedFemSPR = round(median(Fem_SPR_Vals),3)
  EstLow95FemSPR = as.numeric(round(quantile(Fem_SPR_Vals, 0.025),3))
  EstUp95FemSPR = as.numeric(round(quantile(Fem_SPR_Vals, 0.975),3))
  EstFemSPR = data.frame(EstMedFemSPR=EstMedFemSPR,
                         EstLow95FemSPR=EstLow95FemSPR,
                         EstUp95FemSPR=EstUp95FemSPR)

  EstMedMalSPR = round(median(Mal_SPR_Vals),3)
  EstLow95MalSPR = as.numeric(round(quantile(Mal_SPR_Vals, 0.025),3))
  EstUp95MalSPR = as.numeric(round(quantile(Mal_SPR_Vals, 0.975),3))
  EstMalSPR = data.frame(EstMedMalSPR=EstMedMalSPR,
                         EstLow95MalSPR=EstLow95MalSPR,
                         EstUp95MalSPR=EstUp95MalSPR)

  EstMedCombSexSPR = round(median(CombSex_SPR_Vals),3)
  EstLow95CombSexSPR = as.numeric(round(quantile(CombSex_SPR_Vals, 0.025),3))
  EstUp95CombSexSPR = as.numeric(round(quantile(CombSex_SPR_Vals, 0.975),3))
  EstCombSexSPR = data.frame(EstMedCombSexSPR=EstMedCombSexSPR,
                         EstLow95CombSexSPR=EstLow95CombSexSPR,
                         EstUp95CombSexSPR=EstUp95CombSexSPR)

  EstMedEquilRelFemSpBiom = round(median(Equil_RelFemSpBiom_Vals),3)
  Low95EquilRelFemSpBiom = as.numeric(round(quantile(Equil_RelFemSpBiom_Vals, 0.025),3))
  Upp95EquilRelFemSpBiom = as.numeric(round(quantile(Equil_RelFemSpBiom_Vals, 0.975),3))
  EstEquilRelFemSpBiom = data.frame(EstMedEquilRelFemSpBiom=EstMedEquilRelFemSpBiom,
                                    Low95EquilRelFemSpBiom=Low95EquilRelFemSpBiom,
                                    Upp95EquilRelFemSpBiom=Upp95EquilRelFemSpBiom)

  EstMedEquilRelMalSpBiom = round(median(Equil_RelMalSpBiom_Vals),3)
  Low95EquilRelMalSpBiom = as.numeric(round(quantile(Equil_RelMalSpBiom_Vals, 0.025),3))
  Upp95EquilRelMalSpBiom = as.numeric(round(quantile(Equil_RelMalSpBiom_Vals, 0.975),3))
  EstEquilRelMalSpBiom = data.frame(EstMedEquilRelMalSpBiom=EstMedEquilRelMalSpBiom,
                                    Low95EquilRelMalSpBiom=Low95EquilRelMalSpBiom,
                                    Upp95EquilRelMalSpBiom=Upp95EquilRelMalSpBiom)

  EstMedEquilRelCombSexSpBiom = round(median(Equil_RelCombSexSpBiom_Vals),3)
  Low95EquilRelCombSexSpBiom = as.numeric(round(quantile(Equil_RelCombSexSpBiom_Vals, 0.025),3))
  Upp95EquilRelCombSexSpBiom = as.numeric(round(quantile(Equil_RelCombSexSpBiom_Vals, 0.975),3))
  EstEquilRelCombSexSpBiom = data.frame(EstMedEquilRelCombSexSpBiom=EstMedEquilRelCombSexSpBiom,
                                        Low95EquilRelCombSexSpBiom=Low95EquilRelCombSexSpBiom,
                                        Upp95EquilRelCombSexSpBiom=Upp95EquilRelCombSexSpBiom)


  EstMedBMSYratio = round(median(BMSY_Vals),3)
  Low95EstBMSYratio = as.numeric(round(quantile(BMSY_Vals, 0.025),3))
  Upp95EstBMSYratio = as.numeric(round(quantile(BMSY_Vals, 0.975),3))

  ResSummary_with_err <- data.frame(EstMedFemSPR, EstLow95FemSPR, EstUp95FemSPR,
                                    EstMedMalSPR, EstLow95MalSPR, EstUp95MalSPR,
                                    EstMedCombSexSPR, EstLow95CombSexSPR, EstUp95CombSexSPR,
                                    EstMedEquilRelFemSpBiom, Low95EquilRelFemSpBiom, Upp95EquilRelFemSpBiom,
                                    EstMedEquilRelMalSpBiom, Low95EquilRelMalSpBiom, Upp95EquilRelMalSpBiom,
                                    EstMedEquilRelCombSexSpBiom, Low95EquilRelCombSexSpBiom, Upp95EquilRelCombSexSpBiom,
                                    EstMedBMSYratio, Low95EstBMSYratio, Upp95EstBMSYratio)

  colnames(ResSummary_with_err)=c("Fem_SPR", "Fem_LowSPR", "Fem_UppSPR",
                                  "Mal_SPR", "Mal_LowSPR", "Mal_UppSPR",
                                  "CombSex_SPR", "CombSex_LowSPR", "CombSex_UppSPR",
                                  "Fem_EquilSB", "Fem_LowEquilSB", "Fem_UppEquilSB",
                                  "Mal_EquilSB", "Mal_LowEquilSB", "Mal_UppEquilSB",
                                  "CombSex_EquilSB", "CombSex_LowEquilSB", "CombSex_UppEquilSB",
                                  "BMSYratio", "BMSYratio_Low", "BMSYratio_Upp")

  Results = list(PerRec_FValues = PREst$FishMort,
                 Fem_SPR_Vals=Fem_SPR_Vals,
                 Equil_RelFemSpBiom_Vals=Equil_RelFemSpBiom_Vals,
                 BMSY_Vals=BMSY_Vals,
                 Sim_FemSPR=Sim_FemSPR,
                 Sim_Equil_RelFemSpBiom=Sim_Equil_RelFemSpBiom,
                 Sim_Equil_RelMalSpBiom=Sim_Equil_RelMalSpBiom,
                 Sim_Equil_RelCombSexSpBiom=Sim_Equil_RelCombSexSpBiom,
                 EstFemSPR=EstFemSPR,
                 EstMalSPR=EstMalSPR,
                 EstCombSexSPR=EstCombSexSPR,
                 EstEquilRelFemSpBiom=EstEquilRelFemSpBiom,
                 EstEquilRelMalSpBiom=EstEquilRelMalSpBiom,
                 EstEquilRelCombSexSpBiom=EstEquilRelCombSexSpBiom,
                 EstMedBMSYratio=EstMedBMSYratio,
                 ResSummary_with_err=ResSummary_with_err)

  return(Results)

}

#' Plot of SPR and relative equilibrium biomass vs F from age-based per recruit analysis and extended analysis
#' with a stock-recruitment relationship, with error for M, h and F
#'
#' This function provides a plot of SPR and relative equilibrium biomass vs F from age-based per recruit analysis
#' and extended analysis with a stock-recruitment relationship, with specified error for M, h and F
#'
#' @param MaxModelAge maximum age considered in model
#' @param TimeStep model timestep (e.g. 1 = annual, 1/12 = monthly)
#' @param Linf asymptotic length (von Bertalanffy growth parameter, set to NA if inputting lengths at ages directly from any growth model)
#' @param vbK growth coefficient (von Bertalanffy growth parameter, set to NA if inputting lengths at ages directly from any growth model)
#' @param tzero hypothetical age at zero length (von Bertalanffy growth parameter, set to NA if inputting lengths at ages directly from any growth model)
#' @param EstLenAtAge vector or estimated lengths at ages from any growth model, set to NA if von Bertalanffy growth parameters specified
#' @param lenwt_a weight-length parameter (power relationship)
#' @param ln_lenwt_a weight-length parameter (log-log relationship)
#' @param lenwt_b weight-length parameter (power or log-log relationship)
#' @param WLrel_Type 1=power, 2=log-log relationship (set to NA if inputting weights at ages directly)
#' @param EstWtAtAge vector of weights at ages (set to NA if weight-length growth parameters specified)
#' @param ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' @param ReprodPattern reproductive pattern, 1=gonochoristic, 2=protogynous (female to male sex change) hermaphroditism,
#' 3=protandrous hermaphroditism (male to female sex change)
#' @param InitRatioFem proportion of fish that are females at hatching
#' @param FinalSex_Pmax logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_A50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_A95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param mat_A50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_A95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtAge vector of proportion mature at age (set to NA if using age at maturity parameters)
#' @param sel_A50 logistic parameter for gear selectivity curve
#' @param sel_A95 logistic parameter for gear selectivity curve
#' @param EstSelAtAge vector of gear selectivity at age (set to NA if using age at gear selectivity parameters
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_A50 logistic parameter for fish retention curve
#' @param ret_A95 logistic parameter for fish retention curve
#' @param EstRetenAtAge vector of fish retention at age (set to NA if using age at fish retention parameters
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param Steepness_sd specified error for steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param NatMort_sd specified error for natural mortality
#' @param RefPointPlotOpt plotting option for reference points, 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' @param Current_F estimate of current fishing mortality
#' @param Current_F_sd error for estimate of current fishing mortality
#' @param nReps number of random parameter sets from parametric resampling to generate per recruit outputs with error
#' @param MainLabel plot label
#' @param xaxis_lab y axis label
#' @param yaxis_lab x axis label
#' @param xmax maximum x axis value
#' @param xint x axis interval
#' @param ymax maximum y axis value
#' @param yint y axis interval
#'
#' @return plot of SPR and relative equilibrium biomass vs F from per recruit analysis
#' and extended analysis with a stock-recruitment relationship, with specified error for M, h and F
#'
#' @examples
#' # Example. Non-hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 20 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
#' Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.5, 0.5) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstFemLenAtAge=NA, EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- NA # for log-log relationship
#' lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA, EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' EstSelAtAge <- data.frame(EstFemSelAtAge=NA, EstMalSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
#' EstRetenAtAge <- data.frame(EstFemRetenAtAge=rep(1,MaxModelAge+1), EstMalRetenAtAge=rep(1,MaxModelAge+1)) # fish retention at age (from age 0), inputted as values in data frame
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
#' FittedRes=GetPerRecruitResults_AB_with_err(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                                              lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern,
#'                                              InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95, mat_A50, mat_A95,
#'                                              EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
#'                                              EstRetenAtAge, DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd,
#'                                              Current_F, Current_F_sd, nReps)
#' # Plot. Note, can skip above step and set FittedRes=NA (plot function will be slower
#' PlotPerRecruit_Biom_with_err_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern,
#'                                             InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95, mat_A50, mat_A95,
#'                                             EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
#'                                             EstRetenAtAge, DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd,
#'                                             Current_F, Current_F_sd, RefPointPlotOpt, FittedRes, nReps, MainLabel=NA,
#'                                             xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA, ymax=NA, yint=NA)
#' # # Example. Non-hermaphroditic species
#' # InitRecruit <- 1 # Initial recruitment
#' # MaxModelAge <- 20 # maximum age considered by model, years
#' # TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
#' # Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' # vbK <- c(0.5, 0.5) # year-1 - von Bertalanffy growth model parameters - Females, males
#' # tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' # EstLenAtAge <- data.frame(EstFemLenAtAge=NA, EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' # lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' # ln_lenwt_a <- NA # for log-log relationship
#' # lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' # WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' # EstWtAtAge <- data.frame(EstFemWtAtAge=NA, EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' # ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' # ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' # InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' # FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' # FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
#' # FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' # mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' # mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' # EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' # sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' # sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' # EstSelAtAge <- data.frame(EstFemSelAtAge=NA, EstMalSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' # ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' # ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
#' # ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
#' # EstRetenAtAge <- data.frame(EstFemRetenAtAge=rep(1,MaxModelAge+1), EstMalRetenAtAge=rep(1,MaxModelAge+1)) # fish retention at age (from age 0), inputted as values in data frame
#' # DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' # Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' # Steepness_sd <- 0.025
#' # SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' # NatMort <- 0.2 # natural mortality  (year-1)
#' # NatMort_sd <- 0.025
#' # Current_F <- 0.1 # estimate of fishing mortality, e.g. from catch curve analysis
#' # Current_F_sd <- 0.005
#' # RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' # nReps = 50
#' FittedRes=GetPerRecruitResults_AB_with_err(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                                              lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern,
#'                                              InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95, mat_A50, mat_A95,
#'                                              EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
#'                                              EstRetenAtAge, DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd,
#'                                              Current_F, Current_F_sd, nReps)
#' # Plot. Note, can skip above step and set FittedRes=NA (plot function will be slower
#' PlotPerRecruit_Biom_with_err_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern,
#'                                             InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95, mat_A50, mat_A95,
#'                                             EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
#'                                             EstRetenAtAge, DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd,
#'                                             Current_F, Current_F_sd, RefPointPlotOpt, FittedRes, nReps, MainLabel=NA,
#'                                             xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA, ymax=NA, yint=NA)
#' @export
PlotPerRecruit_Biom_with_err_AB <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                            lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern,
                                            InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95, mat_A50, mat_A95,
                                            EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
                                            EstRetenAtAge, DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd,
                                            Current_F, Current_F_sd, RefPointPlotOpt, FittedRes, nReps, MainLabel,
                                            xaxis_lab, yaxis_lab, xmax, xint, ymax, yint) {

  # get BMSY reference points
  res = GetPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                                ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95,
                                mat_A50, mat_A95, EstMatAtAge, sel_A50, sel_A95, EstSelAtAge,
                                ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness,
                                SRrel_Type, NatMort, Current_F)

  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    Res =  FittedRes
  } else {
    Res=GetPerRecruitResults_AB_with_err(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                         lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern,
                                         InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95, mat_A50, mat_A95,
                                         EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
                                         EstRetenAtAge, DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd,
                                         Current_F, Current_F_sd, nReps)
  }

  if (is.na(yaxis_lab)) yaxis_lab = "Relative spawning biomass"
  if (is.na(xaxis_lab)) xaxis_lab = expression(paste(italic("F") ~ (year^{-1})))
  if (is.na(xmax)) xmax = 2
  if (is.na(xint)) xint = 0.5
  if (is.na(ymax)) ymax = 1
  if (is.na(yint)) yint = 0.2

  # Plot per recruit outputs with uncertainty
  EqB_med = apply(Res$Sim_Equil_RelFemSpBiom,2,quantile, probs=c(0.5))
  EqB_lw = apply(Res$Sim_Equil_RelFemSpBiom,2,quantile, probs=c(0.025))
  EqB_hi = apply(Res$Sim_Equil_RelFemSpBiom,2,quantile, probs=c(0.975))
  plot(Res$PerRec_FValues, EqB_med, "l", frame.plot=F, ylim=c(0,ymax), xlim=c(0,xmax),
       main=MainLabel, col="black", yaxt="n", xaxt="n", ylab="", xlab="")
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
  x=which(res$FishMort==xmax)
  polygon(c(res$FishMort[1:x],rev(res$FishMort[1:x])),c(EqB_lw[1:x],rev(EqB_hi[1:x])), col=grey(0.5,0.25),
          border=grey(0.5,0.25))
  axis(1, at=seq(0, xmax, xint), cex.axis=1, lwd=1, lab=F, line=-0.3)
  axis(2, at=seq(0, ymax, yint), cex.axis=1, lwd=1, lab=F, line=-0.3)
  axis(1, at=seq(0, xmax, xint), labels = seq(0, xmax, xint),
       cex.axis=1, line=-0.5, las=1, lwd=1, tick=F)
  axis(2, at=seq(0, ymax, yint), cex.axis=1, line=-0.5, las=1, lwd=1, tick=F)
  mtext(yaxis_lab, las=3, side=2, line=3, cex=1.2, lwd=1.75)
  mtext(xaxis_lab, las=1, side=1, line=3, cex=1.2, lwd=1.75)
  legend("topleft", col="black", pch = 16, legend="Estimate",
         bty="n", cex=1.0, lty=0, inset = 0.05)
}

#' Plot of SPR and relative equilibrium biomass vs F from length-based per recruit analysis and extended analysis
#' with a stock-recruitment relationship, with error for M, h and F
#'
#' This function provides a plot of SPR and relative equilibrium biomass vs F from length-based per recruit analysis
#' and extended analysis with a stock-recruitment relationship, with specified error for M, h and F
#'
#' @param MaxModelAge maximum age considered by model
#' @param TimeStep model time step (in y)
#' @param lbnd lower bounds of length classes
#' @param ubnd upper bounds of length classes
#' @param midpt mid points bounds of length classes
#' @param nLenCl number of length classes
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
#' @param GrowthParams growth parameters of either von Bertalanffy or Schnute model
#' @param RefnceAges reference ages for Schnute model, set to NA if using von Bertalanffy model
#' @param CVSizeAtAge coefficient of variation for size at age
#' @param lenwt_a weight-length parameter
#' @param ln_lenwt_a weight-length parameter
#' @param lenwt_b weight-length parameter
#' @param WLrel_Type 1=power, 2=log-log
#' @param EstWtAtLen user-specified weights at lengths
#' @param ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' @param ReprodPattern reproductive pattern, 1=gonochoristic, 2=protogynous (female to male sex change) hermaphroditism,
#' 3=protandrous hermaphroditism (male to female sex change)
#' @param InitRatioFem proportion of fish that are females at hatching
#' @param FinalSex_Pmax logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param mat_L50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_L95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtLen vector of proportion mature at length (set to NA if using age at maturity parameters)
#' @param sel_L50 logistic parameter for gear selectivity curve
#' @param sel_L95 logistic parameter for gear selectivity curve
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_L50 logistic parameter for fish retention curve
#' @param ret_L95 logistic parameter for fish retention curve
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param Steepness_sd standard deviation for steepness parameter
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param NatMort_sd standard deviation for natural mortality
#' @param Current_F estimated current fishing mortality
#' @param Current_F_sd standard deviation for estimated current fishing mortality
#' @param RefPointPlotOpt # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' @param FittedRes option to input per recruit results into plot function to increase speed, can set to NA
#' @param nReps number of resamping trials
#' @param MainLabel plot label
#' @param xaxis_lab y axis label
#' @param yaxis_lab x axis label
#' @param xmax maximum x axis value
#' @param xint x axis interval
#' @param ymax maximum y axis value
#' @param yint y axis interval
#'
#' @examples
#' MaxModelAge <- 20 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
#' MaxLen = 800
#' LenInc = 20
#' lbnd = seq(0,MaxLen - LenInc, LenInc)
#' ubnd = lbnd + LenInc
#' midpt = lbnd + (LenInc/2)
#' nLenCl = length(midpt)
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf <- c(550, 500) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.2, 0.2) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' GrowthParams = data.frame(Linf=Linf, vbK=vbK, tzero=tzero)
#' RefnceAges = NA
#' # GrowthParams c(Linf, vbK, tzero) single sex von Bertalanffy, or data.frame(Linf=Linf, vbK=vbK, tzero=tzero),
#' #' both sexes von Bertalanffy, or c(y1, y2, a, b) single sex Schnute, or data.frame(y1=y1, y2=y2, a=a, b=b), both sexes Schnute
#' CVSizeAtAge = c(0.05, 0.05)
#' lenwt_a <- 0.000005 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- NA # for log-log relationship
#' lenwt_b <- 3 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 1 # 1=power, 2=log-log relationship
#' EstWtAtLen <- data.frame(EstFemWtAtLen=NA,
#'                          EstMalWtAtLen=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at age zero
#' FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
#' FinalSex_L50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_L95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' mat_L50 <- c(250, 250) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_L95 <- c(300, 300) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
#'                           EstMalMatAtLen=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_L50 <- c(250, 250) # females, males - Logistic age selectivity relationship parameters
#' sel_L95 <- c(300, 300) # females, males - Logistic age selectivity relationship parameters
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_L50 <- c(250, 250) # females, males - Logistic age fish retention at age parameters
#' ret_L95 <- c(300, 300) # females, males - Logistic age fish retention at age parameters
#' DiscMort <- 0.25 # discard mortality (e.g. 50% released fish die = 0.5)
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' Steepness_sd <- 0.025
#' NatMort <- 0.2 # natural mortality  (year-1)
#' NatMort_sd <- 0.025
#' Current_F <- 0.1 # estimate of fishing mortality, e.g. from catch curve analysis
#' Current_F_sd <- 0.005
#' PlotOpt <- 1 # 1=females, 2=males, 3=combined sex
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' nReps = 10 # number of resampling trials. Set to low number to test, then much higher for final analysis.
#' FittedRes=GetPerRecruitResults_LB_with_err(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
#'                                      RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
#'                                      EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
#'                                      FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, ret_Pmax,
#'                                       ret_L50, ret_L95, DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd,
#'                                       Current_F, Current_F_sd, nReps)
#' # Plot. Note, can skip above step and set FittedRes=NA (plot function will be slower
#' PlotPerRecruit_Biom_with_err_LB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern,
#'                                             InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95, mat_A50, mat_A95,
#'                                             EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
#'                                             EstRetenAtAge, DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd,
#'                                             Current_F, Current_F_sd, PlotOpt, RefPointPlotOpt, FittedRes, nReps, MainLabel=NA,
#'                                             xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA, ymax=NA, yint=NA)
#' @export
PlotPerRecruit_Biom_with_err_LB <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                            lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern,
                                            InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95, mat_A50, mat_A95,
                                            EstMatAtAge, sel_A50, sel_A95, EstSelAtAge, ret_Pmax, ret_A50, ret_A95,
                                            EstRetenAtAge, DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd,
                                            Current_F, Current_F_sd, PlotOpt, RefPointPlotOpt, FittedRes, nReps, MainLabel,
                                            xaxis_lab, yaxis_lab, xmax, xint, ymax, yint) {


  # get BMSY reference points
  res=GetPerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                              RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
                              EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
                              FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, ret_Pmax,
                              ret_L50, ret_L95, DiscMort, Steepness, SRrel_Type, NatMort, Current_F)


  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    Res =  FittedRes
  } else {
    Res=GetPerRecruitResults_LB_with_err(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                         RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
                                         EstWtAtLen, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
                                         FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, ret_Pmax,
                                         ret_L50, ret_L95, DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd,
                                         Current_F, Current_F_sd, nReps)
  }

  if (is.na(yaxis_lab)) yaxis_lab = "Relative spawning biomass"
  if (is.na(xaxis_lab)) xaxis_lab = expression(paste(italic("F") ~ (year^{-1})))
  if (is.na(xmax)) xmax = 2
  if (is.na(xint)) xint = 0.5
  if (is.na(ymax)) ymax = 1
  if (is.na(yint)) yint = 0.2

  if (PlotOpt==1) { # plot females
    EqB_med = apply(Res$Sim_Equil_RelFemSpBiom,2,quantile, probs=c(0.5))
    EqB_lw = apply(Res$Sim_Equil_RelFemSpBiom,2,quantile, probs=c(0.025))
    EqB_hi = apply(Res$Sim_Equil_RelFemSpBiom,2,quantile, probs=c(0.975))
    plot(Res$PerRec_FValues, EqB_med, "l", frame.plot=F, ylim=c(0,ymax), xlim=c(0,xmax),
         col="red", yaxt="n", xaxt="n", ylab="", xlab="", main=MainLabel)
    x=which(Res$PerRec_FValues==xmax)
    polygon(c(Res$PerRec_FValues[1:x],rev(Res$PerRec_FValues[1:x])),c(EqB_lw[1:x],rev(EqB_hi[1:x])),
            col="lightpink", border="lightpink")
    lines(Res$PerRec_FValues, EqB_med,col="red")
    points(Current_F, Res$EstEquilRelFemSpBiom[1], cex=1.2, col="red", pch=16)
    lw=as.numeric(Res$EstEquilRelFemSpBiom[2]); up=as.numeric(Res$EstEquilRelFemSpBiom[3])
    arrows(Current_F, lw, Current_F, up,length=0.05, angle=90, code=3,col="red")
    legend("topleft", col="red", pch = 16, legend="Estimate - females",
           bty="n", cex=1,0, lty=0, inset = 0.05)

  }

  if (PlotOpt==2) { # plot males
    EqB_med = apply(Res$Sim_Equil_RelMalSpBiom,2,quantile, probs=c(0.5))
    EqB_lw = apply(Res$Sim_Equil_RelMalSpBiom,2,quantile, probs=c(0.025))
    EqB_hi = apply(Res$Sim_Equil_RelMalSpBiom,2,quantile, probs=c(0.975))
    plot(Res$PerRec_FValues, EqB_med, "l", frame.plot=F, ylim=c(0,ymax), xlim=c(0,xmax),
         col="blue", yaxt="n", xaxt="n", ylab="", xlab="", main=MainLabel)
    x=which(Res$PerRec_FValues==xmax)
    polygon(c(Res$PerRec_FValues[1:x],rev(Res$PerRec_FValues[1:x])),c(EqB_lw[1:x],rev(EqB_hi[1:x])),
            col="lightblue", border="lightblue")
    lines(Res$PerRec_FValues, EqB_med,col="blue")
    lw=as.numeric(Res$EstEquilRelMalSpBiom[2]); up=as.numeric(Res$EstEquilRelMalSpBiom[3])
    arrows(Current_F, lw, Current_F, up,length=0.05, angle=90, code=3,col="blue")
    points(Current_F, Res$EstEquilRelMalSpBiom[1], cex=1.2, col="blue", pch=16)
    legend("topleft", col="blue", pch = 16, legend="Estimate - males",
           bty="n", cex=1,0, lty=0, inset = 0.05)
  }

  if (PlotOpt==3) { # plot combined sex
    EqB_med = apply(Res$Sim_Equil_RelCombSexSpBiom,2,quantile, probs=c(0.5))
    EqB_lw = apply(Res$Sim_Equil_RelCombSexSpBiom,2,quantile, probs=c(0.025))
    EqB_hi = apply(Res$Sim_Equil_RelCombSexSpBiom,2,quantile, probs=c(0.975))
    plot(Res$PerRec_FValues, EqB_med, "l", frame.plot=F, ylim=c(0,ymax), xlim=c(0,xmax),
         col="black", yaxt="n", xaxt="n", ylab="", xlab="", main=MainLabel)
    x=which(Res$PerRec_FValues==xmax)
    polygon(c(Res$PerRec_FValues[1:x],rev(Res$PerRec_FValues[1:x])),c(EqB_lw[1:x],rev(EqB_hi[1:x])),
            col="lightgrey", border="lightgrey")
    lines(Res$PerRec_FValues, EqB_med,col="black")
    lw=as.numeric(Res$EstEquilRelCombSexSpBiom[2]); up=as.numeric(Res$EstEquilRelCombSexSpBiom[3])
    arrows(Current_F, lw, Current_F, up,length=0.05, angle=90, code=3,col="black")
    points(Current_F, Res$EstEquilRelCombSexSpBiom[1], cex=1.2, col="black", pch=16)
    legend("topleft", col="black", pch = 16, legend="Estimate - combined sex",
           bty="n", cex=1,0, lty=0, inset = 0.05)
  }

  axis(1, at=seq(0, xmax, xint), cex.axis=1, lwd=1, lab=F, line=-0.3)
  axis(2, at=seq(0, ymax, yint), cex.axis=1, lwd=1, lab=F, line=-0.3)
  axis(1, at=seq(0, xmax, xint), labels = seq(0, xmax, xint),
       cex.axis=1, line=-0.5, las=1, lwd=1, tick=F)
  axis(2, at=seq(0, ymax, yint), cex.axis=1, line=-0.5, las=1, lwd=1, tick=F)
  mtext(yaxis_lab, las=3, side=2, line=3, cex=1.2, lwd=1.75)
  mtext(xaxis_lab, las=1, side=1, line=3, cex=1.2, lwd=1.75)

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


}


