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

# Alex Hesp, last updated March 2024
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


 #' Generic function to add axes and axes labels to plots
 #'
 #' @keywords internal
 #' @param xmin x axis minimum
 #' @param xmax x axis maximum
 #' @param xint x axis tick interval
 #' @param ymin y axis minimum
 #' @param ymax y axis maximum
 #' @param yint y axis tick interval
 #' @param cexval label size
 #' @param cexaxisval axis size
 #' @param lwdval line widths
 #' @param lineval axis offset
 #' @param lasval axis orientation
 #'
 #' @return adds axes to plots
 AddAxesAndTickLabelsToPlot <- function(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA) {

   if (is.na(xmin)) xmin=0
   if (is.na(ymin)) ymin=0
   if (is.na(cexval)) cexval=1
   if (is.na(cexaxisval)) cexaxisval=1
   if (is.na(lwdval)) lwdval=1
   if (is.na(lasval)) lasval=1
   if (is.na(lineval)) lineval=0

   axis(1, at = seq(xmin, xmax, xint), line = lineval, labels = F)
   axis(2, at = seq(ymin, ymax, yint), line = lineval, labels = F)
   axis(1, at = seq(xmin, xmax, xint), lwd=lwdval, labels=T, line=lineval, cex=cexval, cex.axis=cexaxisval, las=lasval)
   axis(2, at = seq(ymin, ymax, yint), lwd=lwdval, labels=T, line=lineval, cex=cexval, cex.axis=cexaxisval, las=lasval)

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

       tzero = t1 - (y1 ^ b) * (t2 - t1) / (y2 ^ b - y1 ^ b)
       # cat("SchnuteGrowthfunction: 1 Age",Age,"tzero",tzero,'\n')
       if (Age < tzero) {
         y = 0
       } else {
         v = (y1 ^ b + (y2 ^ b - y1 ^ b) * (Age - t1) / (t2 - t1))
         y = v ^ (1 / b)
       }
     }
   } # a == 0

   if (a != 0) {
     if (b == 0) {
       # Eqn(16)
       y = y1 * exp(log(y2 / y1) * (1 - exp(-a * (Age - t1))) / (1 - exp(-a * (t2 - t1))))
       if (y < 0) {
         y = 0 }
     } else {
       # Eqn(15)

       # First. let's work out tzero
       if (t1 == 0 & y1 == 0) { # i.e. if t1 is set to zero, and y1 also set to zero
         tzero = 0
       } else {
         if (1 + (y1 ^ b) * (1 - exp(-a * (t2 - t1))) / (y2 ^ b - y1 ^ b) <= 0) {
           tzero = t1 - log(1E-4) / a
         } else {
           tzero = t1 - log(1 + (y1 ^ b) * (1 - exp(-a * (t2 - t1))) / (y2 ^ b - y1 ^ b)) / a
         }
         if (is.nan(tzero)) {
           cat("SchnuteGrowthfunction: Problem calculating tzero",'\n')
         }
       }
       # cat("SchnuteGrowthfunction: 2 Age",Age,"tzero",tzero,'\n')
       if (Age < tzero) {
         y = 0
       } else {
         v = (y1 ^ b + (y2 ^ b - y1 ^ b)
              * (1 - exp(-a * (Age - t1))) / (1 - exp(-a * (t2 - t1))))
         y = v ^ (1 / b)
       } # else
     } # b == 0
   }  # a != 0
   # cat("SchnuteGrowthfunction: 3",'\n')
   return(y)
 } # end function


 #' Inverse Schnute growth function
 #'
 #' Calculate age given length, from Schnute growth function
 #'
 #' @keywords internal
 #' @param MaxAge maximum age of species to be considered by model
 #' @param FishLen specified length
 #' @param t1 first reference age
 #' @param t2 second reference age
 #' @param y1 length at t1
 #' @param y2 length at t2
 #' @param a growth curve parameter
 #' @param a growth curve parameter
 #'
 #' @return Age at specified length
 InverseSchnuteGrowthfunction <- function (MaxAge, FishLen, t1, t2, y1, y2, a, b) {

   # return age from Length, given Schnute growth parameters
   # (can be used when both a and b are not equal to zero)
   Linf = ((exp(a*t2)*y2^b-exp(a*t1)*y1^b)/(exp(a*t2)-exp(a*t1)))^(1/b)

   #cat("InverseSchnuteGrowthfunction: FishLen",FishLen,"a",a,"b",b,"y1",y1,"y2",y2,'\n')

   if (is.nan(Linf)) {
     Age = MaxAge
   } else if (FishLen >= Linf-1) {
     tempLen = Linf-1
     Age=log(1-(tempLen^b-y1^b)/(y2^b-y1^b)*(1-exp(-a*(t2-t1))))/-a+t1
     # Age = MaxAge - 1
   } else {
     Age=log(1-(FishLen^b-y1^b)/(y2^b-y1^b)*(1-exp(-a*(t2-t1))))/-a+t1
   }
   #cat("Linf",Linf,"MaxAge",MaxAge,"FishLen",FishLen,"Age",Age,'\n')

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
 #' @param midpt MaxAge
 #'
 #' @return Expected lengths (ExpLen)
 CalcLengthAfterGrowthForTimetep <- function(GrowthCurveType, TimeStep, GrowthParamsForSex, RefnceAgesForSex, midpt, MaxAge) {

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
       StartAge[i] = InverseSchnuteGrowthfunction(MaxAge, FishLen, t1, t2, y1, y2, a, b)
       Age <- StartAge[i] + TimeStep
       ExpLen[i] = SchnuteGrowthfunction(Age, t1, t2, y1, y2, a, b)
     }
   }

   return(ExpLen)
 }

#******************************************
# Catch curve analyses - size-based methods
#******************************************

#' Calculate gillnet selectivity using the method of Kirkwood and Walker (1986)
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
      if (MeanLenRec[i] < 5) MeanLenRec[i] = 5 # i.e. assuming mean size at least 5 mm at 1 year of age
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
#' @param DistnType 1 = Multinomial, 2 = Dirichlet multinomial
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
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0 # proportion of fish that die due to natural mortality
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.08
#' GrowthParams = c(Linf, vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # 2 sexes, von Bertalanffy
#' DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = c(700,850)
#' vbK = c(0.3,0.2)
#' CVSizeAtAge = c(0.08,0.08)
#' GrowthParams = data.frame(Linf=Linf, vbK=vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # 1 sex, Schnute
#' DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
#' GrowthCurveType = 2 # 1 = von Bertalanffy, 2 = Schnute
#' t1 = 0.5 # growth - Schnute
#' t2 = 25 # growth - Schnute
#' y1 = 100 # growth - Schnute
#' y2 = 1000 # growth - Schnute
#' a = 0.02 # growth - Schnute
#' b = 3.0 # growth - Schnute
#' GrowthParams = c(y1, y2, a, b)
#' RefnceAges = c(t1,t2)
#' CVSizeAtAge = 0.05
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # 2 sexes, Schnute
#' DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
#' GrowthCurveType = 2 # 1 = von Bertalanffy, 2 = Schnute
#' t1 = c(0.5,0.5) # growth - Schnute
#' t2 = c(25,25) # growth - Schnute
#' y1 = c(100,100) # growth - Schnute
#' y2 = c(1000,1000) # growth - Schnute
#' a = c(0.02,0.02) # growth - Schnute
#' b = c(3,3) # growth - Schnute
#' CVSizeAtAge = c(0.05, 0.05)
#' GrowthParams = data.frame(y1=y1, y2=y2, a=a, b=b)
#' RefnceAges = data.frame(t1=t1,t2=t2)
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen
#' ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
#' PropReleased = NA # proportion of fish released, vector including mean and sd (option probably now obselete)
#' midpt=Res$midpt
#' lbnd=Res$lbnd
#' ubnd=Res$ubnd
#' InitFishMort = 0.25 # specify starting parameters
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
#' InitL50 = 400
#' InitDelta = 100
#' params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
#' # if estimating both selectivity and retention parameters
#' # InitL50Ret = 410
#' # InitDeltaRet = 60
#' # params = c(InitFishMort_logit, log(InitL50), log(InitDelta), log(InitL50Ret), log(InitDeltaRet))
#' FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
#'                                     lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
#' # Example with selectivity specified as a vector
#' # Simulate data
#' SampleSize=5000
#' set.seed(123)
#' MaxAge = 30
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = 0.2
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
#' SelectivityType=1 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' lbnd = seq(0,MaxLen - LenInc, LenInc)
#' midpt = lbnd + (LenInc/2)
#' SelectivityVec = 1 / (1 + exp(-log(19)*(midpt-400)/(500-400)))
#' PropReleased = NA
#' SelParams = c(NA, NA) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0 # proportion of fish that die due to natural mortality
#' # single sex, von Bertalanffy
#' DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.05
#' GrowthParams = c(Linf, vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen
#' ObsDiscCatchFreqAtLen = NA # vector including mean and sd
#' midpt=Res$midpt
#' lbnd=Res$lbnd
#' ubnd=Res$ubnd
#' InitFishMort = 0.25 # specify starting parameters
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
#' params = c(InitFishMort_logit)
#'
#' FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
#'                                     lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
#'
#' # Fit length-based catch curve to length-frequency data generated from Dirchlet multinomial distribution.
#' # First, use SimLenAndAgeFreqData function to calculate expected proportions at length
#' set.seed(123)
#' SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
#' MaxAge = 30
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = 0.2
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(500, 50) # L50, L95-L50 for retention
#' DiscMort = 0 # proportion of fish that die due to natural mortality
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.08
#' GrowthParams = c(Linf, vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # Now simulate length-frequency data from Dirchlet multinomial distribution.
#' set.seed(122)
#' nSampEvents = 50
#' nFishPerSampEvent = 20
#' theta_val = 0.3
#' midpt = Res$midpt
#' ExpPropAtLen = Res$ModelDiag$ExpRetCatchPropAtLen
#' res=SimLenFreqDat_DirMultDistn(nSampEvents, nFishPerSampEvent, theta_val, midpt, ExpPropAtLen)
#' plot(res$midpt, res$simLenFreq, "o")
#' # Fit length-based catch curve, assuming Dirichlet multinomial distribution
#' ObsRetCatchFreqAtLen = res$simLenFreq
#' ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
#' PropReleased = NA # proportion of fish released, vector including mean and sd (option probably now obselete)
#' midpt=Res$midpt
#' lbnd=Res$lbnd
#' ubnd=Res$ubnd
#' InitFishMort = 0.25 # specify starting parameters
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
#' InitL50 = 400
#' InitDelta = 100
#' DistnType=2
#' InitTheta = 0.3 # specify starting parameters
#' InitTheta_logit = log(InitTheta/(1-InitTheta)) # logit transform
#' params = c(InitFishMort_logit, log(InitL50), log(InitDelta),InitTheta_logit)
#' FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
#'                                           lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
#' @export
GetLengthBasedCatchCurveResults <- function (params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                             lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
{
  nlmb <- nlminb(params, CalcObjFunc_LengthBasedCatchCurve,
                 gradient = NULL, hessian = TRUE)
  hess.out = optimHess(nlmb$par, CalcObjFunc_LengthBasedCatchCurve)
  vcov.Params = solve(hess.out)
  ses = sqrt(diag(vcov.Params))

  # compute parameter correlation matrix
  cor.Params = NA
  if (length(params)>1) {
    temp = diag(1/sqrt(diag(vcov.Params)))
    cor.Params=temp %*% vcov.Params %*% temp
  }

  # Get Z estimate and 95 percent confidence limits
  temp = c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1]) # logit space
  EstFMort = 1/(1+exp(-temp)) # inverse logit transformed value
  EstL50_sel = NA; EstDelta_sel = NA
  EstL50_ret = NA; EstDelta_ret = NA
  RetCatch_EffSampleSize = NA; DiscCatch_EffSampleSize=NA

  if (DistnType == 1) { # multinomial distribution
    if (SelectivityType == 1) { # inputted as specified selectivity vector
      ParamEst = t(data.frame(FMort = round(EstFMort, 2)))
    }
    if (SelectivityType == 2) { # include estimates for logistic gear selectivity parameters
      if (length(params)==3) {
        EstL50_sel = exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
        EstDelta_sel = exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
        ParamEst = t(data.frame(FMort = round(EstFMort, 3), L50_sel = round(EstL50_sel, 3), Delta_sel = round(EstDelta_sel, 3)))
      }
      if (length(params)==5) { # include estimates for logistic retention parameters
        EstL50_sel = exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
        EstDelta_sel = exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
        EstL50_ret = exp(c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
        EstDelta_ret = exp(c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5]))
        ParamEst = t(data.frame(FMort = round(EstFMort, 3), L50_sel = round(EstL50_sel, 3), Delta_sel = round(EstDelta_sel, 3),
                                L50_ret = round(EstL50_ret, 3), Delta_ret = round(EstDelta_ret, 3)))
      }
    }
  }

  if (DistnType == 2) { # Dirichlet multinomial distribution
    if (SelectivityType == 1) { # inputted as specified selectivity vector
      temp = c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2])
      EstTheta = 1/(1+exp(-temp));
      ParamEst = t(data.frame(FMort = round(EstFMort, 2), Theta = round(EstTheta, 3)))

    }
    if (SelectivityType == 2) { # include estimates for logistic gear selectivity parameters
      if (length(params)==4) {
        EstL50_sel = exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
        EstDelta_sel = exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
        temp = c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4])
        EstTheta = 1/(1+exp(-temp));
        ParamEst = t(data.frame(FMort = round(EstFMort, 3), L50_sel = round(EstL50_sel, 3),
                                Delta_sel = round(EstDelta_sel, 3), Theta = round(EstTheta, 3)))

      }
      if (length(params)==6) { # include estimates for logistic retention parameters
        EstL50_sel = exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
        EstDelta_sel = exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
        EstL50_ret = exp(c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
        EstDelta_ret = exp(c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5]))
        temp = c(nlmb$par[6], nlmb$par[6] + c(-1.96, 1.96) * ses[6])
        EstTheta = 1/(1+exp(-temp));
        ParamEst = t(data.frame(FMort = round(EstFMort, 3), L50_sel = round(EstL50_sel, 3), Delta_sel = round(EstDelta_sel, 3),
                                L50_ret = round(EstL50_ret, 3), Delta_ret = round(EstDelta_ret, 3), Theta = round(EstTheta, 3)))

        # Dirichlet multinomial effective sample size - discarded catch
        temp = ((1+(EstTheta*sum(ObsDiscCatchFreqAtLen))) / (1 + EstTheta))
        DiscCatch_EffSampleSize <- t(data.frame(DiscCatch_EffSampleSize = temp))
        colnames(DiscCatch_EffSampleSize) = c("Estimate", "lw_95%CL", "up_95%CL")
      }
    }
    # Dirichlet multinomial effective sample size - retained catch
    temp = ((1+(EstTheta*sum(ObsRetCatchFreqAtLen))) / (1 + EstTheta))
    RetCatch_EffSampleSize <- t(data.frame(RetCatch_EffSampleSize = temp))
    colnames(RetCatch_EffSampleSize) = c("Estimate", "lw_95%CL", "up_95%CL")

  }
  colnames(ParamEst) = c("Estimate", "lw_95%CL", "up_95%CL")

  # store some diagnostic outputs from model
  params = nlmb$par
  CatchCurveType=1 #1=length-based, 2=age and length based
  res = AgeAndLengthBasedCatchCurvesCalcs(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge, CatchCurveType,
                                          MLL, SelectivityType, lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)

  # calculate approximate cv for lengths at max age (growth diagnostic), if single timestep
  CV_LenAtMaxAge_Results = GetMaxLenAtAgeCV_AgeAndLengthBasedCatchCurve(TimeStep, res, MaxAge)

  # get proportion and associated sd for randomly-generated data for retained fish
  if (!is.na(ObsDiscCatchFreqAtLen[1])) {
    ObsTotCatchAtLen = ObsRetCatchFreqAtLen + ObsDiscCatchFreqAtLen
  } else {
    ObsTotCatchAtLen = ObsRetCatchFreqAtLen
  }
  ObsPropRetFish = sum(ObsRetCatchFreqAtLen) / sum(ObsTotCatchAtLen)
  ObsPropRetFish_sd = sqrt(ObsPropRetFish * (1 - ObsPropRetFish)) / sum(ObsTotCatchAtLen)

  ModelDiag = list(midpt=res$midpt,
                   MeanSizeAtAge=res$MeanSizeAtAge,
                   MeanEndingLength=res$MeanEndingLength,
                   TimestepGrowthSizeInc=res$TimestepGrowthSizeInc,
                   RecLenDist=res$RecLenDist,
                   SelAtLength = res$SelAtLength,
                   RetAtLength = res$RetAtLength,
                   SelLandAtLength = res$SelLandAtLength,
                   GrowthModelType = res$GrowthModelType,
                   CV_LenAtMaxAge=CV_LenAtMaxAge_Results$CV_LenAtMaxAge,
                   FemCV_LenAtMaxAge=CV_LenAtMaxAge_Results$FemCV_LenAtMaxAge,
                   MalCV_LenAtMaxAge=CV_LenAtMaxAge_Results$MalCV_LenAtMaxAge,
                   ExpRetCatchPropInLenClass = res$ExpRetCatchPropInLenClass,
                   ExpRetCatchPropInLenClass_Fem = res$ExpRetCatchPropInLenClass_Fem,
                   ExpRetCatchPropInLenClass_Mal = res$ExpRetCatchPropInLenClass_Mal,
                   ExpDiscCatchPropInLenClass = res$ExpDiscCatchPropInLenClass,
                   ExpDiscCatchPropInLenClass_Fem = res$ExpDiscCatchPropInLenClass_Fem,
                   ExpDiscCatchPropInLenClass_Mal = res$ExpDiscCatchPropInLenClass_Mal,
                   ExpRetCatchPropLengthGivenIntAge = res$ExpRetCatchPropLengthGivenIntAge,
                   ExpRetCatchPropLengthGivenIntAge_Fem = res$ExpRetCatchPropLengthGivenIntAge_Fem,
                   ExpRetCatchPropLengthGivenIntAge_Mal = res$ExpRetCatchPropLengthGivenIntAge_Mal,
                   RetCatchAtLen = res$RetCatchAtLen,
                   TotCatchAtLen = res$TotCatchAtLen,
                   DiscCatchAtLen = res$DiscCatchAtLen,
                   FAtLen = res$FAtLen,
                   FAtLen_Fem = res$FAtLen_Fem,
                   FAtLen_Mal = res$FAtLen_Mal,
                   ZAtLen = res$ZAtLen,
                   ZAtLen_Fem = res$ZAtLen_Fem,
                   ZAtLen_Mal = res$ZAtLen_Mal,
                   FAtLenReten = res$FAtLenReten,
                   FAtLenReten_Fem = res$FAtLenReten_Fem,
                   FAtLenReten_Mal = res$FAtLenReten_Mal,
                   FAtLenDisc = res$FAtLenDisc,
                   FAtLenDisc_Fem = res$FAtLenDisc_Fem,
                   FAtLenDisc_Mal = res$FAtLenDisc_Mal)

  ResultsSummary = list(ParamEst = ParamEst,
                        convergence = nlmb$convergence,
                        nll = nlmb$objective,
                        RetCatch_SampleSize = sum(ObsRetCatchFreqAtLen),
                        DiscCatch_SampleSize = sum(ObsRetCatchFreqAtLen),
                        RetCatch_EffSampleSize = RetCatch_EffSampleSize,
                        DiscCatch_EffSampleSize = DiscCatch_EffSampleSize,
                        ObsPropRetFish = ObsPropRetFish,
                        ObsPropRetFish_sd = ObsPropRetFish_sd)

  Results = list(ParamEst = ParamEst,
                 EstFMort = EstFMort[1],
                 EstFMort_se = ((EstFMort[3]-EstFMort[2])/2)/1.96,
                 params = nlmb$par,
                 convergence = nlmb$convergence,
                 vcov.Params = vcov.Params,
                 cor.Params = cor.Params,
                 RetCatch_SampleSize = sum(ObsRetCatchFreqAtLen),
                 DiscCatch_SampleSize = sum(ObsDiscCatchFreqAtLen),
                 RetCatch_EffSampleSize = RetCatch_EffSampleSize,
                 DiscCatch_EffSampleSize = DiscCatch_EffSampleSize,
                 ResultsSummary = ResultsSummary,
                 ModelDiag = ModelDiag)

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
  Res = AgeAndLengthBasedCatchCurvesCalcs(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge, CatchCurveType,
                                          MLL, SelectivityType, lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)

  ExpRetCatchPropInLenClass = Res$ExpRetCatchPropInLenClass
  ExpDiscCatchPropInLenClass = Res$ExpDiscCatchPropInLenClass

  NLL_PropReleased = 0
  if (!is.na(PropReleased[1])) {
    EstPropReleased = sum(Res$DiscCatchAtLen) / sum(Res$TotCatchAtLen)
    NLL_PropReleased = -dnorm(EstPropReleased, PropReleased[1], PropReleased[2], log=TRUE)
  }

  if (DistnType==1) { # multinomial distribution
    NLL_RetCatch = CalcMultNLLMargLengthComposition(ObsRetCatchFreqAtLen, ExpRetCatchPropInLenClass)
    NLL_DiscCatch = 0
    if (!is.na(ObsDiscCatchFreqAtLen[1])) {
      NLL_DiscCatch = CalcMultNLLMargLengthComposition(ObsDiscCatchFreqAtLen, ExpDiscCatchPropInLenClass)
    }
  }
  if (DistnType==2) { # Dirichlet multinomial distribution

    if (SelectivityType == 1) { # inputted as specified selectivity vector
      temp = params[2]
      DM_theta = 1/(1+exp(-temp));
    }
    if (SelectivityType == 2) { # include estimates for logistic gear selectivity parameters
      if (length(params)==4) {
        temp = params[4]
        DM_theta = 1/(1+exp(-temp));
      }
      if (length(params)==6) { # include estimates for logistic retention parameters
        temp = params[6]
        DM_theta = 1/(1+exp(-temp));
      }
    }
    cat("CalcObjFunc_LengthBasedCatchCurve: DM_theta",DM_theta,'\n')
    NLL_RetCatch = CalcDirMultNLLMargLengthComposition(ObsRetCatchFreqAtLen, ExpRetCatchPropInLenClass, DM_theta)
    NLL_DiscCatch = 0
    if (!is.na(ObsDiscCatchFreqAtLen[1])) {
      NLL_DiscCatch = CalcDirMultNLLMargLengthComposition(ObsDiscCatchFreqAtLen, ExpDiscCatchPropInLenClass, DM_theta)
    }
  }

  NLL = NLL_RetCatch + NLL_PropReleased + NLL_DiscCatch + Res$L50_Pen + Res$L95_Pen

  cat("NLL", NLL, " NLL_DiscCatch ", NLL_DiscCatch, " NLL_PropReleased ", NLL_PropReleased,
      " FMort", 1/(1+exp(-params[1])), " L50_Pen ", Res$L50_Pen, " L95_Pen " ,Res$L95_Pen, '\n')

  return(NLL)

}


#' Get NLL for age and length-based catch curve
#'
#' @keywords internal
#' @param estimated parameters, including fishing mortality, growth and logistic selectivity parameters
#' @return negative log-likelihood (NLL)
CalcObjFunc_AgeAndLengthBasedCatchCurve <- function(params) {
  # get NLL for length based catch curve, for optimisation

  DistnType=1 # multinomial
  CatchCurveType=2
  GrowthParams = NA
  CVSizeAtAge = NA
  Res = AgeAndLengthBasedCatchCurvesCalcs(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge, CatchCurveType,
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
    Length_NLL = CalcMultNLLMargLengthComposition(ObsRetCatchFreqAtLen, ExpRetCatchPropInLenClass)

    # get NLL for age at length observations
    ExpRetCatchPropAtIntAge = as.matrix(Res$ExpRetCatchPropAtIntAge)
    ExpRetCatchPropLengthGivenIntAge = as.matrix(Res$ExpRetCatchPropLengthGivenIntAge)
    ExpRetCatchPropIntAgeGivenLength = CalcExpCatchPropIntAgeGivenLength_cpp(nLenCl, nAgeCl, ExpRetCatchPropLengthGivenIntAge, ExpRetCatchPropAtIntAge)

    if (TimeStep == 1) {
      ObsRetCatchFreqAtLengthAndIntAge = ObsRetCatchFreqAtLengthAndAge
    } else {
      ObsRetCatchFreqAtLengthAndIntAge = ConvertObsDataFromDecAgesToIntegerAges(TimeStep, MaxAge, nLenCl, ObsRetCatchFreqAtLengthAndAge)
    }
    CondAgeAtLengthNLL = CalcNLLCondAgeAtLength_cpp(nLenCl, nAgeCl, ObsRetCatchFreqAtLengthAndIntAge, ExpRetCatchPropIntAgeGivenLength)
  }
  if (GrowthModelType == 2 | GrowthModelType == 4) {

    # get NLL for marginal length composition
    # females
    ExpRetCatchPropInLenClass_Fem = Res$ExpRetCatchPropInLenClass_Fem
    ObsRetCatchFreqAtLen_Fem = unlist(ObsRetCatchFreqAtLen[1,])
    Length_NLL_Fem = CalcMultNLLMargLengthComposition(ObsRetCatchFreqAtLen_Fem, ExpRetCatchPropInLenClass_Fem)
    # males
    ExpRetCatchPropInLenClass_Mal = Res$ExpRetCatchPropInLenClass_Mal
    ObsRetCatchFreqAtLen_Mal = unlist(ObsRetCatchFreqAtLen[2,])
    Length_NLL_Mal = CalcMultNLLMargLengthComposition(ObsRetCatchFreqAtLen_Mal, ExpRetCatchPropInLenClass_Mal)
    Length_NLL = Length_NLL_Fem + Length_NLL_Mal

    # get NLL for age at length observations
    # females
    ExpRetCatchPropAtIntAge_Fem = as.matrix(Res$ExpRetCatchPropAtIntAge_Fem)
    ExpRetCatchPropLengthGivenIntAge_Fem = as.matrix(Res$ExpRetCatchPropLengthGivenIntAge_Fem)
    ObsRetCatchFreqAtLengthAndAge_Fem = unlist(ObsRetCatchFreqAtLengthAndAge[,,1])
    if (TimeStep == 1) {
      ObsRetCatchFreqAtLengthAndIntAge_Fem = ObsRetCatchFreqAtLengthAndAge_Fem
    } else {
      ObsRetCatchFreqAtLengthAndIntAge_Fem = ConvertObsDataFromDecAgesToIntegerAges(TimeStep, MaxAge, nLenCl, ObsRetCatchFreqAtLengthAndAge_Fem)
    }
    ExpRetCatchPropIntAgeGivenLength_Fem = CalcExpCatchPropIntAgeGivenLength_cpp(nLenCl, nAgeCl, ExpRetCatchPropLengthGivenIntAge_Fem, ExpRetCatchPropAtIntAge_Fem)
    CondAgeAtLengthNLL_Fem = CalcNLLCondAgeAtLength_cpp(nLenCl, nAgeCl, ObsRetCatchFreqAtLengthAndIntAge_Fem, ExpRetCatchPropIntAgeGivenLength_Fem)
    # males
    ExpRetCatchPropAtIntAge_Mal = as.matrix(Res$ExpRetCatchPropAtIntAge_Mal)
    ExpRetCatchPropLengthGivenIntAge_Mal = as.matrix(Res$ExpRetCatchPropLengthGivenIntAge_Mal)
    ObsRetCatchFreqAtLengthAndAge_Mal = unlist(ObsRetCatchFreqAtLengthAndAge[,,2])
    if (TimeStep == 1) {
      ObsRetCatchFreqAtLengthAndIntAge_Mal = ObsRetCatchFreqAtLengthAndAge_Mal
    } else {
      ObsRetCatchFreqAtLengthAndIntAge_Mal = ConvertObsDataFromDecAgesToIntegerAges(TimeStep, MaxAge, nLenCl, ObsRetCatchFreqAtLengthAndAge_Mal)
    }
    ExpRetCatchPropIntAgeGivenLength_Mal = CalcExpCatchPropIntAgeGivenLength_cpp(nLenCl, nAgeCl, ExpRetCatchPropLengthGivenIntAge_Mal, ExpRetCatchPropAtIntAge_Mal)
    CondAgeAtLengthNLL_Mal = CalcNLLCondAgeAtLength_cpp(nLenCl, nAgeCl, ObsRetCatchFreqAtLengthAndIntAge_Mal, ExpRetCatchPropIntAgeGivenLength_Mal)

    CondAgeAtLengthNLL = CondAgeAtLengthNLL_Fem + CondAgeAtLengthNLL_Mal
  }

  NLL = Length_NLL + CondAgeAtLengthNLL + L50_Pen + L95_Pen
  cat("NLL", NLL, "Length_NLL", Length_NLL, "CondAgeAtLengthNLL", CondAgeAtLengthNLL, "L50_Pen", L50_Pen,
      "L95_Pen",L95_Pen, "\n")
  cat("F ", 1/(1+exp(-params[1])), " other params ",exp(params[2:length(params)]), "\n")
  cat("", "\n")

  return(NLL)

}

#' Get NLL for marginal length composition data, assuming Dirichlet multinomial distribution
#'
#' @keywords internal
#'
#' @param DM_theta parameter theta from Dirichlet multinomial distribution
#' @param ObsRetCatchFreqAtLen observed frequencies of fish in length classes
#' @param ExpCatchPropInLenClass observed catch proportions in length classes
#'
#' @return negative log-likelihood,assuming Dirichlet multinomial distribution (Length_NLL)
CalcDirMultNLLMargLengthComposition <- function(ObsRetCatchFreqAtLen, ExpRetCatchPropInLenClass, DM_theta) {

  SampleSize = sum(ObsRetCatchFreqAtLen)
  nLenCl = length(ObsRetCatchFreqAtLen)
  ObsPropAtLen = ObsRetCatchFreqAtLen/SampleSize
  sum1 = 0; sum2 = 0; NLL = 0
  for (t in 1:nLenCl) {
    sum1 = sum1 + lgamma(SampleSize * ObsPropAtLen[t] + 1)
    sum2 = sum2 + (lgamma(SampleSize * ObsPropAtLen[t] + DM_theta * SampleSize * ExpRetCatchPropInLenClass[t])
                   - lgamma(DM_theta * SampleSize * ExpRetCatchPropInLenClass[t]))
  }
  Length_NLL = -(lgamma(SampleSize+1) - sum1 + (lgamma(DM_theta * SampleSize) - lgamma(SampleSize + DM_theta * SampleSize)) + sum2)

  return(Length_NLL)
}


#' Get NLL for marginal length composition data, assuming multinomial distribution
#'
#' @keywords internal
#'
#' @param ObsRetCatchFreqAtLen observed frequencies of fish in length classes
#' @param ExpCatchPropInLenClass observed catch proportions in length classes
#'
#' @return negative log-likelihood, assuming multinomial distribution (Length_NLL)
CalcMultNLLMargLengthComposition <- function(ObsRetCatchFreqAtLen, ExpRetCatchPropInLenClass) {

  Length_NLL = -sum(ObsRetCatchFreqAtLen * log(ExpRetCatchPropInLenClass + 1E-4))

  return(Length_NLL)
}

#' Get selectivity parameter values for age and length based catch curve
#'
#' Get selectivity parameter values for age and length based catch curve given specified catch curve, growth and selectivity model options
#'
#' @keywords internal
#' @param params estimated model parameters (varies, depending on growth curve type, selectivity type and catch curve type)
#' @param CatchCurveType 1=length based, 2=age and length based
#' @param SelectivityType 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#'
#' @return values of growth parameters given specified catch curve, growth and selectivity model option
GetSelectParams_AgeAndLengthBasedCatchCurvesCalcs <- function(Res, params, CatchCurveType, SelectivityType) {

  # get selectivity parameters for length based or length and age based catch curve model
  L50 = NA; L95 = NA
  L50_ret = NA; L95_ret = NA

  if (SelectivityType == 2) { # logistic selectivity
    L50 = exp(params[2])
    L95 = exp(params[3])
    if (CatchCurveType==1) { # length-based catch curve
      if (length(params)==5) { # estimating logistic selectiviy and logistic retention (using discard and retained catch data)
        L50_ret = exp(params[4])
        L95_ret = exp(params[5])
      }
    }
  }

  SelParams = c(L50, L95)
  RetenParams = c(L50_ret, L95_ret)

  result = list(SelParams = SelParams,
                RetenParams = RetenParams)

  return(result)
}


#' Get growth parameter values for age and length based catch curve
#'
#' Get growth parameter values for age and length based catch curve given specified catch curve, growth and selectivity model options
#'
#' @keywords internal
#' @param params estimated model parameters (varies, depending on growth curve type, selectivity type and catch curve type)
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
#' @param SelectivityType 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#'
#' @return list including GrowthParams (von Bertalanffy or Schnute growth parameters), GrowthModelType
#' 1=von Bertalanffy 1 sex, 2= von Bertalanffy 2 sexes, 3=Schnute 1 sex, 4=Schnute 2 sexes, and
#' CVSizeAtAge for 1 or 2 sex model
GetGrowthParams_AgeAndLengthBasedCatchCurvesCalcs <- function(params, GrowthCurveType, SelectivityType) {

  # Determine catch curve model type, selectivity type, growth curve type and number of sexes
  if (GrowthCurveType == 1) { # von Bertalanffy
    if (SelectivityType == 1) {
      if (length(params)==4) { # single sex, selectivity vector
        Linf = exp(params[2])
        vbK = exp(params[3])
        GrowthParams = c(Linf,vbK)
        CVSizeAtAge = exp(params[4])
        GrowthModelType = 1
      }
      if (length(params)==6) { # 2 sexes, selectivity vector
        Linf = exp(params[2:3])
        vbK = exp(params[4:5])
        GrowthParams = data.frame(Linf=Linf,vbK=vbK)
        CVSizeAtAge = c(exp(params[6]),exp(params[6]))
        GrowthModelType = 2
      }
    }
    if (SelectivityType==2) {
      if (length(params)==6) { # single sex, logistic selectivity
        Linf = exp(params[4])
        vbK = exp(params[5])
        GrowthParams = c(Linf,vbK)
        CVSizeAtAge = exp(params[6])
        GrowthModelType = 1
      }

      if (length(params)==8) { # 2 sexes, logistic selectivity
        Linf = exp(params[4:5])
        vbK = exp(params[6:7])
        GrowthParams = data.frame(Linf=Linf,vbK=vbK)
        CVSizeAtAge = c(exp(params[8]),exp(params[8]))
        GrowthModelType = 2
      }
    }
  } # von Bertalanffy

  if (GrowthCurveType == 2) { # Schnute
    if (SelectivityType==1) {
      if (length(params)==5) { # single sex, selectivity vector
        y1 = 0
        y2 = exp(params[2])
        a = params[3]
        b = params[4]
        GrowthParams = c(y1,y2,a,b)
        CVSizeAtAge = exp(params[5])
        GrowthModelType = 3
      }
      if (length(params)==8) { # 2 sexes, selectivity vector
        y1 = c(0,0)
        y2 = exp(params[2:3])
        a = params[4:5]
        b = params[6:7]
        GrowthParams = data.frame(y1=y1,y2=y2,a=a,b=b)
        CVSizeAtAge = c(exp(params[8]),exp(params[8]))
        GrowthModelType = 4
      }
    }
    if (SelectivityType==2) {
      if (length(params)==7) { # single sex, logistic selectivity
        y1 = 0
        y2 = exp(params[4])
        a = params[5]
        b = params[6]
        GrowthParams = c(y1,y2,a,b)
        CVSizeAtAge = exp(params[7])
        GrowthModelType = 3
      }
      if (length(params)==9) { # 2 sexes, logistic selectivity
        y1 = c(0,0)
        y2 = exp(params[4:5])
        a = params[6:7]
        b = params[8:9]
        GrowthParams = data.frame(y1=y1,y2=y2,a=a,b=b)
        CVSizeAtAge = c(exp(params[9]),exp(params[9]))
        GrowthModelType = 4
      }
    }
  } # Schnute



  result = list(GrowthParams = GrowthParams,
                GrowthModelType = GrowthModelType,
                CVSizeAtAge = CVSizeAtAge)

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
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
#' @param GrowthModelType 1=von Bertalanffy - 1 sex, 2=von Bertalanffy - 2 sexes,
#' Schnute - 1 sex, 2=Schnute - 2 sexes
#' @param GrowthParams growth parameters (used for GrowthModelType 2 or 4)
#' @param RefnceAges Schnute reference ages (used for GrowthModelType 4)
#' @param SelectivityType 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#'
#' @return mean size at age from growth curve (MeanSizeAtAge), estimated length after one year for each length class
#' mid point, given growth curve (MeanEndingLength) and associated change in length after one year of growth (TimestepGrowthSizeInc)
GetGrowthInputsForLengthTransitionMatrices <- function(MaxAge, TimeStep, nLenCl, midpt, GrowthCurveType, GrowthModelType, GrowthParams, RefnceAges, SelectivityType) {

  # growth calcs for single sex catch curve (or using combined growth curve for both sexes)
  Ages = seq(TimeStep,MaxAge,TimeStep)
  t1=NA; t2=NA; y1=NA; y2=NA; a=NA; b=NA
  Linf = NA

  if (GrowthModelType == 1) { # von Bertalanffy, single sex
    Linf = GrowthParams[1]
    vbK = GrowthParams[2]
    MeanSizeAtAge = Linf * (1 - exp (-vbK * Ages))
    MeanEndingLength = midpt + (Linf - midpt) * (1 - exp(-vbK*TimeStep))
    TimestepGrowthSizeInc = MeanEndingLength-midpt # amount of annual growth with respect to initial length
  }
  if (GrowthModelType == 2) { # von Bertalanffy, separate sex
    Linf = GrowthParams[,1]
    vbK = GrowthParams[,2]
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
    MeanEndingLength = CalcLengthAfterGrowthForTimetep(GrowthCurveType=2, TimeStep, GrowthParamsForSex, RefnceAgesForSex, midpt, MaxAge)
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
      MeanEndingLength[i,] = CalcLengthAfterGrowthForTimetep(GrowthCurveType=2, TimeStep, GrowthParamsForSex, RefnceAgesForSex, midpt, MaxAge)
      TimestepGrowthSizeInc[i,] = MeanEndingLength[i,] - midpt
    }
  }

  result = list(Ages = Ages,
                MeanSizeAtAge = MeanSizeAtAge,
                MeanEndingLength = MeanEndingLength,
                TimestepGrowthSizeInc = TimestepGrowthSizeInc,
                t1=t1,
                t2=t2,
                y1=y1,
                y2=y2,
                a=a,
                b=b,
                Linf=Linf)

  return(result)
}


#' Allocate number for alternative length based catch curve models
#'
#' Allocate number for alternative single and separate sex length based catch curve models,
#' and growth curve types (1=von Beralanffy, 2=Schnute)
#'
#' @keywords internal
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
#' @param GrowthParams (von Bertalanffy or Schnute growth parameters)
#'
#' @return GrowthModelType 1=length-based von Bertalanffy, 2=length-based Schnute,
#' 3=age and length based von Bertalanffy, 4=age and length based Schnute
GetGrowthModelType <- function (GrowthCurveType, GrowthParams)
{
  # Allocated a model number, given type of catch curve and growth curve
  if (GrowthCurveType==1) { # von Bertalanffy
    if (is.vector(GrowthParams)) {
      GrowthModelType =  1
    } else {
      GrowthModelType = 2
    }
  }
  if (GrowthCurveType==2) { # Schnute
    if (is.vector(GrowthParams)) {
      GrowthModelType = 3
    } else {
      GrowthModelType = 4
    }
  }

  Results = GrowthModelType

  return(Results)

}

#' Calculate asymptotic length from Schnute growth parameters
#'
#' Calculate asymptotic length from Schnute growth parameters, for use in penalty function
#' for selectivity parameters
#'
#' @keywords internal
#' @param GrowthModelType 1=length-based von Bertalanffy, 2=length-based Schnute,
#' 3=age and length based von Bertalanffy, 4=age and length based Schnute
#' @param Res Schnute parameters (output from GetGrowthModelType function)
#' @param GrowthParams Schnute growth parameters for 2 sex catch curve models
#' @param RefnceAges Schnute growth reference ages for 2 sex catch curve models
#'
#' @return Linf asymptotic length parameter
CalcAsymptoticLengthFromSchnuteParams <- function (GrowthModelType, Res, GrowthParams, RefnceAges)
{

  if (GrowthModelType == 3) { # Schnute single sex
    t1=Res$t1; t2=Res$t2; y1=Res$y1; y2=Res$y2; a=Res$a; b=Res$b
    if (a <= 0) a=0.01
    Linf = ((exp(a*t2)*y2^b-exp(a*t1)*y1^b)/(exp(a*t2)-exp(a*t1)))^(1/b)
    # cat("CalcAsymptoticLengthFromSchnuteParams2: a = ",a,"b = ",b, "y1 =",y1, "y2 =",y2,"Linf = ",Linf,'\n')
  }
  if (GrowthModelType == 4) { # Schnute separate sexes
    t1=RefnceAges[,1]; t2=RefnceAges[,2]; y1=GrowthParams[,1]; y2=GrowthParams[,2];
    a=GrowthParams[,3]; b=GrowthParams[,4]
    Linf1 = ((exp(a[1]*t2[1])*y2[1]^b[1]-exp(a[1]*t1[1])*y1[1]^b[1])/(exp(a[1]*t2[1])-exp(a[1]*t1[1])))^(1/b[1])
    Linf2 = ((exp(a[2]*t2[2])*y2[2]^b[2]-exp(a[2]*t1[2])*y1[2]^b[2])/(exp(a[2]*t2[2])-exp(a[2]*t1[2])))^(1/b[2])
    Linf = max(Linf1,Linf2)
  }
  return(Linf)
}

#' Calculate penalties for selectivity and reset parameters values when necessary
#'
#' This function ensures that the length at 95% selectivity does not exceed Linf, and
#' that this length is greater than the length at 50% selectivity by at least 2 mm
#'
#' @keywords internal
#' @param Res list of outputs from GetGrowthInputsForLengthTransitionMatrices function
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
#' @param GrowthModelType # 1=1 sex, von Bertalanffy, 2=2 sexes, von Bertalanffy, 3=1 sex, Schnute, 4=2 sexes, Schnute
#' @param GrowthParams input growth parameters
#' @param Linf input growth reference ages
#' @param SelectivityType 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' @param SelParams selectivity parameters L50, L95-L50
#'
#' @return Linf asymptotic length parameter
CalcSelectivityPenalties <- function (Res, GrowthCurveType, GrowthModelType, GrowthParams, RefnceAges, SelectivityType, SelParams) {

  if (GrowthCurveType==1) { # von Bertalanffy
    if (GrowthModelType == 1) { # single sex, von Bertalanffy
      Linf = Res$Linf
    } else { # 2 sexes, von Bertalanffy
      Linf = max(Res$Linf)
    }
  }

  if (GrowthCurveType==2) { # Schnute
    # calc Linf from Schnute growth curve params
    Linf = CalcAsymptoticLengthFromSchnuteParams(GrowthModelType, Res, GrowthParams, RefnceAges)
  }

  # ensure L95 < Linf
  L95_Pen = 0; L50_Pen = 0
  if (SelectivityType == 2) { # logistic selectivity
    # cat("CalcSelectivityPenalties: GrowthParams = ",GrowthParams,'\n')
    # cat("CalcSelectivityPenalties: SelParams = ",SelParams,'\n')
    # cat("CalcSelectivityPenalties: Linf = ",Linf,'\n')
    if (sum(SelParams) > max(Linf)) {
      L95_Pen = (sum(SelParams) - max(Linf))^2
      SelParams[2] = max(Linf) - SelParams[1]
    }
    # calculate L50 penalty
    if (SelParams[1] > (sum(SelParams) - 2)) {
      L50_Pen = ((sum(SelParams) - 2) - SelParams[1])^2
      SelParams[1] = sum(SelParams) - 2
    }
  }

  Results = list(L95_Pen=L95_Pen,
                 L50_Pen=L50_Pen,
                 SelParams=SelParams)

  return(Results)

}

#' Calculations for length and age and lengths-based catch curves
#'
#' This function undertakes a range of calculations required by the length and age and lengths-based catch curves
#'
#' @keywords internal
#'
#' @param params vector of model parameters in log space (params) to be estimated
#' @param DistnType 1 = Multinomial, 2 = Dirichlet multinomial
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
AgeAndLengthBasedCatchCurvesCalcs <- function (params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge, CatchCurveType,
                                               MLL, SelectivityType, lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
{

  nTimeSteps = length(seq(TimeStep,MaxAge,TimeStep))
  Ages = seq(TimeStep, MaxAge, TimeStep)
  AgeCl = floor(round(Ages,6))
  nAgeCl = length(unique(AgeCl))
  MinAge = min(AgeCl)

  # get GrowthModelType
  if (CatchCurveType==1) { # length based
    GrowthModelType = GetGrowthModelType(GrowthCurveType, GrowthParams)
  }

  if (CatchCurveType==2) { # length and age based
    res2=GetGrowthParams_AgeAndLengthBasedCatchCurvesCalcs(params, GrowthCurveType, SelectivityType)
    GrowthModelType = res2$GrowthModelType
    GrowthParams = res2$GrowthParams
    CVSizeAtAge = res2$CVSizeAtAge
    # cat("AgeAndLengthBasedCatchCurvesCalcs: GrowthParams =", GrowthParams,"\n")
  }

  # get key inputs for length transition matrices
  nLenCl = length(midpt)
  Res = GetGrowthInputsForLengthTransitionMatrices(MaxAge, TimeStep, nLenCl, midpt, GrowthCurveType, GrowthModelType,
                                                   GrowthParams, RefnceAges, SelectivityType)
  MeanSizeAtAge = Res$MeanSizeAtAge
  TimestepGrowthSizeInc = Res$TimestepGrowthSizeInc
  MeanEndingLength = Res$MeanEndingLength

  # get parameters for specified growth curve and catch curve type
  ParamRes = GetSelectParams_AgeAndLengthBasedCatchCurvesCalcs(Res, params, CatchCurveType, SelectivityType)
  SelParams = ParamRes$SelParams; RetenParams = ParamRes$RetenParams

  # ensure L95 < Linf and that L95>L50+2
  ParamPenRes = CalcSelectivityPenalties(Res, GrowthCurveType, GrowthModelType, GrowthParams, RefnceAges, SelectivityType, SelParams)
  L95_Pen = ParamPenRes$L95_Pen; L50_Pen = ParamPenRes$L50_Pen; SelParams = ParamPenRes$SelParams

  # get selectivity and retention
  SelRes = GetSelectivityAndRetention(midpt, SelectivityType, SelParams, SelectivityVec, RetenParams, MLL)
  SelAtLength = SelRes$SelAtLength; RetAtLength = SelRes$RetAtLength

  # get size distribution for juvenile recruits
  RecLenDist = CalcSizeDistOfRecruits(MeanSizeAtAge, CVSizeAtAge, lbnd, ubnd, midpt, nLenCl) # length distribution of 1+ recruits

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
    FAtLenCapt = CatchCurveResults$FAtLenCapt
    FAtLenReten = CatchCurveResults$FAtLenReten
    FAtLenDisc = CatchCurveResults$FAtLenDisc
    FAtLen_Fem = NA; FAtLen_Mal = NA
    ZAtLen_Fem = NA; ZAtLen_Mal = NA
    FAtLenCapt_Fem = NA; FAtLenCapt_Mal = NA
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
    CVSizeAtAge_Fem = CVSizeAtAge[1]
    LTM_Fem = CalcLTM_cpp(TimestepGrowthSizeInc_Fem, CVSizeAtAge_Fem, lbnd, midpt, ubnd, nLenCl) # length-transition matrix - females

    TimestepGrowthSizeInc_Mal = as.vector(unlist(TimestepGrowthSizeInc[2,]))
    CVSizeAtAge_Mal = CVSizeAtAge[2]
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
    FAtLenCapt_Fem = CatchCurveResults_Fem$FAtLenCapt
    FAtLenReten_Fem = CatchCurveResults_Fem$FAtLenReten
    FAtLenDisc_Fem = CatchCurveResults_Fem$FAtLenDisc
    FAtLen_Mal = CatchCurveResults_Mal$FAtLen
    ZAtLen_Mal = CatchCurveResults_Mal$ZAtLen
    FAtLenCapt_Mal = CatchCurveResults_Mal$FAtLenCapt
    FAtLenReten_Mal = CatchCurveResults_Mal$FAtLenReten
    FAtLenDisc_Mal = CatchCurveResults_Mal$FAtLenDisc
    FAtLen = NA
    ZAtLen = NA
    FAtLenCapt = NA
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
                 FAtLenCapt = FAtLenCapt,
                 FAtLenCapt_Fem = FAtLenCapt_Fem,
                 FAtLenCapt_Mal = FAtLenCapt_Mal,
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
#' It provides several statistical outputs including convergence statistics, parameter estimates
#' and associated 95 percent confidence limits and associated variance-covariance matrix, calculated using
#' the MASS package.
#'
#' @param params vector of model parameters (in log space) to be estimated
#' @param RefnceAges Reference ages for Schnute growth curve (set to NA for von Bertalanffy growth curve)
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
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
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified is knife-edged at MLL
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.05
#' RefnceAges = NA
#' GrowthParams = c(Linf, vbK, CVSizeAtAge)
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # 2 sexes, von Bertalanffy
#' # GrowthCurveType = 1 # 1 = von Bertalanffy
#' # Linf = c(700,850)
#' # vbK = c(0.25,0.2)
#' # CVSizeAtAge = c(0.05,0.05)
#' # RefnceAges = NA
#' # GrowthParams = data.frame(Linf=Linf, vbK=vbK, CVSizeAtAge=CVSizeAtAge)
#' # Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#' #                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' lbnd=Res$lbnd
#' midpt=Res$midpt
#' ubnd=Res$ubnd
#' # get data - 1 sex (or combined sexes)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen # 1 sex
#' ObsRetCatchFreqAtLengthAndAge = as.matrix(Res$ObsRetCatchFreqAtLengthAndDecAge) # 1 sex
#' # # get data - 2 sexes
#' # ObsRetCatchFreqAtLen <- data.frame(matrix(nrow = 2, ncol = length(midpt))) # 2 sex
#' # colnames(ObsRetCatchFreqAtLen) <- midpt
#' # ObsRetCatchFreqAtLen[1,] = Res$ObsRetCatchFreqAtLen_Fem
#' # ObsRetCatchFreqAtLen[2,] = Res$ObsRetCatchFreqAtLen_Mal
#' # ObsRetCatchFreqAtLengthAndAge = array(c(unlist(Res$ObsRetCatchFreqAtLengthAndDecAge_Fem), unlist(Res$ObsRetCatchFreqAtLengthAndDecAge_Mal)),
#' #                                       c(nTimeSteps, length(midpt), 2), dimnames=list(rownames(Res$ObsRetCatchFreqAtLengthAndDecAge_Fem),
#' #                                                                                      colnames(Res$ObsRetCatchFreqAtLengthAndDecAge_Fem)))
#' # get params - 1 sex
#' InitFishMort = 0.3 # specify starting parameters
#' InitL50 = 320
#' InitDelta = 50 # L95-L50
#' InitLinf = 800
#' InitvbK = 0.2
#' InitCVSizeAtAge = 0.05
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' params = c(InitFishMort_logit, log(c(InitL50, InitDelta, InitLinf, InitvbK, InitCVSizeAtAge)))
#' # # get params - 2 sexes
#' # InitFishMort = 0.3 # specify starting parameters
#' # InitL50 = 320
#' # InitDelta = 50 # L95-L50
#' # InitLinf = c(800,800)
#' # InitvbK = c(0.25,0.25)
#' # InitCVSizeAtAge = 0.05
#' # InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' # params = c(InitFishMort_logit, log(c(InitL50, InitDelta, InitLinf, InitvbK, InitCVSizeAtAge)))
#' FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
#' # Example with specified selectivity vector
#' # Simulate data
#' SelectivityType=1 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = 1 / (1 + exp(-log(19)*(midpt-400)/(500-400)))
#' SelParams = c(NA, NA) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.05
#' GrowthParams = c(Linf, vbK, CVSizeAtAge)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' lbnd=Res$lbnd
#' midpt=Res$midpt
#' ubnd=Res$ubnd
#' # get data - 1 sex (or combined sexes)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen # 1 sex
#' ObsRetCatchFreqAtLengthAndAge = as.matrix(Res$ObsRetCatchFreqAtLengthAndDecAge) # 1 sex
#' # get params - 1 sex
#' InitFishMort = 0.3 # specify starting parameters
#' InitLinf = 800
#' InitvbK = 0.2
#' InitCVSizeAtAge = 0.05
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' params = c(InitFishMort_logit, log(c(InitLinf, InitvbK, InitCVSizeAtAge)))
#' FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
#' ########################################################################################################
#' # Examples below for when data are available in 'raw' form, i.e. as lengths and ages for individual fish
#' # single (or combined) sexes, using data for individual fish
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
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified is knife-edged at MLL
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.05
#' RefnceAges = NA
#' GrowthParams = c(Linf, vbK, CVSizeAtAge)
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # get individual fish data
#' FishAges = Res$ObsDecAgeRetCatch # ages for individual fish
#' FishLengths = Res$ObsRandLenRetCatch # length for individual fish
#' lbnd = Res$lbnd
#' midpt = Res$midpt
#' ubnd = Res$ubnd
#' LenInc = ubnd[1] - lbnd[1]
#' # single sex - get length frequency
#' lbns = trunc(FishLengths/LenInc)*LenInc
#' ObsRetCatchFreqAtLen  = as.vector(table(factor(lbns, levels=lbnd)))
#' # single sex - get frequencies in each length and age class
#' ObsRetCatchFreqAtLengthAndAge <- data.frame(matrix(nrow = MaxAge, ncol = length(lbnd)))
#' colnames(ObsRetCatchFreqAtLengthAndAge) <- lbnd
#' ObsRetCatchFreqAtLengthAndAge = as.matrix(table(factor(FishAges, levels=1:MaxAge),
#'                                                 factor(lbns, levels=c(lbnd))))
#' # fit model
#' InitFishMort = 0.3 # specify starting parameters
#' InitL50 = 320
#' InitDelta = 50 # L95-L50
#' InitLinf = 800
#' InitvbK = 0.2
#' InitCVSizeAtAge = 0.05
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' params = c(InitFishMort_logit, log(c(InitL50, InitDelta, InitLinf, InitvbK, InitCVSizeAtAge)))
#' FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
#' # 2 sexes, using data for individual fish
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
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified is knife-edged at MLL
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0
#' # 2 sexes, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy
#' Linf = c(700,850)
#' vbK = c(0.25,0.2)
#' CVSizeAtAge = c(0.05,0.05)
#' RefnceAges = NA
#' GrowthParams = data.frame(Linf=Linf, vbK=vbK, CVSizeAtAge=CVSizeAtAge)
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # get individual fish data
#' FishAges = c(Res$ObsDecAgeRetCatch_Fem, Res$ObsDecAgeRetCatch_Mal) # ages for individual fish
#' FishLengths = c(Res$ObsRandLenRetCatch_Fem, Res$ObsRandLenRetCatch_Mal) # length for individual fish
#' lbns = trunc(FishLengths/LenInc)*LenInc
#' FishSex = c(rep(1,length(Res$ObsDecAgeRetCatch_Fem)), rep(2, length(Res$ObsDecAgeRetCatch_Mal)))
#' lbnd = Res$lbnd
#' midpt = Res$midpt
#' ubnd = Res$ubnd
#' LenInc = ubnd[1] - lbnd[1]
#' DecAges=Res$DecAges
#' # 2 sexes - get length frequency
#' ObsRetCatchFreqAtLen <- data.frame(matrix(nrow = 2, ncol = length(midpt))) # 2 sex
#' colnames(ObsRetCatchFreqAtLen) <- midpt
#' ObsRetCatchFreqAtLen = as.matrix(table(factor(FishSex, levels=c(1:2)),
#'                                        factor(lbns, levels=c(lbnd))))
#' # 2 sexes - get frequencies in each length and age class
#' ObsRetCatchFreqAtLengthAndAge <- as.array(table(factor(FishAges, levels=DecAges),
#'                                                 factor(lbns, levels=c(lbnd)),
#'                                                 factor(FishSex, levels=c(1:2))))
#' # fit model
#' # get params - 2 sexes
#' InitFishMort = 0.3 # specify starting parameters
#' InitL50 = 320
#' InitDelta = 50 # L95-L50
#' InitLinf = c(800,800)
#' InitvbK = c(0.25,0.25)
#' InitCVSizeAtAge = 0.05
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' params = c(InitFishMort_logit, log(c(InitL50, InitDelta, InitLinf, InitvbK, InitCVSizeAtAge)))
#' FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
#' ############################################################
#' # Fit single (or combined) sex model using Schnute growth curve
#' # Simulate data
#' set.seed(123)
#' SampleSize=1000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
#' MaxAge = 30
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = 0.2
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0 # proportion of fish that die due to natural mortality
#' # 1 sex, Schnute
#' GrowthCurveType = 2 # 1 = von Bertalanffy, 2 = Schnute
#' t1 = 0 # growth - Schnute
#' t2 = 25 # growth - Schnute
#' y1 = 0 # growth - Schnute
#' y2 = 1000 # growth - Schnute
#' a = 0.1 # growth - Schnute
#' b = 2 # growth - Schnute
#' GrowthParams = c(y1, y2, a, b) # ensure t1 and y1 are fixed at zero (not possible to estimate with this model)
#' RefnceAges = c(t1,t2)
#' CVSizeAtAge = 0.08
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' plot(Res$ObsDecAgeRetCatch, Res$ObsRandLenRetCatch, ylim=c(0,1200))
#' # fit model
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen
#' ObsRetCatchFreqAtLengthAndAge = Res$ObsRetCatchFreqAtLengthAndDecAge
#' lbnd = Res$lbnd
#' midpt = Res$midpt
#' ubnd = Res$ubnd
#' GrowthCurveType=2 # Schnute
#' InitFishMort = 0.2 # specify starting parameters
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' InitL50 = 400
#' InitDelta = 20 # L95-L50
#' RefnceAges = c(0,20)
#' y2=900; a=0.15; b=1.5; InitCVSizeAtAge = 0.1
#' params = c(InitFishMort_logit, log(InitL50), log(InitDelta), log(y2), a, b, log(InitCVSizeAtAge))
#' FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
#' # if wanting to obtain estimated lengths-at-ages over larger age range,
#' # extrapolating beyond MaxAge (using WAFishBiology package)
#' library(WAFishBiology)
#' MaxModelAge <- 50
#' t1=0 # RefnceAges[1]
#' t2=RefnceAges[2]
#' y2=FittedRes$ParamEst[4,1]
#' a=FittedRes$ParamEst[5,1]
#' b=FittedRes$ParamEst[6,1]
#' EstLenAtAge <- rep(0,MaxModelAge)
#' for (Age in 1:MaxModelAge) {
#'   EstLenAtAge[Age]=SchnuteGrowthfunction(Age, t1=0, t2, y1=0, y2, a, b)
#' }
#' plot(1:MaxModelAge, EstLenAtAge, "o")
#' @export
GetAgeAndLengthBasedCatchCurveResults <- function (params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                                   lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
{

  nlmb <- nlminb(params, CalcObjFunc_AgeAndLengthBasedCatchCurve,
                 gradient = NULL, hessian = TRUE, control=list(eval.max=1000, iter.max=1000))


  hess.out = optimHess(nlmb$par, CalcObjFunc_AgeAndLengthBasedCatchCurve)
  vcov.Params = solve(hess.out)
  ses = sqrt(diag(vcov.Params))

  # compute parameter correlation matrix
  temp = diag(1/sqrt(diag(vcov.Params)))
  cor.Params=temp %*% vcov.Params %*% temp


  # Get estimates of parameters and their 95 percent confidence limits
  ParamEst = GetParamRes_AgeAndLengthBasedCatchCurve(GrowthCurveType, SelectivityType, params, nlmb, ses)

  # store some diagnostic outputs from model
  DistnType=1 # multinomial
  CatchCurveType=2
  GrowthParams = NA
  CVSizeAtAge = NA
  params=nlmb$par

  res = AgeAndLengthBasedCatchCurvesCalcs(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge, CatchCurveType,
                                          MLL, SelectivityType, lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)

  # calculate approximate cv for lengths at max age (growth diagnostic), if single timestep
  CV_LenAtMaxAge_Results = GetMaxLenAtAgeCV_AgeAndLengthBasedCatchCurve(TimeStep, res, MaxAge)

  ResultsSummary = list(ParamEst = ParamEst,
                        convergence = nlmb$convergence,
                        nll = nlmb$objective,
                        SampleSize = sum(ObsRetCatchFreqAtLen))

  ModelDiag = list(params = nlmb$par,
                   vcov.Params = vcov.Params,
                   cor.Params = cor.Params,
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
                   ExpRetCatchPropLengthGivenIntAge=res$ExpRetCatchPropLengthGivenIntAge,
                   ExpRetCatchPropLengthGivenIntAge_Fem=res$ExpRetCatchPropLengthGivenIntAge_Fem,
                   ExpRetCatchPropLengthGivenIntAge_Mal=res$ExpRetCatchPropLengthGivenIntAge_Mal,
                   ObsRetCatchFreqAtLen=ObsRetCatchFreqAtLen,
                   ObsRetCatchFreqAtLengthAndAge=ObsRetCatchFreqAtLengthAndAge,
                   CV_LenAtMaxAge=CV_LenAtMaxAge_Results$CV_LenAtMaxAge,
                   FemCV_LenAtMaxAge=CV_LenAtMaxAge_Results$FemCV_LenAtMaxAge,
                   MalCV_LenAtMaxAge=CV_LenAtMaxAge_Results$MalCV_LenAtMaxAge)

  Results = list(ParamEst = ParamEst,
                 EstFMort = ParamEst[1,1],
                 EstFMort_se = ((ParamEst[1,3]-ParamEst[1,2])/2)/1.96,
                 convergence = nlmb$convergence,
                 nll = nlmb$objective,
                 SampleSize = sum(ObsRetCatchFreqAtLen),
                 ResultsSummary = ResultsSummary,
                 ModelDiag = ModelDiag)

  return(Results)
}


#' Calculate CV for length corresponding to max age (model diagnostic)
#'
#' Calculate CV for length corresponding to max age (model diagnostic), for
#' single sex, or both females and males if model is for separate sexes
#'
#' @keywords internal
#'
#' @param Timestep model timestep
#' @param res output from AgeAndLengthBasedCatchCurvesCalcs function
#' @param MaxAge maximum age
#'
#' @return ParamEst
GetMaxLenAtAgeCV_AgeAndLengthBasedCatchCurve <- function (TimeStep, res, MaxAge)
{

  GrowthModelType = res$GrowthModelType # 1=von Bertalaffy, single sex, 2=Schnute
  # single sex, 3=von Bertalanffy 2 sexes, 4=Schnute 2 sexes

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

  CV_LenAtMaxAge_Results = list(CV_LenAtMaxAge=CV_LenAtMaxAge,
                               FemCV_LenAtMaxAge=FemCV_LenAtMaxAge,
                               MalCV_LenAtMaxAge=MalCV_LenAtMaxAge)

  return(CV_LenAtMaxAge_Results)
}


#' Get estimates of model parameters and associated 95 percent confidence intervals
#'
#' Get estimates of model parameters and associated 95 percent confidence intervals. Outputs vary
#' depending on growth model type, selectivity type and number of sexes
#'
#' @keywords internal
#'
#' @param GrowthCurveType 1=von Bertalaffy, 2=Schnute
#' @param SelectivityType 1=selectivity inputted as vector, 2=logistic
#' @param params estimated parameters, varying depending on growth model type, selectivity type and number of sexes
#' @param nlmb output from nlminb optimisation
#' @param ses estimates of standard errors for estimated parameters
#'
#' @return ParamEst
GetParamRes_AgeAndLengthBasedCatchCurve <- function (GrowthCurveType, SelectivityType, params, nlmb, ses)
{

  # estimate of fishing mortality in logit space
  temp = c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1])

  # inverse logit transformed value
  EstFMort = 1/(1+exp(-temp));

  if (GrowthCurveType == 1) { # von Bertalanffy
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
  }

  if (GrowthCurveType == 2) { # Schnute
    if (SelectivityType == 1 & length(params)==5) { # selectivity vector input, single sex input
      Esty2 = exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
      Esta = c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3])
      Estb = c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4])
      EstCV = exp(c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) *  ses[5]))
      ParamEst = t(data.frame(FMort = round(EstFMort, 3), y2 = round(Esty2, 3),
                              a = round(Esta, 3), b = round(Estb, 3), CV = round(EstCV, 3)))
    }
    if (SelectivityType == 2 & length(params)==7) { # logistic selectivity, single sex input
      EstL50 = exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
      EstDelta = exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
      Esty2 = exp(c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
      Esta = c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5])
      Estb = c(nlmb$par[6], nlmb$par[6] + c(-1.96, 1.96) * ses[6])
      EstCV = exp(c(nlmb$par[7], nlmb$par[7] + c(-1.96, 1.96) *  ses[7]))
      ParamEst = t(data.frame(FMort = round(EstFMort, 3), SelL50 = round(EstL50, 3), Delta = round(EstDelta, 3),
                              y2 = round(Esty2, 3), a = round(Esta, 3), b = round(Estb, 3), CV = round(EstCV, 3)))
    }
    if (SelectivityType == 1 & length(params)==8) {   # selectivity vector input, separate sex input
      Esty2_F = exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
      Esty2_M = exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
      Esta_F = c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4])
      Esta_M = c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5])
      Estb_F = c(nlmb$par[6], nlmb$par[6] + c(-1.96, 1.96) * ses[6])
      Estb_M = c(nlmb$par[7], nlmb$par[7] + c(-1.96, 1.96) * ses[7])
      EstCV = exp(c(nlmb$par[8], nlmb$par[8] + c(-1.96, 1.96) *  ses[8]))
      ParamEst = t(data.frame(FMort = round(EstFMort, 3), y2_F = round(Esty2_F, 3), y2_M = round(Esty2_M, 3),
                              a_F = round(Esta_F, 3), a_M = round(Esta_M, 3), b_F = round(Estb_F, 3), b_M = round(Estb_M, 3),
                              CV = round(EstCV, 3)))
    }
    if (SelectivityType == 2 & length(params)==8) { # selectivity vector input, separate sex input # logistic selectivity, separate sex input
      EstL50 = exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
      EstDelta = exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
      Esty2_F = exp(c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
      Esty2_M = exp(c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5]))
      Esta_F = c(nlmb$par[6], nlmb$par[6] + c(-1.96, 1.96) * ses[6])
      Esta_M = c(nlmb$par[7], nlmb$par[7] + c(-1.96, 1.96) * ses[7])
      Estb_F = c(nlmb$par[8], nlmb$par[8] + c(-1.96, 1.96) * ses[8])
      Estb_M = c(nlmb$par[9], nlmb$par[9] + c(-1.96, 1.96) * ses[9])
      EstCV = exp(c(nlmb$par[10], nlmb$par[10] + c(-1.96, 1.96) *  ses[10]))
      ParamEst = t(data.frame(FMort = round(EstFMort, 3), SelL50 = round(EstL50, 3), Delta = round(EstDelta, 3),
                              y2_F = round(Esty2_F, 3), y2_M = round(Esty2_M, 3), a_F = round(Esta_F, 3), a_M = round(Esta_M, 3),
                              b_F = round(Estb_F, 3), b_M = round(Estb_M, 3), CV = round(EstCV, 3)))
    }
  }

  colnames(ParamEst) = c("Estimate", "lw_95%CL", "up_95%CL")
  return(ParamEst)

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
#' @return ObsRetCatchFreqAtLengthAndIntAge
ConvertObsDataFromDecAgesToIntegerAges <- function(TimeStep, MaxAge, nLenCl, ObsRetCatchFreqAtLengthAndAge) {

  Ages = seq(TimeStep,MaxAge,TimeStep)
  nTimeSteps = length(Ages)
  AgeCl = floor(round(Ages,6))
  nAgeCl = length(unique(AgeCl))

  ObsRetCatchFreqAtLengthAndIntAge <- data.frame(matrix(0, nrow = nAgeCl, ncol = nLenCl))
  colnames(ObsRetCatchFreqAtLengthAndIntAge) <- midpt
  ObsRetCatchFreqAtLengthAndIntAge = as.matrix(ObsRetCatchFreqAtLengthAndIntAge)

  for (t in 1:nTimeSteps) {
    if (AgeCl[1]<1) {
      x = AgeCl[t]+1
    } else {
      x = AgeCl[t]
    }

    ObsRetCatchFreqAtLengthAndIntAge[x,] = ObsRetCatchFreqAtLengthAndIntAge[x,] + ObsRetCatchFreqAtLengthAndAge[t,]
  }

  return(ObsRetCatchFreqAtLengthAndIntAge)
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


#' Simulate age frequency data
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


#' Calculate catch proportions within length classes at age
#'
#' Calculate retained and discarded catch proportions within each of the length classes, for each decimal age
#'
#' @keywords internal
#'
#' @param nTimeSteps number of model timesteps
#' @param nLenCl number of length classes
#' @param midpt mid points of length classes
#' @param RetCatchAtDecAgeLen_Fem expected female retained catches at length and decimal age
#' @param RetCatchAtDecAgeLen_Mal expected male retained catches at length and decimal age
#' @param DiscCatchAtDecAgeLen_Fem expected female discarded catches at length and decimal age
#' @param DiscCatchAtDecAgeLen_Mal expected male discarded catches at length and decimal age
#'
#' @return ExpRetCatchPropLengthGivenDecAge_Fem, ExpRetCatchPropLengthGivenDecAge_Mal, ExpDiscCatchPropLengthGivenDecAge_Fem,
#' ExpDiscCatchPropLengthGivenDecAge_Mal
#'
CalcCatchLenPropGivenAge <- function(nTimeSteps, nLenCl, midpt, RetCatchAtDecAgeLen_Fem, RetCatchAtDecAgeLen_Mal,
                                     DiscCatchAtDecAgeLen_Fem, DiscCatchAtDecAgeLen_Mal) {

  EmptyFrame <- data.frame(matrix(nrow = nTimeSteps, ncol = nLenCl))
  colnames(EmptyFrame) <- midpt; EmptyFrame = as.matrix(EmptyFrame)
  ExpRetCatchPropLengthGivenDecAge_Fem = EmptyFrame; ExpRetCatchPropLengthGivenDecAge_Mal = EmptyFrame
  ExpDiscCatchPropLengthGivenDecAge_Fem = EmptyFrame; ExpDiscCatchPropLengthGivenDecAge_Mal = EmptyFrame

  for (j in 1:nTimeSteps) {
    # retained catches
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

    # discarded catches
    if (sum(DiscCatchAtDecAgeLen_Fem[j,])==0) {
      ExpDiscCatchPropLengthGivenDecAge_Fem[j,] = 0
    } else {
      ExpDiscCatchPropLengthGivenDecAge_Fem[j,] = DiscCatchAtDecAgeLen_Fem[j,] / sum(DiscCatchAtDecAgeLen_Fem[j,])
    }
    if (sum(DiscCatchAtDecAgeLen_Mal[j,])==0) {
      ExpDiscCatchPropLengthGivenDecAge_Mal[j,] = 0
    } else {
      ExpDiscCatchPropLengthGivenDecAge_Mal[j,] = DiscCatchAtDecAgeLen_Mal[j,] / sum(DiscCatchAtDecAgeLen_Mal[j,])
    }
  }

  Results = list(ExpRetCatchPropLengthGivenDecAge_Fem = ExpRetCatchPropLengthGivenDecAge_Fem,
                ExpRetCatchPropLengthGivenDecAge_Mal = ExpRetCatchPropLengthGivenDecAge_Mal,
                ExpDiscCatchPropLengthGivenDecAge_Fem = ExpDiscCatchPropLengthGivenDecAge_Fem,
                ExpDiscCatchPropLengthGivenDecAge_Mal = ExpDiscCatchPropLengthGivenDecAge_Mal)

  return(Results)

}


#' Calculate catch proportions within length classes at age
#'
#' Calculate retained and discarded catch proportions within each of the length classes, for each decimal age
#'
#' @keywords internal
#'
#' @param nTimeSteps number of model timesteps
#' @param nLenCl number of length classes
#' @param midpt mid points of length classes
#' @param ExpRetCatchPropLengthGivenDecAge_Fem expected female retained catches at length given decimal age
#' @param ExpRetCatchPropAtDecAge_Fem expected female retained catches at decimal age
#' @param ExpRetCatchPropLengthGivenDecAge_Mal expected male retained catches at length given decimal age
#' @param ExpRetCatchPropAtDecAge_Mal expected male retained catches at decimal age
#' @param ExpDiscCatchPropLengthGivenDecAge_Fem expected female discarded catches at length given decimal age
#' @param ExpDiscCatchPropAtDecAge_Fem expected female discarded catches at decimal age
#' @param ExpDiscCatchPropLengthGivenDecAge_Mal expected male discarded catches at length given decimal age
#' @param ExpDiscCatchPropAtDecAge_Mal expected male discarded catches at decimal age
#'
#' @return ExpRetCatchPropDecAgeGivenLength_Fem, ExpRetCatchPropDecAgeGivenLength_Mal, ExpDiscCatchPropDecAgeGivenLength_Fem,
#' ExpDiscCatchPropDecAgeGivenLength_Mal
CalcCatchAgePropGivenLen <- function(nTimeSteps, nLenCl, midpt, ExpRetCatchPropLengthGivenDecAge_Fem, ExpRetCatchPropAtDecAge_Fem,
                                     ExpRetCatchPropLengthGivenDecAge_Mal, ExpRetCatchPropAtDecAge_Mal,
                                     ExpDiscCatchPropLengthGivenDecAge_Fem, ExpDiscCatchPropAtDecAge_Fem,
                                     ExpDiscCatchPropLengthGivenDecAge_Mal, ExpDiscCatchPropAtDecAge_Mal) {

  EmptyFrame <- data.frame(matrix(nrow = nTimeSteps, ncol = nLenCl))
  colnames(EmptyFrame) <- midpt; EmptyFrame = as.matrix(EmptyFrame)
  ExpRetCatchPropDecAgeGivenLength_Fem = EmptyFrame; ExpRetCatchPropDecAgeGivenLength_Mal = EmptyFrame
  FractDenom_Fem = rep(0,nLenCl); FractDenom_Mal = rep(0,nLenCl)

  # retained catches
  for (i in 1:nLenCl) {
    for (j in 1:nTimeSteps) {
      FractDenom_Fem[i] = FractDenom_Fem[i] + (ExpRetCatchPropLengthGivenDecAge_Fem[j,i] * ExpRetCatchPropAtDecAge_Fem[j])
      if (FractDenom_Fem[i] == 0) {
        FractDenom_Fem[i] = 1E-20
      }
      FractDenom_Mal[i] = FractDenom_Mal[i] + (ExpRetCatchPropLengthGivenDecAge_Mal[j,i] * ExpRetCatchPropAtDecAge_Mal[j])
      if (FractDenom_Mal[i] == 0) {
        FractDenom_Mal[i]= 1E-20
      }
    }
    for (j in 1:nTimeSteps) {
      ExpRetCatchPropDecAgeGivenLength_Fem[j,i] = (ExpRetCatchPropLengthGivenDecAge_Fem[j,i] * ExpRetCatchPropAtDecAge_Fem[j]) / FractDenom_Fem[i]
      ExpRetCatchPropDecAgeGivenLength_Mal[j,i] = (ExpRetCatchPropLengthGivenDecAge_Mal[j,i] * ExpRetCatchPropAtDecAge_Mal[j]) / FractDenom_Mal[i]
    }
  }

  # discarded catches
  ExpDiscCatchPropDecAgeGivenLength_Fem = EmptyFrame; ExpDiscCatchPropDecAgeGivenLength_Mal = EmptyFrame
  FractDenom_Fem = rep(0,nLenCl); FractDenom_Mal = rep(0,nLenCl)
  if (!(is.na(MLL)) | !(is.na(RetenParams[1]))) {
    for (i in 1:nLenCl) {
      for (j in 1:nTimeSteps) {
        FractDenom_Fem[i] = FractDenom_Fem[i] + (ExpDiscCatchPropLengthGivenDecAge_Fem[j,i] * ExpDiscCatchPropAtDecAge_Fem[j])
        if (FractDenom_Fem[i] == 0) {
          FractDenom_Fem[i] = 1E-20
        }
        FractDenom_Mal[i] = FractDenom_Mal[i] + (ExpDiscCatchPropLengthGivenDecAge_Mal[j,i] * ExpDiscCatchPropAtDecAge_Mal[j])
        if (FractDenom_Mal[i] == 0) {
          FractDenom_Mal[i]= 1E-20
        }
      }
      for (j in 1:nTimeSteps) {
        ExpDiscCatchPropDecAgeGivenLength_Fem[j,i] = (ExpDiscCatchPropLengthGivenDecAge_Fem[j,i] * ExpDiscCatchPropAtDecAge_Fem[j]) / FractDenom_Fem[i]
        ExpDiscCatchPropDecAgeGivenLength_Mal[j,i] = (ExpDiscCatchPropLengthGivenDecAge_Mal[j,i] * ExpDiscCatchPropAtDecAge_Mal[j]) / FractDenom_Mal[i]
      }
      if (sum(ExpDiscCatchPropDecAgeGivenLength_Fem[,i]) == 0) {
        ExpDiscCatchPropDecAgeGivenLength_Fem[,i] = 1E-20
      }
      if (sum(ExpDiscCatchPropDecAgeGivenLength_Mal[,i]) == 0) {
        ExpDiscCatchPropDecAgeGivenLength_Mal[,i] = 1E-20
      }
    }
  }

  Results = list(ExpRetCatchPropDecAgeGivenLength_Fem = ExpRetCatchPropDecAgeGivenLength_Fem,
                 ExpRetCatchPropDecAgeGivenLength_Mal = ExpRetCatchPropDecAgeGivenLength_Mal,
                 ExpDiscCatchPropDecAgeGivenLength_Fem = ExpDiscCatchPropDecAgeGivenLength_Fem,
                 ExpDiscCatchPropDecAgeGivenLength_Mal = ExpDiscCatchPropDecAgeGivenLength_Mal)

  return(Results)

}

#' Generate random observed catch frequency data
#'
#' Generate random observed catch frequencies for each length class in each age class, and marginal length and age frequencies
#'
#' @keywords internal
#'
#' @param TimeStep model timestep
#' @param nTimeSteps number of model timesteps
#' @param midpt mid points of length classes
#' @param nLenCl number of length classes
#' @param MLL minimum legal length
#' @param RetenParams fish retention parameters
#' @param SampleSize_Fem specified sample size for female retained fish
#' @param SampleSize_Mal specified sample size for male retained fish
#' @param SampleSize specified sample sizes for retained and discarded fish
#' @param DiscSampleSize_Fem specified sample size for female discarded fish
#' @param DiscSampleSize_Mal specified sample size for male  discarded fish
#' @param DiscSampleSize specified sample size for discarded fish
#' @param ObsRetCatchFreqAtLen_Fem randomly generated female retained catch frequency at length
#' @param ExpRetCatchPropDecAgeGivenLength_Fem expected female retained catch prop at decimal age given length
#' @param ObsRetCatchFreqAtLen_Mal randomly generated male retained catch frequency at length
#' @param ExpRetCatchPropDecAgeGivenLength_Mal expected male retained catch prop at decimal age given length
#' @param ObsDiscCatchFreqAtLen_Fem randomly generated female discarded catch frequency at length
#' @param ExpDiscCatchPropDecAgeGivenLength_Fem expected female discarded catch prop at decimal age given length
#' @param ObsDiscCatchFreqAtLen_Mal randomly generated female discarded catch frequency at length
#' @param ExpDiscCatchPropDecAgeGivenLength_Mal expected male discarded catch prop at decimal age given length
#'
#' @return ObsRetCatchFreqAtLengthAndDecAge_Fem, ObsRetCatchFreqAtLengthAndDecAge_Mal, ObsRetCatchFreqAtLengthAndDecAge,
#' ObsDiscCatchFreqAtLengthAndDecAge_Fem, ObsDiscCatchFreqAtLengthAndDecAge_Mal, ObsDiscCatchFreqAtLengthAndDecAge,
#' ObsDecAgeRetCatch_Fem, ObsDecAgeRetCatch_Mal, ObsDecAgeRetCatch, ObsDecAgeDiscCatch_Fem, ObsDecAgeDiscCatch_Mal,
#' ObsDecAgeDiscCatch, ObsAgeClRetCatch_Fem, ObsAgeClRetCatch_Mal, ObsAgeClRetCatch, ObsAgeClDiscCatch_Fem,
#' ObsAgeClDiscCatch_Mal, ObsAgeClDiscCatch, ObsLenClRetCatchMidPt_Fem, ObsLenClRetCatchMidPt_Mal, ObsLenClRetCatchMidPt,
#' ObsLenClDiscCatchMidPt_Fem, ObsLenClDiscCatchMidPt_Mal, ObsLenClDiscCatchMidPt
GenerateCatchFreqData <- function(TimeStep, nTimeSteps, midpt, nLenCl, MLL, RetenParams, SampleSize_Fem,
                                  SampleSize_Mal, SampleSize, DiscSampleSize_Fem, DiscSampleSize_Mal, DiscSampleSize,
                                  ObsRetCatchFreqAtLen_Fem, ExpRetCatchPropDecAgeGivenLength_Fem,
                                  ObsRetCatchFreqAtLen_Mal, ExpRetCatchPropDecAgeGivenLength_Mal,
                                  ObsDiscCatchFreqAtLen_Fem, ExpDiscCatchPropDecAgeGivenLength_Fem,
                                  ObsDiscCatchFreqAtLen_Mal, ExpDiscCatchPropDecAgeGivenLength_Mal) {

  # generate random observed catch frequencies at length and age
  EmptyFrame <- data.frame(matrix(nrow = nTimeSteps, ncol = nLenCl))
  colnames(EmptyFrame) <- midpt; EmptyFrame = as.matrix(EmptyFrame)
  ObsRetCatchFreqAtLengthAndDecAge_Fem = EmptyFrame; ObsRetCatchFreqAtLengthAndDecAge_Mal = EmptyFrame
  ObsRetCatchFreqAtLengthAndDecAge = EmptyFrame
  ObsDiscCatchFreqAtLengthAndDecAge_Fem = EmptyFrame; ObsDiscCatchFreqAtLengthAndDecAge_Mal = EmptyFrame
  ObsDiscCatchFreqAtLengthAndDecAge = EmptyFrame

  for (i in 1:nLenCl) {

    # retained catches
    ObsRetCatchFreqAtLengthAndDecAge_Fem[,i] = rmultinom(1, ObsRetCatchFreqAtLen_Fem[i], ExpRetCatchPropDecAgeGivenLength_Fem[,i])
    ObsRetCatchFreqAtLengthAndDecAge_Mal[,i] = rmultinom(1, ObsRetCatchFreqAtLen_Mal[i], ExpRetCatchPropDecAgeGivenLength_Mal[,i])
    ObsRetCatchFreqAtLengthAndDecAge[,i] = ObsRetCatchFreqAtLengthAndDecAge_Fem[,i] + ObsRetCatchFreqAtLengthAndDecAge_Mal[,i]

    if (!(is.na(MLL)) | !(is.na(RetenParams[1]))) {
      # discarded catches
      ObsDiscCatchFreqAtLengthAndDecAge_Fem[,i] = rmultinom(1, ObsDiscCatchFreqAtLen_Fem[i], ExpDiscCatchPropDecAgeGivenLength_Fem[,i])
      ObsDiscCatchFreqAtLengthAndDecAge_Mal[,i] = rmultinom(1, ObsDiscCatchFreqAtLen_Mal[i], ExpDiscCatchPropDecAgeGivenLength_Mal[,i])
      ObsDiscCatchFreqAtLengthAndDecAge[,i] = ObsDiscCatchFreqAtLengthAndDecAge_Fem[,i] + ObsDiscCatchFreqAtLengthAndDecAge_Mal[,i]
    }
  }

  # generate data for individual fish (age classes and mid points of length classes)
  ObsDecAgeRetCatch_Fem = rep(NA, SampleSize_Fem); ObsDecAgeRetCatch_Mal = rep(NA, SampleSize_Mal); ObsDecAgeRetCatch = rep(NA, SampleSize[1])
  ObsLenClRetCatchMidPt_Fem = rep(NA, SampleSize_Fem); ObsLenClRetCatchMidPt_Mal = rep(NA, SampleSize_Mal); ObsLenClRetCatchMidPt = rep(NA, SampleSize[1])

  # females - retained catch
  strt=1; fnsh=0
  for (j in 1:nTimeSteps) {
    for (i in 1:nLenCl) {
      x=ObsRetCatchFreqAtLengthAndDecAge_Fem[j,i] # number of fish in current length and age class
      if(x>0) {
        fnsh=strt+x-1
        ObsDecAgeRetCatch_Fem[strt:fnsh]=j*TimeStep
        ObsLenClRetCatchMidPt_Fem[strt:fnsh]=midpt[i]
        strt=strt+x
      }
    }
  }
  # males - retained catch
  strt=1; fnsh=0
  for (j in 1:nTimeSteps) {
    for (i in 1:nLenCl) {
      x=ObsRetCatchFreqAtLengthAndDecAge_Mal[j,i] # number of fish in current length and age class
      if(x>0) {
        fnsh=strt+x-1
        ObsDecAgeRetCatch_Mal[strt:fnsh]=j*TimeStep
        ObsLenClRetCatchMidPt_Mal[strt:fnsh]=midpt[i]
        strt=strt+x
      }
    }
  }

  if (is.na(MLL) & is.na(RetenParams[1])) {
    ObsDecAgeDiscCatch_Fem = NA; ObsDecAgeDiscCatch_Mal = NA; ObsDecAgeDiscCatch = NA
    ObsAgeClDiscCatch_Fem = NA; ObsAgeClDiscCatch_Mal = NA; ObsAgeClDiscCatch = NA
    ObsLenClDiscCatchMidPt_Fem = NA; ObsLenClDiscCatchMidPt_Mal = NA; ObsLenClDiscCatchMidPt = NA

  } else {

    ObsDecAgeDiscCatch_Fem = rep(NA, DiscSampleSize_Fem); ObsDecAgeDiscCatch_Mal = rep(NA, DiscSampleSize_Mal); ObsDecAgeDiscCatch = rep(NA, DiscSampleSize)
    ObsLenClDiscCatchMidPt_Fem = rep(NA, DiscSampleSize_Fem); ObsLenClDiscCatchMidPt_Mal = rep(NA, DiscSampleSize_Mal); ObsLenClDiscCatchMidPt = rep(NA, DiscSampleSize)

    # females - discarded catch
    strt=1; fnsh=0
    for (j in 1:nTimeSteps) {
      for (i in 1:nLenCl) {
        x=ObsDiscCatchFreqAtLengthAndDecAge_Fem[j,i] # number of fish in current length and age class
        if(x>0) {
          fnsh=strt+x-1
          ObsDecAgeDiscCatch_Fem[strt:fnsh]=j*TimeStep
          ObsLenClDiscCatchMidPt_Fem[strt:fnsh]=midpt[i]
          strt=strt+x
        }
      }
    }
    # males - discarded catch
    strt=1; fnsh=0
    for (j in 1:nTimeSteps) {
      for (i in 1:nLenCl) {
        x=ObsDiscCatchFreqAtLengthAndDecAge_Mal[j,i] # number of fish in current length and age class
        if(x>0) {
          fnsh=strt+x-1
          ObsDecAgeDiscCatch_Mal[strt:fnsh]=j*TimeStep
          ObsLenClDiscCatchMidPt_Mal[strt:fnsh]=midpt[i]
          strt=strt+x
        }
      }
    }
  }

  ObsDecAgeRetCatch = c(ObsDecAgeRetCatch_Fem, ObsDecAgeRetCatch_Mal) # decimal ages combined sexes
  ObsLenClRetCatchMidPt = c(ObsLenClRetCatchMidPt_Fem, ObsLenClRetCatchMidPt_Mal)

  # get random age class retained catch frequencies
  ObsAgeClRetCatch_Fem = floor(round(ObsDecAgeRetCatch_Fem,6));
  ObsAgeClRetCatch_Mal = floor(round(ObsDecAgeRetCatch_Mal,6))
  ObsAgeClRetCatch = c(ObsAgeClRetCatch_Fem, ObsAgeClRetCatch_Mal)

  if (is.na(MLL) & is.na(RetenParams[1])) {

    ObsAgeClDiscCatch_Fem = NA; ObsAgeClDiscCatch_Mal = NA
    ObsAgeClDiscCatch = NA # decimal ages combined sexes
    ObsDecAgeDiscCatch = NA
    ObsLenClDiscCatchMidPt = NA

  } else {

    ObsAgeClDiscCatch_Fem = floor(round(ObsDecAgeDiscCatch_Fem,6));
    ObsAgeClDiscCatch_Mal = floor(round(ObsDecAgeDiscCatch_Mal,6))
    ObsAgeClDiscCatch = c(ObsAgeClDiscCatch_Fem, ObsAgeClDiscCatch_Mal) # decimal ages combined sexes
    ObsDecAgeDiscCatch = c(ObsDecAgeDiscCatch_Fem, ObsDecAgeDiscCatch_Mal) # decimal ages combined sexes
    ObsLenClDiscCatchMidPt = c(ObsLenClDiscCatchMidPt_Fem, ObsLenClDiscCatchMidPt_Mal)

  }

  Results = list(ObsRetCatchFreqAtLengthAndDecAge_Fem = ObsRetCatchFreqAtLengthAndDecAge_Fem,
                 ObsRetCatchFreqAtLengthAndDecAge_Mal = ObsRetCatchFreqAtLengthAndDecAge_Mal,
                 ObsRetCatchFreqAtLengthAndDecAge = ObsRetCatchFreqAtLengthAndDecAge,
                 ObsDiscCatchFreqAtLengthAndDecAge_Fem = ObsDiscCatchFreqAtLengthAndDecAge_Fem,
                 ObsDiscCatchFreqAtLengthAndDecAge_Mal = ObsDiscCatchFreqAtLengthAndDecAge_Mal,
                 ObsDiscCatchFreqAtLengthAndDecAge = ObsDiscCatchFreqAtLengthAndDecAge,
                 ObsDecAgeRetCatch_Fem = ObsDecAgeRetCatch_Fem,
                 ObsDecAgeRetCatch_Mal = ObsDecAgeRetCatch_Mal,
                 ObsDecAgeRetCatch = ObsDecAgeRetCatch,
                 ObsDecAgeDiscCatch_Fem = ObsDecAgeDiscCatch_Fem,
                 ObsDecAgeDiscCatch_Mal = ObsDecAgeDiscCatch_Mal,
                 ObsDecAgeDiscCatch = ObsDecAgeDiscCatch,
                 ObsAgeClRetCatch_Fem = ObsAgeClRetCatch_Fem,
                 ObsAgeClRetCatch_Mal = ObsAgeClRetCatch_Mal,
                 ObsAgeClRetCatch = ObsAgeClRetCatch,
                 ObsAgeClDiscCatch_Fem = ObsAgeClDiscCatch_Fem,
                 ObsAgeClDiscCatch_Mal =ObsAgeClDiscCatch_Mal,
                 ObsAgeClDiscCatch = ObsAgeClDiscCatch,
                 ObsLenClRetCatchMidPt_Fem = ObsLenClRetCatchMidPt_Fem,
                 ObsLenClRetCatchMidPt_Mal = ObsLenClRetCatchMidPt_Mal,
                 ObsLenClRetCatchMidPt = ObsLenClRetCatchMidPt,
                 ObsLenClDiscCatchMidPt_Fem = ObsLenClDiscCatchMidPt_Fem,
                 ObsLenClDiscCatchMidPt_Mal = ObsLenClDiscCatchMidPt_Mal,
                 ObsLenClDiscCatchMidPt = ObsLenClDiscCatchMidPt)

  return(Results)

}


#' Get discarded catch summary statistics
#'
#' Calculate prop discarded and mean length of discarded fish
#'
#' @keywords internal
#'
#' @param MLL minimum legal length
#' @param RetenParams parameters for fish retention
#' @param ObsDiscCatchFreqAtLen randomly-generated length frequency data for discarded fish
#' @param ObsTotCatchAtLen randomly-generated length frequency data for discarded and retained fish
#' @param ObsLenClDiscCatchMidPt_Fem randomly-generated length data for discarded females
#' @param ObsLenClDiscCatchMidPt_Mal randomly-generated length data for discarded males
#' @param ObsAgeClDiscCatch_Fem randomly-generated age data for discarded females
#' @param ObsAgeClDiscCatch_Mal randomly-generated age data for discarded males
#'
#' @return ObsPropDiscFish, ObsPropDiscFish_sd, ObsMeanLenDiscFish_Fem, ObsMeanLenDiscFish_Mal, ObsMeanLenDiscFish,
#' ObsLenDiscFish_sd_Fem, ObsLenDiscFish_sd_Mal, ObsLenDiscFish_sd, ObsMeanLenDiscFish_sd_Fem, ObsMeanLenDiscFish_sd_Mal,
#' ObsMeanLenDiscFish_sd, ObsMeanAgeDiscFish_Fem, ObsMeanAgeDiscFish_Mal, ObsMeanAgeDiscFish, ObsAgeDiscFish_sd_Fem,
#' ObsAgeDiscFish_sd_Mal, ObsAgeDiscFish_sd, ObsMeanAgeDiscFish_sd_Fem, ObsMeanAgeDiscFish_sd_Mal, ObsMeanAgeDiscFish_sd
#'
GetDiscardCatchStats <- function(MLL, RetenParams, ObsDiscCatchFreqAtLen, ObsTotCatchAtLen, ObsLenClDiscCatchMidPt_Fem,
                                 ObsLenClDiscCatchMidPt_Mal, ObsAgeClDiscCatch_Fem, ObsAgeClDiscCatch_Mal) {

if (is.na(MLL) & is.na(RetenParams[1])) {

  ObsPropDiscFish = NA; ObsPropDiscFish_sd = NA
  ObsMeanLenDiscFish_Fem = NA; ObsMeanLenDiscFish_Mal = NA; ObsMeanLenDiscFish = NA
  ObsLenDiscFish_sd_Fem = NA; ObsLenDiscFish_sd_Mal = NA; ObsLenDiscFish_sd = NA
  ObsMeanLenDiscFish_sd_Fem = NA; ObsMeanLenDiscFish_sd_Mal = NA; ObsMeanLenDiscFish_sd = NA
  ObsMeanAgeDiscFish_Fem = NA; ObsMeanAgeDiscFish_Mal = NA; ObsMeanAgeDiscFish = NA
  ObsAgeDiscFish_sd_Fem = NA; ObsAgeDiscFish_sd_Mal = NA; ObsAgeDiscFish_sd = NA
  ObsMeanAgeDiscFish_sd_Fem = NA; ObsMeanAgeDiscFish_sd_Mal = NA; ObsMeanAgeDiscFish_sd = NA

} else {

  # get expected proportion and associated sd of fish that are discarded
  ObsPropDiscFish = sum(ObsDiscCatchFreqAtLen) / sum(ObsTotCatchAtLen)
  ObsPropDiscFish_sd = sqrt(ObsPropDiscFish * (1 - ObsPropDiscFish)) / sum(ObsTotCatchAtLen)

  # get expected mean length, sd and se, for discarded fish
  ObsMeanLenDiscFish_Fem = mean(ObsLenClDiscCatchMidPt_Fem)
  ObsMeanLenDiscFish_Mal = mean(ObsLenClDiscCatchMidPt_Mal)
  ObsMeanLenDiscFish = mean(c(ObsLenClDiscCatchMidPt_Fem, ObsLenClDiscCatchMidPt_Mal))

  ObsLenDiscFish_sd_Fem = sd(ObsLenClDiscCatchMidPt_Fem)
  ObsLenDiscFish_sd_Mal = sd(ObsLenClDiscCatchMidPt_Mal)
  ObsLenDiscFish_sd = sd(c(ObsLenClDiscCatchMidPt_Fem, ObsLenClDiscCatchMidPt_Mal))

  ObsMeanLenDiscFish_sd_Fem = ObsLenDiscFish_sd_Fem / sqrt(length(ObsLenClDiscCatchMidPt_Fem)-1)
  ObsMeanLenDiscFish_sd_Mal = ObsLenDiscFish_sd_Mal / sqrt(length(ObsLenClDiscCatchMidPt_Mal)-1)
  ObsMeanLenDiscFish_sd = ObsLenDiscFish_sd / sqrt(length(c(ObsLenClDiscCatchMidPt_Fem, ObsLenClDiscCatchMidPt_Mal))-1)

  # get expected mean age, sd and se, for discarded fish
  ObsMeanAgeDiscFish_Fem = mean(ObsAgeClDiscCatch_Fem)
  ObsMeanAgeDiscFish_Mal = mean(ObsAgeClDiscCatch_Mal)
  ObsMeanAgeDiscFish = mean(c(ObsAgeClDiscCatch_Fem, ObsAgeClDiscCatch_Mal))

  ObsAgeDiscFish_sd_Fem = sd(ObsAgeClDiscCatch_Fem)
  ObsAgeDiscFish_sd_Mal = sd(ObsAgeClDiscCatch_Mal)
  ObsAgeDiscFish_sd = sd(c(ObsAgeClDiscCatch_Fem, ObsAgeClDiscCatch_Mal))

  ObsMeanAgeDiscFish_sd_Fem = ObsAgeDiscFish_sd_Fem / sqrt(length(ObsAgeClDiscCatch_Fem)-1)
  ObsMeanAgeDiscFish_sd_Mal = ObsAgeDiscFish_sd_Mal / sqrt(length(ObsAgeClDiscCatch_Mal)-1)
  ObsMeanAgeDiscFish_sd = ObsAgeDiscFish_sd / sqrt(length(c(ObsAgeClDiscCatch_Fem, ObsAgeClDiscCatch_Mal))-1)

}

  Results = list(ObsPropDiscFish = ObsPropDiscFish,
                 ObsPropDiscFish_sd = ObsPropDiscFish_sd,
                 ObsMeanLenDiscFish_Fem = ObsMeanLenDiscFish_Fem,
                 ObsMeanLenDiscFish_Mal = ObsMeanLenDiscFish_Mal,
                 ObsMeanLenDiscFish = ObsMeanLenDiscFish,
                 ObsLenDiscFish_sd_Fem = ObsLenDiscFish_sd_Fem,
                 ObsLenDiscFish_sd_Mal = ObsLenDiscFish_sd_Mal,
                 ObsLenDiscFish_sd = ObsLenDiscFish_sd,
                 ObsMeanLenDiscFish_sd_Fem = ObsMeanLenDiscFish_sd_Fem,
                 ObsMeanLenDiscFish_sd_Mal = ObsMeanLenDiscFish_sd_Mal,
                 ObsMeanLenDiscFish_sd = ObsMeanLenDiscFish_sd,
                 ObsMeanAgeDiscFish_Fem = ObsMeanAgeDiscFish_Fem,
                 ObsMeanAgeDiscFish_Mal = ObsMeanAgeDiscFish_Mal,
                 ObsMeanAgeDiscFish = ObsMeanAgeDiscFish,
                 ObsAgeDiscFish_sd_Fem = ObsAgeDiscFish_sd_Fem,
                 ObsAgeDiscFish_sd_Mal = ObsAgeDiscFish_sd_Mal,
                 ObsAgeDiscFish_sd = ObsAgeDiscFish_sd,
                 ObsMeanAgeDiscFish_sd_Fem = ObsMeanAgeDiscFish_sd_Fem,
                 ObsMeanAgeDiscFish_sd_Mal = ObsMeanAgeDiscFish_sd_Mal,
                 ObsMeanAgeDiscFish_sd = ObsMeanAgeDiscFish_sd)

  return(Results)

}

test = c(10,30,50)
test2 = c(5, 10, 10)
test3=rep(test,test2)

#' Get retained catch summary statistics
#'
#' Calculate prop retained and mean length of retained fish
#'
#' @keywords internal
#'
#' @param ObsRetCatchFreqAtLen randomly-generated length frequency data for retained fish
#' @param ObsTotCatchAtLen randomly-generated length frequency data for retained and discarded fish
#' @param ObsLenClRetCatchMidPt_Fem randomly-generated length data for females
#' @param ObsLenClRetCatchMidPt_Mal randomly-generated length data for males
#' @param ObsAgeClRetCatch_Fem randomly-generated age class data for females
#' @param ObsAgeClRetCatch_Mal randomly-generated length class data for females
#'
#' @return ObsPropRetFish, ObsPropRetFish_sd, ObsMeanLenRetFish_Fem, ObsMeanLenRetFish_Mal, ObsMeanLenRetFish,
#' ObsLenRetFish_sd_Fem, ObsLenRetFish_sd_Mal, ObsLenRetFish_sd, ObsMeanLenRetFish_sd_Fem, ObsMeanLenRetFish_sd_Mal,
#' ObsMeanLenRetFish_sd, ObsMeanAgeRetFish_Fem, ObsMeanAgeRetFish_Mal, ObsMeanAgeRetFish, ObsAgeRetFish_sd_Fem,
#' ObsAgeRetFish_sd_Mal, ObsAgeRetFish_sd, ObsMeanAgeRetFish_sd_Fem, ObsMeanAgeRetFish_sd_Mal, ObsMeanAgeRetFish_sd
#'
GetRetainedCatchStats <- function(ObsRetCatchFreqAtLen, ObsTotCatchAtLen, ObsLenClRetCatchMidPt_Fem, ObsLenClRetCatchMidPt_Mal,
                                  ObsAgeClRetCatch_Fem, ObsAgeClRetCatch_Mal) {

    # get proportion and associated sd for randomly-generated data for retained fish
    ObsPropRetFish = sum(ObsRetCatchFreqAtLen) / sum(ObsTotCatchAtLen)
    ObsPropRetFish_sd = sqrt(ObsPropRetFish * (1 - ObsPropRetFish)) / sum(ObsTotCatchAtLen)

    # get mean length for randomly-generated data for released fish
    ObsMeanLenRetFish_Fem = mean(ObsLenClRetCatchMidPt_Fem)
    ObsMeanLenRetFish_Mal = mean(ObsLenClRetCatchMidPt_Mal)
    ObsMeanLenRetFish = mean(c(ObsLenClRetCatchMidPt_Fem,ObsLenClRetCatchMidPt_Mal))

    # get sd for length, for randomly-generated data for released fish
    ObsLenRetFish_sd_Fem = sd(ObsLenClRetCatchMidPt_Fem)
    ObsLenRetFish_sd_Mal = sd(ObsLenClRetCatchMidPt_Mal)
    ObsLenRetFish_sd = sd(c(ObsLenClRetCatchMidPt_Fem,ObsLenClRetCatchMidPt_Mal))

    # get standard error for lengths, for randomly-generated data for released fish
    ObsMeanLenRetFish_sd_Fem = ObsLenRetFish_sd_Fem / sqrt(length(ObsLenClRetCatchMidPt_Fem)-1)
    ObsMeanLenRetFish_sd_Mal = ObsLenRetFish_sd_Mal / sqrt(length(ObsLenClRetCatchMidPt_Mal)-1)
    ObsMeanLenRetFish_sd = ObsLenRetFish_sd / sqrt(length(c(ObsLenClRetCatchMidPt_Fem,ObsLenClRetCatchMidPt_Mal))-1)

    # get mean age for randomly-generated data for released fish (based on age class data)
    ObsMeanAgeRetFish_Fem = mean(ObsAgeClRetCatch_Fem)
    ObsMeanAgeRetFish_Mal = mean(ObsAgeClRetCatch_Mal)
    ObsMeanAgeRetFish = mean(c(ObsAgeClRetCatch_Fem, ObsAgeClRetCatch_Mal))

    # note, this is calculated sd for the population, not the sample
    ObsAgeRetFish_sd_Fem = sd(ObsAgeClRetCatch_Fem)
    ObsAgeRetFish_sd_Mal = sd(ObsAgeClRetCatch_Mal)
    ObsAgeRetFish_sd = sd(c(ObsAgeClRetCatch_Fem, ObsAgeClRetCatch_Mal))

    # standard deviation of for mean length of retained fish
    ObsMeanAgeRetFish_sd_Fem = ObsAgeRetFish_sd_Fem / sqrt(length(ObsAgeClRetCatch_Fem)-1)
    ObsMeanAgeRetFish_sd_Mal = ObsAgeRetFish_sd_Mal / sqrt(length(ObsAgeClRetCatch_Mal)-1)
    ObsMeanAgeRetFish_sd = ObsAgeRetFish_sd / sqrt(length(c(ObsAgeClRetCatch_Fem,ObsAgeClRetCatch_Mal))-1)

  Results = list(ObsPropRetFish = ObsPropRetFish,
                 ObsPropRetFish_sd = ObsPropRetFish_sd,
                 ObsMeanLenRetFish_Fem = ObsMeanLenRetFish_Fem,
                 ObsMeanLenRetFish_Mal = ObsMeanLenRetFish_Mal,
                 ObsMeanLenRetFish = ObsMeanLenRetFish,
                 ObsLenRetFish_sd_Fem = ObsLenRetFish_sd_Fem,
                 ObsLenRetFish_sd_Mal = ObsLenRetFish_sd_Mal,
                 ObsLenRetFish_sd = ObsLenRetFish_sd,
                 ObsMeanLenRetFish_sd_Fem = ObsMeanLenRetFish_sd_Fem,
                 ObsMeanLenRetFish_sd_Mal = ObsMeanLenRetFish_sd_Mal,
                 ObsMeanLenRetFish_sd = ObsMeanLenRetFish_sd,
                 ObsMeanAgeRetFish_Fem = ObsMeanAgeRetFish_Fem,
                 ObsMeanAgeRetFish_Mal = ObsMeanAgeRetFish_Mal,
                 ObsMeanAgeRetFish = ObsMeanAgeRetFish,
                 ObsAgeRetFish_sd_Fem = ObsAgeRetFish_sd_Fem,
                 ObsAgeRetFish_sd_Mal = ObsAgeRetFish_sd_Mal,
                 ObsAgeRetFish_sd = ObsAgeRetFish_sd,
                 ObsMeanAgeRetFish_sd_Fem = ObsMeanAgeRetFish_sd_Fem,
                 ObsMeanAgeRetFish_sd_Mal = ObsMeanAgeRetFish_sd_Mal,
                 ObsMeanAgeRetFish_sd = ObsMeanAgeRetFish_sd)

  return(Results)

}


#' Get fish selectivity and retention vectors
#'
#' Get fish selectivity and retention vectors, depending on user inputs
#'
#' @keywords internal
#'
#' @param midpt mid points of length classes
#' @param SelectivityType 1=inputted as vector, 2= inputted as params
#' @param SelParams params for gear selectivity
#' @param SelectivityVec inputted selectivity at length vector
#' @param RetenParams params for fish retention
#' @param MLL minimum legal length
#'
#' @return SelAtLength, RetAtLength
GetSelectivityAndRetention <- function(midpt, SelectivityType, SelParams, SelectivityVec, RetenParams, MLL) {

  # gear selectivity
  if (SelectivityType == 1) { # inputted as vector
    SelAtLength = SelectivityVec
  }
  if (SelectivityType == 2) { # logistic gear selectivity
    if (!is.na(SelParams[1])) { # calculate selectivity curve
      L50=SelParams[1]
      L95=L50 + SelParams[2]
      SelAtLength = CalcLogisticSelOrReten(L50, L95, midpt)
    } else { # selectivity curve is unknown.
      cat("Gear selectivity params or a selectivity vector needs to be specified",'\n')
    }
  }

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

  Results = list(SelAtLength = SelAtLength,
                 RetAtLength = RetAtLength)

  return(Results)

}


#' Get retained catch frequencies at length in each integer age, and in each age class
#'
#' Get retained catch frequencies in each of the length classes at each integer age, and in each age class
#'
#' @keywords internal
#'
#' @param TimeStep model time step
#' @param MaxAge maximum age
#' @param midpt mid points of length classes
#' @param nLenCl number of length classes
#' @param ObsAgeClRetCatch_Fem generated female catch in each integer age class
#' @param ObsAgeClRetCatch_Mal generated male catch in each integer age class
#' @param ObsLenClRetCatchMidPt_Fem generated female catch in each length class
#' @param ObsLenClRetCatchMidPt_Mal generated male catch in each length class
#'
#' @return ObsRetCatchFreqAtLengthAndIntAge_Fem, ObsRetCatchFreqAtLengthAndIntAge_Mal, ObsRetCatchFreqAtLengthAndIntAge,
#' ObsRetCatchFreqAtIntAge_Fem, ObsRetCatchFreqAtIntAge_Mal, ObsRetCatchFreqAtIntAge
#'
GetRetCatchFreqAtLenAndIntAge <- function(TimeStep, MaxAge, midpt, nLenCl, ObsAgeClRetCatch_Fem, ObsAgeClRetCatch_Mal,
                                       ObsLenClRetCatchMidPt_Fem, ObsLenClRetCatchMidPt_Mal) {

  # get catch frequencies at length for each integer age
  MinAge = floor(TimeStep)
  nAgeCl = length(MinAge:MaxAge)
  ObsRetCatchFreqAtLengthAndIntAge_Fem <- data.frame(matrix(0,nrow = nAgeCl, ncol = nLenCl))
  ObsRetCatchFreqAtLengthAndIntAge_Fem = as.matrix(ObsRetCatchFreqAtLengthAndIntAge_Fem)
  colnames(ObsRetCatchFreqAtLengthAndIntAge_Fem)=midpt
  ObsRetCatchFreqAtLengthAndIntAge_Mal=ObsRetCatchFreqAtLengthAndIntAge_Fem
  ObsRetCatchFreqAtLengthAndIntAge=ObsRetCatchFreqAtLengthAndIntAge_Fem
  nTimeSteps = length(seq(TimeStep,MaxAge,TimeStep))
  ObsRetCatchFreqAtIntAge_Fem = rep(0,nAgeCl)
  ObsRetCatchFreqAtIntAge_Mal = ObsRetCatchFreqAtIntAge_Fem
  ObsRetCatchFreqAtIntAge = ObsRetCatchFreqAtIntAge_Fem

  for (i in 1:length(ObsAgeClRetCatch_Fem)) {
    tempAge = ObsAgeClRetCatch_Fem[i]
    tempLen = ObsLenClRetCatchMidPt_Fem[i]
    j = tempAge - MinAge + 1
    k = which(midpt==tempLen)
    ObsRetCatchFreqAtLengthAndIntAge_Fem[j,k] = ObsRetCatchFreqAtLengthAndIntAge_Fem[j,k] + 1
    ObsRetCatchFreqAtIntAge_Fem[j] = ObsRetCatchFreqAtIntAge_Fem[j] + 1
  }
  for (i in 1:length(ObsAgeClRetCatch_Mal)) {
    tempAge = ObsAgeClRetCatch_Mal[i]
    tempLen = ObsLenClRetCatchMidPt_Mal[i]
    j = tempAge - MinAge + 1
    k = which(midpt==tempLen)
    ObsRetCatchFreqAtLengthAndIntAge_Mal[j,k] = ObsRetCatchFreqAtLengthAndIntAge_Mal[j,k] + 1
    ObsRetCatchFreqAtIntAge_Mal[j] = ObsRetCatchFreqAtIntAge_Mal[j] + 1
  }
  ObsRetCatchFreqAtLengthAndIntAge = ObsRetCatchFreqAtLengthAndIntAge_Fem + ObsRetCatchFreqAtLengthAndIntAge_Mal
  ObsRetCatchFreqAtIntAge = ObsRetCatchFreqAtIntAge_Fem + ObsRetCatchFreqAtIntAge_Mal

  Results = list(ObsRetCatchFreqAtLengthAndIntAge_Fem = ObsRetCatchFreqAtLengthAndIntAge_Fem,
                 ObsRetCatchFreqAtLengthAndIntAge_Mal = ObsRetCatchFreqAtLengthAndIntAge_Mal,
                 ObsRetCatchFreqAtLengthAndIntAge = ObsRetCatchFreqAtLengthAndIntAge,
                 ObsRetCatchFreqAtIntAge_Fem = ObsRetCatchFreqAtIntAge_Fem,
                 ObsRetCatchFreqAtIntAge_Mal = ObsRetCatchFreqAtIntAge_Mal,
                 ObsRetCatchFreqAtIntAge = ObsRetCatchFreqAtIntAge)

  return(Results)

}


#' Get discarded catch frequencies at length in each integer age, and in each age class
#'
#' Get discarded catch frequencies in each of the length classes at each integer age, and in each age class
#'
#' @keywords internal
#'
#' @param TimeStep model time step
#' @param MaxAge maximum age
#' @param midpt mid points of length classes
#' @param nLenCl number of length classes
#' @param ObsAgeClDiscCatch_Fem generated female catch in each integer age class
#' @param ObsAgeClDiscCatch_Mal generated male catch in each integer age class
#' @param ObsLenClDiscCatchMidPt_Fem generated female catch in each length class
#' @param ObsLenClDiscCatchMidPt_Mal generated male catch in each length class
#'
#' @return ObsDiscCatchFreqAtLengthAndIntAge_Fem, ObsDiscCatchFreqAtLengthAndIntAge_Mal, ObsDiscCatchFreqAtLengthAndIntAge
GetDiscCatchFreqAtLenAndIntAge <- function(TimeStep, MaxAge, midpt, nLenCl, ObsAgeClDiscCatch_Fem, ObsAgeClDiscCatch_Mal,
                                          ObsLenClDiscCatchMidPt_Fem, ObsLenClDiscCatchMidPt_Mal) {

  # get catch frequencies at length for each integer age
  MinAge = floor(TimeStep)
  nAgeCl = length(MinAge:MaxAge)
  ObsDiscCatchFreqAtLengthAndIntAge_Fem <- data.frame(matrix(0,nrow = nAgeCl, ncol = nLenCl))
  ObsDiscCatchFreqAtLengthAndIntAge_Fem = as.matrix(ObsDiscCatchFreqAtLengthAndIntAge_Fem)
  colnames(ObsDiscCatchFreqAtLengthAndIntAge_Fem)=midpt
  ObsDiscCatchFreqAtLengthAndIntAge_Mal=ObsDiscCatchFreqAtLengthAndIntAge_Fem
  ObsDiscCatchFreqAtLengthAndIntAge=ObsDiscCatchFreqAtLengthAndIntAge_Fem
  nTimeSteps = length(seq(TimeStep,MaxAge,TimeStep))
  ObsDiscCatchFreqAtIntAge_Fem = rep(0,nAgeCl)
  ObsDiscCatchFreqAtIntAge_Mal = ObsDiscCatchFreqAtIntAge_Fem
  ObsDiscCatchFreqAtIntAge = ObsDiscCatchFreqAtIntAge_Fem

  for (i in 1:length(ObsAgeClDiscCatch_Fem)) {
    tempAge = ObsAgeClDiscCatch_Fem[i]
    tempLen = ObsLenClDiscCatchMidPt_Fem[i]
    j = tempAge - MinAge + 1
    k = which(midpt==tempLen)
    ObsDiscCatchFreqAtLengthAndIntAge_Fem[j,k] = ObsDiscCatchFreqAtLengthAndIntAge_Fem[j,k] + 1
    ObsDiscCatchFreqAtIntAge_Fem[j] = ObsDiscCatchFreqAtIntAge_Fem[j] + 1
  }
  for (i in 1:length(ObsAgeClDiscCatch_Mal)) {
    tempAge = ObsAgeClDiscCatch_Mal[i]
    tempLen = ObsLenClDiscCatchMidPt_Mal[i]
    j = tempAge - MinAge + 1
    k = which(midpt==tempLen)
    ObsDiscCatchFreqAtLengthAndIntAge_Mal[j,k] = ObsDiscCatchFreqAtLengthAndIntAge_Mal[j,k] + 1
    ObsDiscCatchFreqAtIntAge_Mal[j] = ObsDiscCatchFreqAtIntAge_Mal[j] + 1
  }
  ObsDiscCatchFreqAtLengthAndIntAge = ObsDiscCatchFreqAtLengthAndIntAge_Fem + ObsDiscCatchFreqAtLengthAndIntAge_Mal
  ObsDiscCatchFreqAtIntAge = ObsDiscCatchFreqAtIntAge_Fem + ObsDiscCatchFreqAtIntAge_Mal

  Results = list(ObsDiscCatchFreqAtLengthAndIntAge_Fem = ObsDiscCatchFreqAtLengthAndIntAge_Fem,
                 ObsDiscCatchFreqAtLengthAndIntAge_Mal = ObsDiscCatchFreqAtLengthAndIntAge_Mal,
                 ObsDiscCatchFreqAtLengthAndIntAge = ObsDiscCatchFreqAtLengthAndIntAge,
                 ObsDiscCatchFreqAtIntAge_Fem = ObsDiscCatchFreqAtIntAge_Fem,
                 ObsDiscCatchFreqAtIntAge_Mal = ObsDiscCatchFreqAtIntAge_Mal,
                 ObsDiscCatchFreqAtIntAge = ObsDiscCatchFreqAtIntAge)

  return(Results)

}



#' Get random retained and discarded fish lengths
#'
#' Get random retained and discarded fish lengths, given random length classes
#'
#' @keywords internal
#'
#' @param MLL minimum legal length
#' @param RetenParams retention parameters
#' @param LenInc size of length classes
#' @param ObsLenClRetCatchMidPt_Fem randomly-generated observed retained catch female length class data
#' @param ObsLenClRetCatchMidPt_Mal randomly-generated observed retained catch male length class data
#' @param ObsLenClDiscCatchMidPt_Fem randomly-generated observed discarded catch female length class data
#' @param ObsLenClDiscCatchMidPt_Mal randomly-generated observed discarded catch male length class data
#' @param SampleSize_Fem specified sample size for female retained catch sample
#' @param SampleSize_Mal specified sample size for male retained catch sample
#' @param DiscSampleSize_Fem specified sample size for female discarded catch sample
#' @param DiscSampleSize_Mal specified sample size for male discarded catch sample
#'
#' @return ObsRandLenRetCatch_Fem, ObsRandLenRetCatch_Mal, ObsRandLenRetCatch, ObsRandLenDiscCatch_Fem,
#' ObsRandLenDiscCatch_Mal, ObsRandLenDiscCatch
GetRandFishLengths <- function(MLL, RetenParams, LenInc, ObsLenClRetCatchMidPt_Fem, ObsLenClRetCatchMidPt_Mal,
                               ObsLenClDiscCatchMidPt_Fem, ObsLenClDiscCatchMidPt_Mal, SampleSize_Fem, SampleSize_Mal,
                               DiscSampleSize_Fem, DiscSampleSize_Mal) {


  # generate random fish lengths, within each length class, for retained catches
  LenInterv = LenInc / 2 # randomising fish lengths, within each length class
  ObsRandLenRetCatch_Fem = round(ObsLenClRetCatchMidPt_Fem + runif(SampleSize_Fem,-LenInterv, LenInterv),0)
  ObsRandLenRetCatch_Mal = round(ObsLenClRetCatchMidPt_Mal + runif(SampleSize_Mal,-LenInterv, LenInterv),0)
  ObsRandLenRetCatch = c(ObsRandLenRetCatch_Fem, ObsRandLenRetCatch_Mal)

  if (is.na(MLL) & is.na(RetenParams[1])) {

    ObsRandLenDiscCatch_Fem = NA
    ObsRandLenDiscCatch_Mal = NA
    ObsRandLenDiscCatch = NA

  } else {

    # random fish lengths, within each length class, for each of the fish in discarded catches
    ObsRandLenDiscCatch_Fem = round(ObsLenClDiscCatchMidPt_Fem + runif(DiscSampleSize_Fem,-LenInterv, LenInterv),0)
    ObsRandLenDiscCatch_Mal = round(ObsLenClDiscCatchMidPt_Mal + runif(DiscSampleSize_Mal,-LenInterv, LenInterv),0)
    ObsRandLenDiscCatch = c(ObsRandLenDiscCatch_Fem, ObsRandLenDiscCatch_Mal)
  }

  Results = list(ObsRandLenRetCatch_Fem = ObsRandLenRetCatch_Fem,
                 ObsRandLenRetCatch_Mal = ObsRandLenRetCatch_Mal,
                 ObsRandLenRetCatch = ObsRandLenRetCatch,
                 ObsRandLenDiscCatch_Fem = ObsRandLenDiscCatch_Fem,
                 ObsRandLenDiscCatch_Mal = ObsRandLenDiscCatch_Mal,
                 ObsRandLenDiscCatch = ObsRandLenDiscCatch)

  return(Results)

}

#' Get required sample sizes for randomly-generated data
#'
#' Get required sample sizes for randomly-generated size and age data for retained and discarded catches
#'
#' @keywords internal
#'
#' @param MLL minimum legal length
#' @param RetenParams retention parameters
#' @param SampleSize size of length classes
#' @param RetCatchAtLen_Fem randomly-generated observed retained catch female length class data
#' @param RetCatchAtLen_Mal randomly-generated observed retained catch male length class data
#' @param DiscCatchAtLen_Fem randomly-generated observed discarded catch female length class data
#' @param DiscCatchAtLen_Mal randomly-generated observed discarded catch male length class data
#'
#' @return SampleSize_Fem, SampleSize_Mal, DiscSampleSize_Fem, DiscSampleSize_Mal, DiscSampleSize
GetReqdSampleSizesForRandData <- function(MLL, RetenParams, SampleSize, RetCatchAtLen_Fem, RetCatchAtLen_Mal,
                                      DiscCatchAtLen_Fem, DiscCatchAtLen_Mal) {

  # generate observed retained catch frequencies at length
  SampleSize_Fem = SampleSize[1] * (sum(RetCatchAtLen_Fem) / (sum(RetCatchAtLen_Fem) + sum(RetCatchAtLen_Mal)))
  SampleSize_Fem = round(SampleSize_Fem,0)
  SampleSize_Mal = SampleSize[1] * (sum(RetCatchAtLen_Mal) / (sum(RetCatchAtLen_Fem) + sum(RetCatchAtLen_Mal)))
  SampleSize_Mal = round(SampleSize_Mal,0)
  if (is.na(MLL) & is.na(RetenParams[1])) { # no discarding
    DiscSampleSize_Fem = 0; DiscSampleSize_Mal = 0; DiscSampleSize = 0
  } else { # discarding
    if (length(SampleSize)==2) {
      DiscSampleSize = SampleSize[2] # specifying sample sizes for both retained and discarded fish
    } else {
      # calculate expected proportion of catch that is discarded, then set total discard sample size to
      # this proportion of specified sample size for retained catch
      DiscSampleSize = SampleSize * ((sum(DiscCatchAtLen_Fem) + sum(DiscCatchAtLen_Mal)) /
                                    ((sum(RetCatchAtLen_Fem) + sum(RetCatchAtLen_Mal)) +
                                    (sum(DiscCatchAtLen_Fem) + sum(DiscCatchAtLen_Mal))))
      DiscSampleSize = round(DiscSampleSize,0)
      # DiscSampleSize = SampleSize # set discarded fish sample size equal to retained fish sample size

    }
    DiscSampleSize_Fem = DiscSampleSize * (sum(DiscCatchAtLen_Fem) / (sum(DiscCatchAtLen_Fem) + sum(DiscCatchAtLen_Mal)))
    DiscSampleSize_Fem = round(DiscSampleSize_Fem,0)
    DiscSampleSize_Mal = DiscSampleSize * (sum(DiscCatchAtLen_Mal) / (sum(DiscCatchAtLen_Fem) + sum(DiscCatchAtLen_Mal)))
    DiscSampleSize_Mal = round(DiscSampleSize_Mal,0)
  }

  Results = list(SampleSize_Fem = SampleSize_Fem,
                 SampleSize_Mal = SampleSize_Mal,
                 DiscSampleSize_Fem = DiscSampleSize_Fem,
                 DiscSampleSize_Mal = DiscSampleSize_Mal,
                 DiscSampleSize = DiscSampleSize)

  return(Results)

}

#' Get length frequency sample data simulated from the Dirichlet multinomial distribution.
#'
#' Get length frequency data simulated from the Dirichlet multinomial distribution. Using this
#' distribution allows for auto-correlation in length frequency data associated with, for example,
#' fish schooling (according to size) behaviours. This routine uses the dirmult package. Expected
#' proportions in each length class may be calculated using the SimLenAndAgeFreqData function.
#'
#' @param nSampEvents minimum legal length
#' @param nFishPerSampEvent retention parameters
#' @param theta_val size of length classes
#' @param midpt randomly-generated observed retained catch female length class data
#' @param ExpPropAtLen randomly-generated observed retained catch male length class data
#'
#' @return midpt, simSampEventLenDat, simLenFreq, TotSampSize
#' @examples
#' # Simulate length-frequency data from Dirchlet multinomial distribution.
#' # First, use SimLenAndAgeFreqData function to calculate expected proportions at length
#' set.seed(123)
#' SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
#' MaxAge = 30
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = 0.2
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(500, 50) # L50, L95-L50 for retention
#' DiscMort = 0 # proportion of fish that die due to natural mortality
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.08
#' GrowthParams = c(Linf, vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # Now simulate length-frequency data from Dirchlet multinomial distribution.
#' set.seed(122)
#' nSampEvents = 50
#' nFishPerSampEvent = 20
#' theta_val = 0.3
#' midpt = Res$midpt
#' ExpPropAtLen = Res$ModelDiag$ExpRetCatchPropAtLen
#' res=SimLenFreqDat_DirMultDistn(nSampEvents, nFishPerSampEvent, theta_val, midpt, ExpPropAtLen)
#' plot(res$midpt, res$simLenFreq, "o")
#' @export
SimLenFreqDat_DirMultDistn <- function(nSampEvents, nFishPerSampEvent, theta_val, midpt, ExpPropAtLen) {

  # Simulate length data from a Dirichlet multinomial distribution
  # J = number of fish sampling events
  # K = number of age classes
  # n = number of fish sampled from each sampling event
  # pi = expected proportion at age
  # theta = amount of autocorrelation between ages of fish within sampling events

  nLenCl=length(Res$midpt)
  simDat = dirmult::simPop(J=nSampEvents, K=nLenCl, n=nFishPerSampEvent, pi=ExpPropAtLen, theta=theta_val)
  simSampEventLenDat = data.frame(simDat$data)
  colnames(simSampEventLenDat)=Res$midpt
  rownames(simSampEventLenDat)=1:nSampEvents
  simLenFreq = as.vector(colSums(simSampEventLenDat))
  TotSampSize = sum(simLenFreq)

  Results=list(midpt=midpt,
               simSampEventLenDat=simSampEventLenDat,
               simLenFreq=simLenFreq,
               TotSampSize=TotSampSize)

  return(Results)

}



#' Simulate length and age data
#'
#' Simulate length and age data, given specified growth, mortality and selectivity
#'
#' @param SampleSize required same size, single value generates this number of retained fish, and routines
#' then determines an expected sample size for discarded fish based on input assumptions for selectivity and
#' retention. Can also specify vector, e.g. SampleSize=c(10000, 2000) produces 10000 retained fish, 2000 discarded fish
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
#' # Simulate data
#' set.seed(123)
#' SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
#' set.seed(123)
#' MaxAge = 30
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = 0.2
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0 # proportion of fish that die due to natural mortality
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.08
#' GrowthParams = c(Linf, vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # 2 sexes, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = c(700,850)
#' vbK = c(0.3,0.2)
#' CVSizeAtAge = c(0.08,0.08)
#' GrowthParams = data.frame(Linf=Linf, vbK=vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # 1 sex, Schnute
#' GrowthCurveType = 2 # 1 = von Bertalanffy, 2 = Schnute
#' t1 = 0.5 # growth - Schnute
#' t2 = 25 # growth - Schnute
#' y1 = 100 # growth - Schnute
#' y2 = 1000 # growth - Schnute
#' a = 0.02 # growth - Schnute
#' b = 3.0 # growth - Schnute
#' GrowthParams = c(y1, y2, a, b)
#' RefnceAges = c(t1,t2)
#' CVSizeAtAge = 0.05
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # 2 sexes, Schnute
#' GrowthCurveType = 2 # 1 = von Bertalanffy, 2 = Schnute
#' t1 = c(0.5,0.5) # growth - Schnute
#' t2 = c(25,25) # growth - Schnute
#' y1 = c(100,100) # growth - Schnute
#' y2 = c(1000,1000) # growth - Schnute
#' a = c(0.02,0.02) # growth - Schnute
#' b = c(3,3) # growth - Schnute
#' CVSizeAtAge = c(0.05, 0.05)
#' GrowthParams = data.frame(y1=y1, y2=y2, a=a, b=b)
#' RefnceAges = data.frame(t1=t1,t2=t2)
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # Example with selectivity specified as a vector
#' # Simulate data
#' SampleSize=5000
#' set.seed(123)
#' MaxAge = 30
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = 0.2
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
#' SelectivityType=1 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' lbnd = seq(0,MaxLen - LenInc, LenInc)
#' midpt = lbnd + (LenInc/2)
#' SelectivityVec = 1 / (1 + exp(-log(19)*(midpt-400)/(500-400)))
#' SelParams = c(NA, NA) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0 # proportion of fish that die due to natural mortality
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.05
#' GrowthParams = c(Linf, vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' @export
SimLenAndAgeFreqData <- function(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                                 SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge) {

  DecAges = seq(TimeStep,MaxAge,TimeStep)
  MinAge = floor(DecAges[1])
  AgeCl = MinAge:MaxAge
  nTimeSteps = length(DecAges)
  lbnd = seq(0,MaxLen - LenInc, LenInc)
  ubnd = lbnd + LenInc
  midpt = lbnd + (LenInc/2)
  nLenCl = length(midpt)
  FishMort_logit = log(FishMort/(1-FishMort)) # logit transform (so F is always between 0 and 1)

  # get selectivity and retention
  SelRes = GetSelectivityAndRetention(midpt, SelectivityType, SelParams, SelectivityVec, RetenParams, MLL)
  SelAtLength = SelRes$SelAtLength
  RetAtLength = SelRes$RetAtLength

  # get growth model type
  GrowthModelType = GetGrowthModelType(GrowthCurveType, GrowthParams)

  # get inputs for length transition matrices
  Res = GetGrowthInputsForLengthTransitionMatrices(MaxAge, TimeStep, nLenCl, midpt, GrowthCurveType, GrowthModelType,
                                                   GrowthParams, RefnceAges, SelectivityType)
  MeanSizeAtAge = Res$MeanSizeAtAge
  TimestepGrowthSizeInc = Res$TimestepGrowthSizeInc

  # get expected size distribution of 1+ recruits
  RecLenDist = CalcSizeDistOfRecruits(MeanSizeAtAge, CVSizeAtAge, lbnd, ubnd, midpt, nLenCl)

  if (GrowthModelType==1 | GrowthModelType==3) { # single sex
    TimestepGrowthSizeInc_Fem = TimestepGrowthSizeInc; TimestepGrowthSizeInc_Mal = TimestepGrowthSizeInc
    CVSizeAtAge_Fem = CVSizeAtAge; CVSizeAtAge_Mal = CVSizeAtAge
    RecLenDist_Fem = RecLenDist; RecLenDist_Mal = RecLenDist
  } else { # separate sex
    TimestepGrowthSizeInc_Fem = as.vector(unlist(TimestepGrowthSizeInc[1,]))
    TimestepGrowthSizeInc_Mal = as.vector(unlist(TimestepGrowthSizeInc[2,]))
    CVSizeAtAge_Fem = CVSizeAtAge[1]; CVSizeAtAge_Mal = CVSizeAtAge[2]
    RecLenDist_Fem = as.vector(unlist(RecLenDist[1,]))
    RecLenDist_Mal = as.vector(unlist(RecLenDist[2,]))
  }

  # calculate length transition matrices for each sex
  LTM_Fem = CalcLTM_cpp(TimestepGrowthSizeInc_Fem, CVSizeAtAge_Fem, lbnd, midpt, ubnd, nLenCl) # length-transition matrix - females
  LTM_Mal = CalcLTM_cpp(TimestepGrowthSizeInc_Mal, CVSizeAtAge_Mal, lbnd, midpt, ubnd, nLenCl) # length-transition matrix - males

  # get expected landed and discarded catches, and selectivity of landings and discards (retention function approach)
  InitRecNumber = 0.5
  CatchCurveResults_Fem = CalcCatches_AgeAndLengthBasedCatchCurves_cpp(FishMort_logit, NatMort, RecLenDist_Fem, InitRecNumber, MaxAge, TimeStep, nTimeSteps,
                                                                       nLenCl, midpt, RetAtLength, SelAtLength, DiscMort, LTM_Fem)

  CatchCurveResults_Mal = CalcCatches_AgeAndLengthBasedCatchCurves_cpp(FishMort_logit, NatMort, RecLenDist_Mal, InitRecNumber, MaxAge, TimeStep, nTimeSteps,
                                                                       nLenCl, midpt, RetAtLength, SelAtLength, DiscMort, LTM_Mal)

  # expected retained catch at length
  RetCatchAtLen_Fem = CatchCurveResults_Fem$RetCatchAtLen; RetCatchAtLen_Mal = CatchCurveResults_Mal$RetCatchAtLen
  RetCatchAtLen = RetCatchAtLen_Fem + RetCatchAtLen_Mal

  # expected discarded catch at length
  DiscCatchAtLen_Fem = CatchCurveResults_Fem$DiscCatchAtLen; DiscCatchAtLen_Mal = CatchCurveResults_Mal$DiscCatchAtLen
  DiscCatchAtLen = DiscCatchAtLen_Fem + DiscCatchAtLen_Mal

  # expected total catch at length (discarded pus retained catches)
  TotCatchAtLen_Fem = CatchCurveResults_Fem$TotCatchAtLen; TotCatchAtLen_Mal = CatchCurveResults_Mal$TotCatchAtLen
  TotCatchAtLen = TotCatchAtLen_Fem + TotCatchAtLen_Mal

  # expected retained catch at decimal age and length
  SelLandAtLength = CatchCurveResults_Fem$SelLandAtLength # selectivity of landings
  RetCatchAtDecAgeLen_Fem = as.matrix(unlist(CatchCurveResults_Fem$RetCatchAtDecAgeLen)) # retained catches, decimal ages and lengths
  RetCatchAtDecAgeLen_Mal = as.matrix(unlist(CatchCurveResults_Mal$RetCatchAtDecAgeLen))
  RetCatchAtDecAgeLen = RetCatchAtDecAgeLen_Fem + RetCatchAtDecAgeLen_Mal # single sex or combined sexes

  # expected discarded catch at decimal age and length
  DiscCatchAtDecAgeLen_Fem = as.matrix(unlist(CatchCurveResults_Fem$DiscCatchAtDecAgeLen)) # retained catches, decimal ages and lengths
  DiscCatchAtDecAgeLen_Mal = as.matrix(unlist(CatchCurveResults_Mal$DiscCatchAtDecAgeLen))
  DiscCatchAtDecAgeLen = DiscCatchAtDecAgeLen_Fem + DiscCatchAtDecAgeLen_Mal # single sex or combined sexes

  # expected total catch at decimal age and length (discards plus retained catches)
  TotCatchAtDecAgeLen_Fem = as.matrix(unlist(CatchCurveResults_Fem$TotCatchAtDecAgeLen)) # retained catches, decimal ages and lengths
  TotCatchAtDecAgeLen_Mal = as.matrix(unlist(CatchCurveResults_Mal$TotCatchAtDecAgeLen))
  TotCatchAtDecAgeLen = TotCatchAtDecAgeLen_Fem + TotCatchAtDecAgeLen_Mal # single sex or combined sexes

  # expected retained catch proportions at length
  ExpRetCatchPropAtLen_Fem = RetCatchAtLen_Fem / sum(RetCatchAtLen_Fem)
  ExpRetCatchPropAtLen_Mal = RetCatchAtLen_Mal / sum(RetCatchAtLen_Mal)
  ExpRetCatchPropAtLen = RetCatchAtLen / sum(RetCatchAtLen)

  # expected discarded catch proportions at length
  ExpDiscCatchPropAtLen_Fem = DiscCatchAtLen_Fem / sum(DiscCatchAtLen_Fem)
  ExpDiscCatchPropAtLen_Mal = DiscCatchAtLen_Mal / sum(DiscCatchAtLen_Mal)
  ExpDiscCatchPropAtLen = DiscCatchAtLen / sum(DiscCatchAtLen)

  # get sizes of samples that need to be randomly-generated for retained and discarded fish
  SampSizeRes = GetReqdSampleSizesForRandData(MLL, RetenParams, SampleSize, RetCatchAtLen_Fem, RetCatchAtLen_Mal,
                                            DiscCatchAtLen_Fem, DiscCatchAtLen_Mal)

  SampleSize_Fem = SampSizeRes$SampleSize_Fem; SampleSize_Mal = SampSizeRes$SampleSize_Mal
  DiscSampleSize_Fem = SampSizeRes$DiscSampleSize_Fem; DiscSampleSize_Mal = SampSizeRes$DiscSampleSize_Mal
  DiscSampleSize = SampSizeRes$DiscSampleSize

  # generate observed retained catch frequencies at length
  ObsRetCatchFreqAtLen_Fem = as.vector(rmultinom(1, SampleSize_Fem, ExpRetCatchPropAtLen_Fem))
  ObsRetCatchFreqAtLen_Mal = as.vector(rmultinom(1, SampleSize_Mal, ExpRetCatchPropAtLen_Mal))
  ObsRetCatchFreqAtLen = ObsRetCatchFreqAtLen_Fem + ObsRetCatchFreqAtLen_Mal # combined sexes

  # generate observed discarded catch frequencies at length
  if (is.na(MLL) & is.na(RetenParams[1])) { # no discarding
    ObsDiscCatchFreqAtLen = NA
    ObsDiscCatchFreqAtLen_Fem = NA
    ObsDiscCatchFreqAtLen_Mal = NA
    ObsTotCatchAtLen = ObsRetCatchFreqAtLen
  } else { # discarding
    ObsDiscCatchFreqAtLen_Fem = as.vector(rmultinom(1, DiscSampleSize_Fem, ExpDiscCatchPropAtLen_Fem))
    ObsDiscCatchFreqAtLen_Mal = as.vector(rmultinom(1, DiscSampleSize_Mal, ExpDiscCatchPropAtLen_Mal))
    ObsDiscCatchFreqAtLen = ObsDiscCatchFreqAtLen_Fem + ObsDiscCatchFreqAtLen_Mal # combined sexes
    ObsTotCatchAtLen = ObsRetCatchFreqAtLen + ObsDiscCatchFreqAtLen
  }

  # expected retained catches at decimal age
  RetCatchAtDecAge_Fem = rowSums(RetCatchAtDecAgeLen_Fem); RetCatchAtDecAge_Mal = rowSums(RetCatchAtDecAgeLen_Mal)
  RetCatchAtDecAge = RetCatchAtDecAge_Fem + RetCatchAtDecAge_Mal

  # expected discarded catches at decimal age
  DiscCatchAtDecAge_Fem = rowSums(DiscCatchAtDecAgeLen_Fem); DiscCatchAtDecAge_Mal = rowSums(DiscCatchAtDecAgeLen_Mal)
  DiscCatchAtDecAge = DiscCatchAtDecAge_Fem + DiscCatchAtDecAge_Mal

  # expected discarded catches at decimal age
  TotCatchAtDecAge_Fem = rowSums(TotCatchAtDecAgeLen_Fem); TotCatchAtDecAge_Mal = rowSums(TotCatchAtDecAgeLen_Mal)
  TotCatchAtDecAge = TotCatchAtDecAge_Fem + TotCatchAtDecAge_Mal

  # expected retained catch proportions at decimal age
  ExpRetCatchPropAtDecAge_Fem = RetCatchAtDecAge_Fem / sum(RetCatchAtDecAge_Fem)
  ExpRetCatchPropAtDecAge_Mal = RetCatchAtDecAge_Mal / sum(RetCatchAtDecAge_Mal)
  ExpRetCatchPropAtDecAge = RetCatchAtDecAge / sum(RetCatchAtDecAge)

  # expected discarded catch proportions at decimal age
  ExpDiscCatchPropAtDecAge_Fem = DiscCatchAtDecAge_Fem / sum(DiscCatchAtDecAge_Fem)
  ExpDiscCatchPropAtDecAge_Mal = DiscCatchAtDecAge_Mal / sum(DiscCatchAtDecAge_Mal)
  ExpDiscCatchPropAtDecAge = DiscCatchAtDecAge / sum(DiscCatchAtDecAge)

  # expected catch proportions at length, in each age class
  CatchLenPropRes = CalcCatchLenPropGivenAge(nTimeSteps, nLenCl, midpt, RetCatchAtDecAgeLen_Fem, RetCatchAtDecAgeLen_Mal,
                                       DiscCatchAtDecAgeLen_Fem, DiscCatchAtDecAgeLen_Mal)
  ExpRetCatchPropLengthGivenDecAge_Fem = CatchLenPropRes$ExpRetCatchPropLengthGivenDecAge_Fem
  ExpRetCatchPropLengthGivenDecAge_Mal = CatchLenPropRes$ExpRetCatchPropLengthGivenDecAge_Mal
  ExpDiscCatchPropLengthGivenDecAge_Fem = CatchLenPropRes$ExpDiscCatchPropLengthGivenDecAge_Fem
  ExpDiscCatchPropLengthGivenDecAge_Mal = CatchLenPropRes$ExpDiscCatchPropLengthGivenDecAge_Mal

  # expected catch proportions at age, given length
  CatchAgePropRes = CalcCatchAgePropGivenLen(nTimeSteps, nLenCl, midpt, ExpRetCatchPropLengthGivenDecAge_Fem, ExpRetCatchPropAtDecAge_Fem,
                                       ExpRetCatchPropLengthGivenDecAge_Mal, ExpRetCatchPropAtDecAge_Mal,
                                       ExpDiscCatchPropLengthGivenDecAge_Fem, ExpDiscCatchPropAtDecAge_Fem,
                                       ExpDiscCatchPropLengthGivenDecAge_Mal, ExpDiscCatchPropAtDecAge_Mal)
  ExpRetCatchPropDecAgeGivenLength_Fem = CatchAgePropRes$ExpRetCatchPropDecAgeGivenLength_Fem
  ExpRetCatchPropDecAgeGivenLength_Mal = CatchAgePropRes$ExpRetCatchPropDecAgeGivenLength_Mal
  ExpDiscCatchPropDecAgeGivenLength_Fem = CatchAgePropRes$ExpDiscCatchPropDecAgeGivenLength_Fem
  ExpDiscCatchPropDecAgeGivenLength_Mal = CatchAgePropRes$ExpDiscCatchPropDecAgeGivenLength_Mal

  # generate random observed catch frequencies in length classes and ages
  CatchFreqRes = GenerateCatchFreqData(TimeStep, nTimeSteps, midpt, nLenCl, MLL, RetenParams, SampleSize_Fem,
                                    SampleSize_Mal, SampleSize, DiscSampleSize_Fem, DiscSampleSize_Mal, DiscSampleSize,
                                    ObsRetCatchFreqAtLen_Fem, ExpRetCatchPropDecAgeGivenLength_Fem,
                                    ObsRetCatchFreqAtLen_Mal, ExpRetCatchPropDecAgeGivenLength_Mal,
                                    ObsDiscCatchFreqAtLen_Fem, ExpDiscCatchPropDecAgeGivenLength_Fem,
                                    ObsDiscCatchFreqAtLen_Mal, ExpDiscCatchPropDecAgeGivenLength_Mal)

  # get randomly generated observed lengths, given random length class data
  ObsLenClRetCatchMidPt_Fem = CatchFreqRes$ObsLenClRetCatchMidPt_Fem
  ObsLenClRetCatchMidPt_Mal = CatchFreqRes$ObsLenClRetCatchMidPt_Mal
  ObsLenClDiscCatchMidPt_Fem = CatchFreqRes$ObsLenClDiscCatchMidPt_Fem
  ObsLenClDiscCatchMidPt_Mal = CatchFreqRes$ObsLenClDiscCatchMidPt_Mal
  RandCatchLenRes = GetRandFishLengths(MLL, RetenParams, LenInc, ObsLenClRetCatchMidPt_Fem, ObsLenClRetCatchMidPt_Mal,
                                 ObsLenClDiscCatchMidPt_Fem, ObsLenClDiscCatchMidPt_Mal, SampleSize_Fem, SampleSize_Mal,
                                 DiscSampleSize_Fem, DiscSampleSize_Mal)

  # get observed retained catch frequencies at length for each integer age
  ObsAgeClRetCatch_Fem = CatchFreqRes$ObsAgeClRetCatch_Fem
  ObsAgeClRetCatch_Mal = CatchFreqRes$ObsAgeClRetCatch_Mal
  RetCatchRes = GetRetCatchFreqAtLenAndIntAge(TimeStep, MaxAge, midpt, nLenCl, ObsAgeClRetCatch_Fem, ObsAgeClRetCatch_Mal,
                                         ObsLenClRetCatchMidPt_Fem, ObsLenClRetCatchMidPt_Mal)

  # get observed discarded catch frequencies at length for each integer age
  ObsAgeClDiscCatch_Fem = CatchFreqRes$ObsAgeClDiscCatch_Fem
  ObsAgeClDiscCatch_Mal = CatchFreqRes$ObsAgeClDiscCatch_Mal
  DiscCatchRes = GetDiscCatchFreqAtLenAndIntAge(TimeStep, MaxAge, midpt, nLenCl, ObsAgeClDiscCatch_Fem, ObsAgeClDiscCatch_Mal,
                                              ObsLenClDiscCatchMidPt_Fem, ObsLenClDiscCatchMidPt_Mal)

  # get summary statistics relating to retained fish
  RetStats = GetRetainedCatchStats(ObsRetCatchFreqAtLen, ObsTotCatchAtLen, ObsLenClRetCatchMidPt_Fem, ObsLenClRetCatchMidPt_Mal,
                                   ObsAgeClRetCatch_Fem, ObsAgeClRetCatch_Mal)

  # get summary statistics relating to discarded fish
  DiscStats = GetDiscardCatchStats(MLL, RetenParams, ObsDiscCatchFreqAtLen, ObsTotCatchAtLen, ObsLenClDiscCatchMidPt_Fem,
                                   ObsLenClDiscCatchMidPt_Mal, ObsAgeClDiscCatch_Fem, ObsAgeClDiscCatch_Mal)

  SummaryStats = list(ObsPropRetFish = RetStats$ObsPropRetFish,
                      ObsPropRetFish_sd = RetStats$ObsPropRetFish_sd,
                      ObsMeanLenRetFish_Fem = RetStats$ObsMeanLenRetFish_Fem,
                      ObsMeanLenRetFish_Mal = RetStats$ObsMeanLenRetFish_Mal,
                      ObsMeanLenRetFish = RetStats$ObsMeanLenRetFish,
                      ObsLenRetFish_sd_Fem = RetStats$ObsLenRetFish_sd_Fem,
                      ObsLenRetFish_sd_Mal = RetStats$ObsLenRetFish_sd_Mal,
                      ObsLenRetFish_sd = RetStats$ObsLenRetFish_sd,
                      ObsMeanLenRetFish_sd_Fem = RetStats$ObsMeanLenRetFish_sd_Fem,
                      ObsMeanLenRetFish_sd_Mal = RetStats$ObsMeanLenRetFish_sd_Mal,
                      ObsMeanLenRetFish_sd = RetStats$ObsMeanLenRetFish_sd,
                      ObsMeanAgeRetFish_Fem = RetStats$ObsMeanAgeRetFish_Fem,
                      ObsMeanAgeRetFish_Mal = RetStats$ObsMeanAgeRetFish_Mal,
                      ObsMeanAgeRetFish = RetStats$ObsMeanAgeRetFish,
                      ObsAgeRetFish_sd_Fem = RetStats$ObsAgeRetFish_sd_Fem,
                      ObsAgeRetFish_sd_Mal = RetStats$ObsAgeRetFish_sd_Mal,
                      ObsAgeRetFish_sd = RetStats$ObsAgeRetFish_sd,
                      ObsMeanAgeRetFish_sd_Fem = RetStats$ObsMeanAgeRetFish_sd_Fem,
                      ObsMeanAgeRetFish_sd_Mal = RetStats$ObsMeanAgeRetFish_sd_Mal,
                      ObsMeanAgeRetFish_sd = RetStats$ObsMeanAgeRetFish_sd,
                      ObsPropDiscFish = DiscStats$ObsPropDiscFish,
                      ObsPropDiscFish_sd = DiscStats$ObsPropDiscFish_sd,
                      ObsMeanLenDiscFish_Fem = DiscStats$ObsMeanLenDiscFish_Fem,
                      ObsMeanLenDiscFish_Mal = DiscStats$ObsMeanLenDiscFish_Mal,
                      ObsMeanLenDiscFish = DiscStats$ObsMeanLenDiscFish,
                      ObsLenDiscFish_sd_Fem = DiscStats$ObsLenDiscFish_sd_Fem,
                      ObsLenDiscFish_sd_Mal = DiscStats$ObsLenDiscFish_sd_Mal,
                      ObsLenDiscFish_sd = DiscStats$ObsLenDiscFish_sd,
                      ObsMeanLenDiscFish_sd_Fem = DiscStats$ObsMeanLenDiscFish_sd_Fem,
                      ObsMeanLenDiscFish_sd_Mal = DiscStats$ObsMeanLenDiscFish_sd_Mal,
                      ObsMeanLenDiscFish_sd = DiscStats$ObsMeanLenDiscFish_sd,
                      ObsMeanAgeDiscFish_Fem = DiscStats$ObsMeanAgeDiscFish_Fem,
                      ObsMeanAgeDiscFish_Mal = DiscStats$ObsMeanAgeDiscFish_Mal,
                      ObsMeanAgeDiscFish = DiscStats$ObsMeanAgeDiscFish,
                      ObsAgeDiscFish_sd_Fem = DiscStats$ObsAgeDiscFish_sd_Fem,
                      ObsAgeDiscFish_sd_Mal = DiscStats$ObsAgeDiscFish_sd_Mal,
                      ObsAgeDiscFish_sd = DiscStats$ObsAgeDiscFish_sd,
                      ObsMeanAgeDiscFish_sd_Fem = DiscStats$ObsMeanAgeDiscFish_sd_Fem,
                      ObsMeanAgeDiscFish_sd_Mal = DiscStats$ObsMeanAgeDiscFish_sd_Mal,
                      ObsMeanAgeDiscFish_sd = DiscStats$ObsMeanAgeDiscFish_sd)

  ModelDiag = list(DecAges = DecAges,
                   AgeCl = AgeCl,
                   RecLenDist = RecLenDist,
                   MeanSizeAtAge = MeanSizeAtAge,
                   TimestepGrowthSizeInc = TimestepGrowthSizeInc,
                   SelAtLength = SelAtLength,
                   RetAtLength = RetAtLength,
                   SelLandAtLength = SelLandAtLength,
                   FAtLen_Fem = CatchCurveResults_Fem$FAtLen,
                   ZAtLen_Fem = CatchCurveResults_Fem$ZAtLen,
                   FAtLenCapt_Fem = CatchCurveResults_Fem$FAtLenCapt,
                   FAtLenReten_Fem = CatchCurveResults_Fem$FAtLenReten,
                   FAtLenDisc_Fem = CatchCurveResults_Fem$FAtLenDisc,
                   FAtLen_Mal = CatchCurveResults_Mal$FAtLen,
                   ZAtLen_Mal = CatchCurveResults_Mal$ZAtLen,
                   FAtLenCapt_Mal = CatchCurveResults_Mal$FAtLenCapt,
                   FAtLenReten_Mal = CatchCurveResults_Mal$FAtLenReten,
                   FAtLenDisc_Mal = CatchCurveResults_Mal$FAtLenDisc,
                   DiscCatchAtLen = DiscCatchAtLen,
                   RetCatchAtLen = RetCatchAtLen,
                   TotCatchAtLen = TotCatchAtLen,
                   SampleSize_Fem = SampleSize_Fem,
                   SampleSize_Mal = SampleSize_Mal,
                   DiscSampleSize_Fem = DiscSampleSize_Fem,
                   DiscSampleSize_Mal = DiscSampleSize_Mal,
                   DiscSampleSize = DiscSampleSize,
                   RetCatchAtDecAge_Fem = RetCatchAtDecAge_Fem,
                   RetCatchAtDecAge_Mal = RetCatchAtDecAge_Mal,
                   RetCatchAtDecAge = RetCatchAtDecAge,
                   DiscCatchAtDecAge_Fem = DiscCatchAtDecAge_Fem,
                   DiscCatchAtDecAge_Mal = DiscCatchAtDecAge_Mal,
                   DiscCatchAtDecAge = DiscCatchAtDecAge,
                   TotCatchAtDecAge_Fem = TotCatchAtDecAge_Fem,
                   TotCatchAtDecAge_Mal = TotCatchAtDecAge_Mal,
                   TotCatchAtDecAge = TotCatchAtDecAge,
                   ExpRetCatchPropAtLen_Fem = ExpRetCatchPropAtLen_Fem,
                   ExpRetCatchPropAtLen_Mal = ExpRetCatchPropAtLen_Mal,
                   ExpRetCatchPropAtLen = ExpRetCatchPropAtLen,
                   ExpDiscCatchPropAtLen_Fem = ExpDiscCatchPropAtLen_Fem,
                   ExpDiscCatchPropAtLen_Mal = ExpDiscCatchPropAtLen_Mal,
                   ExpDiscCatchPropAtLen = ExpDiscCatchPropAtLen,
                   ExpRetCatchPropAtDecAge = ExpRetCatchPropAtDecAge,
                   ExpRetCatchPropAtDecAge_Fem = ExpRetCatchPropAtDecAge_Fem,
                   ExpRetCatchPropAtDecAge_Mal = ExpRetCatchPropAtDecAge_Mal,
                   ExpRetCatchPropLengthGivenDecAge_Fem = ExpRetCatchPropLengthGivenDecAge_Fem,
                   ExpRetCatchPropLengthGivenDecAge_Mal = ExpRetCatchPropLengthGivenDecAge_Mal,
                   ExpRetCatchPropDecAgeGivenLength_Fem = ExpRetCatchPropDecAgeGivenLength_Fem,
                   ExpRetCatchPropDecAgeGivenLength_Mal = ExpRetCatchPropDecAgeGivenLength_Mal,
                   ObsDiscCatchFreqAtLengthAndDecAge_Fem = CatchFreqRes$ObsDiscCatchFreqAtLengthAndDecAge_Fem,
                   ObsDiscCatchFreqAtLengthAndDecAge_Mal = CatchFreqRes$ObsDiscCatchFreqAtLengthAndDecAge_Mal,
                   ObsDiscCatchFreqAtLengthAndDecAge = CatchFreqRes$ObsDiscCatchFreqAtLengthAndDecAge,
                   ObsRetCatchFreqAtLengthAndIntAge_Fem = RetCatchRes$ObsRetCatchFreqAtLengthAndIntAge_Fem,
                   ObsRetCatchFreqAtLengthAndIntAge_Mal = RetCatchRes$ObsRetCatchFreqAtLengthAndIntAge_Mal,
                   ObsRetCatchFreqAtLengthAndIntAge = RetCatchRes$ObsRetCatchFreqAtLengthAndIntAge,
                   ObsRetCatchFreqAtIntAge_Fem = RetCatchRes$ObsRetCatchFreqAtIntAge_Fem,
                   ObsRetCatchFreqAtIntAge_Mal = RetCatchRes$ObsRetCatchFreqAtIntAge_Mal,
                   ObsRetCatchFreqAtIntAge = RetCatchRes$ObsRetCatchFreqAtIntAge,
                   ObsDiscCatchFreqAtLengthAndIntAge_Fem = DiscCatchRes$ObsDiscCatchFreqAtLengthAndIntAge_Fem,
                   ObsDiscCatchFreqAtLengthAndIntAge_Mal = DiscCatchRes$ObsDiscCatchFreqAtLengthAndIntAge_Mal,
                   ObsDiscCatchFreqAtLengthAndIntAge = DiscCatchRes$ObsDiscCatchFreqAtLengthAndIntAge,
                   ObsDiscCatchFreqAtIntAge_Fem = DiscCatchRes$ObsDiscCatchFreqAtIntAge_Fem,
                   ObsDiscCatchFreqAtIntAge_Mal = DiscCatchRes$ObsDiscCatchFreqAtIntAge_Mal,
                   ObsDiscCatchFreqAtIntAge = DiscCatchRes$ObsDiscCatchFreqAtIntAge,
                   ObsAgeClRetCatch = CatchFreqRes$ObsAgeClRetCatch,
                   ObsAgeClRetCatch_Fem = CatchFreqRes$ObsAgeClRetCatch_Fem,
                   ObsAgeClRetCatch_Mal = CatchFreqRes$ObsAgeClRetCatch_Mal,
                   ObsAgeClDiscCatch = CatchFreqRes$ObsAgeClDiscCatch,
                   ObsAgeClDiscCatch_Fem = CatchFreqRes$ObsAgeClDiscCatch_Fem,
                   ObsAgeClDiscCatch_Mal = CatchFreqRes$ObsAgeClDiscCatch_Mal)

  Results = list(DecAges = DecAges,
                 AgeCl = AgeCl,
                 lbnd = lbnd,
                 midpt = midpt,
                 ubnd = ubnd,
                 ObsDecAgeRetCatch = CatchFreqRes$ObsDecAgeRetCatch,
                 ObsDecAgeRetCatch_Fem = CatchFreqRes$ObsDecAgeRetCatch_Fem,
                 ObsDecAgeRetCatch_Mal = CatchFreqRes$ObsDecAgeRetCatch_Mal,
                 ObsDecAgeDiscCatch = CatchFreqRes$ObsDecAgeDiscCatch,
                 ObsDecAgeDiscCatch_Fem = CatchFreqRes$ObsDecAgeDiscCatch_Fem,
                 ObsDecAgeDiscCatch_Mal = CatchFreqRes$ObsDecAgeDiscCatch_Mal,
                 ObsRandLenRetCatch = RandCatchLenRes$ObsRandLenRetCatch,
                 ObsRandLenRetCatch_Fem = RandCatchLenRes$ObsRandLenRetCatch_Fem,
                 ObsRandLenRetCatch_Mal = RandCatchLenRes$ObsRandLenRetCatch_Mal,
                 ObsRandLenDiscCatch = RandCatchLenRes$ObsRandLenDiscCatch,
                 ObsRandLenDiscCatch_Fem = RandCatchLenRes$ObsRandLenDiscCatch_Fem,
                 ObsRandLenDiscCatch_Mal = RandCatchLenRes$ObsRandLenDiscCatch_Mal,
                 ObsRetCatchFreqAtLen = ObsRetCatchFreqAtLen,
                 ObsRetCatchFreqAtLen_Fem = ObsRetCatchFreqAtLen_Fem,
                 ObsRetCatchFreqAtLen_Mal = ObsRetCatchFreqAtLen_Mal,
                 ObsDiscCatchFreqAtLen = ObsDiscCatchFreqAtLen,
                 ObsDiscCatchFreqAtLen_Fem = ObsDiscCatchFreqAtLen_Fem,
                 ObsDiscCatchFreqAtLen_Mal = ObsDiscCatchFreqAtLen_Mal,
                 ObsLenClRetCatchMidPt = CatchFreqRes$ObsLenClRetCatchMidPt,
                 ObsLenClRetCatchMidPt_Fem = CatchFreqRes$ObsLenClRetCatchMidPt_Fem,
                 ObsLenClRetCatchMidPt_Mal = CatchFreqRes$ObsLenClRetCatchMidPt_Mal,
                 ObsRetCatchFreqAtLengthAndDecAge_Fem = CatchFreqRes$ObsRetCatchFreqAtLengthAndDecAge_Fem,
                 ObsRetCatchFreqAtLengthAndDecAge_Mal = CatchFreqRes$ObsRetCatchFreqAtLengthAndDecAge_Mal,
                 ObsRetCatchFreqAtLengthAndDecAge = CatchFreqRes$ObsRetCatchFreqAtLengthAndDecAge,
                 ModelDiag = ModelDiag,
                 SummaryStats = SummaryStats)


  return(Results)

}



#' Produce plots of simulated age and length data
#'
#' @param MaxAge Maximum age
#' @param MaxLen Maximum length
#' @param SimRes Outputs of SimLenAndAgeFreqData function
#' @param PlotOpt # 0=all plots, 1=retained lengths at age, 2=retained plus discarded lengths at age, 3=length frequency, 4=age frequency
#'
#' @return plots of simulated age and length data
#'
#' @examples
#' # Simulate data
#' set.seed(123)
#' SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
#' set.seed(123)
#' MaxAge = 30
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = 0.2
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(400, 50) # L50, L95-L50 for retention
#' DiscMort = 0 # proportion of fish that die due to natural mortality
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.08
#' GrowthParams = c(Linf, vbK)
#' RefnceAges = NA
#' SimRes=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                             SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # plot data
#' PlotOpt=0 # 0=all plots, 1=retained lengths at age, 2=retained plus discarded lengths at age, 3=length frequency, 6
#' # 7=female length frequency, 8=male length frequency, 9=female age frequency, 10=male age frequency,
#' # 11=selectivity/retention, 12=F-at-age reten + disc, 13=F-at-age reten, 14=F-at-age disc
#' PlotSimLenAndAgeFreqData(MaxAge, MaxLen, SimRes, PlotOpt)
#' @export
PlotSimLenAndAgeFreqData <- function(MaxAge, MaxLen, SimRes, PlotOpt) {

  .pardefault <- par(no.readonly = TRUE) # store current par settings
  # combined sexes
  if (PlotOpt==0) {
    par(mfrow=c(2,2), mar=c(5,4,1,1))
  }
  if (PlotOpt!=0) {
    par(mfrow=c(1,1), mar=c(5,4,2,2))
  }

  # plot lengths and ages of retained fish
  if (PlotOpt==0 | PlotOpt == 1) {
    xaxis_lab = "Age (y)"
    yaxis_lab = "Length (mm)"
    xlims = Get_xaxis_scale(c(0,MaxAge))
    xmax = xlims$xmax
    xint = xlims$xint
    ylims = Get_yaxis_scale(c(0,MaxLen))
    ymax = ylims$ymax
    yint = ylims$yint

    plot(SimRes$ObsDecAgeRetCatch, SimRes$ObsRandLenRetCatch, "p", main=NA, cex.main=1.2, pch=16, cex=0.6, xaxt = "n", yaxt = "n",
         xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax))
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
    legend("bottomright", legend="Ret. fish", inset=c(0.13,0),
           lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col="black")
    if (is.vector(SimRes$ModelDiag$MeanSizeAtAge)) {
      sm1 = spline(SimRes$ModelDiag$DecAges, SimRes$ModelDiag$MeanSizeAtAge, n=100, method="natural")
      lines(sm1$x, sm1$y, col="black")
      # legend("topleft", legend=c("Growth - comb. sexes"), inset=c(0.13,0),
      #        lty=1, cex = 0.8, bty="n", seg.len = 2, pch=-16, col="black")
    } else {
      sm1 = spline(SimRes$ModelDiag$DecAges, SimRes$ModelDiag$MeanSizeAtAge[1,], n=100, method="natural")
      lines(sm1$x, sm1$y, col="red")
      sm1 = spline(SimRes$ModelDiag$DecAges, SimRes$ModelDiag$MeanSizeAtAge[2,], n=100, method="natural")
      lines(sm1$x, sm1$y, col="blue")
      legend("topleft", legend=c("Growth - females","Growth - males"), inset=c(0.13,0),
             lty=1, cex = 0.8, bty="n", seg.len = 2, pch=-16, col=c("red","blue"))
    }
  }

  # plot lengths and ages of retained and discarded fish
  if (PlotOpt==0 | PlotOpt == 2) {
    xaxis_lab = "Age (y)"
    yaxis_lab = "Length (mm)"
    xlims = Get_xaxis_scale(c(0,MaxAge))
    xmax = xlims$xmax
    xint = xlims$xint
    ylims = Get_yaxis_scale(c(0,MaxLen))
    ymax = ylims$ymax
    yint = ylims$yint

    plot(SimRes$ObsDecAgeRetCatch, SimRes$ObsRandLenRetCatch, "p", main=NA, cex.main=1.2, pch=16, cex=0.6, xaxt = "n", yaxt = "n",
         xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax))
    points(SimRes$ObsDecAgeDiscCatch, SimRes$ObsRandLenDiscCatch, col='grey', pch=16, cex=0.6, )
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
    legend("bottomright", legend=c("Ret. fish","Disc. fish"), inset=c(0.13,0),
           lty=1, cex = 0.8, bty="n", seg.len = 0, pch=16, col=c("black","grey"))
    if (is.vector(SimRes$ModelDiag$MeanSizeAtAge)) {
      sm1 = spline(SimRes$ModelDiag$DecAges, SimRes$ModelDiag$MeanSizeAtAge, n=100, method="natural")
      lines(sm1$x, sm1$y, col="black")
      # legend("topleft", legend=c("Growth - comb. sexes"), inset=c(0.13,0),
      #        lty=1, cex = 0.8, bty="n", seg.len = 2, pch=-16, col="black")
    } else {
    sm1 = spline(SimRes$ModelDiag$DecAges, SimRes$ModelDiag$MeanSizeAtAge[1,], n=100, method="natural")
    lines(sm1$x, sm1$y, col="red")
    sm1 = spline(SimRes$ModelDiag$DecAges, SimRes$ModelDiag$MeanSizeAtAge[2,], n=100, method="natural")
    lines(sm1$x, sm1$y, col="blue")
    legend("topleft", legend=c("Growth - females","Growth - males"), inset=c(0.13,0),
           lty=1, cex = 0.8, bty="n", seg.len = 2, pch=-16, col=c("red","blue"))
    }
  }

  # plot catch frequency at length
  if (PlotOpt==0 | PlotOpt == 3) {
    xaxis_lab = "Length (mm)"
    yaxis_lab = "Number of fish"
    xlims = Get_xaxis_scale(c(0,MaxLen))
    xmax = xlims$xmax
    xint = xlims$xint
    if(!is.na(SimRes$ObsDiscCatchFreqAtLen[1])) {
      MaxFreq=max(SimRes$ObsRetCatchFreqAtLen+SimRes$ObsDiscCatchFreqAtLen)
    } else {
      MaxFreq=max(SimRes$ObsRetCatchFreqAtLen)
    }
    ylims = Get_yaxis_scale(c(0,MaxFreq))
    ymax = ylims$ymax
    yint = ylims$yint

    plot(SimRes$midpt, SimRes$ObsRetCatchFreqAtLen, "l", main=NA, cex.main=1.2, pch=16, cex=1, xaxt = "n", yaxt = "n",
         xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax), col="red")
    if(!is.na(SimRes$ObsDiscCatchFreqAtLen[1])) lines(SimRes$midpt, SimRes$ObsDiscCatchFreqAtLen, "l", col='blue', pch=16, cex=1)
    if(!is.na(SimRes$ObsDiscCatchFreqAtLen[1])) lines(SimRes$midpt, SimRes$ObsRetCatchFreqAtLen+SimRes$ObsDiscCatchFreqAtLen, "l", col='black', pch=16, cex=1)
    lines(SimRes$midpt, SimRes$ObsRetCatchFreqAtLen, "l", col='red', pch=16, cex=1,lwd=2)
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
    if(!is.na(SimRes$ObsDiscCatchFreqAtLen[1])) {
      legend("topright", legend=c("Ret. fish","Disc. fish","All fish"), inset=c(0.13,0),
             lty=1, cex = 0.8, bty="n", seg.len = 2, lwd=c(2,1,1), pch=-16, col=c("red","blue","black"))
    } else {
      legend("topright", legend="Ret. fish", inset=c(0.13,0),
             lty=1, cex = 0.8, bty="n", seg.len = 2, lwd=2, pch=-16, col="red")
    }
  }

  # plot catch frequency at age
  if (PlotOpt==0 | PlotOpt == 4) {
    xaxis_lab = "Age (y)"
    yaxis_lab = "Number of fish"
    xlims = Get_xaxis_scale(c(0,MaxAge))
    xmax = xlims$xmax
    xint = xlims$xint
    if(sum(SimRes$ModelDiag$ObsDiscCatchFreqAtIntAge)>0) {
      MaxFreq=max(SimRes$ModelDiag$ObsRetCatchFreqAtIntAge+SimRes$ModelDiag$ObsDiscCatchFreqAtIntAge)
    } else {
      MaxFreq=max(SimRes$ModelDiag$ObsRetCatchFreqAtIntAge)
    }
    ylims = Get_yaxis_scale(c(0,MaxFreq))
    ymax = ylims$ymax
    yint = ylims$yint
    plot(SimRes$AgeCl, SimRes$ModelDiag$ObsRetCatchFreqAtIntAge, "l", main=NA, cex.main=1.2, pch=16, cex=1, xaxt = "n", yaxt = "n",
         xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax), col="red")
    if(sum(SimRes$ModelDiag$ObsDiscCatchFreqAtIntAge)>0) lines(SimRes$AgeCl, SimRes$ModelDiag$ObsDiscCatchFreqAtIntAge, "l", col='blue', pch=16, cex=1)
    if(sum(SimRes$ModelDiag$ObsDiscCatchFreqAtIntAge)>0) lines(SimRes$AgeCl, SimRes$ModelDiag$ObsRetCatchFreqAtIntAge+SimRes$ModelDiag$ObsDiscCatchFreqAtIntAge, "l", col='black', pch=16, cex=1)
    lines(SimRes$AgeCl, SimRes$ModelDiag$ObsRetCatchFreqAtIntAge, "l", col='red', pch=16, cex=1,lwd=2)
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
    if(sum(SimRes$ModelDiag$ObsDiscCatchFreqAtIntAge)>0) {
      legend("topright", legend=c("Ret. fish","Disc. fish","All fish"), inset=c(0.13,0),
             lty=1, cex = 0.8, bty="n", seg.len = 2, lwd=c(2,1,1), pch=-16, col=c("red","blue","black"))
    } else {
      legend("topright", legend="Ret. fish", inset=c(0.13,0),
             lty=1, cex = 0.8, bty="n", seg.len = 2, lwd=2, pch=-16, col="red")
    }
  }


  if (is.data.frame(SimRes$ModelDiag$MeanSizeAtAge)) {

    # separate sexes
    if (PlotOpt==0) {
      par(mfrow=c(3,2), mar=c(5,4,1,1))
    }
    if (PlotOpt!=0) {
      par(mfrow=c(1,1), mar=c(5,4,2,2))
    }

    # females
    # plot lengths and ages of retained fish
    if (PlotOpt==0 | PlotOpt == 5) {
      xaxis_lab = "Age (y)"
      yaxis_lab = "Length (mm)"
      xlims = Get_xaxis_scale(c(0,MaxAge))
      xmax = xlims$xmax
      xint = xlims$xint
      ylims = Get_yaxis_scale(c(0,MaxLen))
      ymax = ylims$ymax
      yint = ylims$yint

      plot(SimRes$ObsDecAgeRetCatch_Fem, SimRes$ObsRandLenRetCatch_Fem, "p", main=NA, cex.main=1.2, pch=16, cex=0.6, xaxt = "n", yaxt = "n",
           xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax), col="red")
      points(SimRes$ObsDecAgeRetCatch_Mal, SimRes$ObsRandLenRetCatch_Mal, cex=0.6, pch=16, col="blue")
      points(SimRes$ObsDecAgeRetCatch_Fem, SimRes$ObsRandLenRetCatch_Fem, cex=0.6, pch=16, col="red")
      AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
      legend("bottomright", legend=c("Ret. females","Ret. males"), inset=c(0.13,0),
             lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16), col=c("red","blue"))
      sm1 = spline(SimRes$ModelDiag$DecAges, SimRes$ModelDiag$MeanSizeAtAge[1,], n=100, method="natural")
      lines(sm1$x, sm1$y, col="red")
      sm1 = spline(SimRes$ModelDiag$DecAges, SimRes$ModelDiag$MeanSizeAtAge[2,], n=100, method="natural")
      lines(sm1$x, sm1$y, col="blue")
      legend("topleft", legend=c("Growth - females","Growth - males"), inset=c(0.13,0),
             lty=1, cex = 0.8, bty="n", seg.len = 2, pch=-16, col=c("red","blue"))
    }

    # plot lengths and ages of retained and discarded fish
    if (PlotOpt==0 | PlotOpt == 6) {
      xaxis_lab = "Age (y)"
      yaxis_lab = "Length (mm)"
      xlims = Get_xaxis_scale(c(0,MaxAge))
      xmax = xlims$xmax
      xint = xlims$xint
      ylims = Get_yaxis_scale(c(0,MaxLen))
      ymax = ylims$ymax
      yint = ylims$yint

      plot(SimRes$ObsDecAgeRetCatch_Fem, SimRes$ObsRandLenRetCatch_Fem, "p", main=NA, cex.main=1.2, pch=16, cex=0.6, xaxt = "n", yaxt = "n",
           xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax), col="red")
      points(SimRes$ObsDecAgeRetCatch_Mal, SimRes$ObsRandLenRetCatch_Mal, col='blue', pch=16, cex=0.6)
      points(SimRes$ObsDecAgeRetCatch_Fem, SimRes$ObsRandLenRetCatch_Fem, col='red', pch=16, cex=0.6)
      points(SimRes$ObsDecAgeDiscCatch_Mal, SimRes$ObsRandLenDiscCatch_Mal, col='lightblue', pch=16, cex=0.6)
      points(SimRes$ObsDecAgeDiscCatch_Fem, SimRes$ObsRandLenDiscCatch_Fem, col='pink', pch=16, cex=0.6)
      AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
      legend("bottomright", legend=c("Ret. females","Ret. males","Disc. females","Disc. males"), inset=c(0.13,0),
             lty=1, cex = 0.8, bty="n", seg.len = 0, pch=16, col=c("red","blue","pink","lightblue"))
      sm1 = spline(SimRes$ModelDiag$DecAges, SimRes$ModelDiag$MeanSizeAtAge[1,], n=100, method="natural")
      lines(sm1$x, sm1$y, col="red")
      sm1 = spline(SimRes$ModelDiag$DecAges, SimRes$ModelDiag$MeanSizeAtAge[2,], n=100, method="natural")
      lines(sm1$x, sm1$y, col="blue")
      legend("topleft", legend=c("Growth - females","Growth - males"), inset=c(0.13,0),
             lty=1, cex = 0.8, bty="n", seg.len = 2, pch=-16, col=c("red","blue"))
    }

    # plot catch frequency at length - females
    if (PlotOpt==0 | PlotOpt == 7) {
      xaxis_lab = "Length (mm)"
      yaxis_lab = "Number of fish"
      xlims = Get_xaxis_scale(c(0,MaxLen))
      xmax = xlims$xmax
      xint = xlims$xint
      if(!is.na(SimRes$ObsDiscCatchFreqAtLen[1])) {
        MaxFreq=max(SimRes$ObsRetCatchFreqAtLen_Fem+SimRes$ObsDiscCatchFreqAtLen_Fem)
      } else {
        MaxFreq=max(SimRes$ObsRetCatchFreqAtLen_Fem)
      }
      ylims = Get_yaxis_scale(c(0,MaxFreq))
      ymax = ylims$ymax
      yint = ylims$yint

      plot(SimRes$midpt, SimRes$ObsRetCatchFreqAtLen_Fem, "l", main=NA, cex.main=1.2, pch=16, cex=1, xaxt = "n", yaxt = "n",
           xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax), col="red")
      if(!is.na(SimRes$ObsDiscCatchFreqAtLen_Fem[1])) lines(SimRes$midpt, SimRes$ObsDiscCatchFreqAtLen_Fem, "l", col='pink', pch=16, cex=1)
      lines(SimRes$midpt, SimRes$ObsRetCatchFreqAtLen_Fem, "l", col='red', pch=16, cex=1,lwd=2)
      if(!is.na(SimRes$ObsDiscCatchFreqAtLen_Fem[1])) lines(SimRes$midpt, SimRes$ObsRetCatchFreqAtLen_Fem+SimRes$ObsDiscCatchFreqAtLen_Fem, "l", col='black', pch=16, cex=1)
      AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
      if(!is.na(SimRes$ObsDiscCatchFreqAtLen[1])) {
        legend("topright", legend=c("Ret. females","Disc. females","All females"), inset=c(0.13,0),
               lty=1, cex = 0.8, bty="n", seg.len = 2, lwd=c(2,1,1), pch=-16, col=c("red","pink","black"))
      } else {
        legend("topright", legend="Ret. females", inset=c(0.13,0),
               lty=1, cex = 0.8, bty="n", seg.len = 2, lwd=2, pch=-16, col="red")
      }
    }

    # plot catch frequency at length - males
    if (PlotOpt==0 | PlotOpt == 8) {
      xaxis_lab = "Length (mm)"
      yaxis_lab = "Number of fish"
      xlims = Get_xaxis_scale(c(0,MaxLen))
      xmax = xlims$xmax
      xint = xlims$xint
      if(!is.na(SimRes$ObsDiscCatchFreqAtLen[1])) {
        MaxFreq=max(SimRes$ObsRetCatchFreqAtLen_Mal+SimRes$ObsDiscCatchFreqAtLen_Mal)
      } else {
        MaxFreq=max(SimRes$ObsRetCatchFreqAtLen_Mal)
      }
      ylims = Get_yaxis_scale(c(0,MaxFreq))
      ymax = ylims$ymax
      yint = ylims$yint

      plot(SimRes$midpt, SimRes$ObsRetCatchFreqAtLen_Mal, "l", main=NA, cex.main=1.2, pch=16, cex=1, xaxt = "n", yaxt = "n",
           xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax), col="blue")
      if(!is.na(SimRes$ObsDiscCatchFreqAtLen_Mal[1])) lines(SimRes$midpt, SimRes$ObsDiscCatchFreqAtLen_Mal, "l", col='lightblue', pch=16, cex=1)
      lines(SimRes$midpt, SimRes$ObsRetCatchFreqAtLen_Mal, "l", col='blue', pch=16, cex=1,lwd=2)
      if(!is.na(SimRes$ObsDiscCatchFreqAtLen_Mal[1])) lines(SimRes$midpt, SimRes$ObsRetCatchFreqAtLen_Mal+SimRes$ObsDiscCatchFreqAtLen_Mal, "l", col='black', pch=16, cex=1)
      AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
      if(!is.na(SimRes$ObsDiscCatchFreqAtLen[1])) {
        legend("topright", legend=c("Ret. males","Disc. males","All males"), inset=c(0.13,0),
               lty=1, cex = 0.8, bty="n", seg.len = 2, lwd=c(2,1,1), pch=-16, col=c("blue","lightblue","black"))
      } else {
        legend("topright", legend="Ret. males", inset=c(0.13,0),
               lty=1, cex = 0.8, bty="n", seg.len = 2, lwd=2, pch=-16, col="blue")
      }
    }

    # plot catch frequency at age - females
    if (PlotOpt==0 | PlotOpt == 9) {
      xaxis_lab = "Age (y)"
      yaxis_lab = "Number of fish"
      xlims = Get_xaxis_scale(c(0,MaxAge))
      xmax = xlims$xmax
      xint = xlims$xint
      if(sum(SimRes$ModelDiag$ObsDiscCatchFreqAtIntAge_Fem)>0) {
        MaxFreq=max(SimRes$ModelDiag$ObsRetCatchFreqAtIntAge_Fem+SimRes$ModelDiag$ObsDiscCatchFreqAtIntAge_Fem)
      } else {
        MaxFreq=max(SimRes$ModelDiag$ObsRetCatchFreqAtIntAge_Fem)
      }
      ylims = Get_yaxis_scale(c(0,MaxFreq))
      ymax = ylims$ymax
      yint = ylims$yint
      plot(SimRes$AgeCl, SimRes$ModelDiag$ObsRetCatchFreqAtIntAge_Fem, "l", main=NA, cex.main=1.2, pch=16, cex=1, xaxt = "n", yaxt = "n",
           xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax), col="red")
      if(sum(SimRes$ModelDiag$ObsDiscCatchFreqAtIntAge_Fem)>0) lines(SimRes$AgeCl, SimRes$ModelDiag$ObsDiscCatchFreqAtIntAge_Fem, "l", col='pink', pch=16, cex=1)
      if(sum(SimRes$ModelDiag$ObsDiscCatchFreqAtIntAge_Fem)>0) lines(SimRes$AgeCl, SimRes$ModelDiag$ObsRetCatchFreqAtIntAge_Fem+SimRes$ModelDiag$ObsDiscCatchFreqAtIntAge_Fem, "l", col='black', pch=16, cex=1)
      lines(SimRes$AgeCl, SimRes$ModelDiag$ObsRetCatchFreqAtIntAge_Fem, "l", col='red', pch=16, cex=1,lwd=2)
      AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
      if(sum(SimRes$ModelDiag$ObsDiscCatchFreqAtIntAge_Fem)>0) {
        legend("topright", legend=c("Ret. females","Disc. females","All females"), inset=c(0.13,0),
               lty=1, cex = 0.8, bty="n", seg.len = 2, lwd=c(2,1,1), pch=-16, col=c("red","pink","black"))
      } else {
        legend("topright", legend="Ret. females", inset=c(0.13,0),
               lty=1, cex = 0.8, bty="n", seg.len = 2, lwd=2, pch=-16, col="red")
      }
    }


    # plot catch frequency at age - males
    if (PlotOpt==0 | PlotOpt == 10) {
      xaxis_lab = "Age (y)"
      yaxis_lab = "Number of fish"
      xlims = Get_xaxis_scale(c(0,MaxAge))
      xmax = xlims$xmax
      xint = xlims$xint
      if(sum(SimRes$ModelDiag$ObsDiscCatchFreqAtIntAge_Mal)>0) {
        MaxFreq=max(SimRes$ModelDiag$ObsRetCatchFreqAtIntAge_Mal+SimRes$ModelDiag$ObsDiscCatchFreqAtIntAge_Mal)
      } else {
        MaxFreq=max(SimRes$ModelDiag$ObsRetCatchFreqAtIntAge_Mal)
      }
      ylims = Get_yaxis_scale(c(0,MaxFreq))
      ymax = ylims$ymax
      yint = ylims$yint
      plot(SimRes$AgeCl, SimRes$ModelDiag$ObsRetCatchFreqAtIntAge_Mal, "l", main=NA, cex.main=1.2, pch=16, cex=1, xaxt = "n", yaxt = "n",
           xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax), col="blue")
      if(sum(SimRes$ModelDiag$ObsDiscCatchFreqAtIntAge_Mal)>0) lines(SimRes$AgeCl, SimRes$ModelDiag$ObsDiscCatchFreqAtIntAge_Mal, "l", col='lightblue', pch=16, cex=1)
      if(sum(SimRes$ModelDiag$ObsDiscCatchFreqAtIntAge_Mal)>0) lines(SimRes$AgeCl, SimRes$ModelDiag$ObsRetCatchFreqAtIntAge_Mal+SimRes$ModelDiag$ObsDiscCatchFreqAtIntAge_Mal, "l", col='black', pch=16, cex=1)
      lines(SimRes$AgeCl, SimRes$ModelDiag$ObsRetCatchFreqAtIntAge_Mal, "l", col='blue', pch=16, cex=1,lwd=2)
      AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
      if(sum(SimRes$ModelDiag$ObsDiscCatchFreqAtIntAge_Mal)>0) {
        legend("topright", legend=c("Ret. males","Disc. males","All males"), inset=c(0.13,0),
               lty=1, cex = 0.8, bty="n", seg.len = 2, lwd=c(2,1,1), pch=-16, col=c("blue","lightblue","black"))
      } else {
        legend("topright", legend="Ret. males", inset=c(0.13,0),
               lty=1, cex = 0.8, bty="n", seg.len = 2, lwd=2, pch=-16, col="blue")
      }
    }
  } else {
    if (PlotOpt > 4) cat("This plot option not available when growth inputted as combined sexes")
  }

  if (PlotOpt==0) {
    par(mfrow=c(2,2), mar=c(5,4,1,1))
  }
  if (PlotOpt!=0) {
    par(mfrow=c(1,1), mar=c(5,4,2,2))
  }

  # plot selectivity and retention at length
  if (PlotOpt==0 | PlotOpt == 11) {
    xaxis_lab = "Length (mm)"
    yaxis_lab = "Probability"
    xlims = Get_xaxis_scale(SimRes$midpt)
    xmax = xlims$xmax
    xint = xlims$xint
    ymax = 1.0
    yint = 0.2
    plot(SimRes$midpt, SimRes$ModelDiag$SelAtLength, "l", cex.main=1.2, pch=1, cex=0.6,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax),
         ylim=c(0,ymax), col="red")
    lines(SimRes$midpt, SimRes$ModelDiag$RetAtLength, col="blue")
    lines(SimRes$midpt, SimRes$ModelDiag$SelLandAtLength, col="black")
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
    if (SelectivityType==1) {
      legend("bottomright", pch=-1, legend=c("Sel_gear", "Prob_Reten", "Sel_land"), lty="solid", col=c("red","blue","black"),
             bty='n', cex=0.8,lwd=1, y.intersp=1.1)
    }
    if (SelectivityType==2) { # logistic selectivity
      if(!is.na(SelParams[1])) {
        L50_Est=paste("L50_gear =",round(SelParams[1],0),"mm")
        Delta_Est=paste("Delta_gear =",round(SelParams[2],0),"mm")
      }
      if(!is.na(RetenParams[1])) {
        L50_Est2=paste("L50_ret =",round(RetenParams[1],0),"mm")
        Delta_Est2=paste("Delta_ret =",round(RetenParams[2],0),"mm")
      }

      if(is.na(MLL)) {
        if(!is.na(SelParams[1]) & !is.na(RetenParams[1])) {
          lgnd = c("Sel_gear", "Prob_Reten", "Sel_land", L50_Est, Delta_Est, L50_Est2, Delta_Est2)
          cls = c("red","blue","black","black","black","black","black")
          lwds = c(1,1,1,-1,-1,-1,-1)
          legend("bottomright", pch=-1, legend=lgnd, lty="solid", col=cls, bty='n', cex=0.6, lwd=lwds, y.intersp=1.1)
        }
        if(!is.na(SelParams[1]) & is.na(RetenParams[1])) {
          lgnd = c("Sel_gear", "Prob_Reten", "Sel_land", L50_Est, Delta_Est)
          cls = c("red","blue","black","black","black")
          lwds = c(1,1,1,-1,-1)
          legend("bottomright", pch=-1, legend=lgnd, lty="solid", col=cls, bty='n', cex=0.6, lwd=lwds, y.intersp=1.1)
        }
        if(is.na(SelParams[1]) & !is.na(RetenParams[1])) {
          lgnd = c("Sel_gear", "Prob_Reten", "Sel_land", L50_Est2, Delta_Est2)
          cls = c("red","blue","black","black","black")
          lwds = c(1,1,1,-1,-1)
          legend("bottomright", pch=-1, legend=lgnd, lty="solid", col=cls, bty='n', cex=0.6, lwd=lwds, y.intersp=1.1)
        }
        if(is.na(SelParams[1]) & is.na(RetenParams[1])) {
          lgnd = c("Sel_gear", "Prob_Reten", "Sel_land")
          cls = c("red","blue","black")
          lwds = c(1,1,1)
          legend("bottomright", pch=-1, legend=lgnd, lty="solid", col=cls, bty='n', cex=0.6, lwd=lwds, y.intersp=1.1)
        }
      } # is.na(MLL)

      if (is.numeric(MLL)) {
        if(!is.na(SelParams[1]) & !is.na(RetenParams[1])) {
          lgnd = c("Sel_gear", "Prob_Reten", "Sel_land", L50_Est, Delta_Est, L50_Est2, Delta_Est2, paste("MLL =",MLL))
          cls = c("red","blue","black","black","black","black","black","black")
          lwds = c(1,1,1,-1,-1,-1,-1, -1)
          legend("bottomright", pch=-1, legend=lgnd, lty="solid", col=cls, bty='n', cex=0.6, lwd=lwds, y.intersp=1.1)
        }
        if(!is.na(SelParams[1]) & is.na(RetenParams[1])) {
          lgnd = c("Sel_gear", "Prob_Reten", "Sel_land", L50_Est, Delta_Est, paste("MLL =",MLL))
          cls = c("red","blue","black","black","black","black")
          lwds = c(1,1,1,-1,-1, -1)
          legend("bottomright", pch=-1, legend=lgnd, lty="solid", col=cls, bty='n', cex=0.6, lwd=lwds, y.intersp=1.1)
        }
        if(is.na(SelParams[1]) & !is.na(RetenParams[1])) {
          lgnd = c("Sel_gear", "Prob_Reten", "Sel_land", L50_Est2, Delta_Est2, paste("MLL =",MLL))
          cls = c("red","blue","black","black","black","black")
          lwds = c(1,1,1,-1,-1, -1)
          legend("bottomright", pch=-1, legend=lgnd, lty="solid", col=cls, bty='n', cex=0.6, lwd=lwds, y.intersp=1.1)
        }
        if(is.na(SelParams[1]) & is.na(RetenParams[1])) {
          lgnd = c("Sel_gear", "Prob_Reten", "Sel_land", paste("MLL =",MLL))
          cls = c("red","blue","black","black")
          lwds = c(1,1,1,-1)
          legend("bottomright", pch=-1, legend=lgnd, lty="solid", col=cls, bty='n', cex=0.6, lwd=lwds, y.intersp=1.1)
        }
      }
    }
  }

  # F at age, associated with retention and discard mortality
  if (PlotOpt==0 | PlotOpt == 12) {
    xlims = Get_xaxis_scale(SimRes$midpt)
    xmax = xlims$xmax
    xint = xlims$xint
    ylims = Get_yaxis_scale(c(0,SimRes$ModelDiag$FAtLen_Fem))
    ymax = ylims$ymax
    yint = ylims$yint
    yaxis_lab = expression(paste(italic("F") ~ (year^{-1})))
    xaxis_lab = "Length (mm)"
    plot(SimRes$midpt, SimRes$ModelDiag$FAtLen_Fem, "l", main="F (retention + discard mortality)", cex.main=1.0, pch=1, cex=0.6,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=NA, frame=F, xlim=c(0,xmax),
         ylim=c(0,ymax), col="red")
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
    mtext(expression(paste(italic("F") ~ (year^{-1}))),las=3,side=2,line=2.5,cex=1.0,lwd=1.75)
    lines(SimRes$midpt, SimRes$ModelDiag$FAtLen_Mal, col="blue")
    legend('topright', col=c("red","blue"),lty="solid",legend=c("females","males"),bty='n', cex=0.8,lwd=1.75)
  }

  # F at age, associated with retention
  if (PlotOpt==0 | PlotOpt == 13) {
    xlims = Get_xaxis_scale(SimRes$midpt)
    xmax = xlims$xmax
    xint = xlims$xint
    ylims = Get_yaxis_scale(c(0,SimRes$ModelDiag$FAtLen_Fem))
    ymax = ylims$ymax
    yint = ylims$yint
    yaxis_lab = expression(paste(italic("F") ~ (year^{-1})))
    xaxis_lab = "Length (mm)"
    plot(SimRes$midpt, SimRes$ModelDiag$FAtLenReten_Fem, "l", main="F (retention)", cex.main=1.0, pch=1, cex=0.6,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=NA, frame=F, xlim=c(0,xmax),
         ylim=c(0,ymax), col="blue")
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
    mtext(expression(paste(italic("F") ~ (year^{-1}))),las=3,side=2,line=3,cex=1.0,lwd=1.75)
    lines(SimRes$midpt, SimRes$ModelDiag$FAtLenReten_Mal, col="blue")
    legend('topright', col=c("red","blue"),lty="solid",legend=c("females","males"),bty='n', cex=0.8,lwd=1.75)
  }

  # F at age, associated with fish capture and post-release mortality
  if (PlotOpt==0 | PlotOpt == 14) {
    xlims = Get_xaxis_scale(SimRes$midpt)
    xmax = xlims$xmax
    xint = xlims$xint
    ylims = Get_yaxis_scale(c(0,SimRes$ModelDiag$FAtLen_Fem))
    ymax = ylims$ymax
    yint = ylims$yint
    yaxis_lab = expression(paste(italic("F") ~ (year^{-1})))
    xaxis_lab = "Length (mm)"
    plot(SimRes$midpt, SimRes$ModelDiag$FAtLenDisc_Fem*DiscMort, "l", main="F (discard mortality)", cex.main=1.0, pch=1, cex=0.6,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=NA, frame=F, xlim=c(0,xmax),
         ylim=c(0,ymax), col="brown")
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
    mtext(expression(paste(italic("F") ~ (year^{-1}))),las=3,side=2,line=2.5,cex=1.0,lwd=1.75)
    lines(SimRes$midpt, SimRes$ModelDiag$FAtLenDisc_Mal*DiscMort, col="blue")
    legend('topright', col=c("red","blue"),lty="solid",legend=c("females","males"),bty='n', cex=0.8,lwd=1.75)
  }

  par(.pardefault)

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
#' SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
#' set.seed(123)
#' MaxAge = 30
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = 0.2
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0 # proportion of fish that die due to natural mortality
#' # single sex, von Bertalanffy
#' DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.08
#' GrowthParams = c(Linf, vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # 2 sexes, von Bertalanffy
#' DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = c(700,850)
#' vbK = c(0.3,0.2)
#' CVSizeAtAge = c(0.08,0.08)
#' GrowthParams = data.frame(Linf=Linf, vbK=vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # 1 sex, Schnute
#' DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
#' GrowthCurveType = 2 # 1 = von Bertalanffy, 2 = Schnute
#' t1 = 0.5 # growth - Schnute
#' t2 = 25 # growth - Schnute
#' y1 = 100 # growth - Schnute
#' y2 = 1000 # growth - Schnute
#' a = 0.02 # growth - Schnute
#' b = 3.0 # growth - Schnute
#' GrowthParams = c(y1, y2, a, b)
#' RefnceAges = c(t1,t2)
#' CVSizeAtAge = 0.05
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # 2 sexes, Schnute
#' DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
#' GrowthCurveType = 2 # 1 = von Bertalanffy, 2 = Schnute
#' t1 = c(0.5,0.5) # growth - Schnute
#' t2 = c(25,25) # growth - Schnute
#' y1 = c(100,100) # growth - Schnute
#' y2 = c(1000,1000) # growth - Schnute
#' a = c(0.02,0.02) # growth - Schnute
#' b = c(3,3) # growth - Schnute
#' CVSizeAtAge = c(0.05, 0.05)
#' GrowthParams = data.frame(y1=y1, y2=y2, a=a, b=b)
#' RefnceAges = data.frame(t1=t1,t2=t2)
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen
#' ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
#' PropReleased = NA # proportion of fish released, vector including mean and sd (option probably now obselete)
#' midpt=Res$midpt
#' lbnd=Res$lbnd
#' ubnd=Res$ubnd
#' InitFishMort = 0.25 # specify starting parameters
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
#' InitL50 = 400
#' InitDelta = 100
#' params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
#' FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
#'                                          lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
#' PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
#'                                    SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
#'                                    RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
#'                                    xaxis_lab=NA, yaxis_lab=NA, xmax=1500, xint=500,
#'                                    ymax=0.15, yint=0.05, PlotCLs=TRUE, FittedRes, nReps=200)
#' @export
PlotLengthBasedCatchCurve_RetCatch <- function(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                             SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                             MaxAge, NatMort, TimeStep, MainLabel, xaxis_lab, yaxis_lab, xmax, xint,
                                             ymax, yint, PlotCLs, FittedRes, nReps) {

  .pardefault <- par(no.readonly = TRUE) # store current par settings

  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    res =  FittedRes
  } else {
    res=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                        lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
  }

  params = res$params
  vcov.params = res$vcov.Params
  ExpCatchAtLen = res$ModelDiag$ExpRetCatchPropInLenClass
  ObsRelCatchAtLen = ObsRetCatchFreqAtLen/sum(ObsRetCatchFreqAtLen)

  set.seed(123)
  sims = data.frame(MASS::mvrnorm(n = nReps, params, vcov.params))
  EstPropAtLen.sim = as.matrix(data.frame(matrix(nrow = nReps, ncol = length(midpt))))

  for (j in 1:nReps) {
    params = unlist(sims[j,])
    CatchCurveType=1 #1=length-based, 2=age and length based
    Res=AgeAndLengthBasedCatchCurvesCalcs(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge, CatchCurveType,
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
    legend("bottomright", legend=c("Female","Male"), y.intersp = 1.0, inset=c(0.05,0.05),
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
#' SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
#' set.seed(123)
#' MaxAge = 30
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = 0.2
#' MaxLen = 1200
#' LenInc = 20
#' MLL=500 # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = NA # L50, L95-L50 for retention
#' DiscMort = 0 # proportion of fish that die due to natural mortality
#' # single sex, von Bertalanffy
#' DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.08
#' GrowthParams = c(Linf, vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen
#' ObsDiscCatchFreqAtLen = Res$ObsDiscCatchFreqAtLen
#' PropReleased = NA # proportion of fish released, vector including mean and sd (option probably now obselete)
#' midpt=Res$midpt
#' lbnd=Res$lbnd
#' ubnd=Res$ubnd
#' InitFishMort = 0.25 # specify starting parameters
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
#' InitL50 = 400
#' InitDelta = 100
#' params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
#' FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
#'                                           lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
#' PlotLengthBasedCatchCurve_DiscCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
#'                                     SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
#'                                     xaxis_lab=NA, yaxis_lab=NA, xmax=1500, xint=500,
#'                                     ymax=0.15, yint=0.05, PlotCLs=TRUE, FittedRes, nReps=200)
#' @export
PlotLengthBasedCatchCurve_DiscCatch <- function(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                               SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                               MaxAge, NatMort, TimeStep, MainLabel, xaxis_lab, yaxis_lab, xmax, xint,
                                               ymax, yint, PlotCLs, FittedRes, nReps) {

  if (is.na(ObsDiscCatchFreqAtLen[1])) {
    cat("ObsDiscCatchFreqAtLen = NA: Can't do this plot without discard data",'\n')
  } else {

    .pardefault <- par(no.readonly = TRUE) # store current par settings

    # if model already fitted, can input results rather than refit
    if (is.list(FittedRes)) {
      res =  FittedRes
    } else {
      res=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                          lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
    }

    params = res$params
    vcov.params = res$vcov.Params
    ExpCatchAtLen = res$ModelDiag$ExpDiscCatchPropInLenClass
    ObsRelCatchAtLen = ObsDiscCatchFreqAtLen/sum(ObsDiscCatchFreqAtLen)

    set.seed(123)
    sims = data.frame(MASS::mvrnorm(n = nReps, params, vcov.params))
    EstPropAtLen.sim = as.matrix(data.frame(matrix(nrow = nReps, ncol = length(midpt))))

    for (j in 1:nReps) {
      params = unlist(sims[j,])
      CatchCurveType=1 #1=length-based, 2=age and length based
      Res=AgeAndLengthBasedCatchCurvesCalcs(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge, CatchCurveType,
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
#' SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
#' set.seed(123)
#' MaxAge = 30
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = 0.2
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0 # proportion of fish that die due to natural mortality
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.08
#' GrowthParams = c(Linf, vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # 2 sexes, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = c(700,850)
#' vbK = c(0.3,0.2)
#' CVSizeAtAge = c(0.08,0.08)
#' GrowthParams = data.frame(Linf=Linf, vbK=vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # 1 sex, Schnute
#' GrowthCurveType = 2 # 1 = von Bertalanffy, 2 = Schnute
#' t1 = 0.5 # growth - Schnute
#' t2 = 25 # growth - Schnute
#' y1 = 100 # growth - Schnute
#' y2 = 1000 # growth - Schnute
#' a = 0.02 # growth - Schnute
#' b = 3.0 # growth - Schnute
#' GrowthParams = c(y1, y2, a, b)
#' RefnceAges = c(t1,t2)
#' CVSizeAtAge = 0.05
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # 2 sexes, Schnute
#' GrowthCurveType = 2 # 1 = von Bertalanffy, 2 = Schnute
#' t1 = c(0.5,0.5) # growth - Schnute
#' t2 = c(25,25) # growth - Schnute
#' y1 = c(100,100) # growth - Schnute
#' y2 = c(1000,1000) # growth - Schnute
#' a = c(0.02,0.02) # growth - Schnute
#' b = c(3,3) # growth - Schnute
#' CVSizeAtAge = c(0.05, 0.05)
#' GrowthParams = data.frame(y1=y1, y2=y2, a=a, b=b)
#' RefnceAges = data.frame(t1=t1,t2=t2)
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen
#' ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
#' PropReleased = NA # proportion of fish released, vector including mean and sd (option probably now obselete)
#' midpt=Res$midpt
#' lbnd=Res$lbnd
#' ubnd=Res$ubnd
#' InitFishMort = 0.25 # specify starting parameters
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
#' InitL50 = 400
#' InitDelta = 100
#' params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
#' FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
#'                                     lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
#' PlotLengthBasedCatchCurve_Selectivity(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityVec,
#'                                       PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
#'                                       MaxAge, NatMort, TimeStep, MainLabel=NA, xaxis_lab=NA, yaxis_lab=NA, xmax=1200, xint=200,
#'                                       ymax=1.0, yint=0.2, PlotCLs=TRUE, FittedRes, nReps=200)
#' @export
PlotLengthBasedCatchCurve_Selectivity <- function(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityVec,
                                                  PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                                  MaxAge, NatMort, TimeStep, MainLabel, xaxis_lab, yaxis_lab, xmax, xint,
                                                  ymax, yint, PlotCLs, FittedRes, nReps) {

  .pardefault <- par(no.readonly = TRUE) # store current par settings

  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    res =  FittedRes
  } else {
    res=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                        lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
  }

  FishLen = seq(0,max(ubnd),1)
  nFishLen = length(FishLen)

  SelAtLength.sim = as.matrix(data.frame(matrix(nrow = nReps, ncol = nFishLen)))
  RetAtLength.sim = SelAtLength.sim
  SelLandAtLength.sim = SelAtLength.sim

  if (SelectivityType == 1) { # input selectivity as vector
    EstGearSel.sim.med = SelectivityVec
    EstGearSel.sim.low = SelectivityVec
    EstGearSel.sim.up = SelectivityVec
  } # SelectivityType == 1

  if (SelectivityType == 2) {
    params = res$params
    vcov.params = res$vcov.Params
    set.seed(123)
    sims = data.frame(MASS::mvrnorm(n = nReps, params, vcov.params))

    for (j in 1:nReps) {
      ParamVals = exp(unlist(sims[j,]))

      SelAtLength.sim[j,] = 1 / (1 + exp(-log(19) * (FishLen - ParamVals[2]) / ParamVals[3]))
      if(length(ParamVals==5)) {
        RetAtLength.sim[j,] = 1 / (1 + exp(-log(19) * (FishLen - ParamVals[4]) / ParamVals[5]))
        SelLandAtLength.sim[j,] = SelAtLength.sim[j,] * RetAtLength.sim[j,]
      }
      cat("j",j,"ParamVals",ParamVals,  '\n')
    }

    EstGearSel.sim.med = as.vector(apply(SelAtLength.sim, 2, median, na.rm=T))
    EstGearSel.sim.low = as.vector(apply(SelAtLength.sim, 2, quantile, probs = 0.025, na.rm=T))
    EstGearSel.sim.up = as.vector(apply(SelAtLength.sim, 2, quantile, probs = 0.975, na.rm=T))

    if(length(ParamVals)>=5) {

      EstRet.sim.med = as.vector(apply(RetAtLength.sim, 2, median, na.rm=T))
      EstRet.sim.low = as.vector(apply(RetAtLength.sim, 2, quantile, probs = 0.025, na.rm=T))
      EstRet.sim.up = as.vector(apply(RetAtLength.sim, 2, quantile, probs = 0.975, na.rm=T))

      EstLandSel.sim.med = as.vector(apply(SelLandAtLength.sim, 2, median, na.rm=T))
      EstLandSel.sim.low = as.vector(apply(SelLandAtLength.sim, 2, quantile, probs = 0.025, na.rm=T))
      EstLandSel.sim.up = as.vector(apply(SelLandAtLength.sim, 2, quantile, probs = 0.975, na.rm=T))

    } else {
      if (is.na(MLL)) { # specifying retention as 1, i.e. all fish caught are retained
        EstRet.sim.med = rep(1,length(FishLen))
        EstRet.sim.low = rep(1,length(FishLen))
        EstRet.sim.up = rep(1,length(FishLen))

      } else { # knife edge retention at MLL
        RetAtLength = rep(1E-20,length(FishLen))
        RetAtLength[which(FishLen>=MLL)]=1
        EstRet.sim.med = RetAtLength
        EstRet.sim.low = RetAtLength
        EstRet.sim.up = RetAtLength
      }
      EstLandSel.sim.med = EstGearSel.sim.med * EstRet.sim.med
      EstLandSel.sim.low = EstGearSel.sim.low * EstRet.sim.low
      EstLandSel.sim.up = EstGearSel.sim.up * EstRet.sim.up
    }
  }

  if (is.na(xaxis_lab)) xaxis_lab = "Length (mm)"
  if (is.na(yaxis_lab)) yaxis_lab = "Probability"
  xlims = Get_xaxis_scale(midpt)
  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint
  if (is.na(ymax)) ymax = 1.2
  if (is.na(yint)) yint = 0.2

  plot(FishLen, EstGearSel.sim.med, "l", main=MainLabel, cex.main=1.2, pch=1, cex=0.6,
       xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax),
       ylim=c(0,ymax), col="red")
  lines(FishLen, EstRet.sim.med, col="blue")
  lines(FishLen, EstLandSel.sim.med, col="black")

  if (PlotCLs == TRUE) {
    sm1 = spline(FishLen, EstGearSel.sim.low, n=100, method="natural")
    sm2 = spline(FishLen, EstGearSel.sim.up, n=100, method="natural")
    x = c(sm1$x, rev(sm2$x)) # using shading for 95% CLs
    y = c(sm1$y, rev(sm2$y))
    polygon(x,y, col="pink",border=NA)
    lines(FishLen, EstGearSel.sim.med, col="red")

    if(length(ParamVals)>=5) {

      sm1 = spline(FishLen, EstRet.sim.low, n=100, method="natural")
      sm2 = spline(FishLen, EstRet.sim.up, n=100, method="natural")
      x = c(sm1$x, rev(sm2$x)) # using shading for 95% CLs
      y = c(sm1$y, rev(sm2$y))
      polygon(x,y, col="lightblue",border=NA)
      lines(FishLen, EstRet.sim.med, col="blue")
    }
    sm1 = spline(FishLen, EstLandSel.sim.low, n=100, method="natural")
    sm2 = spline(FishLen, EstLandSel.sim.up, n=100, method="natural")
    x = c(sm1$x, rev(sm2$x)) # using shading for 95% CLs
    y = c(sm1$y, rev(sm2$y))
    polygon(x,y, col="lightgrey",border=NA)
    lines(FishLen, EstLandSel.sim.med, col="black")
    lines(FishLen, EstGearSel.sim.med, col="red", cex=0.6)
    lines(FishLen, EstLandSel.sim.med, col="black", cex=0.6)
  }

  axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
  axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
  axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
  axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)

  if (SelectivityType==1) {
    legend("bottomright", pch=-1, legend=c("Sel_gear", "Prob_Reten", "Sel_land"), lty="solid", col=c("red","blue","black"),
           bty='n', cex=0.8,lwd=1, y.intersp=1.1)
  }

  if (SelectivityType==2) { # logistic selectivity

    L50_Est=paste("L50_gear =",round(exp(params[2]),0),"mm")
    Delta_Est=paste("Delta_gear =",round(exp(params[3]),0),"mm")

    if(length(params)>=5) {
      L502_Est=paste("L50_land =",round(exp(params[4]),0),"mm")
      Delta2_Est=paste("Delta_land =",round(exp(params[5]),0),"mm")
    }

    if(is.na(MLL)) {
      if(length(params)==3 | length(params)==4) {
        lgnd = c("Sel_gear", "Prob_Reten", "Sel_land", L50_Est, Delta_Est)
        cls = c("red","blue","black","black","black")
        lwds = c(1,1,1,-1,-1)
      }
      if(length(params)>=5) {
        lgnd = c("Sel_gear", "Prob_Reten", "Sel_land", L50_Est, Delta_Est, L502_Est, Delta2_Est)
        cls = c("red","blue","black","black","black","black","black")
        lwds = c(1,1,1,-1,-1,-1,-1)
        }
      legend("bottomright", pch=-1, legend=lgnd, lty="solid", col=cls, bty='n', cex=0.8, lwd=lwds, y.intersp=1.1)

    }
  } # selectivity type = 2

  if (is.numeric(MLL)) {
    abline(v=MLL, lty="dotted")
    lines(FishLen, EstRet.sim.med, col="blue")
    lines(FishLen, EstLandSel.sim.med, col="black")
    lgnd = c("Sel_gear", "Prob_Reten", "Sel_land", paste("MLL =",MLL), L50_Est, Delta_Est)
    cls = c("red","blue","black","black","black","black","black")
    ltps = c("solid","solid","solid","dotted","solid","solid","solid")
    lwds = c(1,1,1,1,-1,-1,-1)

    legend("bottomright", pch=-1, legend=lgnd, lty=ltps, col=cls, bty='n', cex=0.8, lwd=lwds, y.intersp=1.1)
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
#' SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
#' set.seed(123)
#' MaxAge = 30
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = 0.2
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0 # proportion of fish that die due to natural mortality
#' # single sex, von Bertalanffy
#' DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.08
#' GrowthParams = c(Linf, vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # 2 sexes, von Bertalanffy
#' DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = c(700,850)
#' vbK = c(0.3,0.2)
#' CVSizeAtAge = c(0.08,0.08)
#' GrowthParams = data.frame(Linf=Linf, vbK=vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # 1 sex, Schnute
#' DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
#' GrowthCurveType = 2 # 1 = von Bertalanffy, 2 = Schnute
#' t1 = 0.5 # growth - Schnute
#' t2 = 25 # growth - Schnute
#' y1 = 100 # growth - Schnute
#' y2 = 1000 # growth - Schnute
#' a = 0.02 # growth - Schnute
#' b = 3.0 # growth - Schnute
#' GrowthParams = c(y1, y2, a, b)
#' RefnceAges = c(t1,t2)
#' CVSizeAtAge = 0.05
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # 2 sexes, Schnute
#' DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
#' GrowthCurveType = 2 # 1 = von Bertalanffy, 2 = Schnute
#' t1 = c(0.5,0.5) # growth - Schnute
#' t2 = c(25,25) # growth - Schnute
#' y1 = c(100,100) # growth - Schnute
#' y2 = c(1000,1000) # growth - Schnute
#' a = c(0.02,0.02) # growth - Schnute
#' b = c(3,3) # growth - Schnute
#' CVSizeAtAge = c(0.05, 0.05)
#' GrowthParams = data.frame(y1=y1, y2=y2, a=a, b=b)
#' RefnceAges = data.frame(t1=t1,t2=t2)
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen
#' ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
#' PropReleased = NA # proportion of fish released, vector including mean and sd (option probably now obselete)
#' midpt=Res$midpt
#' lbnd=Res$lbnd
#' ubnd=Res$ubnd
#' InitFishMort = 0.25 # specify starting parameters
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
#' InitL50 = 400
#' InitDelta = 100
#' params = c(InitFishMort_logit, log(InitL50), log(InitDelta))
#' FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
#'                                     lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
#' PlotLengthBasedCatchCurve_Mortality(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityVec,
#'                                     PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
#'                                     MaxAge, NatMort, TimeStep, xmax=NA, xint=NA, ymax=NA, yint=NA, FittedRes)
#' @export
PlotLengthBasedCatchCurve_Mortality <- function(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityVec,
                                                PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                                MaxAge, NatMort, TimeStep, xmax, xint, ymax, yint, FittedRes) {

  .pardefault <- par(no.readonly = TRUE) # store current par settings

  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    res =  FittedRes
  } else {
    res=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                        lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)

  }


  CatchCurveType=1 #1=length-based, 2=age and length based
  params = res$params
  res = AgeAndLengthBasedCatchCurvesCalcs(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge, CatchCurveType,
                                          MLL, SelectivityType, lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
  GrowthModelType = res$GrowthModelType

  xlims = Get_xaxis_scale(midpt)
  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint
  if (is.na(ymax)) ymax = 0.2
  if (is.na(yint)) yint = 0.05

  par(mfrow=c(2,2), mar=c(5,4,2,2), oma=c(2,2,2,2))
  yaxis_lab = expression(paste(italic("F") ~ (year^{-1})))
  xaxis_lab = "Length (mm)"

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
#' @param RefnceAges Reference ages for Schnute growth curve (set to NA for von Bertalanffy growth curve)
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
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
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified is knife-edged at MLL
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.05
#' RefnceAges = NA
#' GrowthParams = c(Linf, vbK, CVSizeAtAge)
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # 2 sexes, von Bertalanffy
#' # GrowthCurveType = 1 # 1 = von Bertalanffy
#' # Linf = c(700,850)
#' # vbK = c(0.25,0.2)
#' # CVSizeAtAge = c(0.05,0.05)
#' # RefnceAges = NA
#' # GrowthParams = data.frame(Linf=Linf, vbK=vbK, CVSizeAtAge=CVSizeAtAge)
#' # Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#' #                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' lbnd=Res$lbnd
#' midpt=Res$midpt
#' ubnd=Res$ubnd
#' # get data - 1 sex (or combined sexes)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen # 1 sex
#' ObsRetCatchFreqAtLengthAndAge = as.matrix(Res$ObsRetCatchFreqAtLengthAndDecAge) # 1 sex
#' # # get data - 2 sexes
#' # ObsRetCatchFreqAtLen <- data.frame(matrix(nrow = 2, ncol = length(midpt))) # 2 sex
#' # colnames(ObsRetCatchFreqAtLen) <- midpt
#' # ObsRetCatchFreqAtLen[1,] = Res$ObsRetCatchFreqAtLen_Fem
#' # ObsRetCatchFreqAtLen[2,] = Res$ObsRetCatchFreqAtLen_Mal
#' # ObsRetCatchFreqAtLengthAndAge = array(c(unlist(Res$ObsRetCatchFreqAtLengthAndDecAge_Fem), unlist(Res$ObsRetCatchFreqAtLengthAndDecAge_Mal)),
#' #                                       c(nTimeSteps, length(midpt), 2), dimnames=list(rownames(Res$ObsRetCatchFreqAtLengthAndDecAge_Fem),
#' #                                                                                      colnames(Res$ObsRetCatchFreqAtLengthAndDecAge_Fem)))
#' # get params - 1 sex
#' InitFishMort = 0.3 # specify starting parameters
#' InitL50 = 320
#' InitDelta = 50 # L95-L50
#' InitLinf = 800
#' InitvbK = 0.2
#' InitCVSizeAtAge = 0.05
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' params = c(InitFishMort_logit, log(c(InitL50, InitDelta, InitLinf, InitvbK, InitCVSizeAtAge)))
#' # # get params - 2 sexes
#' # InitFishMort = 0.3 # specify starting parameters
#' # InitL50 = 320
#' # InitDelta = 50 # L95-L50
#' # InitLinf = c(800,800)
#' # InitvbK = c(0.25,0.25)
#' # InitCVSizeAtAge = 0.05
#' # InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' # params = c(InitFishMort_logit, log(c(InitL50, InitDelta, InitLinf, InitvbK, InitCVSizeAtAge)))
#' FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
#' # Example with specified selectivity vector
#' # Simulate data
#' SelectivityType=1 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = 1 / (1 + exp(-log(19)*(midpt-400)/(500-400)))
#' SelParams = c(NA, NA) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.05
#' GrowthParams = c(Linf, vbK, CVSizeAtAge)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' lbnd=Res$lbnd
#' midpt=Res$midpt
#' ubnd=Res$ubnd
#' # get data - 1 sex (or combined sexes)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen # 1 sex
#' ObsRetCatchFreqAtLengthAndAge = as.matrix(Res$ObsRetCatchFreqAtLengthAndDecAge) # 1 sex
#' # get params - 1 sex
#' InitFishMort = 0.3 # specify starting parameters
#' InitLinf = 800
#' InitvbK = 0.2
#' InitCVSizeAtAge = 0.05
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' params = c(InitFishMort_logit, log(c(InitLinf, InitvbK, InitCVSizeAtAge)))
#' FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
#'PlotAgeLengthCatchCurve_MargLength(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                   lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep, MainLabel=NA,
#'                                   xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA,
#'                                   ymax=0.10, yint=0.02, PlotCLs=TRUE, FittedRes, nReps=200)
#' @export
PlotAgeLengthCatchCurve_MargLength <- function(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                               lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep, MainLabel,
                                               xaxis_lab, yaxis_lab, xmax, xint,
                                               ymax, yint, PlotCLs, FittedRes, nReps) {

  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    res =  FittedRes
  } else {

    res=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                              lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
  }

  params = res$ModelDiag$params
  vcov.params = res$ModelDiag$vcov.Params
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

  for (j in 1:nReps) {
    params = unlist(sims[j,])
    CatchCurveType=2 #1=length-based, 2=age and length based
    GrowthParams = NA
    CVSizeAtAge = NA
    Res=AgeAndLengthBasedCatchCurvesCalcs(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge, CatchCurveType,
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
      sm1 = spline(Res$midpt, EstProp.sim_low, n=100, method="natural")
      sm2 = spline(Res$midpt, EstProp.sim_up, n=100, method="natural")
      x = c(sm1$x, rev(sm2$x)) # using shading for 95% CLs
      y = c(sm1$y, rev(sm2$y))
      polygon(x,y,col="pink",border=NA)
    }
    points(midpt, ObsRelCatchAtLen, col="black", pch=16, cex=0.8)
    points(midpt, EstProp.sim, col="red", pch=1, cex=0.8)
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
  }

  # separate sexes
  if (is.data.frame(ObsRetCatchFreqAtLen)) {
    if (is.na(yaxis_lab)) yaxis_lab1 = "Proportion - Females"
    plot(midpt, ObsRelCatchAtLen[1,], "p", main=MainLabel, cex.main=1.2, pch=16, cex=0.8,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab1,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax))
    points(midpt, ObsRelCatchAtLen[1,], col="black", pch=16, cex=0.8)
    if (PlotCLs == TRUE) {
      sm1 = spline(Res$midpt, EstPropF.sim_low, n=100, method="natural")
      sm2 = spline(Res$midpt, EstPropF.sim_up, n=100, method="natural")
      x = c(sm1$x, rev(sm2$x)) # using shading for 95% CLs
      y = c(sm1$y, rev(sm2$y))
      polygon(x,y,col="pink",border=NA)
    }
    points(midpt, ObsRelCatchAtLen[1,], col="black", pch=16, cex=0.8)
    points(midpt, EstPropF.sim, col="red", pch=1, cex=0.8)
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
    if (is.na(yaxis_lab)) yaxis_lab2 = "Proportion - Males"
    plot(midpt, ObsRelCatchAtLen[2,], "p", main=MainLabel, cex.main=1.2, pch=16, cex=0.8,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab2,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax))
    points(midpt, ObsRelCatchAtLen[2,], col="black", pch=16, cex=0.8)
    if (PlotCLs == TRUE) {
      sm1 = spline(Res$midpt, EstPropM.sim_low, n=100, method="natural")
      sm2 = spline(Res$midpt, EstPropM.sim_up, n=100, method="natural")
      x = c(sm1$x, rev(sm2$x)) # using shading for 95% CLs
      y = c(sm1$y, rev(sm2$y))
      polygon(x,y,col="light blue",border=NA)
    }
    points(midpt, ObsRelCatchAtLen[2,], col="black", pch=16, cex=0.8)
    points(midpt, EstPropM.sim, col="blue", pch=1, cex=0.8)
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
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
#' @param RefnceAges Reference ages for Schnute growth curve (set to NA for von Bertalanffy growth curve)
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
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
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified is knife-edged at MLL
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.05
#' RefnceAges = NA
#' GrowthParams = c(Linf, vbK, CVSizeAtAge)
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # 2 sexes, von Bertalanffy
#' # GrowthCurveType = 1 # 1 = von Bertalanffy
#' # Linf = c(700,850)
#' # vbK = c(0.25,0.2)
#' # CVSizeAtAge = c(0.05,0.05)
#' # RefnceAges = NA
#' # GrowthParams = data.frame(Linf=Linf, vbK=vbK, CVSizeAtAge=CVSizeAtAge)
#' # Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#' #                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' lbnd=Res$lbnd
#' midpt=Res$midpt
#' ubnd=Res$ubnd
#' # get data - 1 sex (or combined sexes)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen # 1 sex
#' ObsRetCatchFreqAtLengthAndAge = as.matrix(Res$ObsRetCatchFreqAtLengthAndDecAge) # 1 sex
#' # # get data - 2 sexes
#' # ObsRetCatchFreqAtLen <- data.frame(matrix(nrow = 2, ncol = length(midpt))) # 2 sex
#' # colnames(ObsRetCatchFreqAtLen) <- midpt
#' # ObsRetCatchFreqAtLen[1,] = Res$ObsRetCatchFreqAtLen_Fem
#' # ObsRetCatchFreqAtLen[2,] = Res$ObsRetCatchFreqAtLen_Mal
#' # ObsRetCatchFreqAtLengthAndAge = array(c(unlist(Res$ObsRetCatchFreqAtLengthAndDecAge_Fem), unlist(Res$ObsRetCatchFreqAtLengthAndDecAge_Mal)),
#' #                                       c(nTimeSteps, length(midpt), 2), dimnames=list(rownames(Res$ObsRetCatchFreqAtLengthAndDecAge_Fem),
#' #                                                                                      colnames(Res$ObsRetCatchFreqAtLengthAndDecAge_Fem)))
#' # get params - 1 sex
#' InitFishMort = 0.3 # specify starting parameters
#' InitL50 = 320
#' InitDelta = 50 # L95-L50
#' InitLinf = 800
#' InitvbK = 0.2
#' InitCVSizeAtAge = 0.05
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' params = c(InitFishMort_logit, log(c(InitL50, InitDelta, InitLinf, InitvbK, InitCVSizeAtAge)))
#' # # get params - 2 sexes
#' # InitFishMort = 0.3 # specify starting parameters
#' # InitL50 = 320
#' # InitDelta = 50 # L95-L50
#' # InitLinf = c(800,800)
#' # InitvbK = c(0.25,0.25)
#' # InitCVSizeAtAge = 0.05
#' # InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' # params = c(InitFishMort_logit, log(c(InitL50, InitDelta, InitLinf, InitvbK, InitCVSizeAtAge)))
#' FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
#' # Example with specified selectivity vector
#' # Simulate data
#' SelectivityType=1 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = 1 / (1 + exp(-log(19)*(midpt-400)/(500-400)))
#' SelParams = c(NA, NA) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.05
#' GrowthParams = c(Linf, vbK, CVSizeAtAge)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' lbnd=Res$lbnd
#' midpt=Res$midpt
#' ubnd=Res$ubnd
#' # get data - 1 sex (or combined sexes)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen # 1 sex
#' ObsRetCatchFreqAtLengthAndAge = as.matrix(Res$ObsRetCatchFreqAtLengthAndDecAge) # 1 sex
#' # get params - 1 sex
#' InitFishMort = 0.3 # specify starting parameters
#' InitLinf = 800
#' InitvbK = 0.2
#' InitCVSizeAtAge = 0.05
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' params = c(InitFishMort_logit, log(c(InitLinf, InitvbK, InitCVSizeAtAge)))
#' FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
#' PlotAgeLengthCatchCurve_Growth(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep, MainLabel=NA,
#'                                xaxis_lab=NA, yaxis_lab=NA, xmax=40, xint=10,
#'                                ymax=1000, yint=200, PlotCLs=TRUE, FittedRes, nReps=200)
#' @export
PlotAgeLengthCatchCurve_Growth <- function(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                           lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep, MainLabel,
                                           xaxis_lab, yaxis_lab, xmax, xint,
                                           ymax, yint, PlotCLs, FittedRes, nReps) {


  # generate data for individual fish (age classes and mid points of length classes)
  if (is.vector(ObsRetCatchFreqAtLen)) {
    SampleSize = sum(ObsRetCatchFreqAtLen)
    ObsAge = rep(NA, SampleSize)
    ObsLenClRetCatchMidPt = rep(NA, SampleSize)
  }
  if (is.data.frame(ObsRetCatchFreqAtLen)) {
    SampleSize_F = sum(ObsRetCatchFreqAtLen[1,])
    SampleSize_M = sum(ObsRetCatchFreqAtLen[2,])
    ObsAge_F = rep(NA, SampleSize_F)
    ObsAge_M = rep(NA, SampleSize_M)
    ObsLenClRetCatchMidPt_F = rep(NA, SampleSize_F)
    ObsLenClRetCatchMidPt_M = rep(NA, SampleSize_M)
  }

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
          ObsLenClRetCatchMidPt[strt:fnsh]=midpt[j]
          strt=strt+x
        }
      }
      if (is.data.frame(ObsRetCatchFreqAtLen)) { # 2 sexes
        # females
        x=ObsRetCatchFreqAtLengthAndAge[i,j,1] # number of females in current length and age class
        if(x>0) {
          fnshF=strtF+x-1
          ObsAge_F[strtF:fnshF]=i*TimeStep
          ObsLenClRetCatchMidPt_F[strtF:fnshF]=midpt[j]
          strtF=strtF+x
        }
        # males
        x=ObsRetCatchFreqAtLengthAndAge[i,j,2] # number of males in current length and age class
        if(x>0) {
          fnshM=strtM+x-1
          ObsAge_M[strtM:fnshM]=i*TimeStep
          ObsLenClRetCatchMidPt_M[strtM:fnshM]=midpt[j]
          strtM=strtM+x
        }
      }
    }
  }

  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    res =  FittedRes
  } else {
    res=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                              lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
  }
  params = res$ModelDiag$params
  vcov.params = res$ModelDiag$vcov.Params
  set.seed(123)
  sims = data.frame(MASS::mvrnorm(n = nReps, params, vcov.params))
  EstLenAtAge.sim = data.frame(matrix(nrow = nReps, ncol = nTimeSteps))
  EstLenAtAgeF.sim = EstLenAtAge.sim
  EstLenAtAgeM.sim = EstLenAtAge.sim

  if (GrowthCurveType==1) { # von Bertalanffy
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
  }


  if (GrowthCurveType==2) { # Schnute
    t1=RefnceAges[1]; t2=RefnceAges[2]
    for (j in 1:nReps) {
      ParamVals = unlist(sims[j,])
      if (is.vector(ObsRetCatchFreqAtLen)) { # combined sex
        if (SelectivityType == 1) { # input vector
          y1=0; y2=exp(ParamVals[2]); a=ParamVals[3]; b=ParamVals[4]
          for(i in 1:length(DecAges)) {
            Age=DecAges[i]
            EstLenAtAge.sim[j,i] = SchnuteGrowthfunction(Age, t1, t2, y1, y2, a, b)
          }
        }
        if (SelectivityType == 2) { # estimated
          y1=0; y2=exp(ParamVals[4]); a=ParamVals[5]; b=ParamVals[6]
          for(i in 1:length(DecAges)) {
            Age=DecAges[i]
            EstLenAtAge.sim[j,i] = SchnuteGrowthfunction(Age, t1, t2, y1, y2, a, b)
          }
        }
      }
      if (is.data.frame(ObsRetCatchFreqAtLen)) { # separate sexes
        if (SelectivityType == 1) { # input vector
          for(i in 1:length(DecAges)) {
            Age=DecAges[i]
            y1=0; y2=exp(ParamVals[2]); a=ParamVals[4]; b=ParamVals[6]
            EstLenAtAgeF.sim[j,] = SchnuteGrowthfunction(Age, t1, t2, y1, y2, a, b)
            y1=0; y2=exp(ParamVals[3]); a=ParamVals[5]; b=ParamVals[7]
            EstLenAtAgeM.sim[j,] = SchnuteGrowthfunction(Age, t1, t2, y1, y2, a, b)
          }
        }
        if (SelectivityType == 2) { # estimated
          for(i in 1:length(DecAges)) {
            Age=DecAges[i]
            y1=0; y2=exp(ParamVals[4]); a=ParamVals[6]; b=ParamVals[8]
            EstLenAtAgeF.sim[j,] = SchnuteGrowthfunction(Age, t1, t2, y1, y2, a, b)
            y1=0; y2=exp(ParamVals[5]); a=ParamVals[7]; b=ParamVals[9]
            EstLenAtAgeM.sim[j,] = SchnuteGrowthfunction(Age, t1, t2, y1, y2, a, b)
          }
        }
      }
      cat("j",j,'\n')
    }
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
    plot(ObsAge, ObsLenClRetCatchMidPt, "p", main=MainLabel, cex.main=1.2, pch=16, cex=0.6,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax))
    points(ObsAge, ObsLenClRetCatchMidPt, col="black", cex=0.6)
    points(DecAges,EstProp.sim, col="red", cex=0.6)
    if (PlotCLs == TRUE) {
      x = c(DecAges,rev(DecAges)) # using shading for 95% CLs
      y = c(EstProp.sim_low, rev(EstProp.sim_up))
      polygon(x,y,col="pink",border=NA)
    }
    points(ObsAge, ObsLenClRetCatchMidPt, col="black", cex=0.6)
    points(DecAges, EstProp.sim, col="red", pch=1, cex=0.6)
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
  }
  if (is.data.frame(ObsRetCatchFreqAtLen)) { # females
    # females
    ylims = Get_yaxis_scale(midpt)
    if (is.na(ymax)) ymax = ylims$ymax
    if (is.na(yint)) yint = ylims$yint
    if (is.na(yaxis_lab)) yaxis_lab1 = "Length females"
    plot(ObsAge_F, ObsLenClRetCatchMidPt_F, "p", main=MainLabel, cex.main=1.2, pch=16, cex=0.6,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2), ylab=list(yaxis_lab1,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax))
    points(ObsAge_F, ObsLenClRetCatchMidPt_F, col="black", cex=0.6)
    points(DecAges,EstPropF.sim, col="red", cex=0.6)
    if (PlotCLs == TRUE) {
      x = c(DecAges,rev(DecAges)) # using shading for 95% CLs
      y = c(EstPropF.sim_low, rev(EstPropF.sim_up))
      polygon(x,y,col="pink",border=NA)
    }
    points(ObsAge_F, ObsLenClRetCatchMidPt_F, col="black", cex=0.6)
    points(DecAges, EstPropF.sim, col="red", pch=1, cex=0.6)
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
    # males
    if (is.na(yaxis_lab)) yaxis_lab2 = "Length males"
    ylims = Get_yaxis_scale(midpt)
    if (is.na(ymax)) ymax = ylims$ymax
    if (is.na(yint)) yint = ylims$yint
    plot(ObsAge_M, ObsLenClRetCatchMidPt_M, "p", main=MainLabel, cex.main=1.0, pch=16, cex=0.6,
         xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2), ylab=list(yaxis_lab2,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax))
    points(ObsAge_M, ObsLenClRetCatchMidPt_M, col="black", cex=0.6)
    points(DecAges,EstPropM.sim, col="blue", cex=0.6)
    if (PlotCLs == TRUE) {
      x = c(DecAges,rev(DecAges)) # using shading for 95% CLs
      y = c(EstPropM.sim_low, rev(EstPropM.sim_up))
      polygon(x,y,col="light blue",border=NA)
    }
    points(ObsAge_M, ObsLenClRetCatchMidPt_M, col="black", cex=0.6)
    points(DecAges, EstPropM.sim, col="blue", pch=1, cex=0.6)
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
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
#' @param RefnceAges Reference ages for Schnute growth curve (set to NA for von Bertalanffy growth curve)
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
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
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified is knife-edged at MLL
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.05
#' RefnceAges = NA
#' GrowthParams = c(Linf, vbK, CVSizeAtAge)
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # 2 sexes, von Bertalanffy
#' # GrowthCurveType = 1 # 1 = von Bertalanffy
#' # Linf = c(700,850)
#' # vbK = c(0.25,0.2)
#' # CVSizeAtAge = c(0.05,0.05)
#' # RefnceAges = NA
#' # GrowthParams = data.frame(Linf=Linf, vbK=vbK, CVSizeAtAge=CVSizeAtAge)
#' # Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#' #                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' lbnd=Res$lbnd
#' midpt=Res$midpt
#' ubnd=Res$ubnd
#' # get data - 1 sex (or combined sexes)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen # 1 sex
#' ObsRetCatchFreqAtLengthAndAge = as.matrix(Res$ObsRetCatchFreqAtLengthAndDecAge) # 1 sex
#' # # get data - 2 sexes
#' # ObsRetCatchFreqAtLen <- data.frame(matrix(nrow = 2, ncol = length(midpt))) # 2 sex
#' # colnames(ObsRetCatchFreqAtLen) <- midpt
#' # ObsRetCatchFreqAtLen[1,] = Res$ObsRetCatchFreqAtLen_Fem
#' # ObsRetCatchFreqAtLen[2,] = Res$ObsRetCatchFreqAtLen_Mal
#' # ObsRetCatchFreqAtLengthAndAge = array(c(unlist(Res$ObsRetCatchFreqAtLengthAndDecAge_Fem), unlist(Res$ObsRetCatchFreqAtLengthAndDecAge_Mal)),
#' #                                       c(nTimeSteps, length(midpt), 2), dimnames=list(rownames(Res$ObsRetCatchFreqAtLengthAndDecAge_Fem),
#' #                                                                                      colnames(Res$ObsRetCatchFreqAtLengthAndDecAge_Fem)))
#' # get params - 1 sex
#' InitFishMort = 0.3 # specify starting parameters
#' InitL50 = 320
#' InitDelta = 50 # L95-L50
#' InitLinf = 800
#' InitvbK = 0.2
#' InitCVSizeAtAge = 0.05
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' params = c(InitFishMort_logit, log(c(InitL50, InitDelta, InitLinf, InitvbK, InitCVSizeAtAge)))
#' # # get params - 2 sexes
#' # InitFishMort = 0.3 # specify starting parameters
#' # InitL50 = 320
#' # InitDelta = 50 # L95-L50
#' # InitLinf = c(800,800)
#' # InitvbK = c(0.25,0.25)
#' # InitCVSizeAtAge = 0.05
#' # InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' # params = c(InitFishMort_logit, log(c(InitL50, InitDelta, InitLinf, InitvbK, InitCVSizeAtAge)))
#' FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
#' # Example with specified selectivity vector
#' # Simulate data
#' SelectivityType=1 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = 1 / (1 + exp(-log(19)*(midpt-400)/(500-400)))
#' SelParams = c(NA, NA) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.05
#' GrowthParams = c(Linf, vbK, CVSizeAtAge)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' lbnd=Res$lbnd
#' midpt=Res$midpt
#' ubnd=Res$ubnd
#' # get data - 1 sex (or combined sexes)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen # 1 sex
#' ObsRetCatchFreqAtLengthAndAge = as.matrix(Res$ObsRetCatchFreqAtLengthAndDecAge) # 1 sex
#' # get params - 1 sex
#' InitFishMort = 0.3 # specify starting parameters
#' InitLinf = 800
#' InitvbK = 0.2
#' InitCVSizeAtAge = 0.05
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' params = c(InitFishMort_logit, log(c(InitLinf, InitvbK, InitCVSizeAtAge)))
#' FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
#' PlotAgeLengthCatchCurve_Selectivity(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                     lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep, MainLabel=NA,
#'                                     xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA,
#'                                     ymax=NA, yint=NA, PlotCLs=TRUE, FittedRes, nReps=200)
#' @export
PlotAgeLengthCatchCurve_Selectivity <- function(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                                lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep, MainLabel,
                                                xaxis_lab, yaxis_lab, xmax, xint,
                                                ymax, yint, PlotCLs, FittedRes, nReps) {

  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    res =  FittedRes
  } else {
    res=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                              lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
  }

  if (SelectivityType == 1) {
    SelAtLength = res$ModelDiag$SelAtLength
  }
  if (SelectivityType == 2) {
    params = res$ModelDiag$params
    vcov.params = res$ModelDiag$vcov.Params
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
    sm1 = spline(midpt, EstProp.sim_low, n=100, method="natural")
    sm2 = spline(midpt, EstProp.sim_up, n=100, method="natural")
    sm1$y[which(sm1$y<0)]=0; sm1$y[which(sm1$y>1)]=1
    sm2$y[which(sm2$y<0)]=0; sm2$y[which(sm2$y>1)]=1
    x = c(sm1$x, rev(sm2$x)) # using shading for 95% CLs
    y = c(sm1$y, rev(sm2$y))
    polygon(x,y, col="pink",border=NA)

  }
  sm1 = spline(midpt, EstProp.sim, n=100, method="natural")
  sm1$y[which(sm1$y<0)]=0; sm1$y[which(sm1$y>1)]=1
  lines(sm1$x,sm1$y, col="red", cex=0.6)
  points(midpt, EstProp.sim, col="red", pch=1, cex=0.6)
  AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
  if (SelectivityType==2) { # logistic selectivity
    L50est=paste("L50 =",round(exp(params[2]),0),"mm")
    L95est=paste("L95-L50 =",round(exp(params[3]),0),"mm")
    legend("topleft", pch=-1, legend=c(L50est, L95est), lty="solid",col="black",
           bty='n', cex=0.6,lwd=-1, y.intersp=1.0)
  }
}

#' Output observed and expected quantities associated with conditional age-length plots
#'
#' This function outputs several observed and expected quantities associated with conditional age-length plots,
#' including observed vs expected age distribution at length plots, and pearson residual plots
#'
#' @keywords internal
#'
#' @param params vector of model parameters (in log space) to be estimated (if FittedRes=NA)
#' @param RefnceAges Reference ages for Schnute growth curve (set to NA for von Bertalanffy growth curve)
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
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
#' @param FittedRes saved results from GetAgeAndLengthBasedCatchCurveResults function (model will be refitted if set to NA)
#'
#' @return ObsRetCatchFreqAtLengthAndIntAge, ObsCatchPropAgeAtLength, ExpRetCatchPropIntAgeGivenLength (single sex),
#' ObsRetCatchFreqAtLengthAndIntAge_Fem, ObsRetCatchFreqAtLengthAndIntAge_Mal, ObsCatchPropAgeAtLength_Fem,
#' ObsCatchPropAgeAtLength_Mal, ExpRetCatchPropIntAgeGivenLength_Fem, ExpRetCatchPropIntAgeGivenLength_Mal (separaate sex)
GetInputsForPlotting_Cond_AL <- function(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                         lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep, FittedRes) {


  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    res =  FittedRes
  } else {
    res=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                              lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
  }

  MinAge = floor(TimeStep)
  AgeClasses = MinAge:MaxAge
  nAgeCl = length(MinAge:MaxAge)
  nLenCl = length(midpt)

  # calculate expected proportions at age, for each length class
  # single sex
  if (is.vector(ObsRetCatchFreqAtLen)) {
    ExpRetCatchPropAtIntAge=res$ModelDiag$ExpRetCatchPropAtIntAge
    ExpRetCatchPropLengthGivenIntAge=res$ModelDiag$ExpRetCatchPropLengthGivenIntAge
    ExpRetCatchPropIntAgeGivenLength <- data.frame(matrix(nrow = nAgeCl, ncol = nLenCl))
    colnames(ExpRetCatchPropIntAgeGivenLength) <- midpt
    ObsCatchPropAgeAtLength = ExpRetCatchPropIntAgeGivenLength
    ExpRetCatchPropIntAgeGivenLength = CalcExpCatchPropIntAgeGivenLength_cpp(nLenCl, nAgeCl, ExpRetCatchPropLengthGivenIntAge, ExpRetCatchPropAtIntAge)

    # calculate observed proportions at age, for each length class
    if (TimeStep == 1) { # annual time step
      ObsRetCatchFreqAtLengthAndIntAge = ObsRetCatchFreqAtLengthAndAge
    } else { # shorter time step than annual
      ObsRetCatchFreqAtLengthAndIntAge = ConvertObsDataFromDecAgesToIntegerAges(TimeStep, MaxAge, nLenCl, ObsRetCatchFreqAtLengthAndAge)
    }
    for (i in 1:nLenCl) {
      if (sum(ObsRetCatchFreqAtLengthAndIntAge[,i])==0) {
        ObsCatchPropAgeAtLength[,i] = 0
      } else {
        ObsCatchPropAgeAtLength[,i] = ObsRetCatchFreqAtLengthAndIntAge[,i] / sum(ObsRetCatchFreqAtLengthAndIntAge[,i])
      }
    }

    Result = list(ObsRetCatchFreqAtLengthAndIntAge=ObsRetCatchFreqAtLengthAndIntAge,
                  ObsCatchPropAgeAtLength=ObsCatchPropAgeAtLength,
                  ExpRetCatchPropIntAgeGivenLength=ExpRetCatchPropIntAgeGivenLength)
  }
  # 2 sexes
  if (is.data.frame(ObsRetCatchFreqAtLen)) {
    ExpRetCatchPropAtIntAge_Fem=res$ModelDiag$ExpRetCatchPropAtIntAge_Fem
    ExpRetCatchPropAtIntAge_Mal=res$ModelDiag$ExpRetCatchPropAtIntAge_Mal
    ExpRetCatchPropLengthGivenIntAge_Fem=res$ModelDiag$ExpRetCatchPropLengthGivenIntAge_Fem
    ExpRetCatchPropLengthGivenIntAge_Mal=res$ModelDiag$ExpRetCatchPropLengthGivenIntAge_Mal
    ExpRetCatchPropIntAgeGivenLength_Fem <- data.frame(matrix(nrow = nAgeCl, ncol = nLenCl))
    colnames(ExpRetCatchPropIntAgeGivenLength_Fem) <- midpt
    ObsCatchPropAgeAtLength_Fem = ExpRetCatchPropIntAgeGivenLength_Fem
    ObsCatchPropAgeAtLength_Mal = ExpRetCatchPropIntAgeGivenLength_Fem
    ExpRetCatchPropIntAgeGivenLength_Mal = ExpRetCatchPropIntAgeGivenLength_Fem
    ExpRetCatchPropIntAgeGivenLength_Fem = CalcExpCatchPropIntAgeGivenLength_cpp(nLenCl, nAgeCl, ExpRetCatchPropLengthGivenIntAge_Fem, ExpRetCatchPropAtIntAge_Fem)
    ExpRetCatchPropIntAgeGivenLength_Mal = CalcExpCatchPropIntAgeGivenLength_cpp(nLenCl, nAgeCl, ExpRetCatchPropLengthGivenIntAge_Mal, ExpRetCatchPropAtIntAge_Mal)



    if (TimeStep == 1) { # annual time step
      ObsRetCatchFreqAtLengthAndIntAge_Fem = ObsRetCatchFreqAtLengthAndAge[,,1]
      ObsRetCatchFreqAtLengthAndIntAge_Mal = ObsRetCatchFreqAtLengthAndAge[,,2]
    } else { # shorter than annual time step
      ObsRetCatchFreqAtLengthAndAge_Fem = as.matrix(ObsRetCatchFreqAtLengthAndAge[,,1])
      ObsRetCatchFreqAtLengthAndAge_Mal = as.matrix(ObsRetCatchFreqAtLengthAndAge[,,2])
      ObsRetCatchFreqAtLengthAndIntAge_Fem = ConvertObsDataFromDecAgesToIntegerAges(TimeStep, MaxAge, nLenCl, ObsRetCatchFreqAtLengthAndAge_Fem)
      ObsRetCatchFreqAtLengthAndIntAge_Mal = ConvertObsDataFromDecAgesToIntegerAges(TimeStep, MaxAge, nLenCl, ObsRetCatchFreqAtLengthAndAge_Mal)
    }

    for (i in 1:nLenCl) {
      if (sum(ObsRetCatchFreqAtLengthAndIntAge_Fem[,i])==0) {
        ObsCatchPropAgeAtLength_Fem[,i] = 0
      } else {
        ObsCatchPropAgeAtLength_Fem[,i] = ObsRetCatchFreqAtLengthAndIntAge_Fem[,i] / sum(ObsRetCatchFreqAtLengthAndIntAge_Fem[,i])
      }
      if (sum(ObsRetCatchFreqAtLengthAndIntAge_Mal[,i])==0) {
        ObsCatchPropAgeAtLength_Mal[,i] = 0
      } else {
        ObsCatchPropAgeAtLength_Mal[,i] = ObsRetCatchFreqAtLengthAndIntAge_Mal[,i] / sum(ObsRetCatchFreqAtLengthAndIntAge_Mal[,i])
      }
    }

    Result = list(ObsRetCatchFreqAtLengthAndIntAge_Fem=ObsRetCatchFreqAtLengthAndIntAge_Fem,
                  ObsRetCatchFreqAtLengthAndIntAge_Mal=ObsRetCatchFreqAtLengthAndIntAge_Mal,
                  ObsCatchPropAgeAtLength_Fem=ObsCatchPropAgeAtLength_Fem,
                  ObsCatchPropAgeAtLength_Mal=ObsCatchPropAgeAtLength_Mal,
                  ExpRetCatchPropIntAgeGivenLength_Fem=ExpRetCatchPropIntAgeGivenLength_Fem,
                  ExpRetCatchPropIntAgeGivenLength_Mal=ExpRetCatchPropIntAgeGivenLength_Mal)
  }

  return(Result)

}

#' Plot observed and estimated proportions at age in each length class, from age and length catch curve model
#'
#' This function provides a plot of the estimated proportions at age in each length class,
#' as estimated from an age and length-based  catch curve model with length-based selectivity. The model is
#' fitted to a sample of fish length and age data, by minimising the overall negative log-likelihood,
#' including the NLL associated with the marginal length composition and a conditional age at length NLL,
#' given the parameters (selectivity, growth and mortality) and data, using nlminb.
#'
#' @param params vector of model parameters (in log space) to be estimated (if FittedRes=NA)
#' @param RefnceAges Reference ages for Schnute growth curve (set to NA for von Bertalanffy growth curve)
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
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
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified is knife-edged at MLL
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.05
#' RefnceAges = NA
#' GrowthParams = c(Linf, vbK, CVSizeAtAge)
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # 2 sexes, von Bertalanffy
#' # GrowthCurveType = 1 # 1 = von Bertalanffy
#' # Linf = c(700,850)
#' # vbK = c(0.25,0.2)
#' # CVSizeAtAge = c(0.05,0.05)
#' # RefnceAges = NA
#' # GrowthParams = data.frame(Linf=Linf, vbK=vbK, CVSizeAtAge=CVSizeAtAge)
#' # Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#' #                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' lbnd=Res$lbnd
#' midpt=Res$midpt
#' ubnd=Res$ubnd
#' # get data - 1 sex (or combined sexes)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen # 1 sex
#' ObsRetCatchFreqAtLengthAndAge = as.matrix(Res$ObsRetCatchFreqAtLengthAndDecAge) # 1 sex
#' # # get data - 2 sexes
#' # ObsRetCatchFreqAtLen <- data.frame(matrix(nrow = 2, ncol = length(midpt))) # 2 sex
#' # colnames(ObsRetCatchFreqAtLen) <- midpt
#' # ObsRetCatchFreqAtLen[1,] = Res$ObsRetCatchFreqAtLen_Fem
#' # ObsRetCatchFreqAtLen[2,] = Res$ObsRetCatchFreqAtLen_Mal
#' # ObsRetCatchFreqAtLengthAndAge = array(c(unlist(Res$ObsRetCatchFreqAtLengthAndDecAge_Fem), unlist(Res$ObsRetCatchFreqAtLengthAndDecAge_Mal)),
#' #                                       c(nTimeSteps, length(midpt), 2), dimnames=list(rownames(Res$ObsRetCatchFreqAtLengthAndDecAge_Fem),
#' #                                                                                      colnames(Res$ObsRetCatchFreqAtLengthAndDecAge_Fem)))
#' # get params - 1 sex
#' InitFishMort = 0.3 # specify starting parameters
#' InitL50 = 320
#' InitDelta = 50 # L95-L50
#' InitLinf = 800
#' InitvbK = 0.2
#' InitCVSizeAtAge = 0.05
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' params = c(InitFishMort_logit, log(c(InitL50, InitDelta, InitLinf, InitvbK, InitCVSizeAtAge)))
#' # # get params - 2 sexes
#' # InitFishMort = 0.3 # specify starting parameters
#' # InitL50 = 320
#' # InitDelta = 50 # L95-L50
#' # InitLinf = c(800,800)
#' # InitvbK = c(0.25,0.25)
#' # InitCVSizeAtAge = 0.05
#' # InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' # params = c(InitFishMort_logit, log(c(InitL50, InitDelta, InitLinf, InitvbK, InitCVSizeAtAge)))
#' FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
#' # Example with specified selectivity vector
#' # Simulate data
#' SelectivityType=1 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = 1 / (1 + exp(-log(19)*(midpt-400)/(500-400)))
#' SelParams = c(NA, NA) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.05
#' GrowthParams = c(Linf, vbK, CVSizeAtAge)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' lbnd=Res$lbnd
#' midpt=Res$midpt
#' ubnd=Res$ubnd
#' # get data - 1 sex (or combined sexes)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen # 1 sex
#' ObsRetCatchFreqAtLengthAndAge = as.matrix(Res$ObsRetCatchFreqAtLengthAndDecAge) # 1 sex
#' # get params - 1 sex
#' InitFishMort = 0.3 # specify starting parameters
#' InitLinf = 800
#' InitvbK = 0.2
#' InitCVSizeAtAge = 0.05
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' params = c(InitFishMort_logit, log(c(InitLinf, InitvbK, InitCVSizeAtAge)))
#' FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
#' PlotAgeLengthCatchCurve_Cond_AL(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep, MainLabel=NA,
#'                                 xaxis_lab=NA, yaxis_lab=NA, xmax=40, xint=10,
#'                                 ymax=NA, yint=NA, FittedRes)
#' @export
PlotAgeLengthCatchCurve_Cond_AL <- function(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                            lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep, MainLabel,
                                            xaxis_lab, yaxis_lab, xmax, xint, ymax, yint, FittedRes) {

  .pardefault <- par(no.readonly = TRUE) # store current par settings

  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    res =  FittedRes
  } else {
    res=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                              lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
  }

  MinAge = floor(TimeStep)
  AgeClasses = MinAge:MaxAge
  nAgeCl = length(MinAge:MaxAge)
  nLenCl = length(midpt)

  # get required inputs for plotting
  Res = GetInputsForPlotting_Cond_AL(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                     lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep, FittedRes)

  # calculate expected proportions at age, for each length class
  # single sex
  if (is.vector(ObsRetCatchFreqAtLen)) {
    ObsCatchPropAgeAtLength = Res$ObsCatchPropAgeAtLength
    ExpRetCatchPropIntAgeGivenLength = Res$ExpRetCatchPropIntAgeGivenLength
  }
  # 2 sexes
  if (is.data.frame(ObsRetCatchFreqAtLen)) {
    ObsCatchPropAgeAtLength_Fem = Res$ObsCatchPropAgeAtLength_Fem
    ObsCatchPropAgeAtLength_Mal = Res$ObsCatchPropAgeAtLength_Mal
    ExpRetCatchPropIntAgeGivenLength_Fem = Res$ExpRetCatchPropIntAgeGivenLength_Fem
    ExpRetCatchPropIntAgeGivenLength_Mal = Res$ExpRetCatchPropIntAgeGivenLength_Mal
  }

  if (is.na(xaxis_lab)) xaxis_lab = "AgeClass (y)"
  if (is.na(yaxis_lab)) yaxis_lab = "Proportion"
  xlims = Get_xaxis_scale(AgeClasses)
  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint
  if (is.na(ymax)) ymax = 1.0
  if (is.na(yint)) yint = 0.2

  # plot observed and expected age proportions for each length class
  # single sex
  if (is.vector(ObsRetCatchFreqAtLen)) {
    par(mfcol=c(3,3), mar=c(3.5,3.5,1,1), oma=c(1,1,1,0), tck=-0.03)
    k=0
    for (i in 1:nLenCl) {
      if (sum(ObsCatchPropAgeAtLength[,i]>0)) {
        k=k+1
        if (k==10) k=1

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
          MainLab=paste(lbnd[i],"-",ubnd[i],"mm,","n =",ObsRetCatchFreqAtLen[i])
        } else {
          MainLab=MainLabel
        }
        mtext(MainLab,las=1,side=3,line=0,adj=0.5, cex=0.6)
        legend("topright", legend=c("Obs","Exp"), lty="solid", lwd=c(-1,1),
               pch=c(16,-1), col=c("black","red"), bty='n', cex=0.8)
      }
    }
  }

  if (is.data.frame(ObsRetCatchFreqAtLen)) {

    # females
    par(mfcol=c(3,3), mar=c(3.5,3.5,1,1), oma=c(1,1,1,0), tck=-0.03)
    k=0
    for (i in 1:nLenCl) {
      if (sum(ObsCatchPropAgeAtLength_Fem[,i]>0)) {
        k=k+1
        if (k==10) k=1
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
          MainLab=paste(lbnd[i],"-",ubnd[i],"mm,","n =",ObsRetCatchFreqAtLen[1,i])
        } else {
          MainLab=MainLabel
        }
        mtext(MainLab,las=1,side=3,line=0,adj=0.5, cex=0.6)
        legend("topright", legend=c("Fem_Obs","Fem_Exp"), lty="solid", lwd=c(-1,1),
               pch=c(16,-1), col=c("black","red"), bty='n', cex=0.8)
      }
    }
    # males
    par(mfcol=c(3,3), mar=c(3.5,3.5,1,1), oma=c(1,1,1,0), tck=-0.03)
    k=0
    for (i in 1:nLenCl) {
      if (sum(ObsCatchPropAgeAtLength_Mal[,i]>0)) {
        k=k+1
        if (k==10) k=1
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
          MainLab=paste(lbnd[i],"-",ubnd[i],"mm,","n =",ObsRetCatchFreqAtLen[2,i])
        } else {
          MainLab=MainLabel
        }
        mtext(MainLab,las=1,side=3,line=0,adj=0.5, cex=0.6)
        legend("topright", legend=c("Mal_Obs","Mal_Exp"), lty="solid", lwd=c(-1,1),
               pch=c(16,-1), col=c("black","blue"), bty='n', cex=0.8)
      }
    }
  }

  # reset default par options
  par(.pardefault)

}

#' Plot Pearson residuals for proportions at age in each length class, from age and length catch curve model
#'
#' This function provides a plot of Pearson residuals for observed and expected proportions at age in each length class,
#' from the age and length catch curve model. The model is fitted to a sample of fish length and age data, by minimising
#' the overall negative log-likelihood, including the NLL associated with the marginal length composition and a conditional age at length NLL,
#' given the parameters (selectivity, growth and mortality) and data, using nlminb.
#'
#' @param params vector of model parameters (in log space) to be estimated (if FittedRes=NA)
#' @param RefnceAges Reference ages for Schnute growth curve (set to NA for von Bertalanffy growth curve)
#' @param MLL minimum legal length (for setting knife edge retention, set to NA if not assumed)
#' @param GrowthCurveType 1=von Bertalanffy, 2=Schnute
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
#' @param CircleScale size of number increases or decreases sizes of circles for Pearson residuals
#' @param MinLenClFreq minimum number of fish in length class to show Pearson residuals
#' @param FittedRes saved results from GetLengthBasedCatchCurveResults function (model will be refitted if set to NA)
#'
#' @return Plot of Pearson residuals for fit of age and length based catch curve
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
#' MaxLen = 1200
#' LenInc = 20
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified is knife-edged at MLL
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.05
#' RefnceAges = NA
#' GrowthParams = c(Linf, vbK, CVSizeAtAge)
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # 2 sexes, von Bertalanffy
#' # GrowthCurveType = 1 # 1 = von Bertalanffy
#' # Linf = c(700,850)
#' # vbK = c(0.25,0.2)
#' # CVSizeAtAge = c(0.05,0.05)
#' # RefnceAges = NA
#' # GrowthParams = data.frame(Linf=Linf, vbK=vbK, CVSizeAtAge=CVSizeAtAge)
#' # Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#' #                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' lbnd=Res$lbnd
#' midpt=Res$midpt
#' ubnd=Res$ubnd
#' # get data - 1 sex (or combined sexes)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen # 1 sex
#' ObsRetCatchFreqAtLengthAndAge = as.matrix(Res$ObsRetCatchFreqAtLengthAndDecAge) # 1 sex
#' # # get data - 2 sexes
#' # ObsRetCatchFreqAtLen <- data.frame(matrix(nrow = 2, ncol = length(midpt))) # 2 sex
#' # colnames(ObsRetCatchFreqAtLen) <- midpt
#' # ObsRetCatchFreqAtLen[1,] = Res$ObsRetCatchFreqAtLen_Fem
#' # ObsRetCatchFreqAtLen[2,] = Res$ObsRetCatchFreqAtLen_Mal
#' # ObsRetCatchFreqAtLengthAndAge = array(c(unlist(Res$ObsRetCatchFreqAtLengthAndDecAge_Fem), unlist(Res$ObsRetCatchFreqAtLengthAndDecAge_Mal)),
#' #                                       c(nTimeSteps, length(midpt), 2), dimnames=list(rownames(Res$ObsRetCatchFreqAtLengthAndDecAge_Fem),
#' #                                                                                      colnames(Res$ObsRetCatchFreqAtLengthAndDecAge_Fem)))
#' # get params - 1 sex
#' InitFishMort = 0.3 # specify starting parameters
#' InitL50 = 320
#' InitDelta = 50 # L95-L50
#' InitLinf = 800
#' InitvbK = 0.2
#' InitCVSizeAtAge = 0.05
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' params = c(InitFishMort_logit, log(c(InitL50, InitDelta, InitLinf, InitvbK, InitCVSizeAtAge)))
#' # # get params - 2 sexes
#' # InitFishMort = 0.3 # specify starting parameters
#' # InitL50 = 320
#' # InitDelta = 50 # L95-L50
#' # InitLinf = c(800,800)
#' # InitvbK = c(0.25,0.25)
#' # InitCVSizeAtAge = 0.05
#' # InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' # params = c(InitFishMort_logit, log(c(InitL50, InitDelta, InitLinf, InitvbK, InitCVSizeAtAge)))
#' FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
#' # Example with specified selectivity vector
#' # Simulate data
#' SelectivityType=1 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = 1 / (1 + exp(-log(19)*(midpt-400)/(500-400)))
#' SelParams = c(NA, NA) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0
#' # single sex, von Bertalanffy
#' GrowthCurveType = 1 # 1 = von Bertalanffy
#' Linf = 800
#' vbK = 0.2
#' CVSizeAtAge = 0.05
#' GrowthParams = c(Linf, vbK, CVSizeAtAge)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' lbnd=Res$lbnd
#' midpt=Res$midpt
#' ubnd=Res$ubnd
#' # get data - 1 sex (or combined sexes)
#' ObsRetCatchFreqAtLen = Res$ObsRetCatchFreqAtLen # 1 sex
#' ObsRetCatchFreqAtLengthAndAge = as.matrix(Res$ObsRetCatchFreqAtLengthAndDecAge) # 1 sex
#' # get params - 1 sex
#' InitFishMort = 0.3 # specify starting parameters
#' InitLinf = 800
#' InitvbK = 0.2
#' InitCVSizeAtAge = 0.05
#' InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform (so F is always between 0 and 1)
#' params = c(InitFishMort_logit, log(c(InitLinf, InitvbK, InitCVSizeAtAge)))
#' FittedRes=GetAgeAndLengthBasedCatchCurveResults(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
#'
#' PlotAgeLengthCatchCurve_Pears_Resid(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
#'                                                 lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep, MainLabel=NA,
#'                                                 xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA, ymax=NA, yint=NA, CircleScale=5, MinLenClFreq=10, FittedRes)
#'
#' @export
PlotAgeLengthCatchCurve_Pears_Resid <- function(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                                lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep, MainLabel,
                                                xaxis_lab, yaxis_lab, xmax, xint, ymax, yint, CircleScale, MinLenClFreq, FittedRes) {


  MinAge = floor(TimeStep)
  AgeClasses = MinAge:MaxAge
  nAgeCl = length(MinAge:MaxAge)
  nLenCl = length(midpt)

  # get required inputs for plotting
  Res = GetInputsForPlotting_Cond_AL(params, RefnceAges, MLL, GrowthCurveType, SelectivityType, ObsRetCatchFreqAtLen, ObsRetCatchFreqAtLengthAndAge,
                                     lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep, FittedRes)

  # calculate expected proportions at age, for each length class
  # single sex
  if (is.vector(ObsRetCatchFreqAtLen)) {
    ObsCatchPropAgeAtLength = Res$ObsCatchPropAgeAtLength
    ObsRetCatchFreqAtLengthAndIntAge = Res$ObsRetCatchFreqAtLengthAndIntAge
    ExpRetCatchPropIntAgeGivenLength = Res$ExpRetCatchPropIntAgeGivenLength
  }
  # 2 sexes
  if (is.data.frame(ObsRetCatchFreqAtLen)) {
    ObsCatchPropAgeAtLength_Fem = Res$ObsCatchPropAgeAtLength_Fem
    ObsCatchPropAgeAtLength_Mal = Res$ObsCatchPropAgeAtLength_Mal
    ObsRetCatchFreqAtLengthAndIntAge_Fem = Res$ObsRetCatchFreqAtLengthAndIntAge_Fem
    ObsRetCatchFreqAtLengthAndIntAge_Mal = Res$ObsRetCatchFreqAtLengthAndIntAge_Mal
    ExpRetCatchPropIntAgeGivenLength_Fem = Res$ExpRetCatchPropIntAgeGivenLength_Fem
    ExpRetCatchPropIntAgeGivenLength_Mal = Res$ExpRetCatchPropIntAgeGivenLength_Mal
  }


  # calculate the Pearson residuals (for a multinomial distribution)
  # single sex - Pearson residuals

  if (is.vector(ObsRetCatchFreqAtLen)) {
    PearResid <- data.frame(matrix(nrow = nAgeCl, ncol = nLenCl))
    colnames(PearResid) <- midpt
    for (i in 1:nLenCl) {
      PearResid[,i]=0
      if (sum(ObsRetCatchFreqAtLengthAndIntAge[,i]) > 1) {
        Prop_Resid = (ObsCatchPropAgeAtLength[,i] - ExpRetCatchPropIntAgeGivenLength[,i])
        cat('\n')
        Prop_Var = sqrt(ExpRetCatchPropIntAgeGivenLength[,i]*(1-ExpRetCatchPropIntAgeGivenLength[,i]) /
                          sum(ObsRetCatchFreqAtLengthAndIntAge[,i]))
        Prop_Var[which(Prop_Var<=0.01)]=0.01
        PearResid[,i] = Prop_Resid / Prop_Var
      }
    }
  }

  # 2 sexes - Pearson residuals
  if (is.data.frame(ObsRetCatchFreqAtLen)) {
    PearResid_Fem <- data.frame(matrix(nrow = nAgeCl, ncol = nLenCl))
    colnames(PearResid_Fem) <- midpt
    PearResid_Mal <- PearResid_Fem
    for (i in 1:nLenCl) {
      PearResid_Fem[,i] = 0; PearResid_Mal[,i] = 0;
      if (sum(ObsRetCatchFreqAtLengthAndIntAge_Fem[,i]) > 1) {
        Prop_Resid = (ObsCatchPropAgeAtLength_Fem[,i] - ExpRetCatchPropIntAgeGivenLength_Fem[,i])
        Prop_Var = sqrt(ExpRetCatchPropIntAgeGivenLength_Fem[,i]*(1-ExpRetCatchPropIntAgeGivenLength_Fem[,i]) /
                          sum(ObsRetCatchFreqAtLengthAndIntAge_Fem[,i]))
        Prop_Var[which(Prop_Var<=0.01)]=0.01
        PearResid_Fem[,i] = Prop_Resid / Prop_Var
      }
      if (sum(ObsRetCatchFreqAtLengthAndIntAge_Mal[,i]) > 1) {
        Prop_Resid = (ObsCatchPropAgeAtLength_Mal[,i] - ExpRetCatchPropIntAgeGivenLength_Mal[,i])
        Prop_Var = sqrt(ExpRetCatchPropIntAgeGivenLength_Mal[,i]*(1-ExpRetCatchPropIntAgeGivenLength_Mal[,i]) /
                          sum(ObsRetCatchFreqAtLengthAndIntAge_Mal[,i]))
        Prop_Var[which(Prop_Var<=0.01)]=0.01
        PearResid_Mal[,i] = Prop_Resid / Prop_Var
      }

    }
  }

  xlims = Get_xaxis_scale(ubnd)
  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint
  ylims = Get_yaxis_scale(1:MaxAge)
  if (is.na(ymax)) ymax = ylims$ymax
  if (is.na(yint)) yint = ylims$yint
  if (is.na(xaxis_lab)) xaxis_lab = "Length, mm"
  if (is.na(yaxis_lab)) yaxis_lab = "Age class, y"
  par(mfrow=c(1,1), mar=c(5,4,2,2))

  # combined sexes
  if (is.vector(ObsRetCatchFreqAtLen)) {
    plot(1,1, xlim=c(0,xmax), ylim=c(0,ymax), col=0,
         ylab='',xlab='',xaxt='n',yaxt='n')
    for (i in 1:nLenCl) {
      if (sum(ObsRetCatchFreqAtLengthAndIntAge[,i])>MinLenClFreq) {
        for (k in 1:MaxAge) {
          CircleSize = CircleScale * abs(PearResid[k,i])
          # positive residuals
          if (PearResid[k,i] > 0) {
            symbols(midpt[i],k,circles=CircleSize,bg="blue",lwd=1.5,inches=F, add=T)
          }
          # negative residuals
          if (PearResid[k,i] < 0) {
            symbols(midpt[i],k,circles=CircleSize,bg="red",lwd=1.5,inches=F, add=T)
          }
        }
      }
    }
    mtext(yaxis_lab,side=2,cex=1.2,line=2.5)
    mtext(xaxis_lab,side=1,cex=1.2,line=2.5)
    AddAxesAndTickLabelsToPlot(xmin=0, xmax, xint, ymin=0, ymax, yint,
                               cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
    legend("bottomright", legend=c("+ve resid.","-ve resid."), inset=c(0.13,0),
           cex=0.8, bty="n", seg.len = 0, pch=16, border=T, col=c("blue","red"))
  }

  # separate sexes
  if (is.data.frame(ObsRetCatchFreqAtLen)) {
    # females
    plot(1,1, xlim=c(0,xmax), ylim=c(0,ymax), col=0,
         ylab='',xlab='',xaxt='n',yaxt='n')
    for (i in 1:nLenCl) {
      if (sum(ObsRetCatchFreqAtLengthAndIntAge_Fem[,i])>MinLenClFreq) {
        for (k in 1:MaxAge) {
          CircleSize = CircleScale * abs(PearResid_Fem[k,i])
          # positive residuals
          if (PearResid_Fem[k,i] > 0) {
            symbols(midpt[i],k,circles=CircleSize,bg="orange",lwd=1.5,inches=F, add=T)
          }
          # negative residuals
          if (PearResid_Fem[k,i] < 0) {
            symbols(midpt[i],k,circles=CircleSize,bg="red",lwd=1.5,inches=F, add=T)
          }
        }
      }
    }
    mtext(yaxis_lab,side=2,cex=1.2,line=2.5)
    mtext(xaxis_lab,side=1,cex=1.2,line=2.5)
    AddAxesAndTickLabelsToPlot(xmin=0, xmax, xint, ymin=0, ymax, yint,
                               cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
    legend("bottomright", legend=c("Fem. +ve resid.","Fem. -ve resid."), inset=c(0.13,0),
           cex=0.8, bty="n", seg.len = 0, pch=16, border=T, col=c("orange","red"))

    # males
    plot(1,1, xlim=c(0,xmax), ylim=c(0,ymax), col=0,
         ylab='',xlab='',xaxt='n',yaxt='n')
    for (i in 1:nLenCl) {
      if (sum(ObsRetCatchFreqAtLengthAndIntAge_Mal[,i])>MinLenClFreq) {
        for (k in 1:MaxAge) {
          CircleSize = CircleScale * abs(PearResid_Mal[k,i])
          # positive residuals
          if (PearResid_Mal[k,i] > 0) {
            symbols(midpt[i],k,circles=CircleSize,bg="orange",lwd=1.5,inches=F, add=T)
          }
          # negative residuals
          if (PearResid_Mal[k,i] < 0) {
            symbols(midpt[i],k,circles=CircleSize,bg="blue",lwd=1.5,inches=F, add=T)
          }
        }
      }
    }
    mtext(yaxis_lab,side=2,cex=1.2,line=2.5)
    mtext(xaxis_lab,side=1,cex=1.2,line=2.5)
    AddAxesAndTickLabelsToPlot(xmin=0, xmax, xint, ymin=0, ymax, yint,
                               cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA)
    legend("bottomright", legend=c("Mal. +ve resid.","Mal. -ve resid."), inset=c(0.13,0),
           cex=0.8, bty="n", seg.len = 0, pch=16, border=T, col=c("orange","blue"))

  }
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
#' @param MaxAge maximum age of species to be considered by the model
#'
#' @return length class with peak frequency (PeakLencl), last length class included in analysis (LastLenCl),
#' observed data to which catch curve is fitted (ObsCatchFreqAtLen2), relatives age at midpoints of length classes
#' (Age_midptlencl), natural logarithms of changes in numbers with respect to time taken to grow through
#' each length class (ln_n_dt)
#'
Calcs_LenCovertCatchCurve_Schnute <- function(GrowthParams, RefnceAges, ObsRetCatchFreqAtLen,
                                                  MinFreq, lbnd, midpt, ubnd, MaxAge) {


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
    FishLen = lbnd[j]
    Age_lbndlencl[i] = InverseSchnuteGrowthfunction(MaxAge, FishLen, t1, t2, y1, y2, a, b)
    FishLen = ubnd[j]
    Age_ubndlencl[i] = InverseSchnuteGrowthfunction(MaxAge, FishLen, t1, t2, y1, y2, a, b)
    FishLen = midpt[j]
    Age_midptlencl[i] = InverseSchnuteGrowthfunction(MaxAge, FishLen, t1, t2, y1, y2, a, b)
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
#' @param MaxAge maximum age to be considered by model
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
#'                                 MinFreq, lbnd, midpt, ubnd, MaxAge)
#' @export
GetLenConvCatchCurveResults <- function(ModelType, GrowthParams, RefnceAges, ObsRetCatchFreqAtLen,
                                        MinFreq, lbnd, midpt, ubnd, MaxAge) {

  # Pauly's length-converted catch curve
  if (ModelType == 1) {
    res=Calcs_PaulyLenConvertCatchCurve(GrowthParams, ObsRetCatchFreqAtLen,
                                        MinFreq, lbnd, midpt, ubnd)
  }
  if (ModelType == 2) { # length conv catch curve using Schnute function
    res=Calcs_LenCovertCatchCurve_Schnute(GrowthParams, RefnceAges, ObsRetCatchFreqAtLen,
                                          MinFreq, lbnd, midpt, ubnd, MaxAge)
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
                 EstZMort = ZMort,
                 EstZMort_se = ZMort_se,
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
#' SelParams = c(300, 50) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
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
#'                                MinFreq, lbnd, midpt, ubnd, MaxAge)
#' PlotLenConvCatchCurveResults(MaxAge, ModelType, GrowthParams, RefnceAges, ObsRetCatchFreqAtLen, MinFreq,
#'                             lbnd, midpt, ubnd)
#' @export
#'
#'
PlotLenConvCatchCurveResults <- function(MaxAge, ModelType, GrowthParams, RefnceAges, ObsRetCatchFreqAtLen, MinFreq,
                                         lbnd, midpt, ubnd) {

  # .pardefault <- par(no.readonly = TRUE) # store default par settings

    res=GetLenConvCatchCurveResults(ModelType, GrowthParams, RefnceAges, ObsRetCatchFreqAtLen,
                                    MinFreq, lbnd, midpt, ubnd, MaxAge)
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
                 EstZMort = round(Est_ZMort,3),
                 EstZMort_se = round(ZMort_se,3),
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
#' @return total mortality and approximate lower and upper 95 percent confidence
#' limits, (EstZMort), approximate standard error for total mortality (EstZMort_se),
#' median estimate and lower and upper 95 percent confidence limits for Z from resampling (EstZMort_resamp),
#' variables for plotting, including expected frequencies at age with associated confidence limits, calculated from resampling
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
  EstZ = c(ZMort, ZMort + c(-1.96, 1.96) * ZMort_se)
  EstZMort = t(data.frame(EstZ = round(EstZ, 3)))
  colnames(EstZMort) = c("Estimate", "lw_95%CL", "up_95%CL")

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

  EstZ_resamp = c(ZMort_resamp, ZMort_lowresamp, ZMort_upresamp)
  EstZMort_resamp = t(data.frame(EstZ_resamp = round(EstZ_resamp, 3)))
  colnames(EstZMort_resamp) = c("Estimate", "lw_95%CL", "up_95%CL")

  results = list(EstZMort = EstZMort,
                 EstZMort_se = round(ZMort_se, 3),
                 EstZMort_resamp = EstZMort_resamp,
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
#' This function returns the multinomial or Dirichlet multinomial negative log-likelihood, given age
#' frequency data, a value total mortality, and age-based selectivity parameters, for a catch curve with
#' age-based, logistic selectivity. Function requires an estimate of natural mortality (NatMort), a value
#' for maximum age (MaxAge) and age frequency data (stored in memory in R)
#'
#' @keywords internal
#'
#' @param params model parameters log(c(FMort, SelA50, SelDelta) for multinomial NLL or
#' log(c(FMort, SelA50, SelDelta, theta) for Dirichlet multinomial NLL
#'
#' @return negative log-likelihood (NLL)
Calculate_NLL_LogisticCatchCurve <- function(params) {

  FMort = exp(params[1])
  SelA50 = exp(params[2])
  SelDelta = exp(params[3])
  SelA95 = SelA50 + SelDelta
  # inverse logit transformed value
  if (length(params)==4) {
    temp = params[4]
    DM_theta = 1/(1+exp(-temp));
  }

  SelAtAge = rep(0,length(Ages))
  SelAtAge = 1 / (1 + exp(-log(19) * (Ages - SelA50) / (SelA95 - SelA50)))
  FAtAge = SelAtAge * FMort
  ZAtAge = NatMort + FAtAge

  N = numeric(length(Ages))
  N[1] = 1
  k=1
  MinAge = min(Ages)
  MaxAge = max(Ages)
  for (i in seq(MinAge+1,MaxAge,1)) {
    k=k+1
    if (i < MaxAge) {
      N[k] = N[k-1] * exp(-ZAtAge[k-1])
    } else {
      N[k] = N[k-1] * exp(-ZAtAge[k-1]) / (1 - exp(-ZAtAge[k]))
    }
  }
  CatchAtAge = N * (FAtAge / ZAtAge) * (1 - exp(-ZAtAge))
  ExpPropAtAge = CatchAtAge / sum(CatchAtAge)

  # calculate F penalty
  F_Pen = 0
  if (FMort > 2.0) {
    F_Pen = 1000 * (FMort - 2.0)^2
  }

  # calculate multinomial negative log-likelihood
  if (length(params==3)) {
    NLL = -sum((ObsAgeFreq * log(ExpPropAtAge + 1E-4))) + F_Pen
  }

  # calculate Dirichlet multinomial negative log-likelihood
  if (length(params)==4) {
    SampleSize = sum(ObsAgeFreq)
    nAges = length(ObsAgeFreq)
    ObsPropAtAge = ObsAgeFreq/SampleSize
    sum1 = 0; sum2 = 0; NLL = 0
    for (t in 1:nAges) {
      sum1 = sum1 + lgamma(SampleSize * ObsPropAtAge[t] + 1)
      sum2 = sum2 + (lgamma(SampleSize * ObsPropAtAge[t] + DM_theta * SampleSize * ExpPropAtAge[t])
                     - lgamma(DM_theta * SampleSize * ExpPropAtAge[t]))
    }
    NLL = -(lgamma(SampleSize+1) - sum1 + (lgamma(DM_theta * SampleSize) - lgamma(SampleSize + DM_theta * SampleSize)) + sum2)
    NLL = NLL + F_Pen

  }

  if (length(params)==3) {
    cat("NLL",NLL,"F_Pen",F_Pen,"params",exp(params),'\n')
  } else {
    cat("NLL",NLL,"F_Pen",F_Pen,"params",c(FMort, SelA50, SelDelta, DM_theta),'\n')
  }

  return(NLL)

}

#' Get statistical outputs from a fitted catch curve with age-based, logistic selectivity
#'
#' This function fits a catch curve with an asymptotic, age-based logistic selectivity curve,
#' to a sample of fish age frequency data, by minimising either a multinomial or Dirichlet
#' multinomial the negative log-likelihood associated with the parameters and data, using nlminb.
#' It provides various statistical outputs including convergence statistics, parameter estimates
#' and associated 95 percent confidence limits and associated variance-covariance matrix, calculated
#' using the MASS package
#'
#' @param params model parameters log(c(FMort, SelA50, SelDelta) for multinomial NLL or
#' log(c(FMort, SelA50, SelDelta, theta) for Dirichlet multinomial NLL
#' @param NatMort natural mortality
#' @param Ages ages in observed data
#' @param ObsAgeFreq observed age frequency data
#'
#' @return negative log-likelihood (nll), nlminb convergence diagnostic (convergence)
#' sample size (SampleSize), parameter estimates with lower and upper 95 percent
#' confidence limits (ParamEst), point estimates for parameters estimated in log space (Estlnparams),
#' variance-covariance matrix for estimated parameters (vcov.Params), parameter correlation matrix (cor.Params),
#' standard errors for estimated parameters (lnEstFMort_se, lnEstSelA50_se, lnEstFSelA95_se), selectivity at age (SelAtAge), fishing
#' mortality at age (FAtAge), estimated frequencies at age with associated 95 percent confidence limits
#' (EstFreq, EstFreq_Zlow, EstFreq_Zup), random values of parameters in log space from parametric resampling,
#' using rmultinom function (lnparams.sims), and associated median and lower 2.5 and upper 97.5 percentiles
#' of estimates for frequency at age (EstFreq.sim), selectivity parameters in normal space (SelA50.sim, SelA95.sim),
#' fishing mortality (FMort.sim) and total mortality (EstZMort.sim), estimated parameters from saved par
#' object (Estparams), Dirichlet multinomial effective sample size estimate (EffSampSize), if using this
#' objective function
#' @examples
#' # fit logistic catch cuve model using multinomial likelihood
#' # simulate data
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
#' # fit model
#' Init_FMort = 0.2
#' Init_SelA50 = 3
#' Init_SelDelta = 2
#' params = log(c(Init_FMort, Init_SelA50, Init_SelDelta))
#' res=GetLogisticCatchCurveResults(params, NatMort, Ages, ObsAgeFreq)
#' # fit logistic catch curve model using Dirichlet multinomial likelihood
#' # get expected catch proportions at age, for simulating data
#' MinAge = 1
#' MaxAge = 40
#' Ages = MinAge:MaxAge
#' NatMort <- exp(1.46 - (1.01 * (log(MaxAge)))) # i.e. Hoenig's (1983) eqn for fish
#' FMort = 0.1
#' SelA50 = 5
#' SelA95 = 7
#' SelAtAge = rep(0, length(Ages))
#' N = rep(0, length(Ages))
#' SelAtAge = 1/(1 + exp(-log(19) * (Ages - SelA50)/(SelA95 - SelA50)))
#' FAtAge = SelAtAge * FMort
#' ZAtAge = NatMort + FAtAge
#' k = 1
#' N[1] = 1
#' for (i in seq(MinAge + 1, MaxAge, 1)) {
#'   k = k + 1
#'   if (i < MaxAge) {
#'     N[k] = N[k-1] * exp(-ZAtAge[k-1])
#'   }
#'   else {
#'     N[k] = N[k-1] * exp(-ZAtAge[k-1])/(1 - exp(-ZAtAge[k]))
#'   }
#' }
#' CatchAtAge = N * (FAtAge/ZAtAge) * (1 - exp(-ZAtAge))
#' PropAtAge = CatchAtAge/sum(CatchAtAge)
#' library(dirmult)
#' # Simulate data from a Dirichlet multinomial distribution
#' # J = number of fish sampling events
#' # K = number of age classes
#' # n = number of fish sampled from each sampling event
#' # pi = expected proportion at age
#' # theta = amount of autocorrelation between ages of fish within sampling events
#' set.seed(123)
#' theta_val = 0.3
#' simDat = simPop(J=50, K=nAges, n=10, pi=PropAtAge, theta=theta_val)
#' simAges = data.frame(simDat$data)
#' colnames(simAges)=Ages
#' simAgeFreq = colSums(simAges)
#' ObsAgeFreq = as.vector(simAgeFreq)
#' # fit catch curve
#' Init_FMort = 0.2
#' Init_SelA50 = 3
#' Init_SelDelta = 2
#' Init_theta = 0.5
#' Init_theta_logit = log(Init_theta/(1-Init_theta)) # logit transform (so theta is always between 0 and 1)
#' params = c(log(Init_FMort), log(Init_SelA50), log(Init_SelDelta), Init_theta_logit)
#' res=GetLogisticCatchCurveResults(params, NatMort, Ages, ObsAgeFreq)
#' @export
GetLogisticCatchCurveResults <- function (params, NatMort, Ages, ObsAgeFreq)
{

  nlmb = nlminb(params, Calculate_NLL_LogisticCatchCurve)
  hess.out = optimHess(nlmb$par, Calculate_NLL_LogisticCatchCurve)
  vcov.Params = solve(hess.out)
  ses = sqrt(diag(vcov.Params))

  # compute parameter correlation matrix
  temp = diag(1/sqrt(diag(vcov.Params)))
  cor.Params=temp %*% vcov.Params %*% temp

  EstFMort = c(exp(nlmb$par[1]), exp(nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
  EstSelA50 = c(exp(nlmb$par[2]), exp(nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
  EstSelDelta = c(exp(nlmb$par[3]), exp(nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
  EstSelA95 = EstSelA50[1] + EstSelDelta[1]

  # inverse logit transformed value for Dirichlet multinomial theta parameter
  if (length(params)==4) {
    temp = c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4])
    EstTheta = 1/(1+exp(-temp));
  }

  SelAtAge = rep(0, length(Ages))
  SelAtAge = 1/(1 + exp(-log(19) * (Ages - EstSelA50[1])/(EstSelA95[1] - EstSelA50[1])))
  FAtAge = SelAtAge * EstFMort[1]
  ZAtAge = NatMort + FAtAge

  if (length(params)==3) {
    ParamEst = t(data.frame(FMort = round(EstFMort, 3), SelA50 = round(EstSelA50,3),
                            EstSelDelta = round(EstSelDelta, 3)))
  }
  if (length(params)==4) {
    ParamEst = t(data.frame(FMort = round(EstFMort, 3), SelA50 = round(EstSelA50,3),
                            EstSelDelta = round(EstSelDelta, 3), Theta = round(EstTheta, 3)))
  }
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
  SelDelta.sim = rep(0,500)
  SelA95.sim = rep(0,500)
  FMort.sim = rep(0,500)
  EstZMort.sim = rep(0,500)
  if (length(params)==4) EstTheta.sim = rep(0,500)

  for (j in 1:500) {
    FMort.sim[j] = exp(sims[j, 1])
    SelA50.sim[j] = exp(sims[j, 2])
    SelDelta.sim[j] = exp(sims[j, 3])
    SelA95.sim[j] = SelA50.sim[j] + SelDelta.sim[j]

    if (length(params)==4) {
      temp = sims[j, 4]
      EstTheta.sim= 1/(1+exp(-temp))
    }

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
        N.sim[k] = N.sim[k-1] * exp(-ZAtAge.sim[k-1])
      }
      else {
        N.sim[k] = N.sim[k-1] * exp(-ZAtAge.sim[k-1])/(1 - exp(-ZAtAge.sim[k]))
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

  ModelDiag = list(lnEstFMort.se = ses[1],
                   lnEstSelA50.se = ses[2],
                   lnEstFSelA95.se = ses[3],
                   SelAtAge = SelAtAge,
                   FAtAge = FAtAge,
                   EstFreq = EstFreq,
                   EstFreq_Zlow = EstFreq_Zlow,
                   EstFreq_Zup = EstFreq_Zup,
                   params.sims = sims,
                   EstFreq.sim = EstFreq.sim,
                   SelA50.sim = SelA50.sim,
                   SelDelta.sim = SelDelta.sim,
                   SelA95.sim = SelA95.sim,
                   FMort.sim = FMort.sim,
                   EstZMort.sim = EstZMort.sim,
                   Estparams = nlmb$par)

  results = list(nll = nll,
                 convergence = convergence,
                 SampleSize = SampleSize,
                 ParamEst = ParamEst,
                 EstFMort = EstFMort,
                 EstFMort_se = ses[1],
                 vcov.Params = vcov.Params,
                 cor.Params = cor.Params,
                 ModelDiag = ModelDiag)

  # Dirichlet multinomial effective sample size
  if (length(params)==4) {
    temp = ((1+ (EstTheta*sum(ObsAgeFreq))) / (1 + EstTheta))
    EffSampSize <- t(data.frame(EffSampSize = temp))
    colnames(EffSampSize) = c("Estimate", "lw_95%CL", "up_95%CL")
    results$EffSampSize = EffSampSize
  }


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
    Z_value = round(Res$EstZMort[1],digits=3)
    x=which(Ages==Res$RecAge) # RecAge position
    xx=length(which(Res$EstFreq_Zup>0))+x # position of last age
    xxx=length(which(Res$EstFreq_Zlow>0))+x # position of last age
    xxxx=length(which(Res$EstFreq>0))+x # position of last age
  }
  # Linear
  if (CatchCurveModel == 2) {
    Z_value = round(Res$EstZMort[1],digits=3)
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
    if (CatchCurveModel == 3) {
      sm1 = spline(Ages[j], Res$ModelDiag$EstFreq_Zlow[1:length(jj)], n=100, method="natural")
      sm2 = spline(Ages[j], Res$ModelDiag$EstFreq_Zup[1:length(j)], n=100, method="natural")
    } else {
      sm1 = spline(Ages[jj], Res$EstFreq_Zlow[1:length(jj)], n=100, method="natural")
      sm2 = spline(Ages[j], Res$EstFreq_Zup[1:length(j)], n=100, method="natural")
    }
    if (length(which(sm1$y<0))>0) sm1$y[1:max(which(sm1$y<0))]=0
    if (length(which(sm2$y<0))>0) sm2$y[1:max(which(sm2$y<0))]=0
    x = c(sm1$x, rev(sm2$x)) # using shading for 95% CLs
    y = c(sm1$y, rev(sm2$y))
    polygon(x,y, col="pink",border=NA)
  }
  # lines(Ages[jjj], Res$EstFreq[1:length(jjj)], col="black")
  points(Ages, ObsAgeFreq, pch=16, cex=0.8)
  if (CatchCurveModel == 3) {
    points(Ages[j], Res$ModelDiag$EstFreq[1:length(j)], col="red", pch=1, cex=0.8)
  } else {
    points(Ages[j], Res$EstFreq[1:length(j)], col="red", pch=1, cex=0.8)
  }
  if (PlotCLs == FALSE) { # if not plotting confidence intervals, can include last age
    xx=which(Ages==max(Ages)) # position of last age
    points(Ages[x:xx], Res$EstFreq[1:length(x:xx)], col="red", pch=1, cex=0.8)
  }

  legend("topright", legend=bquote(paste("Z = ", .(Z_value), " ",y^-1)), y.intersp = 1.5, inset=c(0.13,0),
         lty=1, cex = 1, bty="n",seg.len = 0)
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
    Z_value = round(Res$EstZMort[1],digits=3)
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
    Z_value = round(Res$EstZMort[1],digits=3)
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
    x=min(which(log(Res$ModelDiag$EstFreq_Zlow)>0)) # RecAge position
    xx=length(which(log(Res$ModelDiag$EstFreq_Zlow) > -1)) # position of last age
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
    points(Ages[jjj], log(Res$EstFreq[1:length(jjj)]), pch=1, col="red", cex=0.6)
  }
  if (CatchCurveModel == 3) {
    if (PlotCLs == TRUE) {
      sm1 = spline(Ages[j], log(Res$ModelDiag$EstFreq_Zup[j]), n=100, method="natural")
      sm2 = spline(Ages[j], log(Res$ModelDiag$EstFreq_Zlow[j]), n=100, method="natural")
      x = c(sm1$x, rev(sm2$x)) # using shading for 95% CLs
      y = c(sm1$y, rev(sm2$y))
      polygon(x,y, col="pink",border=NA)
    }
    points(Ages[j], log(Res$ModelDiag$EstFreq[j]), pch=1, col="red", cex=0.6)
  }

  points(Ages, log(ObsAgeFreq), pch=16, cex=0.8)
  legend("topright", legend=bquote(paste("Z = ", .(Z_value), " ",y^-1)), y.intersp = 1.5, inset=c(0.13,0),
         lty=1, cex = 1, bty="n",seg.len = 0)
}


#*********************
# Per recruit analyses
#*********************


#' Calculates weight at age for each sex for age-based per recruit analysis
#'
#' This function calculates mean weight at age for each sex, based on power curve,
#' log-log relationship, or allows for these to be used as vectors inputted by user,
#' for age-based per recruit analysis
#'
#' @keywords internal
#'
#' @param lenwt_a weight-length parameter
#' @param ln_lenwt_a weight-length parameter
#' @param lenwt_b weight-length parameter
#' @param WLrel_Type 1=power, 2=log-log
#' @param EstWtAtAge user-specified weights at ages
#' @param FemLenAtAge female lengths at each age (from growth curve)
#' @param MalWtAtAge male lengths at each age (from growth curve)
#'
#' @return FemWtAtAge
#' @return MalWtAtAge
CalcWeightAtAge <- function(lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, FemLenAtAge, MalLenAtAge) {

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

  results = list(FemWtAtAge=FemWtAtAge,
                 MalWtAtAge=MalWtAtAge)

  return(results)

}


#' Calculates length at age for each sex for age-based per recruit analysis
#'
#' This function calculates mean length at age for each sex, based on von Bertalanffy growth curve or
#' allows for these to be used as vectors inputted by user, for age-based per recruit analysis
#'
#' @keywords internal
#'
#' @param Linf von Bertalanffy growth parameter
#' @param vbK von Bertalanffy growth parameter
#' @param tzero von Bertalanffy growth parameter
#' @param EstLenAtAge lengths at age, inputted as vector
#' @param Ages
#'
#' @return FemLenAtAge
#' @return MalLenAtAge
CalcLengthAtAge <- function(Linf, vbK, tzero, EstLenAtAge, Ages) {


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

  results = list(FemLenAtAge=FemLenAtAge,
                 MalLenAtAge=MalLenAtAge)

  return(results)

}


#' Calculate proportion of fish at age that are females
#'
#' This function calculates proportion of fish at age that are females,
#' based on reproductive pattern (gonochoristic or hermaphroditic), and
#' sex change parameters
#'
#' @keywords internal
#'
#' @param FinalSex_A50 logistic sex change parameter
#' @param FinalSex_A95 logistic sex change parameter
#' @param EstSexRatioAtAge input vecto sex change at length
#' @param ReprodPattern 1=gonochoristic, 2=protogynous, 3=protandrous
#' @param Ages specified ages
#'
#' @return PropFemAtAge
CalcPropFemAtAge <- function(FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, ReprodPattern, Ages) {

  # for sex-changing species, calculate proportion female at age
  if (ReprodPattern == 1) { # gonochoristic (separate sexes)
    PropFemAtAge = NA
  }
  if (ReprodPattern == 2) { # protogynous hermaphroditism (female to male sex change)
    if (is.na(EstSexRatioAtAge[1])) {
      PropFemAtAge = 1 - (1 / (1 + exp(-log(19) * (Ages-FinalSex_A50) / (FinalSex_A95 - FinalSex_A50))))
    }
    if (!is.na(EstSexRatioAtAge[1])) {
      PropFemAtAge <- EstSexRatioAtAge
    }
  }
  if (ReprodPattern == 3) { # protandrous hermaphroditism (male to female sex change)
    if (is.na(EstSexRatioAtAge[1])) {
      PropFemAtAge = 1 / (1 + exp(-log(19) * (Ages-FinalSex_A50) / (FinalSex_A95 - FinalSex_A50)))
    }
    if (!is.na(EstSexRatioAtAge[1])) {
      PropFemAtAge <- EstSexRatioAtAge
    }
  }

  results = PropFemAtAge

  return(results)

}

#' Calculate proportion of fish at age that are mature
#'
#' This function calculates proportion of fish at age that are mature,
#' based on input maturity parameters, or from set input maturity schedules
#'
#' @keywords internal
#'
#' @param mat_A50 logistic maturity parameter
#' @param mat_A95 logistic maturity parameter
#' @param EstMatAtAge input maturity schedules (could be set to NA)
#' @param Ages ages
#'
#' @return FemPropMatAtAge; MalPropMatAtAge
CalcPropMatureAtAge <- function(mat_A50, mat_A95, EstMatAtAge, Ages) {

  if (is.na(EstMatAtAge[1,1])) {
    FemPropMatAtAge <- 1 / (1 + exp(-log(19) * (Ages-mat_A50[1]) / (mat_A95[1] - mat_A50[1])))
    MalPropMatAtAge <- 1 / (1 + exp(-log(19) * (Ages-mat_A50[2]) / (mat_A95[2] - mat_A50[2])))
  }
  if (!is.na(EstMatAtAge[1,1])) {
    FemPropMatAtAge <- EstMatAtAge[,1]
    MalPropMatAtAge <- EstMatAtAge[,2]
  }

  results = list(FemPropMatAtAge=FemPropMatAtAge,
                 MalPropMatAtAge=MalPropMatAtAge)

  return(results)

}

#' Calculate population sex ratio and fertilisation rate (length-based per recruit analysis)
#'
#' This function calculates the population sex ratio, by number, based on mature fish,
#' and expected egg fertilisation rate, based on value (NA or between 0.2-1) for the
#' parameter, EggFertParam. This function is of the same form as a Beverton-Holt stock
#' recruitment relationship, with a value of 0.2 for EggFertParam resulting in a direct
#' effect of current sex ratio, relative to the unfished level, on egg fertilisation rate,
#' vs no effect if EggFertParam = 1. The value of EggFertParam is set based on understanding
#' of fish spawning behaviour, using a higher value of EggFertParam for group spawning vs lower
#' value for pair spawning.
#'
#' @keywords internal
#'
#' @param Unfish_MalNPerRecLen unfished male survival at length (in numbers)
#' @param Fish_MalNPerRecLen fished male survival at length (in numbers)
#' @param Unfish_FemNPerRecLen unfished female survival at length (in numbers)
#' @param Fish_FemNPerRecLen fished female survival at length (in numbers)
#' @param MalPropMatAtLen proportion of males mature at length
#' @param FemPropMatAtLen proportion of females mature at length
#' @param EggFertParam egg fertilisation parameters (NA, or ~0.2-1)
#'
#' @return UnfishMalToFemRatio, FishMalToFemRatio, MalDeplRatio, Eq_FertRate
CalcPopnSexRatioAndFertRate_LB <- function(Unfish_MalNPerRecLen, Fish_MalNPerRecLen, Unfish_FemNPerRecLen, Fish_FemNPerRecLen,
                                           MalPropMatAtLen, FemPropMatAtLen, EggFertParam) {

  # calculate ratio of mature males to mature females, by number (for calculating egg fertilisation rates)
  UnfishMalToFemProp <- sum(Unfish_MalNPerRecLen * MalPropMatAtLen) /
    sum((Unfish_FemNPerRecLen * FemPropMatAtLen) + (Unfish_MalNPerRecLen * MalPropMatAtLen))
  FishMalToFemProp <- sum(Fish_MalNPerRecLen * MalPropMatAtLen) /
    sum((Fish_FemNPerRecLen * FemPropMatAtLen) + (Fish_MalNPerRecLen * MalPropMatAtLen))
  MalDeplRatio <- FishMalToFemProp / UnfishMalToFemProp

  # calculate egg fertilisation rate
  if (!is.na(EggFertParam)) {
    # Brooks et al. (2008)
    if (EggFertParam >= 0.2 & EggFertParam <= 1) {
      Eq_FertRate = (4 * EggFertParam * MalDeplRatio) /
        ((1 - EggFertParam) + (5 * EggFertParam - 1) * MalDeplRatio)
    } else {
      cat("Problem: Eq_FertRate must be >= 0.2 and <= 1.0")
    }
  }
  if (is.na(EggFertParam)) {
    Eq_FertRate = 1
  }


  results = list(UnfishMalToFemProp=UnfishMalToFemProp,
                 FishMalToFemProp=FishMalToFemProp,
                 MalDeplRatio=MalDeplRatio,
                 Eq_FertRate=Eq_FertRate)

  return(results)

}


#' Calculate population sex ratio and fertilisation rate (age-based per recruit analysis)
#'
#' This function calculates the population sex ratio, by number, based on mature fish,
#' and expected egg fertilisation rate, based on value (NA or between 0.2-1) for the
#' parameter, EggFertParam. This function is of the same form as a Beverton-Holt stock
#' recruitment relationship, with a value of 0.2 for EggFertParam resulting in a direct
#' effect of current sex ratio, relative to the unfished level, on egg fertilisation rate,
#' vs no effect if EggFertParam = 1. The value of EggFertParam is set based on understanding
#' of fish spawning behaviour, using a higher value of EggFertParam for group spawning vs lower
#' value for pair spawning.
#'
#' @keywords internal
#'
#' @param UnfishMalSurvAtAge unfished male survival at age (in numbers)
#' @param FishMalSurvAtAge fished male survival at age (in numbers)
#' @param UnfishFemSurvAtAge unfished female survival at age (in numbers)
#' @param FishFemSurvAtAge fished female survival at age (in numbers)
#' @param MalPropMatAtAge proportion of males mature at age
#' @param FemPropMatAtAge proportion of females mature at age
#' @param EggFertParam egg fertilisation parameters (NA, or ~0.2-1)
#'
#' @return UnfishMalToFemRatio, FishMalToFemRatio, MalDeplRatio, Eq_FertRate
CalcPopnSexRatioAndFertRate_AB <- function(UnfishMalSurvAtAge, FishMalSurvAtAge, UnfishFemSurvAtAge, FishFemSurvAtAge,
                                        MalPropMatAtAge, FemPropMatAtAge, EggFertParam) {

  # calculate ratio of mature males to mature females, by number (for calculating egg fertilisation rates)
  UnfishMalToFemProp <- sum(UnfishMalSurvAtAge * MalPropMatAtAge) /
    sum((UnfishFemSurvAtAge * FemPropMatAtAge) + (UnfishMalSurvAtAge * MalPropMatAtAge))
  FishMalToFemProp <- sum(FishMalSurvAtAge * MalPropMatAtAge) /
    sum((FishFemSurvAtAge * FemPropMatAtAge) + (FishMalSurvAtAge * MalPropMatAtAge))
  MalDeplRatio <- FishMalToFemProp / UnfishMalToFemProp

  # calculate egg fertilisation rate
  if (!is.na(EggFertParam)) {
    # Brooks et al. (2008)
    if (EggFertParam >= 0.2 & EggFertParam <= 1) {
      Eq_FertRate = (4 * EggFertParam * MalDeplRatio) /
        ((1 - EggFertParam) + (5 * EggFertParam - 1) * MalDeplRatio)
    } else {
      cat("Problem: Eq_FertRate must be >= 0.2 and <= 1.0")
    }
  }
  if (is.na(EggFertParam)) {
    Eq_FertRate = 1
  }


  results = list(UnfishMalToFemProp=UnfishMalToFemProp,
                 FishMalToFemProp=FishMalToFemProp,
                 MalDeplRatio=MalDeplRatio,
                 Eq_FertRate=Eq_FertRate)

  return(results)

}


#' Calculate unfished and fished survival of fish at age
#'
#' Calculate unfished and fished survival of fish at age. Rountine also allows for sex change
#' for  hermaphroditic species.
#'
#' @keywords internal
#'
#' @param nTimeSteps number of model timesteps
#' @param InitRatioFem sex ratio at birth
#' @param NatMort natural mortality
#' @param FemZAtAge female total mortality at age
#' @param MalZAtAge male total mortality at age
#' @param ReprodPattern 1=gonochoristic, 2=protogynous, 3=protandrous
#' @param Ages PropFemAtAge proportion of fish at age that are female
#'
#' @return UnfishFemSurvAtAge, UnfishMalSurvAtAge, FishFemSurvAtAge, FishMalSurvAtAge
CalcSurvivalAtAge <- function(nTimeSteps, InitRatioFem, NatMort, FMort, FemZAtAge, MalZAtAge,
                              ReprodPattern, PropFemAtAge) {

  # calculate relative unfished female and male survival at age
  UnfishFemSurvAtAge <- rep(0,nTimeSteps); UnfishMalSurvAtAge <- rep(0,nTimeSteps)
  tempUnfishFemSurvAtAge <- rep(0,nTimeSteps); tempUnfishMalSurvAtAge <- rep(0,nTimeSteps)
  FishFemSurvAtAge <- rep(0,nTimeSteps); FishMalSurvAtAge <- rep(0,nTimeSteps)
  tempFishFemSurvAtAge <- rep(0,nTimeSteps); tempFishMalSurvAtAge <- rep(0,nTimeSteps)

  for (j in 1:nTimeSteps) {
    if (j==1) { # age=0
      UnfishFemSurvAtAge[j] <- InitRatioFem
      UnfishMalSurvAtAge[j] <- (1 - InitRatioFem)
      FishFemSurvAtAge[j] <- InitRatioFem
      FishMalSurvAtAge[j] <- (1 - InitRatioFem)
    }
    else if (j < nTimeSteps) {
      UnfishFemSurvAtAge[j] <- UnfishFemSurvAtAge[j-1] * exp(-NatMort * TimeStep)
      UnfishMalSurvAtAge[j] <- UnfishMalSurvAtAge[j-1] * exp(-NatMort * TimeStep)
      FishFemSurvAtAge[j] <- FishFemSurvAtAge[j-1] * exp(-FemZAtAge[j-1] * TimeStep)
      FishMalSurvAtAge[j] <- FishMalSurvAtAge[j-1] * exp(-MalZAtAge[j-1] * TimeStep)
    }
    else if (j==nTimeSteps) { # maximum model age, plus group
      UnfishFemSurvAtAge[j] <- (UnfishFemSurvAtAge[j-1] * exp(-NatMort * TimeStep)) / (1 - exp(-NatMort))
      UnfishMalSurvAtAge[j] <- (UnfishMalSurvAtAge[j-1] * exp(-NatMort * TimeStep)) / (1 - exp(-NatMort))
      FishFemSurvAtAge[j] <- FishFemSurvAtAge[j-1] * exp(-FemZAtAge[j-1] * TimeStep) / (1 - exp(-(FMort+NatMort)))
      FishMalSurvAtAge[j] <- FishMalSurvAtAge[j-1] * exp(-MalZAtAge[j-1] * TimeStep) / (1 - exp(-(FMort+NatMort)))
    }

    if (ReprodPattern > 1) { # hermaphroditic species
      # assuming that, regardless of fishing pattern, the sex ratio at age is determined by
      # by the logistic curve describing probability of sex change at age
      tempUnfishFemSurvAtAge[j] = PropFemAtAge[j] * (UnfishFemSurvAtAge[j] + UnfishMalSurvAtAge[j])
      tempUnfishMalSurvAtAge[j] = (1 - PropFemAtAge[j]) * (UnfishFemSurvAtAge[j] + UnfishMalSurvAtAge[j])
      UnfishFemSurvAtAge[j] = tempUnfishFemSurvAtAge[j]
      UnfishMalSurvAtAge[j] = tempUnfishMalSurvAtAge[j]

      tempFishFemSurvAtAge[j] = PropFemAtAge[j] * (FishFemSurvAtAge[j] + FishMalSurvAtAge[j])
      tempFishMalSurvAtAge[j] = (1 - PropFemAtAge[j]) * (FishFemSurvAtAge[j] + FishMalSurvAtAge[j])
      FishFemSurvAtAge[j] = tempFishFemSurvAtAge[j]
      FishMalSurvAtAge[j] = tempFishMalSurvAtAge[j]
    }
  }

  results = list(UnfishFemSurvAtAge = UnfishFemSurvAtAge,
                 UnfishMalSurvAtAge = UnfishMalSurvAtAge,
                 FishFemSurvAtAge = FishFemSurvAtAge,
                 FishMalSurvAtAge = FishMalSurvAtAge)
  return(results)

}

#' Calculate fishing mortality at age (discards and landings)
#'
#' Calculate fishing mortality at age associated with discarding, landings, and
#' combined.
#'
#' @keywords internal
#'
#' @param FMort fully-selected fishing mortality
#' @param DiscMort proportion of discarded fish that die
#' @param FemSelDiscAtAge selectivity of discarding at age for females
#' @param MalSelDiscAtAge selectivity of discarding at age for males
#' @param FemSelLandAtAge selectivity of landings at age for females
#' @param MalSelLandAtAge selectivity of landings at age for males
#'
#' @return FemDiscFAtAge, MalDiscFAtAge, FemLandFAtAge, MalLandFAtAge, FemFAtAge, MalFAtAge
CalcFishingMortalityAtAge <- function(FMort, DiscMort, FemSelDiscAtAge, MalSelDiscAtAge,
                                      FemSelLandAtAge, MalSelLandAtAge) {

  # calculate female and male fishing mortality at age associated with discarding of undersize fish
  FemDiscFAtAge <- FMort * DiscMort * FemSelDiscAtAge
  MalDiscFAtAge <- FMort * DiscMort * MalSelDiscAtAge

  # calculate female and male fishing mortality at age associated with landings
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

  results = list(FemDiscFAtAge = FemDiscFAtAge,
                 MalDiscFAtAge = MalDiscFAtAge,
                 FemLandFAtAge = FemLandFAtAge,
                 MalLandFAtAge = MalLandFAtAge,
                 FemFAtAge = FemFAtAge,
                 MalFAtAge = MalFAtAge)

  return(results)

}

#' Calculate gear selectivity and retention at age
#'
#' Calculate gear selectivity and retention at age, given input parameters or vectors
#'
#' @keywords internal
#'
#' @param EstGearSelAtAge input vector for gear selectivity at age (may be set to NA)
#' @param EstRetenAtAge input vector for fish retention at age (may be set to NA)
#' @param Ages ages for analysis calculated from max age and model timestep
#' @param Gear_sel_A50 logistic gear selectivity curve parameter
#' @param Gear_sel_A95 logistic gear selectivity curve parameter
#' @param Land_sel_A50 logistic parameter for selectivity of landings curve
#' @param Land_sel_A95 logistic parameter for selectivity of landings curve
#' @param EstLandSelAtAge vector for for selectivity of landings at age (set to NA if using age at selectivity of landings parameters)
#' @param ret_A50 logistic gear retention curve parameter
#' @param ret_A95 logistic gear retention curve parameter
#' @param ret_Pmax logistic gear retention curve parameter
#' @param DiscMort proportion of fish that die following capture and release
#'
#' @return FemGearSelAtAge, MalGearSelAtAge, FemRetProbAtAge, MalRetProbAtAge, FemSelLandAtAge, MalSelLandAtAge, FemSelDiscAtAge, MalSelDiscAtAge
CalcSelectivityAndRetentionAtAge <- function(EstGearSelAtAge, EstRetenAtAge, Ages, Gear_sel_A50, Gear_sel_A95,
                                             Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_A50, ret_A95, ret_Pmax, DiscMort) {

  FemGearSelAtAge <- NA; MalGearSelAtAge <- NA
  if (!is.na(EstGearSelAtAge[1,1])) {
    FemGearSelAtAge <- EstGearSelAtAge[,1]
    MalGearSelAtAge <- EstGearSelAtAge[,2]
  } else if (!is.na(Gear_sel_A50[1])) {
    FemGearSelAtAge <- 1/(1+exp(-log(19) * (Ages - Gear_sel_A50[1]) / (Gear_sel_A95[1] - Gear_sel_A50[1])))
    MalGearSelAtAge <- 1/(1+exp(-log(19) * (Ages - Gear_sel_A50[2]) / (Gear_sel_A95[2] - Gear_sel_A50[2])))
  }
  # note, gear selectivity may be unknown, and need to be calculated from fish retention and selectivity of landings

  # Calculate probability of fish retention at age
  FemRetProbAtAge <- NA; MalRetProbAtAge <- NA
  if (!is.na(EstRetenAtAge[1,1])) {
    FemRetProbAtAge <- EstRetenAtAge[,1]
    MalRetProbAtAge <- EstRetenAtAge[,2]
  } else if (!is.na(ret_A50[1])) {
    FemRetProbAtAge <- ret_Pmax[1] / (1+exp(-log(19) * (Ages - ret_A50[1]) / (ret_A95[1] - ret_A50[1])))
    MalRetProbAtAge <- ret_Pmax[2] / (1+exp(-log(19) * (Ages - ret_A50[2]) / (ret_A95[2] - ret_A50[2])))
  }
  # note, probability of fish retention may be unknown, and need to be calculated from selectivity of landings and gear selectivity


  # Calculate selectivity of landings at age
  FemSelLandAtAge <- NA; MalSelLandAtAge <- NA
  if (!is.na(EstLandSelAtAge[1,1])) {
    FemSelLandAtAge <- EstLandSelAtAge[,1]
    MalSelLandAtAge <- EstLandSelAtAge[,2]
  } else if (!is.na(Land_sel_A50[1])) {
    FemSelLandAtAge <- 1/(1+exp(-log(19) * (Ages - Land_sel_A50[1]) / (Land_sel_A95[1] - Land_sel_A50[1])))
    MalSelLandAtAge <- 1/(1+exp(-log(19) * (Ages - Land_sel_A50[2]) / (Land_sel_A95[2] - Land_sel_A50[2])))
  } else if (!is.na(FemGearSelAtAge[1]) & !is.na(FemRetProbAtAge[1])) {
    FemSelLandAtAge <- FemGearSelAtAge * FemRetProbAtAge
    MalSelLandAtAge <- MalGearSelAtAge * MalRetProbAtAge
  }
  # note, it is assumed that selectivity of landings is inputted in some way

  # Calculate selectivity of discards at age
  if (DiscMort==0) {
    FemSelDiscAtAge <- 0
    MalSelDiscAtAge <- 0
  }

  if (DiscMort>0) {

    # if gear selectivity unknown, but selectivity of landings and and probability of retention inputted.
    # probability of retention may sometimes be assumed from MLL.
    if (is.na(FemGearSelAtAge[1])) {
      if (!is.na(FemSelLandAtAge[1]) & !is.na(FemRetProbAtAge[1])) {
        FemGearSelAtAge <- FemSelLandAtAge / FemRetProbAtAge
        MalGearSelAtAge <- MalSelLandAtAge / MalRetProbAtAge
      }
    }

    # if probability of retention is unknown, but selectivity of landings and gear selectivity are known.
    # probability of retention may sometimes be assumed from MLL.
    if (is.na(FemRetProbAtAge[1])) {
      if (!is.na(FemSelLandAtAge[1]) & !is.na(FemGearSelAtAge[1])) {
        FemRetProbAtAge <- FemSelLandAtAge / FemGearSelAtAge
        MalRetProbAtAge <- MalSelLandAtAge / MalGearSelAtAge
      }
    }

    if (is.na(FemGearSelAtAge[1])) cat("Problem, FemGearSelAtAge[1]=NA. Need to input more info.",'\n')
    if (is.na(FemRetProbAtAge[1])) cat("Problem, FemRetProbAtAge[1]=NA. Need to input more info.",'\n')
    if (is.na(FemSelLandAtAge[1])) cat("Problem, FemSelLandAtAge[1]=NA. Need to input more info.",'\n')

    FemSelDiscAtAge <- FemGearSelAtAge * (1 - FemRetProbAtAge)
    MalSelDiscAtAge <- MalGearSelAtAge * (1 - MalRetProbAtAge)
  }

  results = list(FemGearSelAtAge = FemGearSelAtAge,
                 MalGearSelAtAge = MalGearSelAtAge,
                 FemRetProbAtAge = FemRetProbAtAge,
                 MalRetProbAtAge = MalRetProbAtAge,
                 FemSelLandAtAge = FemSelLandAtAge,
                 MalSelLandAtAge = MalSelLandAtAge,
                 FemSelDiscAtAge = FemSelDiscAtAge,
                 MalSelDiscAtAge = MalSelDiscAtAge)

  return(results)

}

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
#' @param FinalSex_A50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_A95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param EstSexRatioAtAge NA  # sex ratio at age (from age 0) inputted as vector
#' @param EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' @param mat_A50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_A95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtAge vector of proportion mature at age (set to NA if using age at maturity parameters)
#' @param Gear_sel_A50 logistic parameter for gear selectivity curve
#' @param Gear_sel_A95 logistic parameter for gear selectivity curve
#' @param EstGearSelAtAge vector for gear selectivity at age (set to NA if using age at gear selectivity parameters)
#' @param Land_sel_A50 logistic parameter for selectivity of landings curve
#' @param Land_sel_A95 logistic parameter for selectivity of landings curve
#' @param EstLandSelAtAge vector for for selectivity of landings at age (set to NA if using age at selectivity of landings parameters)
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
#' (UnfishFemSpBiomAtAge, UnfishMalSpBiomAtAge), fished female survival at age (FishFemSurvAtAge, FishMalSurvAtAge),
#' fished female and male biomass mature biomass at age (FishFemSpBiomAtAge, FishMalSpBiomAtAge), female and male
#' retained catch at age in numbers (FemCatchAtAgeNum, MalCatchAtAgeNum), female and male
#' release catch at age in numbers (FemRelCatchAtAgeNum, MalRelCatchAtAgeNum), female and male retained catch at age in biomass
#' (FemCatchAtAge, MalCatchAtAge), equilibrium recruitment for either Beverton-Holt or Ricker relationship (Eq_Rec), equilibrium catch (Eq_Catch), equilibrium
#' female and male and spawning biomass (Eq_FemSpBiom, Eq_MalSpBiom), equilibrium relative female,
#' male and combined sex spawning biomass (Eq_SPR, Eq_MalRelSpBiom, Eq_CombSexRelSpBiom),
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
#' FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtAge <- NA  # sex ratio at age (from age 0) inputted as vector
#' EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, ~0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' Gear_sel_A50 <- NA # females, males - Logistic age selectivity relationship parameters
#' Gear_sel_A95 <- NA # females, males - Logistic age selectivity relationship parameters
#' EstGearSelAtAge <- data.frame(FemGearSelAtAge=NA, MalGearSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' Land_sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' Land_sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' EstLandSelAtAge <- data.frame(FemSelLandAtAge=NA, MalSelLandAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
#' EstRetenAtAge <- data.frame(FemRetProbAtAge=NA, MalRetProbAtAge=NA) # fish retention at age (from age 0), inputted as values in data frame
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 0.2 # natural mortality  (year-1)
#' FMort <- 0.4 # estimate of fishing mortality, e.g. from catch curve analysis
#' Res=CalcYPRAndSPRForFMort_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
#'                             ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge,
#'                             EggFertParam, mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95,
#'                             EstGearSelAtAge, Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax,
#'                             ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, FMort)
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
#' FinalSex_A50 <- 35 # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- 60 # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtAge <- NA  # sex ratio at age (from age 0) inputted as vector
#' EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_A50 <- c(20, 20) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(30, 30) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' Gear_sel_A50 <- c(20, 20) # females, males - Logistic age selectivity relationship parameters
#' Gear_sel_A95 <- c(30, 30) # females, males - Logistic age selectivity relationship parameters
#' EstGearSelAtAge <- data.frame(FemGearSelAtAge=NA, MalGearSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' Land_sel_A50 <- c(25, 25) # females, males - Logistic age selectivity relationship parameters
#' Land_sel_A95 <- c(35, 35) # females, males - Logistic age selectivity relationship parameters
#' EstLandSelAtAge <- data.frame(FemSelLandAtAge=NA, MalSelLandAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
#' EstRetenAtAge <- data.frame(FemRetProbAtAge=NA, MalRetProbAtAge=NA) # fish retention at age (from age 0), inputted as values in data frame
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 0.07 # natural mortality  (year-1)
#' FMort <- 0.07 # estimate of fishing mortality, e.g. from catch curve analysis
#' Res=CalcYPRAndSPRForFMort_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
#'                             ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge,
#'                             EggFertParam, mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95,
#'                             EstGearSelAtAge, Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax,
#'                             ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, FMort)
#' @export
CalcYPRAndSPRForFMort_AB <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                  lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                                  ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge,
                                  EggFertParam, mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95,
                                  EstGearSelAtAge, Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax,
                                  ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, FMort) {

  # Determine ages to be modelled, specified on maximum age and timestep
  Ages <- seq(0,MaxModelAge,TimeStep)

  # number of model time steps
  nTimeSteps <- length(Ages)

  # calculate length at age - von Bertalanffy growth curve
  res=CalcLengthAtAge(Linf, vbK, tzero, EstLenAtAge, Ages)
  FemLenAtAge = res$FemLenAtAge; MalLenAtAge = res$MalLenAtAge

  # calculate weight (g) at age
  res=CalcWeightAtAge(lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, FemLenAtAge, MalLenAtAge)
  FemWtAtAge <- res$FemWtAtAge; MalWtAtAge <- res$MalWtAtAge

  # for sex-changing species, calculate proportion female at age
  PropFemAtAge=CalcPropFemAtAge(FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, ReprodPattern, Ages)

  # calculate proportion mature at age
  res=CalcPropMatureAtAge(mat_A50, mat_A95, EstMatAtAge, Ages)
  FemPropMatAtAge = res$FemPropMatAtAge; MalPropMatAtAge = res$MalPropMatAtAge

  # Calculate gear selectivity and retention at age, and selectivity of landings and discards
  res=CalcSelectivityAndRetentionAtAge(EstGearSelAtAge, EstRetenAtAge, Ages, Gear_sel_A50, Gear_sel_A95,
                                       Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_A50, ret_A95, ret_Pmax, DiscMort)
  FemGearSelAtAge = res$FemGearSelAtAge; MalGearSelAtAge = res$MalGearSelAtAge
  FemRetProbAtAge = res$FemRetProbAtAge; MalRetProbAtAge = res$MalRetProbAtAge
  FemSelLandAtAge = res$FemSelLandAtAge; MalSelLandAtAge = res$MalSelLandAtAge
  FemSelDiscAtAge = res$FemSelDiscAtAge; MalSelDiscAtAge = res$MalSelDiscAtAge

  # Calculate fishing mortality at age associated with discarding, landings, and combined.
  res=CalcFishingMortalityAtAge(FMort, DiscMort, FemSelDiscAtAge, MalSelDiscAtAge, FemSelLandAtAge, MalSelLandAtAge)
  FemDiscFAtAge = res$FemDiscFAtAge; MalDiscFAtAge = res$MalDiscFAtAge
  FemLandFAtAge = res$FemLandFAtAge; MalLandFAtAge = res$MalLandFAtAge
  FemFAtAge = res$FemFAtAge; MalFAtAge = res$MalFAtAge

  # calculate female and male total mortality at age
  FemZAtAge <- FemFAtAge + NatMort
  MalZAtAge <- MalFAtAge + NatMort

  # calculate relative unfished female and male survival at age
  res=CalcSurvivalAtAge(nTimeSteps, InitRatioFem, NatMort, FMort, FemZAtAge, MalZAtAge, ReprodPattern, PropFemAtAge)
  UnfishFemSurvAtAge = res$UnfishFemSurvAtAge; UnfishMalSurvAtAge = res$UnfishMalSurvAtAge
  FishFemSurvAtAge = res$FishFemSurvAtAge; FishMalSurvAtAge = res$FishMalSurvAtAge


  # calculate ratio of mature males to mature females, by number and egg fertilisation rates
  res=CalcPopnSexRatioAndFertRate_AB(UnfishMalSurvAtAge, FishMalSurvAtAge, UnfishFemSurvAtAge, FishFemSurvAtAge,
                                          MalPropMatAtAge, FemPropMatAtAge, EggFertParam)
  UnfishMalToFemProp = res$UnfishMalToFemProp
  FishMalToFemProp = res$FishMalToFemProp
  MalDeplRatio = res$MalDeplRatio
  Eq_FertRate = res$Eq_FertRate

  # calculate female and male unfished spawning biomass at age
  UnfishFemSpBiomAtAge <- UnfishFemSurvAtAge * FemPropMatAtAge * (((FemWtAtAge * 1000) ^ ReprodScale) / 1000)
  UnfishMalSpBiomAtAge <- UnfishMalSurvAtAge * MalPropMatAtAge * (((MalWtAtAge * 1000) ^ ReprodScale) / 1000)

  # Total unfished spawning biomass in kg
  UnfishFemSpBiom <- sum(UnfishFemSpBiomAtAge)
  UnfishMalSpBiom <- sum(UnfishMalSpBiomAtAge)
  UnfishCombSexSpBiom <- UnfishFemSpBiom + UnfishMalSpBiom

  # calculate female and male mature biomass at age for fished population
  FishFemSpBiomAtAge <- FishFemSurvAtAge * FemPropMatAtAge * (((FemWtAtAge * 1000) ^ ReprodScale) / 1000)
  FishMalSpBiomAtAge <- FishMalSurvAtAge * MalPropMatAtAge * (((MalWtAtAge * 1000) ^ ReprodScale) / 1000)

  # calculate retained female and male catches at age (in numbers) - Baranov catch equation
  FemCatchAtAgeNum <- FishFemSurvAtAge * (FemLandFAtAge/FemZAtAge) * (1 - exp(-(FemZAtAge * TimeStep)))
  MalCatchAtAgeNum <- FishMalSurvAtAge * (MalLandFAtAge/MalZAtAge) * (1 - exp(-(MalZAtAge * TimeStep)))

  # calculate female and male discarded catch at age (in numbers) - Baranov catch equation
  # (including fish that either survive or die after release)
  FemRelCatchAtAgeNum <- FishFemSurvAtAge * ((FMort * FemSelDiscAtAge)/FemZAtAge) * (1 - exp(-(FemZAtAge * TimeStep)))
  MalRelCatchAtAgeNum <- FishMalSurvAtAge * ((FMort * MalSelDiscAtAge)/MalZAtAge) * (1 - exp(-(MalZAtAge * TimeStep)))

  # calculate retained female and male catch at age (in biomass)
  FemCatchAtAge <- FemCatchAtAgeNum * FemWtAtAge
  MalCatchAtAge <- MalCatchAtAgeNum * MalWtAtAge

  # calculate yield per recruit in kg
  YPR <- max(0,sum(FemCatchAtAge) + sum(MalCatchAtAge))

  # calculate female and male spawning biomass per recruit in kg
  FishFemSpBiom <- sum(FishFemSpBiomAtAge)
  FishMalSpBiom <- sum(FishMalSpBiomAtAge)
  FishCombSexSpBiom <- FishFemSpBiom + FishMalSpBiom

  # calculate spawning potential ratio (SPR)
  Fem_SPR <- max(0,FishFemSpBiom / UnfishFemSpBiom)
  Mal_SPR <- max(0,FishMalSpBiom / UnfishMalSpBiom)
  CombSex_SPR <- max(0,((FishFemSpBiom + FishMalSpBiom)  / (UnfishFemSpBiom + UnfishMalSpBiom)))

  if (ReprodPattern == 1) { # gonochoristic species
    UnfishSpBiom = UnfishFemSpBiom
    FishSpBiom = FishFemSpBiom
  }
  if (ReprodPattern > 1) { # hermaphroditic species
    if (!is.na(EggFertParam)) {
      # if accounting for egg fertilisation rate, then recruitment is calculated based on
      # number of eggs fertilised (assumed to be proportional to female spawning biomass measure)
      UnfishSpBiom = UnfishFemSpBiom
      FishSpBiom = FishFemSpBiom
    } else {
      # based recruitment on combined sex spawning biomass for hermaphroditic species, as default
      # if fertilisation rate not being considered
      UnfishSpBiom = UnfishCombSexSpBiom
      FishSpBiom = FishCombSexSpBiom
    }
  }

  # calculate equilibrium recruitment from stock recruitment relationship
  res = CalcEquilibriumRecruitment(SRrel_Type, Steepness, FishSpBiom, UnfishSpBiom, Eq_FertRate)
  Eq_Rec = res$Eq_Rec
  Eq_Rec_AllEggFert = res$Eq_Rec_AllEggFert # assuming all eggs fertilised, i.e. Eq_FertRate=1

  # calculate equilibrium catch
  Eq_Catch <- max(0,Eq_Rec * YPR)

  # calculate equilibrium female and male spawning biomass
  Eq_FemSpBiom <- Eq_Rec * FishFemSpBiom
  Eq_FemSpBiom_AllEggFert <- Eq_Rec_AllEggFert * FishFemSpBiom
  Eq_MalSpBiom <- Eq_Rec * FishMalSpBiom

  # calculate equilibrium model SPR
  Eq_FemRelSpBiom <- max(0,Eq_FemSpBiom / UnfishFemSpBiom)
  Eq_FemRelSpBiom_AllEggFert <- max(0,Eq_FemSpBiom_AllEggFert / UnfishFemSpBiom)
  Eq_MalRelSpBiom <- max(0,Eq_MalSpBiom / UnfishMalSpBiom)
  Eq_CombSexRelSpBiom <- max(0,(Eq_FemSpBiom + Eq_MalSpBiom) / (UnfishFemSpBiom + UnfishMalSpBiom))

  ModelDiag = data.frame(Ages = Ages,
                         UnfishFemSurvAtAge = UnfishFemSurvAtAge,
                         UnfishMalSurvAtAge = UnfishMalSurvAtAge,
                         UnfishFemSpBiomAtAge = UnfishFemSpBiomAtAge,
                         UnfishMalSpBiomAtAge = UnfishMalSpBiomAtAge,
                         FishFemSurvAtAge = FishFemSurvAtAge,
                         FishMalSurvAtAge = FishMalSurvAtAge,
                         FishFemSpBiomAtAge = FishFemSpBiomAtAge,
                         FishMalSpBiomAtAge = FishMalSpBiomAtAge,
                         FemCatchAtAgeNum = FemCatchAtAgeNum,
                         MalCatchAtAgeNum = MalCatchAtAgeNum,
                         FemCatchAtAge = FemCatchAtAge,
                         MalCatchAtAge = MalCatchAtAge,
                         FemRelCatchAtAgeNum = FemRelCatchAtAgeNum,
                         MalRelCatchAtAgeNum = MalRelCatchAtAgeNum,
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

  Results = list(YPR = YPR,
                 Fem_SPR = Fem_SPR,
                 Mal_SPR = Mal_SPR,
                 CombSex_SPR = CombSex_SPR,
                 Eq_Rec = Eq_Rec,
                 Eq_Rec_AllEggFert = Eq_Rec_AllEggFert,
                 Eq_Catch = Eq_Catch,
                 Eq_FemRelSpBiom = Eq_FemRelSpBiom,
                 Eq_FemRelSpBiom_AllEggFert = Eq_FemRelSpBiom_AllEggFert,
                 Eq_MalRelSpBiom = Eq_MalRelSpBiom,
                 Eq_CombSexRelSpBiom = Eq_CombSexRelSpBiom,
                 UnfishMalToFemProp = UnfishMalToFemProp,
                 FishMalToFemProp = FishMalToFemProp,
                 MalDeplRatio = MalDeplRatio,
                 Eq_FertRate = Eq_FertRate,
                 ModelDiag = ModelDiag)

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
CalcMeanSizeAtAge <- function(GrowthCurveType, MaxModelAge, Ages, GrowthParams, RefnceAges) {

  MeanSizeAtAge <- data.frame(matrix(nrow = 2, ncol = length(Ages)))
  colnames(MeanSizeAtAge) <- round(Ages,2)
  MeanSizeAtAge <- as.matrix(MeanSizeAtAge)

  # calculate length at age - von Bertalanffy growth curve
  if (GrowthCurveType == 1) {
    if (is.vector(GrowthParams)) { # growth params inputted as single sex
      FemLenAtAge <- GrowthParams[1] * (1 - exp(-GrowthParams[2] * (Ages - GrowthParams[3])))
      MalLenAtAge <- GrowthParams[1] * (1 - exp(-GrowthParams[2] * (Ages - GrowthParams[3])))
    }
    if (is.data.frame(GrowthParams)) { # growth params inputted as 2 sexes
      FemLenAtAge <- GrowthParams[1,1] * (1 - exp(-GrowthParams[1,2] * (Ages - GrowthParams[1,3])))
      MalLenAtAge <- GrowthParams[2,1] * (1 - exp(-GrowthParams[2,2] * (Ages - GrowthParams[2,3])))
    }
  }

  # Schnute growth curve
  if (GrowthCurveType == 2) {
    FemLenAtAge <- rep(0,length(Ages))
    MalLenAtAge <- rep(0,length(Ages))
    for (i in 1:2) {

      if (is.vector(RefnceAges)) { # growth params inputted as single sex
        t1=RefnceAges[1]; t2=RefnceAges[2]
        y1=GrowthParams[1]; y2=GrowthParams[2]
        a=GrowthParams[3]; b=GrowthParams[4]
      }
      if (is.data.frame(RefnceAges)) { # growth params inputted as 2 sexes
        t1=RefnceAges[i,1]; t2=RefnceAges[i,2]
        y1=GrowthParams[i,1]; y2=GrowthParams[i,2]
        a=GrowthParams[i,3]; b=GrowthParams[i,4]
      }

      k=0
      for (t in Ages) {
        k=k+1
        if (i==1) {
          FemLenAtAge[k] = SchnuteGrowthfunction(t, t1, t2, y1, y2, a, b)
        }
        if (i==2) {
          MalLenAtAge[k] = SchnuteGrowthfunction(t, t1, t2, y1, y2, a, b)
        }
      } # t
    } # i
  } # if

  # set negative lengths to zero
  FemLenAtAge[which(FemLenAtAge<0)] = 0
  MalLenAtAge[which(MalLenAtAge<0)] = 0
  MeanSizeAtAge[1,] <- FemLenAtAge
  MeanSizeAtAge[2,] <- MalLenAtAge

  return(MeanSizeAtAge)

}


#' Calculates weight at length for each sex for length-based per recruit analysis
#'
#' This function calculates weight at length for each sex, based on a power curve,
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
CalcWeightAtLength <- function(lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen, midpt) {


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

#' Calculate proportion of fish at length that are females
#'
#' This function calculates proportion of fish at length that are females,
#' based on reproductive pattern (gonochoristic or hermaphroditic), and
#' sex change parameters
#'
#' @keywords internal
#'
#' @param FinalSex_L50 logistic sex change parameter
#' @param FinalSex_L95 logistic sex change parameter
#' @param EstSexRatioAtLen input vector sex change at length
#' @param ReprodPattern 1=gonochoristic, 2=protogynous, 3=protandrous
#' @param midpt length class mid-points
#'
#' @return PropFemAtLen
CalcPropFemAtLength <- function(FinalSex_L50, FinalSex_L95, EstSexRatioAtLen, ReprodPattern, midpt) {


  # for sex-changing species, calculate proportion female at age
  if (ReprodPattern == 1) { # gonochoristic (separate sexes)
    PropFemAtLen = NA
  }
  if (ReprodPattern == 2) { # protogynous hermaphroditism (female to male sex change)
    if (is.na(EstSexRatioAtLen[1])) {
      PropFemAtLen = 1 - (1 / (1 + exp(-log(19) * (midpt - FinalSex_L50) / (FinalSex_L95 - FinalSex_L50))))
    }
    if (!is.na(EstSexRatioAtLen[1])) {
      PropFemAtLen <- EstSexRatioAtLen
    }
  }
  if (ReprodPattern == 3) { # protandrous hermaphroditism (male to female sex change)
    if (is.na(EstSexRatioAtLen[1])) {
      PropFemAtLen = 1 / (1 + exp(-log(19) * (midpt - FinalSex_L50) / (FinalSex_L95 - FinalSex_L50)))
    }
    if (!is.na(EstSexRatioAtAge[1])) {
      PropFemAtLen <- EstSexRatioAtLen
    }
  }

  results = PropFemAtLen

  return(results)

}


#' Calculate proportion of fish at length that are mature
#'
#' This function calculates proportion of fish at length that are mature,
#' based on input maturity parameters, or from set input maturity schedules
#'
#' @keywords internal
#'
#' @param mat_L50 logistic maturity parameter
#' @param mat_L95 logistic maturity parameter
#' @param EstMatAtLen input maturity schedules (could be set to NA)
#' @param midpt length class mid-points
#'
#' @return FemPropMatAtLen
#' @return MalPropMatAtLen
CalcPropMatureAtLength <- function(mat_L50, mat_L95, EstMatAtLen, midpt) {

  if (is.na(EstMatAtLen[1,1])) {
    FemPropMatAtLen <- 1 / (1 + exp(-log(19) * (midpt - mat_L50[1]) / (mat_L95[1] - mat_L50[1])))
    MalPropMatAtLen <- 1 / (1 + exp(-log(19) * (midpt - mat_L50[2]) / (mat_L95[2] - mat_L50[2])))
  }
  if (!is.na(EstMatAtLen[1,1])) {
    FemPropMatAtLen <- EstMatAtLen[,1]
    MalPropMatAtLen <- EstMatAtLen[,2]
  }

  results = list(FemPropMatAtLen=FemPropMatAtLen,
                 MalPropMatAtLen=MalPropMatAtLen)

  return(results)

}


#' Calculate relative equilibrium recruitment
#'
#' This function calculates expected relative equilibrium recruitment from stock-recruitment
#' relationship, given steepness, and fished and unfished spawning biomass per recruitment
#'
#' @keywords internal
#'
#' @param SRrel_Type stock recruitment relationship 1=Beverton-Holt, 2=Ricker
#' @param Steepness steepness of the stock-recruitmnet relationship
#' @param FishSpBiom fished spawning biomass per recruit (females for gonochoristic species, combined sex for hermaphroditic species)
#' @param UnfishSpBiom fished spawning biomass per recruit (females for gonochoristic species, combined sex for hermaphroditic species)
#' @param Eq_FertRate (NA or from 0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#'
#' @return Eq_Rec
CalcEquilibriumRecruitment <- function(SRrel_Type, Steepness, FishSpBiom, UnfishSpBiom, Eq_FertRate) {


  # Calculate equilbrium recruitment
  if (SRrel_Type == 1) { # Beverton-Holt
    Eq_Rec_AllEggFert = ((4 * Steepness * FishSpBiom) - (1 - Steepness) * UnfishSpBiom) /
      (5 * (Steepness - 0.2) * FishSpBiom)
  }
  if (SRrel_Type == 2) { # Ricker
    # Ricker
    Eq_Rec_AllEggFert <- (UnfishSpBiom * (1- ((4 * log(UnfishSpBiom/FishSpBiom)) /
                                           (5 * log(5 * Steepness))))) / FishSpBiom
  }
  #check = (Unfish_FemBiomPerRecSpawnSeas-SR_alpha)/(SR_beta*Unfish_FemBiomPerRecSpawnSeas)

  # Alternative parameterisation for Beverton-Holt relationship
  # if (SRrel_Type == 1) { # Beverton-Holt
  #   BH_SRRa <- (UnfishSpBiom / 1.0) * ((1-Steepness) / (4*Steepness))
  #   BH_SRRb <- (Steepness - 0.2) / (0.8 * Steepness * 1.0)
  # Check that the specified initial recruitment can be recovered,
  # given the S_R parameters and unfished female biomass
  # (Check <- UnfishSpBiom/(BH_SRRa + BH_SRRb * UnfishSpBiom))

  # calculate equilibrium recruitment
  # BH_Eq_Rec <- (FishSpBiom-BH_SRRa) / (BH_SRRb*FishSpBiom)

  # Equations from Norm Hall - Beverton and Holt and Ricker relationship,
  # when both parameterised using steepness
  # }
  # if (SRrel_Type == 2) { # Ricker
  #   BH_SRRa = NA
  #   BH_SRRb = NA
  #   # BH_Eq_Rec = NA
  # }

  # account for egg fertilisation rate, based on current male-female sex ratio for mature
  # fish (in numbers), relative to unfished level/
  Eq_Rec = Eq_Rec_AllEggFert * Eq_FertRate

  results = list(Eq_Rec=Eq_Rec,
                 Eq_Rec_AllEggFert=Eq_Rec_AllEggFert)

  return(results)

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
GetLTMInputsForPerRecruitAnalysis <- function(GrowthCurveType, TimeStep, MaxAge, GrowthParams, RefnceAges, midpt) {

  nLenCl = length(midpt)
  MeanEndingLength <- data.frame(matrix(nrow = 2, ncol = nLenCl))
  colnames(MeanEndingLength) <- midpt
  MeanEndingLength <- as.matrix(MeanEndingLength)
  TimestepGrowthSizeInc <- MeanEndingLength

  # get key inputs for length transition matrices
  if (GrowthCurveType == 1) { # von Bertalanffy
    for (i in 1:2) {
      if (is.vector(GrowthParams)) { # params inputted as separate sex
        MeanEndingLength[i,] = midpt + (GrowthParams[1] - midpt) * (1 - exp(-GrowthParams[2]*TimeStep))
      }
      if (is.data.frame(GrowthParams)) { # params inputted as 2 sexes
        MeanEndingLength[i,] = midpt + (GrowthParams[i,1] - midpt) * (1 - exp(-GrowthParams[i,2]*TimeStep))
      }
      TimestepGrowthSizeInc[i,] = MeanEndingLength[i,] - midpt # amount of annual growth with respect to initial length
    }
  }
  if (GrowthCurveType == 2) { # Schnute
    for (i in 1:2) {
      if (is.vector(GrowthParams)) { # params inputted as separate sex
        y1=GrowthParams[1]
        y2=GrowthParams[2]
        a=GrowthParams[3]
        b=GrowthParams[4]
        t1=RefnceAges[1]
        t2=RefnceAges[2]
      }
      if (is.data.frame(GrowthParams)) { # params inputted as 2 sexes
        y1=GrowthParams[i,1]
        y2=GrowthParams[i,2]
        a=GrowthParams[i,3]
        b=GrowthParams[i,4]
        t1=RefnceAges[i,1]
        t2=RefnceAges[i,2]
      }

      GrowthParamsForSex = c(y1, y2, a, b)
      RefnceAgesForSex = c(t1, t2)
      MeanEndingLength[i,] = CalcLengthAfterGrowthForTimetep(GrowthCurveType=2, TimeStep, GrowthParamsForSex, RefnceAgesForSex, midpt, MaxAge)
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
#' prediction intervals for mean length (CatchLen.60qntl, CatchLen.95qntl). Values are calculated
#' for females, males and combined sexes.
#'
#' @keywords internal
#'
#' @param midpt length class mid-points
#' @param FemCatchNumLen female per recruit catch numbers at length
#' @param MalCatchNumLen male per recruit catch numbers at length
#' @param CombSexCatchNumLen
#'
#' @return PerRecRandFishLen
#' @return PerRecMeanCatchLen
#' @return PerRecLenFreq
#' @return MeanCatchLen.95ci
#' @return CatchLen.60qntl
#' @return CatchLen.95qntl
#' @return PerRecRandFishWt
#' @return PerRecMeanCatchWt
#' @return MeanCatchWt.95ci
#' @return CatchWt.60qntl
#' @return CatchWt.95qntl
#'
CalcPerRecruitFishLenAndFishWtStats<- function(lbnd, midpt, ubnd, FemCatchNumLen, MalCatchNumLen, CombSexCatchNumLen,
                                               lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen) {


  # get approx. 95 and 60% confidence intervals for mean catch length
  LenInterv = (ubnd[1] - lbnd[1]) / 2 # randomising fish lengths, within each length class
  FemProbs = FemCatchNumLen/sum(FemCatchNumLen)
  FemLenFreq = as.vector(rmultinom(1,10000,FemProbs))
  FemFishLen = rep(midpt, FemLenFreq) + runif(10000,-LenInterv, LenInterv)
  MalProbs = MalCatchNumLen/sum(MalCatchNumLen)
  MalLenFreq = as.vector(rmultinom(1,10000,MalProbs))
  MalFishLen = rep(midpt, MalLenFreq) + runif(10000,-LenInterv, LenInterv)
  CombSexProbs = CombSexCatchNumLen/sum(CombSexCatchNumLen)
  CombSexLenFreq = as.vector(rmultinom(1,10000,CombSexProbs))
  CombSexFishLen = rep(midpt, CombSexLenFreq) + runif(10000,-LenInterv, LenInterv)
  res=CalcWeightAtLength(lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen, FemFishLen)
  FemFishWt = res$FemWtAtLen
  res=CalcWeightAtLength(lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen, MalFishLen)
  MalFishWt = res$MalWtAtLen
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
  FemMeanCatchLen = mean(FemFishLen)
  MalMeanCatchLen = mean(MalFishLen)
  CombSexMeanCatchLen = mean(CombSexFishLen)
  PerRecMeanCatchLen = data.frame(FemMeanCatchLen=FemMeanCatchLen,
                                  MalMeanCatchLen=MalMeanCatchLen,
                                  CombSexMeanCatchLen=CombSexMeanCatchLen)

  # calculate mean catch weight
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

  # 60 and 95 percent quantile ranges for catch length
  FemCatchLen.95qntl = c(quantile(FemFishLen, 0.025),quantile(FemFishLen, 0.975))
  MalCatchLen.95qntl = c(quantile(MalFishLen, 0.025),quantile(MalFishLen, 0.975))
  CombSexCatchLen.95qntl = c(quantile(CombSexFishLen, 0.025),quantile(CombSexFishLen, 0.975))
  CatchLen.95qntl = data.frame(FemCatchLen.95qntl=FemCatchLen.95qntl,
                                 MalCatchLen.95qntl=MalCatchLen.95qntl,
                                 CombSexCatchLen.95qntl=CombSexCatchLen.95qntl)

  FemCatchLen.60qntl = c(quantile(FemFishLen, 0.2),quantile(FemFishLen, 0.8))
  MalCatchLen.60qntl = c(quantile(MalFishLen, 0.2),quantile(MalFishLen, 0.8))
  CombSexCatchLen.60qntl = c(quantile(CombSexFishLen, 0.2),quantile(CombSexFishLen, 0.8))
  CatchLen.60qntl = data.frame(FemCatchLen.60qntl=FemCatchLen.60qntl,
                                 MalCatchLen.60qntl=MalCatchLen.60qntl,
                                 CombSexCatchLen.60qntl=CombSexCatchLen.60qntl)


  # 60 and 95 percent quantile ranges for catch weight
  FemCatchWt.95qntl = c(quantile(FemFishWt, 0.025),quantile(FemFishWt, 0.975))
  MalCatchWt.95qntl = c(quantile(MalFishWt, 0.025),quantile(MalFishWt, 0.975))
  CombSexCatchWt.95qntl = c(quantile(CombSexFishWt, 0.025),quantile(CombSexFishWt, 0.975))
  CatchWt.95qntl = data.frame(FemCatchWt.95qntl=FemCatchWt.95qntl,
                                 MalCatchWt.95qntl=MalCatchWt.95qntl,
                                 CombSexCatchWt.95qntl=CombSexCatchWt.95qntl)

  FemCatchWt.60qntl = c(quantile(FemFishWt, 0.2),quantile(FemFishWt, 0.8))
  MalCatchWt.60qntl = c(quantile(MalFishWt, 0.2),quantile(MalFishWt, 0.8))
  CombSexCatchWt.60qntl = c(quantile(CombSexFishWt, 0.2),quantile(CombSexFishWt, 0.8))
  CatchWt.60qntl = data.frame(FemCatchWt.60qntl=FemCatchWt.60qntl,
                                 MalCatchWt.60qntl=MalCatchWt.60qntl,
                                 CombSexCatchWt.60qntl=CombSexCatchWt.60qntl)

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
                 CatchLen.60qntl=CatchLen.60qntl,
                 CatchLen.95qntl=CatchLen.95qntl,
                 PerRecRandFishWt=PerRecRandFishWt,
                 PerRecMeanCatchWt=PerRecMeanCatchWt,
                 MeanCatchWt.95ci=MeanCatchWt.95ci,
                 CatchWt.60qntl=CatchWt.60qntl,
                 CatchWt.95qntl=CatchWt.95qntl,
                 RandFishLen = RandFishLen,
                 RandFishWt = RandFishWt)

  return(Results)

}


#' Calculate gear selectivity and retention at length
#'
#' Calculate gear selectivity and retention at the mid point for eac length class,
#' given input parameters
#'
#' @keywords internal
#'
#' @param Ages ages for analysis calculated from max age and model timestep
#' @param sel_L50 logistic gear selectivity curve parameter
#' @param sel_L95 logistic gear selectivity curve parameter
#' @param EstGearSelAtLen gear selectivity inputted as vector
#' @param ret_L50 logistic gear retention curve parameter
#' @param ret_L95 logistic gear retention curve parameter
#' @param ret_Pmax logistic gear retention curve parameter
#' @param EstRetenAtLen lretention curve inputted as vector
#'
#' @return FemGearSelAtAge, MalGearSelAtAge, FemRetProbAtAge, MalRetProbAtAge
CalcSelectivityAndRetentionAtLen <- function(midpt, sel_L50, sel_L95, EstGearSelAtLen,
                                             ret_L50, ret_L95, ret_Pmax, EstRetenAtLen) {

  # Calculate gear selectivity at length
  FemGearSelAtLen <- NA; MalGearSelAtLen <- NA
  if (!is.na(EstGearSelAtLen[1,1])) {
    FemGearSelAtLen <- EstGearSelAtLen[,1]
    MalGearSelAtLen <- EstGearSelAtLen[,2]
  } else if (!is.na(sel_L50[1])){
    FemGearSelAtLen <- 1/(1+exp(-log(19) * (midpt - sel_L50[1]) / (sel_L95[1] - sel_L50[1])))
    MalGearSelAtLen <- 1/(1+exp(-log(19) * (midpt - sel_L50[2]) / (sel_L95[2] - sel_L50[2])))
  }

  # Calculate fish retention at length
  FemRetProbAtLen <- NA; MalRetProbAtLen <- NA
  if (!is.na(EstRetenAtLen[1,1])) {
    FemRetProbAtLen <- EstRetenAtLen[,1]
    MalRetProbAtLen <- EstRetenAtLen[,2]
  } else if (!is.na(ret_L50[1])) {
    FemRetProbAtLen <- ret_Pmax[1] / (1+exp(-log(19) * (midpt - ret_L50[1]) / (ret_L95[1] - ret_L50[1])))
    MalRetProbAtLen <- ret_Pmax[2] / (1+exp(-log(19) * (midpt - ret_L50[2]) / (ret_L95[2] - ret_L50[2])))
  }

  # Calculate selectivity of landings at age
  FemSelLandAtLen <- FemGearSelAtLen * FemRetProbAtLen
  MalSelLandAtLen <- MalGearSelAtLen * MalRetProbAtLen

  # Calculate selectivity of discards at age
  FemSelDiscAtLen <- FemGearSelAtLen * (1 - FemRetProbAtLen)
  MalSelDiscAtLen <- MalGearSelAtLen * (1 - MalRetProbAtLen)

  results = list(FemGearSelAtLen = FemGearSelAtLen,
                 MalGearSelAtLen = MalGearSelAtLen,
                 FemRetProbAtLen = FemRetProbAtLen,
                 MalRetProbAtLen = MalRetProbAtLen,
                 FemSelLandAtLen = FemSelLandAtLen,
                 MalSelLandAtLen = MalSelLandAtLen,
                 FemSelDiscAtLen = FemSelDiscAtLen,
                 MalSelDiscAtLen = MalSelDiscAtLen)

  return(results)

}

#' Calculate fishing mortality at length (discards and landings)
#'
#' Calculate fishing mortality at the mid point of each length class, associated with discarding, landings, and
#' combined.
#'
#' @keywords internal
#'
#' @param FMort fully-selected fishing mortality
#' @param DiscMort proportion of discarded fish that die
#' @param FemSelDiscAtLen selectivity of discarding at length for females
#' @param MalSelDiscAtLen selectivity of discarding at length for males
#' @param FemSelLandAtLen selectivity of landings at length for females
#' @param MalSelLandAtLen selectivity of landings at length for males
#'
#' @return FemDiscFAtLen, MalDiscFAtLen, FemLandFAtLen, MalLandFAtLen, FemFAtLen, MalFAtLen
CalcFishingMortalityAtLen <- function(FMort, DiscMort, FemSelDiscAtLen, MalSelDiscAtLen,
                                      FemSelLandAtLen, MalSelLandAtLen) {

  # calculate female and male fishing mortality at length associated with discarding of undersize fish
  FemDiscFAtLen <- FMort * DiscMort * FemSelDiscAtLen
  MalDiscFAtLen <- FMort * DiscMort * MalSelDiscAtLen

  # calculate female and male fishing mortality at length associated with landings
  FemLandFAtLen <- FMort * FemSelLandAtLen
  MalLandFAtLen <- FMort * MalSelLandAtLen

  # calculate (total) female and male fishing mortality at length
  if (DiscMort >= 0.001) {
    FemFAtLen <- FMort * (FemSelLandAtLen + (DiscMort * FemSelDiscAtLen))
    MalFAtLen <- FMort * (MalSelLandAtLen + (DiscMort * MalSelDiscAtLen))
  } else {
    FemFAtLen <- FMort * FemSelLandAtLen
    MalFAtLen <- FMort * MalSelLandAtLen
  }

  results = list(FemDiscFAtLen = FemDiscFAtLen,
                 MalDiscFAtLen = MalDiscFAtLen,
                 FemLandFAtLen = FemLandFAtLen,
                 MalLandFAtLen = MalLandFAtLen,
                 FemFAtLen = FemFAtLen,
                 MalFAtLen = MalFAtLen)

  return(results)

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
#' @param FinalSex_L50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param EstSexRatioAtLen NA  # sex ratio at length, inputted as vector
#' @param EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' @param mat_L50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_L95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtLen vector of proportion mature at length (set to NA if using age at maturity parameters)
#' @param sel_L50 logistic parameter for gear selectivity curve
#' @param sel_L95 logistic parameter for gear selectivity curve
#' @param EstGearSelAtLen gear selectivity curve inputted as vector
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_L50 logistic parameter for fish retention curve
#' @param ret_L95 logistic parameter for fish retention curve
#' @param EstRetenAtLen retention curve inputted as vector
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
#' equilibrium recruitment for either Beverton-Holt or Ricker relationship (Eq_Rec),
#' equilibrium catch (Eq_Catch), female and male and spawning biomass (Eq_FemSpBiom, Eq_MalSpBiom), equilibrium relative
#' female, male and combined sex spawning biomass (Eq_FemRelSpBiom, Eq_MalRelSpBiom, Eq_CombSexRelSpBiom), mean size at
#' age of females and males (MeanSizeAtAge), female and male weight at length (FemWtAtLen, MalWtAtLen), female and male proportion
#' mature at length (FemPropMatAtLen, MalPropMatAtLen), female and male gear selectivity at length (FemGearSelAtLen, MalGearSelAtLen),
#' female and male retention at length (FemRetProbAtLen, MalRetProbAtLen), female and male selectivity of landings at length
#' (FemSelLandAtLen, MalSelLandAtLen), female and male selectivity of discards at length (FemSelDiscAtLen, MalSelDiscAtLen),
#' female and male discard mortality at length (FemDiscFAtLen, MalDiscFAtLen), female and male mortality associated with landings
#' at length (FemLandFAtLen, MalLandFAtLen), female and male fishing mortality at length (FemFAtLen, MalFAtLen),
#' female and male total mortality at length (FemZAtLen, MalZAtLen)
#'
#' @examples
#' # non-hermaphroditic species
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
#'                          EstMalWtAtLen=NA) # weight at age, inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at recruitment age/length
#' FinalSex_L50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_L95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtLen <- NA  # sex ratio at length inputted as vector
#' EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, ~0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_L50 <- c(250, 250) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_L95 <- c(300, 300) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
#'                           EstMalMatAtLen=NA) # maturity at length, inputted as values in data frame
#' sel_L50 <- c(250, 250) # females, males - Logistic length selectivity relationship parameters
#' sel_L95 <- c(300, 300) # females, males - Logistic length selectivity relationship parameters
#' EstGearSelAtLen <- data.frame(EstFemGearSelAtLen=NA,
#'                               EstMalGearSelAtLen=NA)
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_L50 <- c(250, 250) # females, males - Logistic fish retention at length parameters
#' ret_L95 <- c(300, 300) # females, males - Logistic fish retention at length parameters
#' EstRetenAtLen <- data.frame(EsFemtRetenAtLen=NA,
#'                             EstMalRetenAtLen=NA)
#' DiscMort <- 0.25 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 4.22 / MaxModelAge # natural mortality  (year-1)
#' FMort = 0.2
#' Res=CalcYPRAndSPRForFMort_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
#'                              RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
#'                              ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
#'                              EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
#'                              ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, FMort)
#'
#' # Hermaphroditic species
#' MaxModelAge <- 100 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
#' MaxLen = 1400
#' LenInc = 20
#' lbnd = seq(0,MaxLen - LenInc, LenInc)
#' ubnd = lbnd + LenInc
#' midpt = lbnd + (LenInc/2)
#' nLenCl = length(midpt)
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf <- c(682, 982) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.14, 0.08) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0.06, -0.48) # years - von Bertalanffy growth model parameters - Females, males
#' GrowthParams = data.frame(Linf=Linf, vbK=vbK, tzero=tzero)
#' RefnceAges = NA
#' # GrowthParams c(Linf, vbK, tzero) single sex von Bertalanffy, or data.frame(Linf=Linf, vbK=vbK, tzero=tzero),
#' #' both sexes von Bertalanffy, or c(y1, y2, a, b) single sex Schnute, or data.frame(y1=y1, y2=y2, a=a, b=b), both sexes Schnute
#' CVSizeAtAge = c(0.05,0.05)
#' lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- -11.017 # for log-log relationship
#' lenwt_b <- 3.041 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 2 # 1=power, 2=log-log relationship
#' EstWtAtLen <- data.frame(EstFemWtAtLen=NA,
#'                          EstMalWtAtLen=NA) # weight at length, inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 2 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 1.0 # Ratio of females to males at recruitment age
#' FinalSex_L50 <- 821 # Logistic sex change relationship parameters (inflection point)
#' FinalSex_L95 <- 930 # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtLen <- NA  # sex ratio at length inputted as vector
#' EggFertParam <- 0.5 # (NA or from ~0.2-1) NA = no effect, ~0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_L50 <- c(653, 0) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_L95 <- c(950, 1) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
#'                           EstMalMatAtLen=NA) # maturity at length, inputted as values in data frame
#' sel_L50 <- c(500, 500) # females, males - Logistic length selectivity relationship parameters
#' sel_L95 <- c(600, 600) # females, males - Logistic length selectivity relationship parameters
#' EstGearSelAtLen <- data.frame(EstFemGearSelAtLen=NA,
#'                               EstMalGearSelAtLen=NA)
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_L50 <- c(499, 499) # females, males - Logistic fish retention at length parameters
#' ret_L95 <- c(500, 500) # females, males - Logistic fish retention at length parameters
#' EstRetenAtLen <- data.frame(EsFemtRetenAtLen=NA,
#'                             EstMalRetenAtLen=NA)
#' DiscMort <- 0.25 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort <- 4.22 / 70 # natural mortality  (year-1)
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' PlotOpt <- 0 # 0=all plots, 1=len at-age, 2=wt at length, 3=fem mat/sel/ret at length, 4=mal mat/sel/ret at length,
#' # 5=fem F at length, 6=mal F at length, 7=fem rel surv, 8=mal rel surv, 9=fem biom at age, 10=fem biom at age,
#' # 11=ypr/eq catch, 12=fem SPR/Brel, 13=mal SPR/Brel, 14=eq recruit
#' Current_F = 0.05
#' FMort = 0.05
#' Res=CalcYPRAndSPRForFMort_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
#'                              RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
#'                              ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
#'                              EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
#'                              ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, FMort)
#' @export
CalcYPRAndSPRForFMort_LB<- function(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                    RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                                    ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                                    EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                                    ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, FMort) {

  if (FMort < 0.0001) FMort = 0.0001 # set negligible fishing mortality for unfished stock,
  # to allow calculations associated with expected catch size distribtion, i.e. from a survey

  # number of model time steps
  # nTimeSteps <- 1 + (MaxModelAge / TimeStep)
  Ages <- seq(TimeStep,MaxModelAge, TimeStep)
  nTimeSteps <- length(Ages)

  # calculate mean size at age
  MeanSizeAtAge = CalcMeanSizeAtAge(GrowthCurveType, MaxModelAge, Ages, GrowthParams, RefnceAges)

  # calculate weight at length
  res=CalcWeightAtLength(lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen, midpt)
  FemWtAtLen = res$FemWtAtLen; MalWtAtLen = res$MalWtAtLen

  # calculate proportion female at length
  PropFemAtLen=CalcPropFemAtLength(FinalSex_L50, FinalSex_L95, EstSexRatioAtLen, ReprodPattern, midpt)

  # calculate proportion mature at length
  res=CalcPropMatureAtLength(mat_L50, mat_L95, EstMatAtLen, midpt)
  FemPropMatAtLen=res$FemPropMatAtLen; MalPropMatAtLen=res$MalPropMatAtLen

  # calculate selectivity and retention at length, including for landings and discards
  # using retention-function approach
  res=CalcSelectivityAndRetentionAtLen(midpt, sel_L50, sel_L95, EstGearSelAtLen,
                                       ret_L50, ret_L95, ret_Pmax, EstRetenAtLen)
  FemGearSelAtLen = res$FemGearSelAtLen; MalGearSelAtLen = res$MalGearSelAtLen
  FemRetProbAtLen = res$FemRetProbAtLen; MalRetProbAtLen = res$MalRetProbAtLen
  FemSelLandAtLen = res$FemSelLandAtLen; MalSelLandAtLen = res$MalSelLandAtLen
  FemSelDiscAtLen = res$FemSelDiscAtLen; MalSelDiscAtLen = res$MalSelDiscAtLen

  # Calculate fishing mortality at length associated with discarding, landings, and combined.
  res=CalcFishingMortalityAtLen(FMort, DiscMort, FemSelDiscAtLen, MalSelDiscAtLen, FemSelLandAtLen, MalSelLandAtLen)
  FemDiscFAtLen = res$FemDiscFAtLen; MalDiscFAtLen = res$MalDiscFAtLen
  FemLandFAtLen = res$FemLandFAtLen; MalLandFAtLen = res$MalLandFAtLen
  FemFAtLen = res$FemFAtLen; MalFAtLen = res$MalFAtLen

  # calculate female and male total mortality at age
  FemZAtLen <- FemFAtLen + NatMort
  MalZAtLen <- MalFAtLen + NatMort

  # calculate initial size distribution of female and male recruits
  RecLenDist = CalcSizeDistOfRecruits(MeanSizeAtAge, CVSizeAtAge, lbnd, ubnd, midpt, nLenCl)

  # get required inputs to calculate length transition matrices
  MaxAge=MaxModelAge
  results = GetLTMInputsForPerRecruitAnalysis(GrowthCurveType, TimeStep, MaxAge, GrowthParams, RefnceAges, midpt)
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

  Unfish_FemNPerRecAtAge=PopnRes$Unfish_FemNPerRecAtAge; Fish_FemNPerRecAtAge=PopnRes$Fish_FemNPerRecAtAge
  Unfish_MalNPerRecAtAge=PopnRes$Unfish_MalNPerRecAtAge; Fish_MalNPerRecAtAge=PopnRes$Fish_MalNPerRecAtAge
  Unfish_FemBiomPerRecAtAge = PopnRes$Unfish_FemBiomPerRecAtAge; Fish_FemBiomPerRecAtAge = PopnRes$Fish_FemBiomPerRecAtAge
  Unfish_MalBiomPerRecAtAge = PopnRes$Unfish_MalBiomPerRecAtAge; Fish_MalBiomPerRecAtAge = PopnRes$Fish_MalBiomPerRecAtAge
  Unfish_FemBiomAtAge = PopnRes$Unfish_FemBiomAtAge; Fish_FemBiomAtAge = PopnRes$Fish_FemBiomAtAge
  Unfish_MalBiomAtAge = PopnRes$Unfish_MalBiomAtAge; Fish_MalBiomAtAge = PopnRes$Fish_MalBiomAtAge
  Unfish_FemNPerRecLen=PopnRes$Unfish_FemNPerRec; Fish_FemNPerRecLen=PopnRes$Fish_FemNPerRec
  Unfish_MalNPerRecLen=PopnRes$Unfish_MalNPerRec; Fish_MalNPerRecLen=PopnRes$Fish_MalNPerRec

  # calculate ratio of mature males to mature females, by number and egg fertilisation rates
  res=CalcPopnSexRatioAndFertRate_LB(Unfish_MalNPerRecLen, Fish_MalNPerRecLen, Unfish_FemNPerRecLen, Fish_FemNPerRecLen,
                                     MalPropMatAtLen, FemPropMatAtLen, EggFertParam)
  UnfishMalToFemProp = res$UnfishMalToFemProp
  FishMalToFemProp = res$FishMalToFemProp
  MalDeplRatio = res$MalDeplRatio
  Eq_FertRate = res$Eq_FertRate

  # Unfished spawning biomass in kg
  UnfishFemSpBiom <- sum(Unfish_FemBiomAtAge)
  UnfishMalSpBiom <- sum(Unfish_MalBiomAtAge)
  UnfishCombSexSpBiom <- UnfishFemSpBiom + UnfishMalSpBiom

  # Fished spawning biomass in kg
  FishFemSpBiom <- sum(Fish_FemBiomAtAge)
  FishMalSpBiom <- sum(Fish_MalBiomAtAge)
  FishCombSexSpBiom <- FishFemSpBiom + FishMalSpBiom

  # calculate catch yield (biomass)
  FemCatchBiomLen = (FemLandFAtLen/FemZAtLen) * (1-exp(-FemZAtLen)) * Fish_FemNPerRecLen * FemWtAtLen
  MalCatchBiomLen = (MalLandFAtLen/MalZAtLen) * (1-exp(-MalZAtLen)) * Fish_MalNPerRecLen * MalWtAtLen
  CombSexCatchBiomLen = FemCatchBiomLen + MalCatchBiomLen

  # calculate catch yield (numbers)
  FemCatchNumLen = (FemLandFAtLen/FemZAtLen) * (1-exp(-FemZAtLen)) * Fish_FemNPerRecLen
  MalCatchNumLen = (MalLandFAtLen/MalZAtLen) * (1-exp(-MalZAtLen)) * Fish_MalNPerRecLen
  CombSexCatchNumLen = FemCatchNumLen + MalCatchNumLen

  # calculate yield per recruit in kg
  YPR <- max(0,sum(FemCatchBiomLen) + sum(MalCatchBiomLen))

  # calculate spawning potential ratio (SPR)
  Fem_SPR <- max(0,FishFemSpBiom / UnfishFemSpBiom)
  Mal_SPR <- max(0,FishMalSpBiom / UnfishMalSpBiom)
  CombSex_SPR <- max(0,((FishFemSpBiom + FishMalSpBiom)  / (UnfishFemSpBiom + UnfishMalSpBiom)))

  if (ReprodPattern == 1) { # gonochoristic species
    UnfishSpBiom = UnfishFemSpBiom
    FishSpBiom = FishFemSpBiom
  }
  if (ReprodPattern > 1) { # hermaphroditic species
    UnfishSpBiom = UnfishCombSexSpBiom
    FishSpBiom = FishCombSexSpBiom
  }

  # Calculate equilibrium recruitment from stock recruitment relationship
  res = CalcEquilibriumRecruitment(SRrel_Type, Steepness, FishSpBiom, UnfishSpBiom, Eq_FertRate)
  Eq_Rec = res$Eq_Rec
  Eq_Rec_AllEggFert = res$Eq_Rec_AllEggFert # assuming all eggs fertilised, i.e. Eq_FertRate=1

  # calculate equilibrium female and male spawning biomass
  Eq_FemSpBiom <- Eq_Rec * FishFemSpBiom
  Eq_FemSpBiom_AllEggFert <- Eq_Rec_AllEggFert * FishFemSpBiom
  Eq_MalSpBiom <- Eq_Rec * FishMalSpBiom

  # calculate equilibrium catch
  Eq_Catch <- max(0,Eq_Rec * YPR)

  # calculate equilibrium model SPR
  Eq_FemRelSpBiom <- max(0,Eq_FemSpBiom / UnfishFemSpBiom)
  Eq_FemRelSpBiom_AllEggFert <- max(0,Eq_FemSpBiom_AllEggFert / UnfishFemSpBiom)
  Eq_MalRelSpBiom <- max(0,Eq_MalSpBiom / UnfishMalSpBiom)
  Eq_CombSexRelSpBiom <- max(0,(Eq_FemSpBiom + Eq_MalSpBiom) / (UnfishFemSpBiom + UnfishMalSpBiom))

  ModelDiag <- list(FemWtAtLen = FemWtAtLen,
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
                          MalZAtLen=MalZAtLen,
                          RecLenDist=RecLenDist,
                          PropFemAtLen=PropFemAtLen,
                          Unfish_FemNPerRecLen=Unfish_FemNPerRecLen,
                          Unfish_MalNPerRecLen=Unfish_MalNPerRecLen,
                          Fish_FemNPerRecLen=Fish_FemNPerRecLen,
                          Fish_MalNPerRecLen=Fish_MalNPerRecLen,
                          FemCatchBiomLen=FemCatchBiomLen,
                          MalCatchBiomLen=MalCatchBiomLen,
                          CombSexCatchBiomLen=CombSexCatchBiomLen,
                          FemCatchNumLen=FemCatchNumLen,
                          MalCatchNumLen=MalCatchNumLen,
                          CombSexCatchNumLen=CombSexCatchNumLen,
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
                          MeanSizeAtAge = MeanSizeAtAge)

  Results = list(YPR = YPR,
                 Fem_SPR = Fem_SPR,
                 Mal_SPR = Mal_SPR,
                 CombSex_SPR = CombSex_SPR,
                 Eq_Rec = Eq_Rec,
                 Eq_Rec_AllEggFert = Eq_Rec_AllEggFert,
                 Eq_Catch = Eq_Catch,
                 Eq_FemRelSpBiom = Eq_FemRelSpBiom,
                 Eq_FemRelSpBiom_AllEggFert = Eq_FemRelSpBiom_AllEggFert,
                 Eq_MalRelSpBiom = Eq_MalRelSpBiom,
                 Eq_CombSexRelSpBiom = Eq_CombSexRelSpBiom,
                 UnfishMalToFemProp = UnfishMalToFemProp,
                 FishMalToFemProp = FishMalToFemProp,
                 MalDeplRatio = MalDeplRatio,
                 Eq_FertRate = Eq_FertRate,
                 ModelDiag = ModelDiag)

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
#' @param FinalSex_A50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_A95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param EstSexRatioAtAge NA  # sex ratio at age (from age 0) inputted as vector
#' @param EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' @param mat_A50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_A95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtAge vector of proportion mature at age (set to NA if using age at maturity parameters)
#' @param Gear_sel_A50 logistic parameter for gear selectivity curve
#' @param Gear_sel_A95 logistic parameter for gear selectivity curve
#' @param EstGearSelAtAge vector for gear selectivity at age (set to NA if using age at gear selectivity parameters)
#' @param Land_sel_A50 logistic parameter for selectivity of landings curve
#' @param Land_sel_A95 logistic parameter for selectivity of landings curve
#' @param EstLandSelAtAge vector for for selectivity of landings at age (set to NA if using age at selectivity of landings parameters)
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
#' unfished female and male mature biomass at age (UnfishFemSpBiomAtAge, UnfishMalSpBiomAtAge), fished female and
#' male spawning potential ratio (Fem_SPR, Mal_SPR), fished female and male per recruit survival at age
#' (FishFemSurvAtAge, FishMalSurvAtAge), fished female and male mature biomass at age (FishFemSpBiomAtAge,
#' FishMalSpBiomAtAge), female and male per recruit retained catches at age in numbers (FemCatchAtAgeNum, MalCatchAtAgeNum),
#' female and male per recruit released catches at age in numbers, (FemRelCatchAtAgeNum, MalRelCatchAtAgeNum),
#' female and male per recruit retained catches at age in biomass (FemCatchAtAge, MalCatchAtAge), and for the extended model, equilibrium recruitment (Eq_Rec),
#' equilibrium catch (Eq_Catch), equilibrium female and male spawning biomass (Eq_FemSpBiom, Eq_MalSpBiom) and
#' relative female, male and combined sex spawning biomass (Eq_FemRelSpBiom, Eq_MalRelSpBiom, Eq_CombSexRelSpBiom),
#' female and male gear selectivity at age (FemGearSelAtAge, MalGearSelAtAge), female and male retention at age probabilities,
#' (FemRetProbAtAge, MalRetProbAtAge), selectivity of female and male fish landings (FemSelLandAtAge, MalSelLandAtAge),
#' female and male selectivity of discards (FemSelDiscAtAge, MalSelDiscAtAge), female and male fishing mortality
#' associated with discarding (FemDiscFAtAge, MalDiscFAtAge), female and male fishing mortality associated with
#' fish landings (FemLandFAtAge, MalLandFAtAge), female and male total fishing mortality at age (FemFAtAge, MalFAtAge),
#' female and male total mortality at age (FemZAtAge, MalZAtAge)range of fishing mortality values applied for which per quantities are calculated
#' (FishMort = seq(0,2,0.01)), equilibrium recruitment vs FMort (Eq_Rec), equilibrium catch vs FMort (Eq_Catch),
#' equilibrium female spawning biomass vs FMort (Eq_FemSpBiom), equilibrium spawning potential ratio vs FMort
#' (Eq_SPR), maximum yield per recruit (YPRmax), maximum equilibrium catch (maxeqCatch), fishing mortality
#' associated with YPRmax (Fmax), fishing mortality associated with maxeqCatch (F_MSY or FmaxeqCatch), biomass target at 1.2maxeqCatch,
#' (BMSY_Targ), biomass threshold at maxeqCatch (BMSY_Thresh), biomass threshold at 0.5maxeqCatch (BMSY_Lim),
#' YPR vs FishMort (YPRResults), Eq_Catch vs FishMort (Eq_CatchResults), Eq_FemRelSpBiom vs FishMort (Fem_SPRResults),
#' Eq_MalRelSpBiom vs FishMort (Mal_SPRResults), Eq_FemRelSpBiom vs FishMort (Eq_FemRelSpBiomResults),
#' Eq_MalRelSpBiom vs FishMort (Eq_MalRelSpBiomResults), Eq_CombSexRelSpBiom vs FishMort (CombSex_SPRResults),
#' Eq_Rec vs FishMort (Eq_RecResults)
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
#' FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtAge <- NA  # sex ratio at age (from age 0) inputted as vector
#' EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' Gear_sel_A50 <- NA # females, males - Logistic age selectivity relationship parameters
#' Gear_sel_A95 <- NA # females, males - Logistic age selectivity relationship parameters
#' EstGearSelAtAge <- data.frame(FemGearSelAtAge=NA, MalGearSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' Land_sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' Land_sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' EstLandSelAtAge <- data.frame(FemSelLandAtAge=NA, MalSelLandAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
#' EstRetenAtAge <- data.frame(EstFemRetenAtAge=rep(1,MaxModelAge+1), EstMalRetenAtAge=rep(1,MaxModelAge+1)) # fish retention at age (from age 0), inputted as values in data frame
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 0.2 # natural mortality  (year-1)
#' Current_F <- 0.07 # estimate of fishing mortality, e.g. from catch curve analysis
#' Res = GetPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge, lenwt_a, ln_lenwt_a,
#'                               lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern, InitRatioFem,
#'                               FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam, mat_A50, mat_A95,
#'                               EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge, Land_sel_A50, Land_sel_A95,
#'                               EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort,
#'                               Steepness, SRrel_Type, NatMort, Current_F)
#' # Example 2: hermaphroditic species
#' InitRecruit <- 1 # Initial recruitment
#' MaxModelAge <- 100 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
#' Linf <- c(1000, 1000) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.1, 0.1) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0, 0) # years - von Bertalanffy growth model parameters - Females, males
#' EstLenAtAge <- data.frame(EstFemLenAtAge=NA, EstMalLenAtAge=NA) # length at age (from age 0), inputted as values in data frame
#' lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- -11.0 # for log-log relationship
#' lenwt_b <- 3.0 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 2 # 1=power, 2=log-log relationship
#' EstWtAtAge <- data.frame(EstFemWtAtAge=NA, EstMalWtAtAge=NA) # weight at age (from age 0), inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 2 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 1 # Ratio of females to males at age zero
#' FinalSex_A50 <- 35 # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- 60 # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtAge <- NA  # sex ratio at age (from age 0) inputted as vector
#' EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_A50 <- c(20, 20) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(30, 30) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' Gear_sel_A50 <- c(20, 20) # females, males - Logistic age selectivity relationship parameters
#' Gear_sel_A95 <- c(30, 30) # females, males - Logistic age selectivity relationship parameters
#' EstGearSelAtAge <- data.frame(FemGearSelAtAge=NA, MalGearSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' Land_sel_A50 <- c(25, 25) # females, males - Logistic age selectivity relationship parameters
#' Land_sel_A95 <- c(35, 35) # females, males - Logistic age selectivity relationship parameters
#' EstLandSelAtAge <- data.frame(FemSelLandAtAge=NA, MalSelLandAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
#' EstRetenAtAge <- data.frame(EstFemRetenAtAge=rep(1,MaxModelAge+1), EstMalRetenAtAge=rep(1,MaxModelAge+1)) # fish retention at age (from age 0), inputted as values in data frame
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 0.07 # natural mortality  (year-1)
#' Current_F <- 0.07 # estimate of fishing mortality, e.g. from catch curve analysis
#' Res = GetPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge, lenwt_a, ln_lenwt_a,
#'                               lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern, InitRatioFem,
#'                               FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam, mat_A50, mat_A95,
#'                               EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge, Land_sel_A50, Land_sel_A95,
#'                               EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort,
#'                               Steepness, SRrel_Type, NatMort, Current_F)
#' @export
GetPerRecruitResults_AB <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge, lenwt_a, ln_lenwt_a,
                                    lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern, InitRatioFem,
                                    FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam, mat_A50, mat_A95,
                                    EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge, Land_sel_A50, Land_sel_A95,
                                    EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort,
                                    Steepness, SRrel_Type, NatMort, Current_F) {

  FishMort <- seq(0,2,0.01)
  nFVals <- length(FishMort) # fishing mortality
  YPRResults <- rep(0,nFVals)
  Fem_SPRResults <- rep(0,nFVals)
  Mal_SPRResults <- rep(0,nFVals)
  CombSex_SPRResults <- rep(0,nFVals)

  Eq_RecResults <- rep(0,nFVals)
  Eq_Rec_AllEggFertResults = rep(0,nFVals)
  Eq_CatchResults <- rep(0,nFVals)
  Eq_FemRelSpBiomResults <- rep(0,nFVals)
  Eq_FemRelSpBiom_AllEggFert_Results <- rep(0,nFVals)

  Eq_MalRelSpBiomResults <- rep(0,nFVals)
  Eq_CombSexRelSpBiomResults <- rep(0,nFVals)
  Eq_FertRateResults <- rep(0,nFVals)
  UnfishMalToFemPropResults <- rep(0,nFVals)
  FishMalToFemPropResults <- rep(0,nFVals)
  Eq_MalDeplRatioResults <- rep(0,nFVals)
  Eq_FertRateResults <- rep(0,nFVals)

  for (k in 1:nFVals) {
    FMort = FishMort[k]
    Res = CalcYPRAndSPRForFMort_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                   lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                                   ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge,
                                   EggFertParam, mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95,
                                   EstGearSelAtAge, Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax,
                                   ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, FMort)
    # per recruit results
    YPRResults[k] <- Res$YPR
    Fem_SPRResults[k] = Res$Fem_SPR
    Mal_SPRResults[k] = Res$Mal_SPR
    CombSex_SPRResults[k] = Res$CombSex_SPR

    # extended model results
    Eq_RecResults[k] = Res$Eq_Rec
    Eq_Rec_AllEggFertResults[k] = Res$Eq_Rec_AllEggFert
    Eq_CatchResults[k] = Res$Eq_Catch
    Eq_FemRelSpBiomResults[k] = Res$Eq_FemRelSpBiom
    Eq_FemRelSpBiom_AllEggFert_Results[k] = Res$Eq_FemRelSpBiom_AllEggFert
    Eq_MalRelSpBiomResults[k] = Res$Eq_MalRelSpBiom
    Eq_CombSexRelSpBiomResults[k] = Res$Eq_CombSexRelSpBiom
    Eq_MalDeplRatioResults[k] <- Res$MalDeplRatio
    Eq_FertRateResults[k] <- Res$Eq_FertRate
    UnfishMalToFemPropResults[k] <- Res$UnfishMalToFemProp
    FishMalToFemPropResults[k] <- Res$FishMalToFemProp

  }
  YPRmax <- max(YPRResults) # maximum yield per recruit
  maxeqCatch <- max(Eq_CatchResults) # maximum equilbrium catch
  Fmax <- FishMort[which(YPRResults==YPRmax)] # fishing mortality at ypr maximum
  FmaxeqCatch <- FishMort[which(Eq_CatchResults==maxeqCatch)] # fishing mortality at maxeqCatch

  # Calculate biological reference points
  # assuming B_thresh=BMSY, B_lim=0.5BMSY and B_targ=1.2BMSY
  # and F_thresh, F_lim and F_targ are values of F resulting in B_thresh, B_lim and B_targ, respectively

  # gonochoristic species or considering sex ratio effect on egg fertilisation rate
  if (ReprodPattern == 1 | !is.na(EggFertParam)) {
    F_MSY=FmaxeqCatch
    x=which(FishMort==F_MSY)
    BMSY_Thresh=Eq_FemRelSpBiomResults[x]
    BMSY_Lim=0.5*BMSY_Thresh
    x=which.min(abs(Eq_FemRelSpBiomResults - BMSY_Lim))
    F_Lim = FishMort[x]
    BMSY_Targ=1.2*BMSY_Thresh
    x=which.min(abs(Eq_FemRelSpBiomResults - BMSY_Targ))
    F_Targ = FishMort[x]
  } else { # hermaphroditic species and not considering sex ratio effect on egg fertilisation rate
    F_MSY=FmaxeqCatch
    x=which(FishMort==F_MSY)
    BMSY_Thresh=Eq_CombSexRelSpBiomResults[x]
    BMSY_Lim=0.5*BMSY_Thresh
    x=which.min(abs(Eq_CombSexRelSpBiomResults - BMSY_Lim))
    F_Lim = FishMort[x]
    BMSY_Targ=1.2*BMSY_Thresh
    x=which.min(abs(Eq_CombSexRelSpBiomResults - BMSY_Targ))
    F_Targ = FishMort[x]
  }

  # get results for current F
  FMort = Current_F
  Res2 = CalcYPRAndSPRForFMort_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                  lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                                  ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge,
                                  EggFertParam, mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95,
                                  EstGearSelAtAge, Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax,
                                  ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, FMort)

  Diagnotistcs = data.frame(Ages = Res2$ModelDiag$Ages,
                         FemLenAtAge = Res2$ModelDiag$FemLenAtAge,
                         MalLenAtAge = Res2$ModelDiag$MalLenAtAge,
                         FemWtAtAge = Res2$ModelDiag$FemWtAtAge,
                         MalWtAtAge = Res2$ModelDiag$MalWtAtAge,
                         PropFemAtAge = Res2$ModelDiag$PropFemAtAge,
                         FemPropMatAtAge = Res2$ModelDiag$FemPropMatAtAge,
                         MalPropMatAtAge = Res2$ModelDiag$MalPropMatAtAge,
                         UnfishFemSurvAtAge = Res2$ModelDiag$UnfishFemSurvAtAge,
                         UnfishMalSurvAtAge = Res2$ModelDiag$UnfishMalSurvAtAge,
                         UnfishFemSpBiomAtAge = Res2$ModelDiag$UnfishFemSpBiomAtAge,
                         UnfishMalSpBiomAtAge = Res2$ModelDiag$UnfishMalSpBiomAtAge,
                         FishFemSurvAtAge = Res2$ModelDiag$FishFemSurvAtAge,
                         FishMalSurvAtAge = Res2$ModelDiag$FishMalSurvAtAge,
                         FishFemSpBiomAtAge = Res2$ModelDiag$FishFemSpBiomAtAge,
                         FishMalSpBiomAtAge = Res2$ModelDiag$FishMalSpBiomAtAge,
                         FemCatchAtAgeNum = Res2$ModelDiag$FemCatchAtAgeNum,
                         MalCatchAtAgeNum = Res2$ModelDiag$MalCatchAtAgeNum,
                         FemRelCatchAtAgeNum = Res2$ModelDiag$FemRelCatchAtAgeNum,
                         MalRelCatchAtAgeNum = Res2$ModelDiag$MalRelCatchAtAgeNum,
                         FemCatchAtAge = Res2$ModelDiag$FemCatchAtAge,
                         MalCatchAtAge = Res2$ModelDiag$MalCatchAtAge,
                         FemGearSelAtAge = Res2$ModelDiag$FemGearSelAtAge,
                         MalGearSelAtAge = Res2$ModelDiag$MalGearSelAtAge,
                         FemRetProbAtAge = Res2$ModelDiag$FemRetProbAtAge,
                         MalRetProbAtAge = Res2$ModelDiag$MalRetProbAtAge,
                         FemSelLandAtAge = Res2$ModelDiag$FemSelLandAtAge,
                         MalSelLandAtAge = Res2$ModelDiag$MalSelLandAtAge,
                         FemSelDiscAtAge = Res2$ModelDiag$FemSelDiscAtAge,
                         MalSelDiscAtAge = Res2$ModelDiag$MalSelDiscAtAge,
                         FemDiscFAtAge = Res2$ModelDiag$FemDiscFAtAge,
                         MalDiscFAtAge = Res2$ModelDiag$MalDiscFAtAge,
                         FemLandFAtAge = Res2$ModelDiag$FemLandFAtAge,
                         MalLandFAtAge = Res2$ModelDiag$MalLandFAtAge,
                         FemFAtAge = Res2$ModelDiag$FemFAtAge,
                         MalFAtAge = Res2$ModelDiag$MalFAtAge,
                         FemZAtAge = Res2$ModelDiag$FemZAtAge,
                         MalZAtAge = Res2$ModelDiag$MalZAtAge)

  ModelDiag = Diagnotistcs

  Results = list(YPR = Res2$YPR,
                 Fem_SPR = Res2$Fem_SPR,
                 Mal_SPR = Res2$Mal_SPR,
                 CombSex_SPR = Res2$CombSex_SPR,
                 Eq_Rec = Res2$Eq_Rec,
                 MalDeplRatio = Res2$MalDeplRatio,
                 Eq_FertRate = Res2$Eq_FertRate,
                 Eq_Rec_AllEggFert = Res2$Eq_Rec_AllEggFert,
                 Eq_Catch = Res2$Eq_Catch,
                 Eq_FemRelSpBiom = Res2$Eq_FemRelSpBiom,
                 Eq_FemRelSpBiom_AllEggFert = Res2$Eq_FemRelSpBiom_AllEggFert,
                 Eq_MalRelSpBiom = Res2$Eq_MalRelSpBiom,
                 Eq_CombSexRelSpBiom = Res2$Eq_CombSexRelSpBiom,
                 YPRmax = YPRmax,
                 FmaxeqCatch = FmaxeqCatch,
                 Fmax = Fmax,
                 F_MSY = F_MSY,
                 BMSY = BMSY_Thresh,
                 F_Targ = F_Targ,
                 F_Thresh = F_MSY,
                 F_Lim = F_Lim,
                 BMSY_Targ = BMSY_Targ,
                 BMSY_Thresh = BMSY_Thresh,
                 BMSY_Lim = BMSY_Lim,
                 FishMort = FishMort,
                 YPRResults = YPRResults,
                 Eq_CatchResults = Eq_CatchResults,
                 Fem_SPRResults = Fem_SPRResults,
                 Mal_SPRResults = Mal_SPRResults,
                 CombSex_SPRResults = CombSex_SPRResults,
                 Eq_FemRelSpBiomResults = Eq_FemRelSpBiomResults,
                 Eq_FemRelSpBiom_AllEggFert_Results = Eq_FemRelSpBiom_AllEggFert_Results,
                 Eq_MalRelSpBiomResults = Eq_MalRelSpBiomResults,
                 Eq_CombSexRelSpBiomResults = Eq_CombSexRelSpBiomResults,
                 Eq_RecResults = Eq_RecResults,
                 Eq_Rec_AllEggFertResults = Eq_Rec_AllEggFertResults,
                 UnfishMalToFemPropResults = UnfishMalToFemPropResults,
                 FishMalToFemPropResults = FishMalToFemPropResults,
                 Eq_MalDeplRatioResults = Eq_MalDeplRatioResults,
                 Eq_FertRateResults = Eq_FertRateResults,
                 ModelDiag = ModelDiag)

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
#' @param FinalSex_L50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param EstSexRatioAtLen NA  # sex ratio at length, inputted as vector
#' @param EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' @param mat_L50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at length)
#' @param mat_L95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at length)
#' @param EstMatAtLen vector of proportion mature at length (set to NA if using length at maturity parameters)
#' @param sel_L50 logistic parameter for gear selectivity curve
#' @param sel_L95 logistic parameter for gear selectivity curve
#' @param EstGearSelAtLen gear selectivity curve inputted as vector
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_L50 logistic parameter for fish retention curve
#' @param ret_L95 logistic parameter for fish retention curve
#' @param EstRetenAtLen retention curve inputted as vector
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param Current_F estimated current fishing mortality
#' @param Output_Opt 1 = standard, 2=include expected length and weight distributions at F (slower)
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
#' (FemCatchBiom, MalCatchBiom), equilibrium recruitment for either Beverton-Holt or Ricker relationship (Eq_Rec),
#' equilibrium catch (Eq_Catch), female and male and spawning biomass (Eq_FemSpBiom, Eq_MalSpBiom), equilibrium relative
#' female, male and combined sex spawning biomass (Eq_FemRelSpBiom, Eq_MalRelSpBiom, Eq_CombSexRelSpBiom),
#' female and male gear selectivity at length (FemGearSelAtLen, MalGearSelAtLen), female and male retention at length
#' (FemRetProbAtLen, MalRetProbAtLen), female and male selectivity of landings at length (FemSelLandAtLen, MalSelLandAtLen),
#' female and male selectivity of discards at length (FemSelDiscAtLen, MalSelDiscAtLen),
#' female and male discard mortality at length (FemDiscFAtLen, MalDiscFAtLen), female and male mortality associated with landings
#' at length (FemLandFAtLen, MalLandFAtLen), female and male fishing mortality at length (FemFAtLen, MalFAtLen),
#' female and male total mortality at length (FemZAtLen, MalZAtLen), fishing mortality values evaluated (FishMort = seq(0,2,0.01)),
#' maximum yield per recruit (YPRmax), fishing mortality associated with maxeqCatch (F_MSY or FmaxeqCatch),
#' biomass target at 1.2maxeqCatch, (BMSY_Targ), biomass threshold at maxeqCatch (BMSY_Thresh), biomass threshold at 0.5maxeqCatch (BMSY_Lim),
#' YPR vs FishMort (YPRResults), Eq_Catch vs FishMort (Eq_CatchResults), Eq_FemRelSpBiom vs FishMort (Fem_SPRResults),
#' Eq_MalRelSpBiom vs FishMort (Mal_SPRResults), Eq_FemRelSpBiom vs FishMort (Eq_FemRelSpBiomResults),
#' Eq_MalRelSpBiom vs FishMort (Eq_MalRelSpBiomResults), Eq_CombSexRelSpBiom vs FishMort (Eq_CombSexRelSpBiomResults),
#' Eq_Rec vs FishMort (Eq_RecResults)
#'
#' @examples
#' # Non-hermaphroditic species
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
#'                          EstMalWtAtLen=NA) # weight at length, inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at recruitment age/length
#' FinalSex_L50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_L95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtLen <- NA  # sex ratio at length inputted as vector
#' EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, ~0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_L50 <- c(250, 250) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_L95 <- c(300, 300) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
#'                           EstMalMatAtLen=NA) # maturity at length, inputted as values in data frame
#' sel_L50 <- c(250, 250) # females, males - Logistic length selectivity relationship parameters
#' sel_L95 <- c(300, 300) # females, males - Logistic length selectivity relationship parameters
#' EstGearSelAtLen <- data.frame(EstFemGearSelAtLen=NA,
#'                               EstMalGearSelAtLen=NA)
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_L50 <- c(250, 250) # females, males - Logistic fish retention at length parameters
#' ret_L95 <- c(300, 300) # females, males - Logistic fish retention at length parameters
#' EstRetenAtLen <- data.frame(EsFemtRetenAtLen=NA,
#'                             EstMalRetenAtLen=NA)
#' DiscMort <- 0.25 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 4.22 / MaxModelAge # natural mortality  (year-1)
#' Current_F = 0.2
#' Output_Opt = 1 # 1=standard output, 2=with added length and weight outputs (slower)
#' Res=GetPerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
#'                             RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
#'                             ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
#'                             EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
#'                             ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, Current_F, Output_Opt)
#' # Hermaphroditic species
#' MaxModelAge <- 100 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
#' MaxLen = 1400
#' LenInc = 20
#' lbnd = seq(0,MaxLen - LenInc, LenInc)
#' ubnd = lbnd + LenInc
#' midpt = lbnd + (LenInc/2)
#' nLenCl = length(midpt)
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf <- c(682, 982) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.14, 0.08) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0.06, -0.48) # years - von Bertalanffy growth model parameters - Females, males
#' GrowthParams = data.frame(Linf=Linf, vbK=vbK, tzero=tzero)
#' RefnceAges = NA
#' # GrowthParams c(Linf, vbK, tzero) single sex von Bertalanffy, or data.frame(Linf=Linf, vbK=vbK, tzero=tzero),
#' #' both sexes von Bertalanffy, or c(y1, y2, a, b) single sex Schnute, or data.frame(y1=y1, y2=y2, a=a, b=b), both sexes Schnute
#' CVSizeAtAge = c(0.05,0.05)
#' lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- -11.017 # for log-log relationship
#' lenwt_b <- 3.041 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 2 # 1=power, 2=log-log relationship
#' EstWtAtLen <- data.frame(EstFemWtAtLen=NA,
#'                          EstMalWtAtLen=NA) # weight at length, inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 2 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 1.0 # Ratio of females to males at recruitment age
#' FinalSex_L50 <- 821 # Logistic sex change relationship parameters (inflection point)
#' FinalSex_L95 <- 930 # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtLen <- NA  # sex ratio at length inputted as vector
#' EggFertParam <- 0.5 # (NA or from ~0.2-1) NA = no effect, ~0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_L50 <- c(653, 0) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_L95 <- c(950, 1) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
#'                           EstMalMatAtLen=NA) # maturity at length, inputted as values in data frame
#' sel_L50 <- c(500, 500) # females, males - Logistic length selectivity relationship parameters
#' sel_L95 <- c(600, 600) # females, males - Logistic length selectivity relationship parameters
#' EstGearSelAtLen <- data.frame(EstFemGearSelAtLen=NA,
#'                               EstMalGearSelAtLen=NA)
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_L50 <- c(499, 499) # females, males - Logistic fish retention at length parameters
#' ret_L95 <- c(500, 500) # females, males - Logistic fish retention at length parameters
#' EstRetenAtLen <- data.frame(EsFemtRetenAtLen=NA,
#'                             EstMalRetenAtLen=NA)
#' DiscMort <- 0.25 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort <- 4.22 / 70 # natural mortality  (year-1)
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' PlotOpt <- 0 # 0=all plots, 1=len at-age, 2=wt at length, 3=fem mat/sel/ret at length, 4=mal mat/sel/ret at length,
#' # 5=fem F at length, 6=mal F at length, 7=fem rel surv, 8=mal rel surv, 9=fem biom at age, 10=fem biom at age,
#' # 11=ypr/eq catch, 12=fem SPR/Brel, 13=mal SPR/Brel, 14=eq recruit
#' Current_F = 0.05
#' Output_Opt = 1 # 1=standard output, 2=with added length and weight outputs (slower)
#' Res=GetPerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
#'                             RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
#'                             ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
#'                             EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
#'                             ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, Current_F, Output_Opt)
#' @export
GetPerRecruitResults_LB <- function(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                    RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                                    ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                                    EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                                    ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, Current_F, Output_Opt) {
  FishMort <- seq(0,2,0.01)
  nFVals <- length(FishMort) # fishing mortality

  if (Output_Opt==2) { # include predicted length and weight distn. vs F, for output
    # mean length of catch
    FemMeanCatchLen <- rep(0,nFVals); MalMeanCatchLen <- rep(0,nFVals); CombSexMeanCatchLen <- rep(0,nFVals)
    # 95 percent confidence limits
    FemMeanCatchLen.lw95ci <- rep(0,nFVals); FemMeanCatchLen.up95ci <- rep(0,nFVals)
    MalMeanCatchLen.lw95ci <- rep(0,nFVals); MalMeanCatchLen.up95ci <- rep(0,nFVals)
    CombSexMeanCatchLen.lw95ci <- rep(0,nFVals); CombSexMeanCatchLen.up95ci <- rep(0,nFVals)
    # 60 quantile range
    FemCatchLen.20qntl <- rep(0,nFVals); FemCatchLen.80qntl <- rep(0,nFVals)
    MalCatchLen.20qntl <- rep(0,nFVals); MalCatchLen.80qntl <- rep(0,nFVals)
    CombSexCatchLen.20qntl <- rep(0,nFVals); CombSexCatchLen.80qntl <- rep(0,nFVals)
    # 95 quantile range
    FemCatchLen.2.5qntl <- rep(0,nFVals); FemCatchLen.97.5qntl <- rep(0,nFVals)
    MalCatchLen.2.5qntl <- rep(0,nFVals); MalCatchLen.97.5qntl <- rep(0,nFVals)
    CombSexCatchLen.2.5qntl <- rep(0,nFVals); CombSexCatchLen.97.5qntl <- rep(0,nFVals)
    # mean weight of catch
    FemMeanCatchWt <- rep(0,nFVals); MalMeanCatchWt <- rep(0,nFVals); CombSexMeanCatchWt <- rep(0,nFVals)
    # 95 percent confidence limits
    FemMeanCatchWt.lw95ci <- rep(0,nFVals); FemMeanCatchWt.up95ci <- rep(0,nFVals)
    MalMeanCatchWt.lw95ci <- rep(0,nFVals); MalMeanCatchWt.up95ci <- rep(0,nFVals)
    CombSexMeanCatchWt.lw95ci <- rep(0,nFVals); CombSexMeanCatchWt.up95ci <- rep(0,nFVals)
    # 60 quantile range
    FemCatchWt.20qntl <- rep(0,nFVals); FemCatchWt.80qntl <- rep(0,nFVals)
    MalCatchWt.20qntl <- rep(0,nFVals); MalCatchWt.80qntl <- rep(0,nFVals)
    CombSexCatchWt.20qntl <- rep(0,nFVals); CombSexCatchWt.80qntl <- rep(0,nFVals)
    # 95 quantile range
    FemCatchWt.2.5qntl <- rep(0,nFVals); FemCatchWt.97.5qntl <- rep(0,nFVals)
    MalCatchWt.2.5qntl <- rep(0,nFVals); MalCatchWt.97.5qntl <- rep(0,nFVals)
    CombSexCatchWt.2.5qntl <- rep(0,nFVals); CombSexCatchWt.97.5qntl <- rep(0,nFVals)
  }

  YPRResults <- rep(0,nFVals)
  Fem_SPRResults <- rep(0,nFVals)
  Mal_SPRResults <- rep(0,nFVals)
  CombSex_SPRResults <- rep(0,nFVals)

  Eq_RecResults <- rep(0,nFVals)
  Eq_Rec_AllEggFertResults = rep(0,nFVals)
  Eq_CatchResults <- rep(0,nFVals)
  Eq_FemRelSpBiomResults <- rep(0,nFVals)
  Eq_FemRelSpBiom_AllEggFert_Results <- rep(0,nFVals)
  Eq_MalRelSpBiomResults <- rep(0,nFVals)
  Eq_CombSexRelSpBiomResults <- rep(0,nFVals)
  Eq_FertRateResults <- rep(0,nFVals)
  UnfishMalToFemPropResults <- rep(0,nFVals)
  FishMalToFemPropResults <- rep(0,nFVals)
  Eq_MalDeplRatioResults <- rep(0,nFVals)
  Eq_FertRateResults <- rep(0,nFVals)

  for (k in 1:nFVals) {
    FMort = FishMort[k]
    # cat("k",k,"FMort",FMort,'\n')
    Res = CalcYPRAndSPRForFMort_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                   RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                                   ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                                   EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                                   ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, FMort)



    # per recruit results
    YPRResults[k] <- Res$YPR
    Fem_SPRResults[k] = Res$Fem_SPR
    Mal_SPRResults[k] = Res$Mal_SPR
    CombSex_SPRResults[k] = Res$CombSex_SPR

    # extended model results
    Eq_RecResults[k] = Res$Eq_Rec
    Eq_Rec_AllEggFertResults[k] = Res$Eq_Rec_AllEggFert
    Eq_CatchResults[k] = Res$Eq_Catch
    Eq_FemRelSpBiomResults[k] = Res$Eq_FemRelSpBiom
    Eq_FemRelSpBiom_AllEggFert_Results[k] = Res$Eq_FemRelSpBiom_AllEggFert
    Eq_MalRelSpBiomResults[k] = Res$Eq_MalRelSpBiom
    Eq_CombSexRelSpBiomResults[k] = Res$Eq_CombSexRelSpBiom
    UnfishMalToFemPropResults[k] = Res$UnfishMalToFemProp
    FishMalToFemPropResults[k] = Res$FishMalToFemProp
    Eq_MalDeplRatioResults[k] = Res$MalDeplRatio
    Eq_FertRateResults[k] = Res$Eq_FertRate

    if (Output_Opt==2) {
      # Random data and stats for expected lengths and weights in catches
      FemCatchNumLen = Res$ModelDiag$FemCatchNumLen
      MalCatchNumLen = Res$ModelDiag$MalCatchNumLen
      CombSexCatchNumLen = Res$ModelDiag$CombSexCatchNumLen
      FemWtAtLen = Res$ModelDiag$FemWtAtLen
      MalWtAtLen = Res$ModelDiag$MalWtAtLen
      CatchLenWtRes = CalcPerRecruitFishLenAndFishWtStats(lbnd, midpt, ubnd, FemCatchNumLen, MalCatchNumLen, CombSexCatchNumLen,
                                                          lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen)


      FemMeanCatchLen[k] <- CatchLenWtRes$PerRecMeanCatchLen$FemMeanCatchLen
      MalMeanCatchLen[k] <- CatchLenWtRes$PerRecMeanCatchLen$MalMeanCatchLen
      CombSexMeanCatchLen[k] <- CatchLenWtRes$PerRecMeanCatchLen$CombSexMeanCatchLen
      FemMeanCatchLen.lw95ci[k] <- CatchLenWtRes$MeanCatchLen.95ci$FemMeanCatchLen.95ci[1]   # 95 percent confidence limits
      FemMeanCatchLen.up95ci[k] <- CatchLenWtRes$MeanCatchLen.95ci$FemMeanCatchLen.95ci[2]
      MalMeanCatchLen.lw95ci[k] <- CatchLenWtRes$MeanCatchLen.95ci$MalMeanCatchLen.95ci[1]
      MalMeanCatchLen.up95ci[k] <- CatchLenWtRes$MeanCatchLen.95ci$MalMeanCatchLen.95ci[2]
      CombSexMeanCatchLen.lw95ci[k] <- CatchLenWtRes$MeanCatchLen.95ci$CombSexMeanCatchLen.95ci[1]
      CombSexMeanCatchLen.up95ci[k] <- CatchLenWtRes$MeanCatchLen.95ci$CombSexMeanCatchLen.95ci[2]
      FemCatchLen.20qntl[k] <- CatchLenWtRes$CatchLen.60qntl$FemCatchLen.60qntl[1] # 60 percent confidence prediction limits
      FemCatchLen.80qntl[k] <- CatchLenWtRes$CatchLen.60qntl$FemCatchLen.60qntl[2]
      MalCatchLen.20qntl[k] <- CatchLenWtRes$CatchLen.60qntl$MalCatchLen.60qntl[1]
      MalCatchLen.80qntl[k] <- CatchLenWtRes$CatchLen.60qntl$MalCatchLen.60qntl[2]
      CombSexCatchLen.20qntl[k] <- CatchLenWtRes$CatchLen.60qntl$CombSexCatchLen.60qntl[1]
      CombSexCatchLen.80qntl[k] <- CatchLenWtRes$CatchLen.60qntl$CombSexCatchLen.60qntl[2]
      FemCatchLen.2.5qntl[k] <- CatchLenWtRes$CatchLen.95qntl$FemCatchLen.95qntl[1] # 95 percent confidence prediction limits
      FemCatchLen.97.5qntl[k] <- CatchLenWtRes$CatchLen.95qntl$FemCatchLen.95qntl[2]
      MalCatchLen.2.5qntl[k] <- CatchLenWtRes$CatchLen.95qntl$MalCatchLen.95qntl[1]
      MalCatchLen.97.5qntl[k] <- CatchLenWtRes$CatchLen.95qntl$MalCatchLen.95qntl[2]
      CombSexCatchLen.2.5qntl[k] <- CatchLenWtRes$CatchLen.95qntl$CombSexCatchLen.95qntl[1]
      CombSexCatchLen.97.5qntl[k] <- CatchLenWtRes$CatchLen.95qntl$CombSexCatchLen.95qntl[2]

      FemMeanCatchWt[k] <- CatchLenWtRes$PerRecMeanCatchWt$FemMeanCatchWt # mean weight of catch
      MalMeanCatchWt[k] <- CatchLenWtRes$PerRecMeanCatchWt$MalMeanCatchWt
      CombSexMeanCatchWt[k] <- CatchLenWtRes$PerRecMeanCatchWt$CombSexMeanCatchWt
      FemMeanCatchWt.lw95ci[k] <- CatchLenWtRes$MeanCatchWt.95ci$FemMeanCatchWt.95ci[1] # 95 percent confidence limits
      FemMeanCatchWt.up95ci[k] <- CatchLenWtRes$MeanCatchWt.95ci$FemMeanCatchWt.95ci[2]
      MalMeanCatchWt.lw95ci[k] <- CatchLenWtRes$MeanCatchWt.95ci$MalMeanCatchWt.95ci[1]
      MalMeanCatchWt.up95ci[k] <- CatchLenWtRes$MeanCatchWt.95ci$MalMeanCatchWt.95ci[2]
      CombSexMeanCatchWt.lw95ci[k] <- CatchLenWtRes$MeanCatchWt.95ci$CombSexMeanCatchWt.95ci[1]
      CombSexMeanCatchWt.up95ci[k] <- CatchLenWtRes$MeanCatchWt.95ci$CombSexMeanCatchWt.95ci[2]
      FemCatchWt.20qntl[k] <- CatchLenWtRes$CatchWt.60qntl$FemCatchWt.60qntl[1] # 60 percent confidence prediction limits
      FemCatchWt.80qntl[k] <- CatchLenWtRes$CatchWt.60qntl$FemCatchWt.60qntl[2]
      MalCatchWt.20qntl[k] <- CatchLenWtRes$CatchWt.60qntl$MalCatchWt.60qntl[1]
      MalCatchWt.80qntl[k] <- CatchLenWtRes$CatchWt.60qntl$MalCatchWt.60qntl[2]
      CombSexCatchWt.20qntl[k] <- CatchLenWtRes$CatchWt.60qntl$CombSexCatchWt.60qntl[1]
      CombSexCatchWt.80qntl[k] <- CatchLenWtRes$CatchWt.60qntl$CombSexCatchWt.60qntl[2]
      FemCatchWt.2.5qntl[k] <- CatchLenWtRes$CatchWt.95qntl$FemCatchWt.95qntl[1] # 95 percent confidence prediction limits
      FemCatchWt.97.5qntl[k] <- CatchLenWtRes$CatchWt.95qntl$FemCatchWt.95qntl[2]
      MalCatchWt.2.5qntl[k] <- CatchLenWtRes$CatchWt.95qntl$MalCatchWt.95qntl[1]
      MalCatchWt.97.5qntl[k] <- CatchLenWtRes$CatchWt.95qntl$MalCatchWt.95qntl[2]
      CombSexCatchWt.2.5qntl[k] <- CatchLenWtRes$CatchWt.95qntl$CombSexCatchWt.95qntl[1]
      CombSexCatchWt.97.5qntl[k] <- CatchLenWtRes$CatchWt.95qntl$CombSexCatchWt.95qntl[2]
    }

  }

  # Length results
  if (Output_Opt==1) {

    MeanCatchLenResults=NA
    MeanCatchWtResults=NA
    MeanCatchLenResults_CurrF=NA
    MeanCatchWtResults_CurrF=NA

  } else {

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
                                     FemCatchLen.20qntl=FemCatchLen.20qntl,
                                     MalCatchLen.20qntl=MalCatchLen.20qntl,
                                     CombSexCatchLen.20qntl=CombSexCatchLen.20qntl,
                                     FemCatchLen.80qntl=FemCatchLen.80qntl,
                                     MalCatchLen.80qntl=MalCatchLen.80qntl,
                                     CombSexCatchLen.80qntl=CombSexCatchLen.80qntl,
                                     FemCatchLen.2.5qntl=FemCatchLen.2.5qntl,
                                     MalCatchLen.2.5qntl=MalCatchLen.2.5qntl,
                                     CombSexCatchLen.2.5qntl=CombSexCatchLen.2.5qntl,
                                     FemCatchLen.97.5qntl=FemCatchLen.97.5qntl,
                                     MalCatchLen.97.5qntl=MalCatchLen.97.5qntl,
                                     CombSexCatchLen.97.5qntl=CombSexCatchLen.97.5qntl)


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
                                    FemCatchWt.20qntl=FemCatchWt.20qntl,
                                    MalCatchWt.20qntl=MalCatchWt.20qntl,
                                    CombSexCatchWt.20qntl=CombSexCatchWt.20qntl,
                                    FemCatchWt.80qntl=FemCatchWt.80qntl,
                                    MalCatchWt.80qntl=MalCatchWt.80qntl,
                                    CombSexCatchWt.80qntl=CombSexCatchWt.80qntl,
                                    FemCatchWt.2.5qntl=FemCatchWt.2.5qntl,
                                    MalCatchWt.2.5qntl=MalCatchWt.2.5qntl,
                                    CombSexCatchWt.2.5qntl=CombSexCatchWt.2.5qntl,
                                    FemCatchWt.97.5qntl=FemCatchWt.97.5qntl,
                                    MalCatchWt.97.5qntl=MalCatchWt.97.5qntl,
                                    CombSexCatchWt.97.5qntl=CombSexCatchWt.97.5qntl)
  }

  YPRmax <- max(YPRResults) # maximum yield per recruit
  maxeqCatch <- max(Eq_CatchResults) # maximum equilibrium catch
  Fmax <- FishMort[which(YPRResults==YPRmax)] # fishing mortality at ypr maximum
  FmaxeqCatch <- FishMort[which(Eq_CatchResults==maxeqCatch)] # fishing mortality at maxeqCatch

  # Calculate biological reference points
  # assuming B_thresh=BMSY, B_lim=0.5BMSY and B_targ=1.2BMSY
  # and F_thresh, F_lim and F_targ are values of F resulting in B_thresh, B_lim and B_targ, respectively

  # gonochoristic species or considering sex ratio effect on egg fertilisation rate
  if (ReprodPattern == 1 | !is.na(EggFertParam)) {
    F_MSY=FmaxeqCatch
    x=which(FishMort==F_MSY)
    BMSY_Thresh=Eq_FemRelSpBiomResults[x]
    BMSY_Lim=0.5*BMSY_Thresh
    x=which.min(abs(Eq_FemRelSpBiomResults - BMSY_Lim))
    F_Lim = FishMort[x]
    BMSY_Targ=1.2*BMSY_Thresh
    x=which.min(abs(Eq_FemRelSpBiomResults - BMSY_Targ))
    F_Targ = FishMort[x]
  } else { # hermaphroditic species and not considering sex ratio effect on egg fertilisation rate
    F_MSY=FmaxeqCatch
    x=which(FishMort==F_MSY)
    BMSY_Thresh=Eq_CombSexRelSpBiomResults[x]
    BMSY_Lim=0.5*BMSY_Thresh
    x=which.min(abs(Eq_CombSexRelSpBiomResults - BMSY_Lim))
    F_Lim = FishMort[x]
    BMSY_Targ=1.2*BMSY_Thresh
    x=which.min(abs(Eq_CombSexRelSpBiomResults - BMSY_Targ))
    F_Targ = FishMort[x]
  }

  # get results for current F
  FMort = Current_F
  Res2 = CalcYPRAndSPRForFMort_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                  RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                                  ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                                  EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                                  ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, FMort)


  # get results associated with predicted length and weight distributions
  if (Output_Opt==2) {
    FemCatchNumLen = Res2$ModelDiag$FemCatchNumLen
    MalCatchNumLen = Res2$ModelDiag$MalCatchNumLen
    CombSexCatchNumLen = Res2$ModelDiag$CombSexCatchNumLen
    FemWtAtLen = Res2$ModelDiag$FemWtAtLen
    MalWtAtLen = Res2$ModelDiag$MalWtAtLen
    CatchLenWtRes = CalcPerRecruitFishLenAndFishWtStats(lbnd, midpt, ubnd, FemCatchNumLen, MalCatchNumLen, CombSexCatchNumLen,
                                                        lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen)

    MeanCatchLenResults_CurrF <- list(PerRecMeanCatchLen = CatchLenWtRes$PerRecMeanCatchLen,
                                      MeanCatchLen.95ci = CatchLenWtRes$MeanCatchLen.95ci,
                                      CatchLen.60qntl = CatchLenWtRes$CatchLen.60qntl,
                                      CatchLen.95qntl = CatchLenWtRes$CatchLen.95qntl,
                                      MeanCatchLenResults = MeanCatchLenResults,
                                      PerRecLenFreq = CatchLenWtRes$PerRecLenFreq,
                                      PerRecRandFishLen = CatchLenWtRes$PerRecRandFishLen)

    MeanCatchWtResults_CurrF <- list(PerRecMeanCatchWt = CatchLenWtRes$PerRecMeanCatchWt,
                                     MeanCatchWt.95ci = CatchLenWtRes$MeanCatchWt.95ci,
                                     CatchWt.60qntl = CatchLenWtRes$CatchWt.60qntl,
                                     CatchWt.95qntl = CatchLenWtRes$CatchWt.95qntl,
                                     MeanCatchWtResults = MeanCatchWtResults,
                                     PerRecRandFishWt = CatchLenWtRes$PerRecRandFishWt)
  }

  Results = list(YPR = Res2$YPR,
                 YPRmax = YPRmax,
                 Fmax = Fmax,
                 Fem_SPR = Res2$Fem_SPR,
                 Mal_SPR = Res2$Mal_SPR,
                 CombSex_SPR = Res2$CombSex_SPR,
                 Eq_Rec = Res2$Eq_Rec,
                 MalDeplRatio = Res2$MalDeplRatio,
                 Eq_FertRate = Res2$Eq_FertRate,
                 Eq_Rec_AllEggFert = Res2$Eq_Rec_AllEggFert,
                 Eq_Catch = Res2$Eq_Catch,
                 Eq_FemRelSpBiom = Res2$Eq_FemRelSpBiom,
                 Eq_FemRelSpBiom_AllEggFert = Res2$Eq_FemRelSpBiom_AllEggFert,
                 Eq_MalRelSpBiom = Res2$Eq_MalRelSpBiom,
                 Eq_CombSexRelSpBiom = Res2$Eq_CombSexRelSpBiom,
                 F_MSY = F_MSY,
                 BMSY = BMSY_Thresh,
                 F_Targ = F_Targ,
                 F_Thresh = F_MSY,
                 F_Lim = F_Lim,
                 BMSY_Targ = BMSY_Targ,
                 BMSY_Thresh = BMSY_Thresh,
                 BMSY_Lim = BMSY_Lim,
                 maxeqCatch = maxeqCatch,
                 FmaxeqCatch = FmaxeqCatch,
                 FishMort = FishMort,
                 YPRResults = YPRResults,
                 Eq_CatchResults = Eq_CatchResults,
                 Fem_SPRResults = Fem_SPRResults,
                 Mal_SPRResults = Mal_SPRResults,
                 CombSex_SPRResults = CombSex_SPRResults,
                 Eq_FemRelSpBiomResults = Eq_FemRelSpBiomResults,
                 Eq_FemRelSpBiom_AllEggFert_Results = Eq_FemRelSpBiom_AllEggFert_Results,
                 Eq_MalRelSpBiomResults = Eq_MalRelSpBiomResults,
                 Eq_CombSexRelSpBiomResults = Eq_CombSexRelSpBiomResults,
                 Eq_RecResults = Eq_RecResults,
                 Eq_Rec_AllEggFertResults = Eq_Rec_AllEggFertResults,
                 UnfishMalToFemPropResults = UnfishMalToFemPropResults,
                 FishMalToFemPropResults = FishMalToFemPropResults,
                 Eq_MalDeplRatioResults = Eq_MalDeplRatioResults,
                 Eq_FertRateResults = Eq_FertRateResults,
                 MeanCatchLenResults = MeanCatchLenResults,
                 MeanCatchLenResults_CurrF = MeanCatchLenResults_CurrF,
                 MeanCatchWtResults = MeanCatchWtResults,
                 MeanCatchWtResults_CurrF = MeanCatchWtResults_CurrF,
                 ModelDiag = Res2$ModelDiag)

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
#' @param FinalSex_A50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_A95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param EstSexRatioAtAge NA  # sex ratio at age (from age 0) inputted as vector
#' @param EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' @param mat_A50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_A95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtAge vector of proportion mature at age (set to NA if using age at maturity parameters)
#' @param Gear_sel_A50 logistic parameter for gear selectivity curve
#' @param Gear_sel_A95 logistic parameter for gear selectivity curve
#' @param EstGearSelAtAge vector for gear selectivity at age (set to NA if using age at gear selectivity parameters)
#' @param Land_sel_A50 logistic parameter for selectivity of landings curve
#' @param Land_sel_A95 logistic parameter for selectivity of landings curve
#' @param EstLandSelAtAge vector for for selectivity of landings at age (set to NA if using age at selectivity of landings parameters)
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_A50 logistic parameter for fish retention curve
#' @param ret_A95 logistic parameter for fish retention curve
#' @param EstRetenAtAge vector of fish retention at age (set to NA if using age at fish retention parameters)
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param PlotOpt 0=all plots, 1=len at-age, 2=wt at age, 3=fem mat/sel/ret at age, 4=mal mat/sel/ret at age,
#' 5=fem F at age, 6=mal F at age, 7=fem rel surv, 8=mal rel surv, 9=fem biom at age, 10=fem biom at age,
#' 11=catch at age, 12=ypr/eq catch, 13=fem SPR/Brel, 14=mal SPR/Brel, 15=comb sex SPR/Brel, 16=eq recruit
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
#' FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtAge <- NA  # sex ratio at age (from age 0) inputted as vector
#' EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' Gear_sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' Gear_sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' EstGearSelAtAge <- data.frame(FemGearSelAtAge=NA, MalGearSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' Land_sel_A50 <- NA # females, males - Logistic age selectivity relationship parameters
#' Land_sel_A95 <- NA # females, males - Logistic age selectivity relationship parameters
#' EstLandSelAtAge <- data.frame(FemSelLandAtAge=NA, MalSelLandAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
#' EstRetenAtAge <- data.frame(EstFemRetenAtAge=rep(1,MaxModelAge+1), EstMalRetenAtAge=rep(1,MaxModelAge+1)) # fish retention at age (from age 0), inputted as values in data frame
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort <- 0.2 # natural mortality  (year-1)
#' Current_F <- 0.07 # estimate of fishing mortality, e.g. from catch curve analysis
#' PlotOpt <- 0 # 0=all plots, 1=len at-age, 2=wt at age, 3=fem mat/sel_land at age, 4=mal mat/sel_land at age,
#' # 5=fem sel_land/sel_disc/prob reten/gear sel, 6=fem sel_land/sel_disc/prob reten/gear sel, 7=fem F at age,
#' # 8=mal F at age, 9=fem rel surv, 10=mal rel surv, 11=fem biom at age, 12=fem biom at age,
#' # 13=catch at age, 14=ypr/eq catch, 15=fem SPR/Brel, 16=mal SPR/Brel, 17=comb sex SPR/Brel, 18=eq recruit
#' # 19=male depletion vs F, 20=plot prop male vs F, 21=plot male depletion vs egg. fert rate,
#' # 22=plot male depletion vs eq. recruitment (plots 1-22 plotted if !is.na(EggFertParam))
#' PlotPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                          lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
#'                          ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam,
#'                          mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge,
#'                          Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95,
#'                          EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, PlotOpt, Current_F)
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
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 2 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 1 # Ratio of females to males at age zero
#' FinalSex_A50 <- 35 # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- 60 # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtAge <- NA  # sex ratio at age (from age 0) inputted as vector
#' EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_A50 <- c(20, 20) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(30, 30) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' Gear_sel_A50 <- c(20, 20) # females, males - Logistic age selectivity relationship parameters
#' Gear_sel_A95 <- c(30, 30) # females, males - Logistic age selectivity relationship parameters
#' EstGearSelAtAge <- data.frame(FemGearSelAtAge=NA, MalGearSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' Land_sel_A50 <- c(25, 25) # females, males - Logistic age selectivity relationship parameters
#' Land_sel_A95 <- c(35, 35) # females, males - Logistic age selectivity relationship parameters
#' EstLandSelAtAge <- data.frame(FemSelLandAtAge=NA, MalSelLandAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' ret_Pmax <- c(1,1)  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- c(15, 15)  # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- c(25, 25)  # females, males - Logistic age fish retention at age parameters
#' EstRetenAtAge <- data.frame(EstFemRetenAtAge=NA, EstMalRetenAtAge=NA) # fish retention at age (from age 0), inputted as values in data frame
#' DiscMort <- 0.5 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort <- 0.08 # natural mortality  (year-1)
#' Current_F <- 0.08 # estimate of fishing mortality, e.g. from catch curve analysis
#' PlotOpt <- 0 # 0=all plots, 1=len at-age, 2=wt at age, 3=fem mat/sel_land at age, 4=mal mat/sel_land at age,
#' # 5=fem sel_land/sel_disc/prob reten/gear sel, 6=fem sel_land/sel_disc/prob reten/gear sel, 7=fem F at age,
#' # 8=mal F at age, 9=fem rel surv, 10=mal rel surv, 11=fem biom at age, 12=fem biom at age,
#' # 13=catch at age, 14=ypr/eq catch, 15=fem SPR/Brel, 16=mal SPR/Brel, 17=comb sex SPR/Brel, 18=eq recruit
#' # 19=male depletion vs F, 20=plot prop male vs F, 21=plot male depletion vs egg. fert rate,
#' # 22=plot male depletion vs eq. recruitment (plots 1-22 plotted if !is.na(EggFertParam))
#' PlotPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                          lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
#'                          ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam,
#'                          mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge,
#'                          Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95,
#'                          EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, PlotOpt, Current_F)
#' @export
PlotPerRecruitResults_AB <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                  lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale,
                                  ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam,
                                  mat_A50, mat_A95, EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge,
                                  Land_sel_A50, Land_sel_A95, EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95,
                                  EstRetenAtAge, DiscMort, Steepness, SRrel_Type, NatMort, PlotOpt, Current_F) {

  .pardefault <- par(no.readonly = TRUE) # store current par settings

  Res = GetPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge, lenwt_a, ln_lenwt_a,
                                lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern, InitRatioFem,
                                FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam, mat_A50, mat_A95,
                                EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge, Land_sel_A50, Land_sel_A95,
                                EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort,
                                Steepness, SRrel_Type, NatMort, Current_F)

  #Plot 1:
  if (PlotOpt==0) {
    par(mfrow = c(3,2), mar=c(3.5,4,2,2),
        oma=c(1,1,1,1),tck=-0.03,mgp = c(3, 0.5, 0))
  } else {
    # don't change user settings
  }

  # PlotOpt <- 0 # 0=all plots, 1=len at-age, 2=wt at age, 3=fem mat/sel_land at age, 4=mal mat/sel_land at age,
  # 5=fem sel_land/sel_disc/prob reten/gear sel, 6=mal sel_land/sel_disc/prob reten/gear sel, 7=fem F at age,
  # 8=mal F at age, 9=fem rel surv, 10=mal rel surv, 11=fem biom at age, 12=mal biom at age,
  # 13=catch at age, 14=ypr/eq catch, 15=fem SPR/Brel, 16=mal SPR/Brel, 17=comb sex SPR/Brel, 18=eq recruit
  # 19=male depletion vs F, 20=plot prop male vs F, 21=plot male depletion vs egg. fert rate,
  # 22=plot male depletion vs eq. recruitment

  # plot growth curve
  if (PlotOpt==0 | PlotOpt==1) {
    y1=max(Res$ModelDiag$FemLenAtAge)
    y2=max(Res$ModelDiag$MalLenAtAge)
    if (y1 > y2) {
      ylims = Get_yaxis_scale(c(0,Res$ModelDiag$FemLenAtAge))
    } else {
      ylims = Get_yaxis_scale(c(0,Res$ModelDiag$MalLenAtAge))
    }
    ymax = ylims$ymax; yint = ylims$yint
    xlims = Get_xaxis_scale(Res$ModelDiag$Ages)
    xmax = xlims$xmax; xint = xlims$xint
    plot(Res$ModelDiag$Ages,Res$ModelDiag$FemLenAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),
         col="red",yaxt="n",xaxt="n",ylab="",xlab="")
    lines(Res$ModelDiag$Ages, Res$ModelDiag$MalLenAtAge,col="blue")
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Fish length (mm"),plain(")"))),las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    legend('bottomright', col=c("red","blue"),legend=c("Female","Male"),bty='n', cex=0.8,lwd=1.75)
  }

  # plot weight at age
  if (PlotOpt==0 | PlotOpt==2) {
    y1=max(Res$ModelDiag$FemWtAtAge)
    y2=max(Res$ModelDiag$MalWtAtAge)
    if (y1 > y2) {
      ylims = Get_yaxis_scale(Res$ModelDiag$FemWtAtAge)
    } else {
      ylims = Get_yaxis_scale(Res$ModelDiag$MalWtAtAge)
    }
    ymax = ylims$ymax; yint = ylims$yint
    xlims = Get_xaxis_scale(Res$ModelDiag$Ages)
    xmax = xlims$xmax; xint = xlims$xint
    plot(Res$ModelDiag$Ages,Res$ModelDiag$FemWtAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),col="red",yaxt="n",xaxt="n",
         ylab="",xlab="")
    lines(Res$ModelDiag$Ages,Res$ModelDiag$MalWtAtAge,col="blue")
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Total weight (kg"),plain(")"))),las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    legend('bottomright', col=c("red","blue"),legend=c("Female","Male"),bty='n', cex=0.8,lwd=1.75)
  }

  # plot female maturity and selectivity of landings at age (and prop female for hermaphroditic species)
  if (PlotOpt==0 | PlotOpt==3) {
    xlims = Get_xaxis_scale(Res$ModelDiag$Ages)
    xmax = xlims$xmax; xint = xlims$xint
    plot(Res$ModelDiag$Ages,Res$ModelDiag$FemPropMatAtAge,"l", pch=16,frame.plot=F,ylim=c(0,1),xlim=c(0,xmax),col="red",yaxt="n",xaxt="n",
         ylab="",xlab="",cex=0.8, lwd=2)
    lines(Res$ModelDiag$Ages, Res$ModelDiag$FemSelLandAtAge, "l", col="black",lty="solid", cex=0.8)
    legend('bottomright', col=c("red","black"),
           legend=c("Fem. prop_mat","Fem. sel_land"),
           lty=c("solid","solid"),bty='n', cex=0.8,lwd=1.75)
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax=1, yint=0.5, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Proportion"))),las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    if (ReprodPattern > 1) {
      lines(Res$ModelDiag$Ages,Res$ModelDiag$PropFemAtAge, "l", col="red",lty="dashed",cex=1)
      legend('topright', col="red",legend="Prop. Fem.",
             lty="dashed",bty='n', cex=0.8,lwd=1)
    }
  }

  # plot male maturity and selectivity of landings at age (and prop male for hermaphroditic species)
  if (PlotOpt==0 | PlotOpt==4) {
    xlims = Get_xaxis_scale(Res$ModelDiag$Ages)
    xmax = xlims$xmax; xint = xlims$xint
    plot(Res$ModelDiag$Ages, Res$ModelDiag$MalPropMatAtAge,"l", pch=16, frame.plot=F,ylim=c(0,1),xlim=c(0,xmax),col="blue",yaxt="n",xaxt="n",
         ylab="",xlab="", cex=0.8, lwd=2)
    lines(Res$ModelDiag$Ages, Res$ModelDiag$MalSelLandAtAge, "l", col="black",lty="solid", cex=0.8)

    legend('bottomright', col=c("blue","black"),
           legend=c("Mal. prop_mat","Mal. sel_land"),
           lty=c("solid","solid"),bty='n', cex=0.8,lwd=1.75)
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax=1, yint=0.5, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Proportion"))),las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    if (ReprodPattern > 1) {
      lines(Res$ModelDiag$Ages,1-Res$ModelDiag$PropFemAtAge, "l", col="blue",lty="dashed",cex=1)
      legend('topright', col="blue",legend="Prop. Male",
             lty="dashed",bty='n', cex=0.8,lwd=1)
    }
  }

  # plot female selectivity of landings, selectivity of discards, prob retention, and gear selectivity at age
  if (PlotOpt==0 | PlotOpt==5) {
    xlims = Get_xaxis_scale(Res$ModelDiag$Ages)
    xmax = xlims$xmax; xint = xlims$xint
    plot(Res$ModelDiag$Ages,Res$ModelDiag$FemSelLandAtAge,"l", pch=16,frame.plot=F,ylim=c(0,1),xlim=c(0,xmax),col="black",yaxt="n",xaxt="n",
         ylab="",xlab="",cex=0.8, lwd=2)
    lines(Res$ModelDiag$Ages,Res$ModelDiag$FemSelDiscAtAge, "l", col="purple",lty="dotted",cex=0.8)
    lines(Res$ModelDiag$Ages, Res$ModelDiag$FemRetProbAtAge, "l", col="brown",lty="dashed", cex=0.8)
    lines(Res$ModelDiag$Ages, Res$ModelDiag$FemGearSelAtAge, "l", col="orange",lty="dotted", cex=0.8)
    legend('bottomright', col=c("black","purple","brown","orange"),
           legend=c("Fem. sel_land","Fem. sel_disc","Fem. prob_ret","Fem. gear_sel"),
           lty=c("solid","dotted","dashed","dotted"),bty='n', cex=0.8,lwd=c(2,1,1,1))
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax=1, yint=0.5, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Proportion"))),las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=c(2,1,1,1))
  }

  # plot male selectivity of landings, selectivity of discards, prob retention, and gear selectivity at age
  if (PlotOpt==0 | PlotOpt==6) {
    xlims = Get_xaxis_scale(Res$ModelDiag$Ages)
    xmax = xlims$xmax; xint = xlims$xint
    plot(Res$ModelDiag$Ages, Res$ModelDiag$MalSelLandAtAge,"l", pch=16, frame.plot=F,ylim=c(0,1),xlim=c(0,xmax),col="black",yaxt="n",xaxt="n",
         ylab="",xlab="", cex=0.8, lwd=2)
    lines(Res$ModelDiag$Ages,Res$ModelDiag$MalSelDiscAtAge, "l", col="purple",lty="dotted",cex=0.8)
    lines(Res$ModelDiag$Ages, Res$ModelDiag$MalRetProbAtAge, "l", col="brown",lty="dashed", cex=0.8)
    lines(Res$ModelDiag$Ages, Res$ModelDiag$MalGearSelAtAge, "l", col="orange",lty="dotted", cex=0.8)

    legend('bottomright', col=c("black","purple","brown","orange"),
           legend=c("Mal. sel_land","Mal. sel_disc","Mal. prob_ret","Mal. gear_sel"),
           lty=c("solid","dotted","dashed","dotted"),bty='n', cex=0.8,lwd=c(2,1,1,1))
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax=1, yint=0.5, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    axis(2,at=seq(0,1,0.5), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
    mtext(expression(paste(plain("Proportion"))),las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  }


  #plot 2:
  if (PlotOpt==0) {
    par(mfrow = c(3,2), mar=c(3.5,4,2,2),
        oma=c(1,1,1,1),tck=-0.03,mgp = c(3, 0.5, 0))
  } else {
    # don't change user settings
  }

  # plot female mortality at age
  if (PlotOpt==0 | PlotOpt==7) {
    ylims = Get_yaxis_scale(Res$ModelDiag$FemFAtAge)
    ymax = ylims$ymax; yint = ylims$yint
    xlims = Get_xaxis_scale(Res$ModelDiag$Ages)
    xmax = xlims$xmax; xint = xlims$xint
    plot(Res$ModelDiag$Ages, Res$ModelDiag$FemFAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),col="red",yaxt="n",xaxt="n",
         ylab="",xlab="")
    lines(Res$ModelDiag$Ages,Res$ModelDiag$FemDiscFAtAge,lty="dotted",col="brown")
    lines(Res$ModelDiag$Ages,Res$ModelDiag$FemLandFAtAge,lty="dotted",col="purple")
    lines(Res$ModelDiag$Ages,rep(NatMort,length(Res$ModelDiag$Ages)),lty="dashed")
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Mortality") ~ (year^{-1}))),las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    legend('topright', col=c("red","brown","purple","black"),lty=c("solid","dotted","dotted","dashed"),
           legend=c("Fem. F","Fem. DiscF","Fem. LandF","Fem. M"),bty='n', cex=1.0,lwd=1.75)
  }

  # plot male mortality at age
  if (PlotOpt==0 | PlotOpt==8) {
    ylims = Get_yaxis_scale(Res$ModelDiag$MalFAtAge)
    ymax = ylims$ymax; yint = ylims$yint
    xlims = Get_xaxis_scale(Res$ModelDiag$Ages)
    xmax = xlims$xmax; xint = xlims$xint
    plot(Res$ModelDiag$Ages, Res$ModelDiag$MalFAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),col="blue",yaxt="n",xaxt="n",
         ylab="",xlab="")
    lines(Res$ModelDiag$Ages,Res$ModelDiag$MalDiscFAtAge,lty="dotted",col="brown")
    lines(Res$ModelDiag$Ages,Res$ModelDiag$MalLandFAtAge,lty="dotted",col="purple")
    lines(Res$ModelDiag$Ages,rep(NatMort,length(Res$ModelDiag$Ages)),lty="dashed")
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Mortality") ~ (year^{-1}))),las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    legend('topright', col=c("blue","brown","purple","black"),lty=c("solid","dotted","dotted","dashed"),
           legend=c("Mal. F","Mal. DiscF","Mal.LandF","Mal. M"),bty='n', cex=1.0,lwd=1.75)
  }

  # plot fished and unfished female survival in terms of numbers given specified current fully-selected fishing mortality
  if (PlotOpt==0 | PlotOpt==9) {
    ylims = Get_yaxis_scale(Res$ModelDiag$UnfishFemSurvAtAge)
    ymax = ylims$ymax; yint = ylims$yint
    if (ymax > 1) {
      ymax=1
      yint=0.2
    }
    xlims = Get_xaxis_scale(Res$ModelDiag$Ages)
    xmax = xlims$xmax; xint = xlims$xint
    plot(Res$ModelDiag$Ages, Res$ModelDiag$UnfishFemSurvAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),col="red",yaxt="n",xaxt="n",
         ylab="",xlab="")
    lines(Res$ModelDiag$Ages, Res$ModelDiag$FishFemSurvAtAge,col="red",lty="dotted")
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Rel. survival"))),las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    legend('topright', col="red",lty=c("solid","dotted"),
           legend=c("Fem. unfish","Fem. fish"),bty='n', cex=1.0,lwd=1.75)
  }

  # plot fished and unfished male survival in terms of numbers given specified current fully-selected fishing mortality
  if (PlotOpt==0 | PlotOpt==10) {
    ylims = Get_yaxis_scale(Res$ModelDiag$UnfishMalSurvAtAge)
    ymax = ylims$ymax; yint = ylims$yint
    if (ymax > 1) {
      ymax=1
      yint=0.2
    }
    xlims = Get_xaxis_scale(Res$ModelDiag$Ages)
    xmax = xlims$xmax; xint = xlims$xint
    plot(Res$ModelDiag$Ages, Res$ModelDiag$UnfishMalSurvAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),col="blue",yaxt="n",xaxt="n",
         ylab="",xlab="")
    lines(Res$ModelDiag$Ages, Res$ModelDiag$FishMalSurvAtAge,col="blue",lty="dotted")
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Rel. survival"))),las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    legend('topright', col="blue",lty=c("solid","dotted"),
           legend=c("Mal. unfish","Mal. fish"),bty='n', cex=1.0,lwd=1.75)
  }

  # plot fished and unfished mature female biomass at age given specified current fully-selected fishing mortality
  if (PlotOpt==0 | PlotOpt==11) {
    ylims = Get_yaxis_scale(Res$ModelDiag$UnfishFemSpBiomAtAge)
    ymax = ylims$ymax; yint = ylims$yint
    if (ymax == 0) {
      ymax = round(1.4 * max(Res$ModelDiag$UnfishFemSpBiomAtAge),1)
      yint = ymax/4
    }
    xlims = Get_xaxis_scale(Res$ModelDiag$Ages)
    xmax = xlims$xmax; xint = xlims$xint
    plot(Res$ModelDiag$Ages, Res$ModelDiag$UnfishFemSpBiomAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),
         col="red",yaxt="n",xaxt="n",ylab="",xlab="")
    lines(Res$ModelDiag$Ages, Res$ModelDiag$FishFemSpBiomAtAge,col="red",lty="dotted")
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Mat. biom. at age (kg"),plain(")"))),las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    legend('topright', col=c("red","red"),lty=c("solid","dotted"),legend=c("Fem. unfish","Fem. fish"),bty='n', cex=1.0,lwd=1.75)
  }

  # plot fished and unfished mature male biomass at age given specified current fully-selected fishing mortality
  if (PlotOpt==0 | PlotOpt==12) {
    ylims = Get_yaxis_scale(Res$ModelDiag$UnfishMalSpBiomAtAge)
    ymax = ylims$ymax; yint = ylims$yint
    xlims = Get_xaxis_scale(Res$ModelDiag$Ages)
    xmax = xlims$xmax; xint = xlims$xint
    if (ymax == 0) {
      ymax = round(1.4 * max(Res$ModelDiag$UnfishFemSpBiomAtAge),1)
      yint = ymax/4
    }
    plot(Res$ModelDiag$Ages, Res$ModelDiag$UnfishMalSpBiomAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),
         col="blue",yaxt="n",xaxt="n",ylab="",xlab="")
    lines(Res$ModelDiag$Ages, Res$ModelDiag$FishMalSpBiomAtAge,col="blue",lty="dotted")
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Mat. biom. at age (kg"),plain(")"))),las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    legend('topright', col=c("blue","blue"),lty=c("solid","dotted"),legend=c("Mal. unfish","Mal. fish"),bty='n', cex=1.0,lwd=1.75)
  }

  #Plot 3:
  # plot female and male catch at age, given specified current fully-selected fishing mortality
    if (PlotOpt==0) {
      par(mfrow = c(3,2), mar=c(3.5,4,2,2),
          oma=c(1,1,1,1),tck=-0.03,mgp = c(3, 0.5, 0))
    } else {
      # don't change user settings
    }

    if (PlotOpt==0 | PlotOpt==13) {
      FemCatchNumAtAgeProp <- Res$ModelDiag$FemCatchAtAgeNum / sum(Res$ModelDiag$FemCatchAtAgeNum)
      MalCatchNumAtAgeProp <- Res$ModelDiag$MalCatchAtAgeNum / sum(Res$ModelDiag$MalCatchAtAgeNum)
      ylims = Get_yaxis_scale(FemCatchNumAtAgeProp)
      ymax = ylims$ymax; yint = ylims$yint
      xlims = Get_xaxis_scale(Res$ModelDiag$Ages)
      xmax = xlims$xmax; xint = xlims$xint
      plot(Res$ModelDiag$Ages, FemCatchNumAtAgeProp,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),
           col="red",yaxt="n",xaxt="n",ylab="",xlab="")
      lines(Res$ModelDiag$Ages, MalCatchNumAtAgeProp,col="blue","l")
      AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
      mtext("Catch at age Prop.",las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
      mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
      legend('topright', col=c("red","blue"),lty="solid",legend=c("females","males"),bty='n', cex=0.8,lwd=1.75)
    }

    # plot yield per recruit (per recruit analysis) and equilibrium catch (equilibrium age-structured model)
    # given specified current fully-selected fishing mortality
    if (PlotOpt==0 | PlotOpt==14) {
      ylims = Get_yaxis_scale(Res$YPRResults)
      ymax = ylims$ymax; yint = ylims$yint

      plot(Res$FishMort, Res$YPRResults,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,max(Res$FishMort)),
           col="black",yaxt="n",xaxt="n",ylab="",xlab="")
      points(Current_F, Res$YPR,cex=1.2,col="black",pch=16)
      lines(Res$FishMort, Res$Eq_CatchResults,col="blue")
      points(Current_F, Res$Eq_Catch, cex=1.2,col="blue",pch=16)
      AddAxesAndTickLabelsToPlot(xmin=NA, xmax=max(Res$FishMort), xint=0.5, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
      mtext(expression(paste(plain("YPR / Eq.Catch (kg"),plain(")"))),las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
      mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
      legend('topright', col=c("black","blue"),lty=c("solid","solid"),
             legend=c("YPR","Eq.catch"),bty='n', cex=0.8,lwd=1.75)
    }

    # plot spawning potential ratio (per recruit analysis) and relative spawning biomass (equilbrium age-structured model)
    # given specified current fully-selected fishing mortality, for female
    if (PlotOpt==0 | PlotOpt==15) {
      plot(Res$FishMort, Res$Fem_SPRResults,"l",frame.plot=F,ylim=c(0,1.0),xlim=c(0,max(Res$FishMort)),
           col="red",yaxt="n",xaxt="n",ylab="",xlab="", lty="dotted")
      lines(Res$FishMort, Res$Eq_FemRelSpBiomResults,col="red",lty="solid")
      points(Current_F, Res$Fem_SPR,cex=1.2,col="red",pch=16)
      points(Current_F, Res$Eq_FemRelSpBiom,cex=1.2,col="red",pch=1)
      AddAxesAndTickLabelsToPlot(xmin=NA, xmax=max(Res$FishMort), xint=0.5, ymin=NA, ymax=1, yint=0.2, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
      mtext(expression(paste(plain("Biom. ratio"))),las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
      mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)

      if (!is.na(EggFertParam) &  EggFertParam < 1) {
        lines(Res$FishMort, Res$Eq_FemRelSpBiom_AllEggFert_Results, lty="dotted")
        points(Current_F, Res$Eq_FemRelSpBiom_AllEggFert, cex=1.2,col="black",pch=16)
        legend('topright', col=c("red","black","red"),lty=c("dotted","dotted","solid"),
        legend=c("Fem SPR","Fem. rel biom. (100% fert.)",
                        paste0("Fem. rel biom (calc. fert. rate)")),bty='n', cex=0.8,lwd=1.75)
      } else {
        legend('topright', col=c("red","red"),lty=c("dotted","solid"),
               legend=c("Fem SPR","Fem. Rel biom"),bty='n', cex=0.8,lwd=1.75)
      }
    }



    # plot spawning potential ratio (per recruit analysis) and relative spawning biomass (equilbrium age-structured model)
    # given specified current fully-selected fishing mortality, for each sex
    if (PlotOpt==0 | PlotOpt==16) {
      plot(Res$FishMort, Res$Mal_SPRResults,"l",frame.plot=F,ylim=c(0,1.0),xlim=c(0,max(Res$FishMort)),
           col="blue",yaxt="n",xaxt="n",ylab="",xlab="", lty="dotted")
      lines(Res$FishMort, Res$Eq_MalRelSpBiomResults,col="blue",lty="solid")
      points(Current_F, Res$Mal_SPR,cex=1.2,col="blue",pch=16)
      points(Current_F, Res$Eq_MalRelSpBiom,cex=1.2,col="blue",pch=1)
      AddAxesAndTickLabelsToPlot(xmin=NA, xmax=max(Res$FishMort), xint=0.5, ymin=NA, ymax=1, yint=0.2, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
      mtext(expression(paste(plain("Biom. ratio"))),las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
      mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
      legend('topright', col=c("blue","blue"),lty=c("dotted","solid"),
             legend=c("Mal SPR","Mal. Rel biom"),bty='n', cex=0.8,lwd=1.75)
    }

    # plot spawning potential ratio (per recruit analysis) and relative spawning biomass (equilbrium age-structured model)
    # given specified current fully-selected fishing mortality, for combined sexes
    if (PlotOpt==0 | PlotOpt==17) {
      plot(Res$FishMort, Res$CombSex_SPRResults,"l",frame.plot=F,ylim=c(0,1.0),xlim=c(0,max(Res$FishMort)),
           col="black",yaxt="n",xaxt="n",ylab="",xlab="", lty="dotted")
      lines(Res$FishMort, Res$Eq_CombSexRelSpBiomResults,col="black",lty="solid")
      points(Current_F, Res$CombSex_SPR,cex=1.2,col="black",pch=16)
      points(Current_F, Res$Eq_CombSexRelSpBiom,cex=1.2,col="black",pch=1)
      AddAxesAndTickLabelsToPlot(xmin=NA, xmax=max(Res$FishMort), xint=0.5, ymin=NA, ymax=1, yint=0.2, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
      mtext(expression(paste(plain("Biom. ratio"))),las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
      mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
      legend('topright', col=c("black","black"),lty=c("dotted","solid"),
             legend=c("CombSex SPR","CombSex Rel biom"),bty='n', cex=0.8,lwd=1.75)
    }

    # plot equilibrium recruitment vs F
    if (PlotOpt==0 | PlotOpt==18) {
      ymax = 1.0
      yint = 0.2
      plot(Res$FishMort, Res$Eq_RecResults,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,max(Res$FishMort)),
           col="red",yaxt="n",xaxt="n",ylab="",xlab="")
      points(Current_F, Res$Eq_Rec, cex=1.2,col="red",pch=16)
      AddAxesAndTickLabelsToPlot(xmin=NA, xmax=max(Res$FishMort), xint=0.5, ymin=NA, ymax=1, yint=0.2, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
      mtext(expression(paste(plain("Equil. Recruitment"))),las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
      mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    }

  # plot 4
  if (!is.na(EggFertParam)) {
    if (PlotOpt==0) {
      par(mfrow = c(2,2), mar=c(3.5,4,2,2),
          oma=c(1,1,1,1),tck=-0.03,mgp = c(3, 0.5, 0))
    } else {
      # don't change user settings
    }

    # plot male depletion vs F
    # Male depletion measure is the ratio of the proportion of mature males:mature females in the population (in numbers)
    # at current fishing pressure relative to the proportion mature males: mature females at the unfished level
    if (PlotOpt==0 | PlotOpt==19) {
      ymax = 1; yint = 0.2
      xmax = max(Res$FishMort)
      plot(Res$FishMort,Res$Eq_MalDeplRatioResults,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),
           col="blue",yaxt="n",xaxt="n",ylab="",xlab="")
      points(Current_F, Res$MalDeplRatio, col="blue", pch=16)
      AddAxesAndTickLabelsToPlot(xmin=NA, xmax=max(Res$FishMort), xint=0.5, ymin=NA, ymax=1, yint=0.2, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
      mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
      mtext(expression("Male depletion (init. vs fin. sex ratio)"),las=3,side=2,line=2,cex=0.7,lwd=1.75)
      legend('topright', col="blue",pch=16, legend="Current F",bty='n', cex=0.8,lwd=-1)
    }

    # plot prop male vs F
    if (PlotOpt==0 | PlotOpt==20) {
      ylims = Get_yaxis_scale(Res$FishMalToFemPropResults)
      ymax = ylims$ymax; yint = ylims$yint
      xmax = max(Res$FishMort)
      x=which(Res$FishMort==Current_F)
      plot(Res$FishMort,Res$FishMalToFemPropResults,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),
           col="blue",yaxt="n",xaxt="n",ylab="",xlab="")
      points(Current_F, Res$FishMalToFemPropResults[x], col="blue", pch=16)
      AddAxesAndTickLabelsToPlot(xmin=NA, xmax=max(Res$FishMort), xint=0.5, ymin=NA, ymax=ymax, yint=yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
      mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
      mtext(expression("Prop. male to mat. fem."),las=3,side=2,line=2,cex=0.7,lwd=1.75)
      legend('topright', col="blue",pch=16, legend="Current F",bty='n', cex=0.8,lwd=-1)
    }

    # plot male depletion vs egg. fert rate
    if (PlotOpt==0 | PlotOpt==21) {
      ymax = 1; yint = 0.2
      xmax = 1; xint = 0.2
      plot(Res$FishMort, Res$Eq_FertRateResults,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),
           col="red",yaxt="n",xaxt="n",ylab="",xlab="")
      points(Res$FishMort[x], Res$Eq_FertRateResults[x], col="red", pch=16)
      AddAxesAndTickLabelsToPlot(xmin=NA, xmax=xmax, xint=xint, ymin=NA, ymax=ymax, yint=yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
      mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
      mtext(expression("Egg fert. rate"),las=3,side=2,line=2,cex=0.7,lwd=1.75)
      legend('topright', col="red",pch=16, legend="Current F",bty='n', cex=0.8,lwd=-1)
    }

    # plot male depletion vs eq. recruitment
    if (PlotOpt==0 | PlotOpt==22) {
      ymax = 1.0
      yint = 0.2
      plot(Res$FishMort, Res$Eq_RecResults,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,max(Res$FishMort)),
           col="red",yaxt="n",xaxt="n",ylab="",xlab="")
      points(Current_F, Res$Eq_Rec, cex=1.2,col="red",pch=16)
      AddAxesAndTickLabelsToPlot(xmin=NA, xmax=max(Res$FishMort), xint=0.5, ymin=NA, ymax=1, yint=0.2, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
      mtext(expression(paste(plain("Equil. Recruitment"))),las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
      mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
      if (!is.na(EggFertParam) &  EggFertParam < 1) {
        lines(Res$FishMort, Res$Eq_Rec_AllEggFertResults, lty="dotted")
        points(Current_F, Res$Eq_Rec_AllEggFert, cex=1.2,col="black",pch=16)
        legend('topright', col=c("black","red"),lty=c("dotted","solid"),
               legend=c("100% egg fert. rate","calc. egg fert. rate"),bty='n', cex=0.8,lwd=1.75)
      }
    }
  }

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
#' @param FinalSex_L50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param EstSexRatioAtLen NA  # sex ratio at length, inputted as vector
#' @param EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' @param mat_L50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at length)
#' @param mat_L95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at length)
#' @param EstMatAtLen vector of proportion mature at length (set to NA if using length at maturity parameters)
#' @param sel_L50 logistic parameter for gear selectivity curve
#' @param sel_L95 logistic parameter for gear selectivity curve
#' @param EstGearSelAtLen gear selectivity curve inputted as vector
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_L50 logistic parameter for fish retention curve
#' @param ret_L95 logistic parameter for fish retention curve
#' @param EstRetenAtLen retention curve inputted as vector
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param Current_F estimated current fishing mortality
#' @param PlotOpt 0=all plots, 1=len at-age, 2=wt at age, 3=fem mat/sel/ret at age, 4=mal mat/sel/ret at age,
#'  5=fem F at age, 6=mal F at age, 7=fem rel surv, 8=mal rel surv, 9=fem biom at age, 10=fem biom at age,
#'  11=catch at age, 12=ypr/eq catch, 13=fem SPR/Brel, 14=mal SPR/Brel, 15=comb sex SPR/Brel, 16=eq recruit
#'
#' @examples
#' # Non-hermaphroditic species
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
#'                          EstMalWtAtLen=NA) # weight at length, inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at recruitment age
#' FinalSex_L50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_L95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtLen <- NA  # sex ratio at length inputted as vector
#' EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, ~0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_L50 <- c(250, 250) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_L95 <- c(300, 300) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
#'                           EstMalMatAtLen=NA) # maturity at length, inputted as values in data frame
#' sel_L50 <- c(250, 250) # females, males - Logistic length selectivity relationship parameters
#' sel_L95 <- c(300, 300) # females, males - Logistic length selectivity relationship parameters
#' EstGearSelAtLen <- data.frame(EstFemGearSelAtLen=NA,
#'                               EstMalGearSelAtLen=NA)
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_L50 <- c(250, 250) # females, males - Logistic fish retention at length parameters
#' ret_L95 <- c(300, 300) # females, males - Logistic fish retention at length parameters
#' EstRetenAtLen <- data.frame(EsFemtRetenAtLen=NA,
#'                             EstMalRetenAtLen=NA)
#' DiscMort <- 0.25 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort <- 4.22 / MaxModelAge # natural mortality  (year-1)
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' PlotOpt <- 0 # 0=all plots, 1=len at-age, 2=wt at length, 3=fem mat/sel_land at length, 4=mal mat/sel_land at length,
#' # 5=fem sel_land/sel_disc/prob reten/gear sel, 6=fem sel_land/sel_disc/prob reten/gear sel, 7=fem F at length,
#' # 8=mal F at length, 9=fem rel surv, 10=mal rel surv, 11=fem biom at age, 12=fem biom at age,
#' # 13=catch at length, 14=ypr/eq catch, 15=fem SPR/Brel, 16=mal SPR/Brel, 17=comb sex SPR/Brel, 18=eq recruit
#' # 19=male depletion vs F, 20=plot prop male vs F, 21=plot male depletion vs egg. fert rate,
#' # 22=plot male depletion vs eq. recruitment (plots 1-22 plotted if !is.na(EggFertParam))
#' Current_F = 0.2
#' PlotPerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
#'                          RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
#'                          ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
#'                          EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
#'                          ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, Current_F, PlotOpt)
#' # Non-hermaphroditic species
#' MaxModelAge <- 100 # maximum age considered by model, years
#' TimeStep <- 1 # Model time step (y) (for shorter-lived species, might be appropriate to use a smaller time step)
#' MaxLen = 1400
#' LenInc = 20
#' lbnd = seq(0,MaxLen - LenInc, LenInc)
#' ubnd = lbnd + LenInc
#' midpt = lbnd + (LenInc/2)
#' nLenCl = length(midpt)
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf <- c(682, 982) # mm - von Bertalanffy growth model parameters - Females, males
#' vbK <- c(0.14, 0.08) # year-1 - von Bertalanffy growth model parameters - Females, males
#' tzero <- c(0.06, -0.48) # years - von Bertalanffy growth model parameters - Females, males
#' GrowthParams = data.frame(Linf=Linf, vbK=vbK, tzero=tzero)
#' RefnceAges = NA
#' # GrowthParams c(Linf, vbK, tzero) single sex von Bertalanffy, or data.frame(Linf=Linf, vbK=vbK, tzero=tzero),
#' #' both sexes von Bertalanffy, or c(y1, y2, a, b) single sex Schnute, or data.frame(y1=y1, y2=y2, a=a, b=b), both sexes Schnute
#' CVSizeAtAge = c(0.05,0.05)
#' lenwt_a <- NA # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' ln_lenwt_a <- -11.017 # for log-log relationship
#' lenwt_b <- 3.041 # combined sexes - weight (g) vs length (mm, TL) relationship parameters
#' WLrel_Type <- 2 # 1=power, 2=log-log relationship
#' EstWtAtLen <- data.frame(EstFemWtAtLen=NA,
#'                          EstMalWtAtLen=NA) # weight at length, inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 2 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 1.0 # Ratio of females to males at recruitment age
#' FinalSex_L50 <- 821 # Logistic sex change relationship parameters (inflection point)
#' FinalSex_L95 <- 930 # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtLen <- NA  # sex ratio at length inputted as vector
#' EggFertParam <- 0.5 # (NA or from ~0.2-1) NA = no effect, ~0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_L50 <- c(653, 0) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_L95 <- c(950, 1) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
#'                           EstMalMatAtLen=NA) # maturity at length, inputted as values in data frame
#' sel_L50 <- c(500, 500) # females, males - Logistic length selectivity relationship parameters
#' sel_L95 <- c(600, 600) # females, males - Logistic length selectivity relationship parameters
#' EstGearSelAtLen <- data.frame(EstFemGearSelAtLen=NA,
#'                               EstMalGearSelAtLen=NA)
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_L50 <- c(499, 499) # females, males - Logistic fish retention at length parameters
#' ret_L95 <- c(500, 500) # females, males - Logistic fish retention at length parameters
#' EstRetenAtLen <- data.frame(EsFemtRetenAtLen=NA,
#'                             EstMalRetenAtLen=NA)
#' DiscMort <- 0.25 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort <- 4.22 / 70 # natural mortality  (year-1)
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' PlotOpt <- 0 # 0=all plots, 1=len at-age, 2=wt at length, 3=fem mat/sel_land at length, 4=mal mat/sel_land at length,
#' # 5=fem sel_land/sel_disc/prob reten/gear sel, 6=fem sel_land/sel_disc/prob reten/gear sel, 7=fem F at length,
#' # 8=mal F at length, 9=fem rel surv, 10=mal rel surv, 11=fem biom at age, 12=fem biom at age,
#' # 13=catch at length, 14=ypr/eq catch, 15=fem SPR/Brel, 16=mal SPR/Brel, 17=comb sex SPR/Brel, 18=eq recruit
#' # 19=male depletion vs F, 20=plot prop male vs F, 21=plot male depletion vs egg. fert rate,
#' # 22=plot male depletion vs eq. recruitment (plots 1-22 plotted if !is.na(EggFertParam))
#' Current_F = 0.05
#' FMort = 0.05
#' PlotPerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
#'                          RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
#'                          ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
#'                          EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
#'                          ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, Current_F, PlotOpt)
#' @export
PlotPerRecruitResults_LB <- function(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                     RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                                     ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                                     EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                                     ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, Current_F, PlotOpt) {

  .pardefault <- par(no.readonly = TRUE) # store current par settings

  Ages <- seq(TimeStep,MaxModelAge,TimeStep)
  nTimeSteps <- length(Ages)

  # get prediction intervals for growth curves
  FMort = Current_F
  Res=CalcYPRAndSPRForFMort_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                               RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                               ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                               EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                               ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, FMort)

  Growth_95PLs = GetPerRecruitGrowthPredIntervals_LB(nTimeSteps, nLenCl, midpt, lbnd, ubnd, Res)


  Output_Opt = 1 # 1=standard output, 2=with added length and weight outputs (slower)
  Res = GetPerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                                ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                                EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                                ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, Current_F, Output_Opt)


  #Plot 1:
  if (PlotOpt==0) {
    par(mfrow = c(3,2), mar=c(3.5,4,2,2),
        oma=c(1,1,1,1),tck=-0.03,mgp = c(3, 0.5, 0))
  }

  # plot growth curve
  if (PlotOpt==0 | PlotOpt==1) {
      x=which(!is.na(Growth_95PLs$FemLenAtAge_hi))
      y1=max(Growth_95PLs$FemLenAtAge_hi[x])
      x=which(!is.na(Growth_95PLs$MalLenAtAge_hi))
      y2=max(Growth_95PLs$MalLenAtAge_hi[x])
      if (y1 > y2) {
        x=which(!is.na(Growth_95PLs$FemLenAtAge_hi))
        ylims = Get_yaxis_scale(Growth_95PLs$FemLenAtAge_hi[x])
      } else {
        x=which(!is.na(Growth_95PLs$MalLenAtAge_hi))
        ylims = Get_yaxis_scale(Growth_95PLs$MalLenAtAge_hi[x])
      }
      ymax = ylims$ymax; yint = ylims$yint
      xlims = Get_xaxis_scale(Ages)
      xmax = xlims$xmax; xint = xlims$xint
      plot(Ages,Res$ModelDiag$MeanSizeAtAge[1,],"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),
           col="red",yaxt="n",xaxt="n",ylab="",xlab="")
      x = c(Ages,rev(Ages))
      y = c(Growth_95PLs[,1],rev(Growth_95PLs[,2]))
      polygon(x,y, col="pink",border=NA)
      y = c(Growth_95PLs[,3],rev(Growth_95PLs[,4]))
      polygon(x,y, col="light blue",border=NA)
      lines(Ages, Res$ModelDiag$MeanSizeAtAge[1,],col="red", lwd=2)
      lines(Ages, Res$ModelDiag$MeanSizeAtAge[2,],col="blue", lwd=2)
      lines(Ages, Growth_95PLs[,1],col="red",lty="dotted")
      lines(Ages, Growth_95PLs[,2],col="red",lty="dotted")
      lines(Ages, Growth_95PLs[,3],col="blue",lty="dotted")
      lines(Ages, Growth_95PLs[,4],col="blue",lty="dotted")

      AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
      mtext(expression(paste(plain("Fish length (mm"),plain(")"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
      mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
      legend('bottomright', col=c("red","blue"),legend=c("Female","Male"),bty='n', cex=0.8,lwd=1.75)
  }

  # plot weight at length
  if (PlotOpt==0 | PlotOpt==2) {
    y1=max(Res$ModelDiag$FemWtAtLen)
    y2=max(Res$ModelDiag$MalWtAtLen)
    if (y1 > y2) {
      ylims = Get_yaxis_scale(c(0,Res$ModelDiag$FemWtAtLen))
    } else {
      ylims = Get_yaxis_scale(c(0,Res$ModelDiag$MalWtAtLen))
    }
    ymax = ylims$ymax; yint = ylims$yint
    xlims=Get_xaxis_scale(0:MaxLen)
    xmax = xlims$xmax; xint = xlims$xint

    plot(midpt,Res$ModelDiag$FemWtAtLen,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),col="red",yaxt="n",xaxt="n",
         ylab="",xlab="")
    lines(midpt,Res$ModelDiag$MalWtAtLen,col="blue")
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Total weight (kg"),plain(")"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Length (mm"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    legend('topleft', col=c("red","blue"),legend=c("Female","Male"),bty='n', cex=0.8,lwd=1.75)
  }

  # plot female maturity and selectivity at length
  if (PlotOpt==0 | PlotOpt==3) {
    xlims=Get_xaxis_scale(0:MaxLen)
    xmax = xlims$xmax; xint = xlims$xint
    plot(midpt,Res$ModelDiag$FemPropMatAtLen,"l", pch=16,frame.plot=F,ylim=c(0,1),xlim=c(0,xmax),col="red",yaxt="n",xaxt="n",
         ylab="",xlab="",cex=0.8)
    if (DiscMort == 0) {
      lines(midpt, Res$ModelDiag$FemSelLandAtLen, "l", col="red",lty="dotted", cex=0.8)
      legend('topright', col=c("red","red"),legend=c("Fem. mature","Fem. sel. land."),
             lty=c("solid","dotted"),bty='n', cex=0.8,lwd=1.75)
    } else {
      lines(midpt,Res$ModelDiag$FemSelLandAtLen, "l", col="red",lty="dotted",cex=0.8)
      lines(midpt,Res$ModelDiag$FemSelDiscAtLen, "l", col="brown",lty="dotted",cex=0.8)
      legend('topright', col=c("red","red","brown"),legend=c("Fem. mature","Fem. sel. land.","Fem. sel. disc."),
             lty=c("solid","dotted","dotted"),bty='n', cex=0.8,lwd=1.75)
    }
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax=1, yint=0.5, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Proportion"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Length (mm"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    if (ReprodPattern > 1) {
      lines(midpt,Res$ModelDiag$PropFemAtLen, "l", col="black",lty="solid",cex=0.8)
      legend('bottomright', col="black",legend="Prop. Fem.",
             lty="solid",bty='n', cex=0.8,lwd=1.75)
    }
  }

  # plot male maturity and selectivity at length
  if (PlotOpt==0 | PlotOpt==4) {
    xlims=Get_xaxis_scale(0:MaxLen)
    xmax = xlims$xmax; xint = xlims$xint
    plot(midpt, Res$ModelDiag$MalPropMatAtLen,"l", pch=16, frame.plot=F,ylim=c(0,1),xlim=c(0,xmax),col="blue",yaxt="n",xaxt="n",
         ylab="",xlab="", cex=0.8)
    if (DiscMort == 0) {
      lines(midpt, Res$ModelDiag$FemSelLandAtLen, "l", col="blue",lty="dotted", cex=0.8)
      legend('topright', col=c("blue","blue"),legend=c("Male mature","Male sel. land."),
             lty=c("solid","dotted"),bty='n', cex=0.8,lwd=1.75)
    } else {
      lines(midpt, Res$ModelDiag$FemSelLandAtLen, "l", col="blue",lty="dotted", cex=0.8)
      lines(midpt, Res$ModelDiag$FemSelDiscAtLen, "l", col="purple",lty="dotted", cex=0.8)
      legend('topright', col=c("blue","blue","purple"),legend=c("Male mature","Male sel. land.", "Male disc."),
             lty=c("solid","dotted","dotted"),bty='n', cex=0.8,lwd=1.75)
    }
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax=1, yint=0.5, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Proportion"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Length (mm"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    if (ReprodPattern > 1) {
      lines(midpt,1-Res$ModelDiag$PropFemAtLen, "l", col="black",lty="solid",cex=0.8)
      legend('bottomright', col="black",legend="Prop. Male",
             lty="solid",bty='n', cex=0.8,lwd=1.75)
    }
  }


  # plot female selectivity of landings, selectivity of discards, prob retention, and gear selectivity at age
  if (PlotOpt==0 | PlotOpt==5) {
    xlims = Get_xaxis_scale(midpt)
    xmax = xlims$xmax; xint = xlims$xint
    plot(midpt,Res$ModelDiag$FemSelLandAtLen,"l", pch=16,frame.plot=F,ylim=c(0,1),xlim=c(0,xmax),col="black",yaxt="n",xaxt="n",
         ylab="",xlab="",cex=0.8, lwd=2)
    lines(midpt,Res$ModelDiag$FemSelDiscAtLen, "l", col="purple",lty="dotted",cex=0.8)
    lines(midpt, Res$ModelDiag$FemRetProbAtLen, "l", col="brown",lty="dashed", cex=0.8)
    lines(midpt, Res$ModelDiag$FemGearSelAtLen, "l", col="orange",lty="dotted", cex=0.8)
    legend('bottomright', col=c("black","purple","brown","orange"),
           legend=c("Fem. sel_land","Fem. sel_disc","Fem. prob_ret","Fem. gear_sel"),
           lty=c("solid","dotted","dashed","dotted"),bty='n', cex=0.8,lwd=c(2,1,1,1))
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax=1, yint=0.5, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Proportion"))),las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Length (mm"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  }

  # plot male selectivity of landings, selectivity of discards, prob retention, and gear selectivity at age
  if (PlotOpt==0 | PlotOpt==6) {
    xlims = Get_xaxis_scale(midpt)
    xmax = xlims$xmax; xint = xlims$xint
    plot(midpt, Res$ModelDiag$MalSelLandAtLen,"l", pch=16, frame.plot=F,ylim=c(0,1),xlim=c(0,xmax),col="black",yaxt="n",xaxt="n",
         ylab="",xlab="", cex=0.8, lwd=2)
    lines(midpt,Res$ModelDiag$MalSelDiscAtLen, "l", col="purple",lty="dotted",cex=0.8)
    lines(midpt, Res$ModelDiag$MalRetProbAtLen, "l", col="brown",lty="dashed", cex=0.8)
    lines(midpt, Res$ModelDiag$MalGearSelAtLen, "l", col="orange",lty="dotted", cex=0.8)

    legend('bottomright', col=c("black","purple","brown","orange"),
           legend=c("Mal. sel_land","Mal. sel_disc","Mal. prob_ret","Mal. gear_sel"),
           lty=c("solid","dotted","dashed","dotted"),bty='n', cex=0.8,lwd=c(2,1,1,1))
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax=1, yint=0.5, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    axis(2,at=seq(0,1,0.5), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
    mtext(expression(paste(plain("Proportion"))),las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Length (mm"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  }


  #plot 2:
  if (PlotOpt==0) {
    par(mfrow = c(3,2), mar=c(3.5,4,2,2),
        oma=c(1,1,1,1),tck=-0.03,mgp = c(3, 0.5, 0))
  }

  # plot female mortality at length
  if (PlotOpt==0 | PlotOpt==7) {
    ylims = Get_yaxis_scale(Res$ModelDiag$FemFAtLen)
    ymax = ylims$ymax; yint = ylims$yint
    xlims=Get_xaxis_scale(0:MaxLen)
    xmax = xlims$xmax; xint = xlims$xint
    plot(midpt, Res$ModelDiag$FemFAtLen,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),col="red",yaxt="n",xaxt="n",
         ylab="",xlab="")
    # lines(Ages,Res$ModelDiag$FemZAtAge,lty="dotted",col="red")
    lines(midpt,Res$ModelDiag$FemDiscFAtLen,lty="dotted",col="brown")
    lines(midpt,Res$ModelDiag$FemLandFAtLen,lty="dotted",col="purple")
    lines(midpt,rep(NatMort,length(midpt)),lty="dashed")
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Mortality") ~ (year^{-1}))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Length (mm"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    legend('topright', col=c("red","brown","purple","black"),lty=c("solid","dotted","dotted","dashed"),
           legend=c("Fem. F","Fem. DiscF","Fem. LandF","Fem. M"),bty='n', cex=1.0,lwd=1.75)
  }

  # plot male mortality at age
  if (PlotOpt==0 | PlotOpt==8) {
    ylims = Get_yaxis_scale(Res$ModelDiag$MalFAtLen)
    ymax = ylims$ymax; yint = ylims$yint
    xlims=Get_xaxis_scale(0:MaxLen)
    xmax = xlims$xmax; xint = xlims$xint
    plot(midpt, Res$ModelDiag$MalFAtLen,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),col="blue",yaxt="n",xaxt="n",
         ylab="",xlab="")
    # lines(Ages, Res$ModelDiag$MalZAtAge,lty="dotted",col="blue")
    lines(midpt,Res$ModelDiag$MalDiscFAtLen,lty="dotted",col="brown")
    lines(midpt,Res$ModelDiag$MalLandFAtLen,lty="dotted",col="purple")
    lines(midpt,rep(NatMort,length(midpt)),lty="dashed")
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Mortality") ~ (year^{-1}))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Length (mm"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    legend('topright', col=c("blue","brown","purple","black"),lty=c("solid","dotted","dotted","dashed"),
           legend=c("Mal. F","Mal. DiscF","Mal.LandF","Mal. M"),bty='n', cex=1.0,lwd=1.75)
  }

  # plot fished and unfished female survival in terms of numbers given specified current fully-selected fishing mortality
  if (PlotOpt==0 | PlotOpt==9) {
    ylims = Get_yaxis_scale(Res$ModelDiag$Unfish_FemNPerRecLen)
    ymax = ylims$ymax; yint = ylims$yint
    if (ymax <= 0) {
      ymax = round(max(1.2*Res$ModelDiag$Unfish_FemNPerRecLen),1)
      yint = ymax
    }
    xlims=Get_xaxis_scale(0:MaxLen)
    xmax = xlims$xmax; xint = xlims$xint
    plot(midpt, Res$ModelDiag$Unfish_FemNPerRecLen,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),col="red",yaxt="n",xaxt="n",
         ylab="",xlab="")
    lines(midpt, Res$ModelDiag$Fish_FemNPerRecLen,col="red",lty="dotted")
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Rel. survival"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Length (mm"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    legend('topright', col="red",lty=c("solid","dotted"),
           legend=c("Fem. unfish","Fem. fish"),bty='n', cex=1.0,lwd=1.75)
  }

  # plot fished and unfished male survival in terms of numbers given specified current fully-selected fishing mortality
  if (PlotOpt==0 | PlotOpt==10) {
    ylims = Get_yaxis_scale(Res$ModelDiag$Unfish_MalNPerRecLen)
    ymax = ylims$ymax; yint = ylims$yint
    if (ymax <= 0) {
      ymax = round(max(1.2*Res$ModelDiag$Unfish_FemNPerRecLen),1)
      yint = ymax
    }
    xlims=Get_xaxis_scale(0:MaxLen)
    xmax = xlims$xmax; xint = xlims$xint
    plot(midpt, Res$ModelDiag$Unfish_MalNPerRecLen,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),col="blue",yaxt="n",xaxt="n",
         ylab="",xlab="")
    lines(midpt, Res$ModelDiag$Fish_MalNPerRecLen,col="blue",lty="dotted")
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Rel. survival"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Length (mm"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    legend('topright', col="blue",lty=c("solid","dotted"),
           legend=c("Mal. unfish","Mal. fish"),bty='n', cex=1.0,lwd=1.75)
  }

  # plot fished and unfished mature female biomass at age given specified current fully-selected fishing mortality
  if (PlotOpt==0 | PlotOpt==11) {
    ylims = Get_yaxis_scale(Res$ModelDiag$Unfish_FemBiomAtAge)
    ymax = ylims$ymax; yint = ylims$yint
    if (max(Res$ModelDiag$Unfish_FemBiomAtAge) < 0.2 * ymax) {
      ymax = ymax/5
      yint = yint/5
    }
    xlims=Get_xaxis_scale(c(0,Ages))
    xmax = xlims$xmax; xint = xlims$xint
    plot(Ages, Res$ModelDiag$Unfish_FemBiomAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),
         col="red",yaxt="n",xaxt="n",ylab="",xlab="")
    lines(Ages, Res$ModelDiag$Fish_FemBiomAtAge,col="red",lty="dotted")
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Biom. at age (kg"),plain(")"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    legend('topright', col=c("red","red"),lty=c("solid","dotted"),legend=c("Fem. unfish","Fem fish"),bty='n', cex=1.0,lwd=1.75)
  }

  # plot fished and unfished mature male biomass at age given specified current fully-selected fishing mortality
  if (PlotOpt==0 | PlotOpt==12) {
    ylims = Get_yaxis_scale(Res$ModelDiag$Unfish_MalBiomAtAge)
    ymax = ylims$ymax; yint = ylims$yint
    if (max(Res$ModelDiag$Unfish_MalBiomAtAge) < 0.2 * ymax) {
      ymax = ymax/5
      yint = yint/5
    }
    xlims=Get_xaxis_scale(c(0,Ages))
    xmax = xlims$xmax; xint = xlims$xint
    plot(Ages, Res$ModelDiag$Unfish_MalBiomAtAge,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),
         col="blue",yaxt="n",xaxt="n",ylab="",xlab="")
    lines(Ages, Res$ModelDiag$Fish_MalBiomAtAge,col="blue",lty="dotted")
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Biom. at age (kg"),plain(")"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Age (y"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    legend('topright', col=c("blue","blue"),lty=c("solid","dotted"),legend=c("Mal unfish","Mal. fish"),bty='n', cex=1.0,lwd=1.75)
  }

  #Plot 3:
  # plot female and male catch at age, given specified current fully-selected fishing mortality
  if (PlotOpt==0) {
  par(mfrow = c(3,2), mar=c(3.5,4,2,2),
      oma=c(1,1,1,1),tck=-0.03,mgp = c(3, 0.5, 0))
  }

  if (PlotOpt==0 | PlotOpt==13) {
    FemCatchNumAtLenProp <- Res$ModelDiag$FemCatchNumLen / sum(Res$ModelDiag$FemCatchNumLen)
    MalCatchNumAtLenProp <- Res$ModelDiag$MalCatchNumLen / sum(Res$ModelDiag$MalCatchNumLen)
    ylims = Get_yaxis_scale(c(FemCatchNumAtLenProp, MalCatchNumAtLenProp))
    ymax = ylims$ymax; yint = ylims$yint
    xlims = Get_xaxis_scale(midpt)
    xmax = xlims$xmax; xint = xlims$xint
    plot(midpt, FemCatchNumAtLenProp,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),
         col="red",yaxt="n",xaxt="n",ylab="",xlab="")
    lines(midpt, MalCatchNumAtLenProp,col="blue","l")
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext("Catch at length prop.",las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
    mtext(expression(paste(plain("Length (mm"),plain(")"))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    legend('topright', col=c("red","blue"),lty="solid",legend=c("females","males"),bty='n', cex=0.8,lwd=1.75)
  }

  # plot yield per recruit (per recruit analysis) and equilibrium catch (equilibrium age-structured model)
  # given specified current fully-selected fishing mortality
  if (PlotOpt==0 | PlotOpt==14) {
    ylims = Get_yaxis_scale(Res$YPRResults)
    ymax = ylims$ymax; yint = ylims$yint
    if (max(Res$YPRResults) < 0.2 * ymax) {
      ymax = ymax/5
      yint = yint/5
    }
    plot(Res$FishMort, Res$YPRResults,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,max(Res$FishMort)),
         col="black",yaxt="n",xaxt="n",ylab="",xlab="")
    points(Current_F, Res$YPR,cex=1.2,col="black",pch=16)
    lines(Res$FishMort, Res$Eq_CatchResults,col="blue")
    points(Current_F, Res$Eq_Catch, cex=1.2,col="blue",pch=16)
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax=max(Res$FishMort), xint=0.5, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    axis(2,at=seq(0,ymax,yint), cex.axis=0.8,line=0,las = 1,lwd=1.5,tick=F) #add y labels
    mtext(expression(paste(plain("YPR / Eq.Catch (kg"),plain(")"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
    mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    legend('topright', col=c("black","blue"),lty=c("solid","solid"),
           legend=c("YPR","Eq.catch"),bty='n', cex=0.8,lwd=1.75)
  }

  # plot spawning potential ratio (per recruit analysis) and relative spawning biomass (equilbrium age-structured model)
  # given specified current fully-selected fishing mortality, for female
  if (PlotOpt==0 | PlotOpt==15) {
    plot(Res$FishMort, Res$Fem_SPRResults,"l",frame.plot=F,ylim=c(0,1.0),xlim=c(0,max(Res$FishMort)),
         col="red",yaxt="n",xaxt="n",ylab="",xlab="", lty="dotted")
    lines(Res$FishMort, Res$Eq_FemRelSpBiomResults,col="red",lty="solid")
    points(Current_F, Res$Fem_SPR,cex=1.2,col="red",pch=16)
    points(Current_F, Res$Eq_FemRelSpBiom,cex=1.2,col="red",pch=1)
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax=max(Res$FishMort), xint=0.5, ymin=NA, ymax=1, yint=0.2, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Biom. ratio"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
    mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    # legend('topright', col=c("red","red"),lty=c("dotted","solid"),
    #        legend=c("Fem SPR","Fem. Rel biom"),bty='n', cex=0.8,lwd=1.75)
    if (!is.na(EggFertParam) &  EggFertParam < 1) {
      lines(Res$FishMort, Res$Eq_FemRelSpBiom_AllEggFert_Results, lty="dotted")
      points(Current_F, Res$Eq_FemRelSpBiom_AllEggFert, cex=1.2,col="black",pch=16)
      legend('topright', col=c("red","black","red"),lty=c("dotted","dotted","solid"),
             legend=c("Fem SPR","Fem. rel biom. (100% fert.)",
                      paste0("Fem. rel biom (calc. fert. rate)")),bty='n', cex=0.8,lwd=1.75)
    } else {
      legend('topright', col=c("red","red"),lty=c("dotted","solid"),
             legend=c("Fem SPR","Fem. Rel biom"),bty='n', cex=0.8,lwd=1.75)
    }
  }


  # plot spawning potential ratio (per recruit analysis) and relative spawning biomass (equilbrium age-structured model)
  # given specified current fully-selected fishing mortality, for each sex
  if (PlotOpt==0 | PlotOpt==16) {
    plot(Res$FishMort, Res$Mal_SPRResults,"l",frame.plot=F,ylim=c(0,1.0),xlim=c(0,max(Res$FishMort)),
         col="blue",yaxt="n",xaxt="n",ylab="",xlab="", lty="dotted")
    lines(Res$FishMort, Res$Eq_MalRelSpBiomResults,col="blue",lty="solid")
    points(Current_F, Res$Mal_SPR,cex=1.2,col="blue",pch=16)
    points(Current_F, Res$Eq_MalRelSpBiom,cex=1.2,col="blue",pch=1)
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax=max(Res$FishMort), xint=0.5, ymin=NA, ymax=1, yint=0.2, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Biom. ratio"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
    mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    legend('topright', col=c("blue","blue"),lty=c("dotted","solid"),
           legend=c("Mal SPR","Mal. Rel biom"),bty='n', cex=0.8,lwd=1.75)
  }

  if (PlotOpt==0 | PlotOpt==17) {
    plot(Res$FishMort, Res$CombSex_SPRResults,"l",frame.plot=F,ylim=c(0,1.0),xlim=c(0,max(Res$FishMort)),
         col="black",yaxt="n",xaxt="n",ylab="",xlab="", lty="dotted")
    lines(Res$FishMort, Res$Eq_CombSexRelSpBiomResults,col="black",lty="solid")
    points(Current_F, Res$CombSex_SPR,cex=1.2,col="black",pch=16)
    points(Current_F, Res$Eq_CombSexRelSpBiom,cex=1.2,col="black",pch=1)
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax=max(Res$FishMort), xint=0.5, ymin=NA, ymax=1, yint=0.2, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Biom. ratio"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
    mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
    legend('topright', col=c("black","black"),lty=c("dotted","solid"),
           legend=c("Comb. sex SPR","Comb. sex Rel biom"),bty='n', cex=0.8,lwd=1.75)
  }

  # plot equilibrium recruitment vs F
  if (PlotOpt==0 | PlotOpt==18) {
    ymax = 1.0
    yint = 0.2
    plot(Res$FishMort, Res$Eq_RecResults,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,max(Res$FishMort)),
         col="red",yaxt="n",xaxt="n",ylab="",xlab="")
    points(Current_F, Res$Eq_Rec, cex=1.2,col="red",pch=16)
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax=max(Res$FishMort), xint=0.5, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    mtext(expression(paste(plain("Equil. Recruitment"))),las=3,side=2,line=2,cex=0.7,lwd=1.75)
    mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
  }

  # plot 4
  if (!is.na(EggFertParam)) {
    if (PlotOpt==0) {
      par(mfrow = c(2,2), mar=c(3.5,4,2,2),
          oma=c(1,1,1,1),tck=-0.03,mgp = c(3, 0.5, 0))
    } else {
      # don't change user settings
    }

    # plot male depletion vs F
    # Male depletion measure is the ratio of the proportion of mature males:mature females in the population (in numbers)
    # at current fishing pressure relative to the proportion mature males: mature females at the unfished level
    if (PlotOpt==0 | PlotOpt==19) {
      ymax = 1; yint = 0.2
      xmax = max(Res$FishMort)
      plot(Res$FishMort,Res$Eq_MalDeplRatioResults,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),
           col="blue",yaxt="n",xaxt="n",ylab="",xlab="")
      points(Current_F, Res$MalDeplRatio, col="blue", pch=16)
      AddAxesAndTickLabelsToPlot(xmin=NA, xmax=max(Res$FishMort), xint=0.5, ymin=NA, ymax=1, yint=0.2, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
      mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
      mtext(expression("Male depletion (init. vs fin. sex ratio)"),las=3,side=2,line=2,cex=0.7,lwd=1.75)
      legend('topright', col="blue",pch=16, legend="Current F",bty='n', cex=0.8,lwd=-1)
    }

    # plot prop male vs F
    if (PlotOpt==0 | PlotOpt==20) {
      ylims = Get_yaxis_scale(Res$FishMalToFemPropResults)
      ymax = ylims$ymax; yint = ylims$yint
      xmax = max(Res$FishMort)
      x=which(Res$FishMort==Current_F)
      plot(Res$FishMort,Res$FishMalToFemPropResults,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),
           col="blue",yaxt="n",xaxt="n",ylab="",xlab="")
      points(Current_F, Res$FishMalToFemPropResults[x], col="blue", pch=16)
      AddAxesAndTickLabelsToPlot(xmin=NA, xmax=max(Res$FishMort), xint=0.5, ymin=NA, ymax=ymax, yint=yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
      mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
      mtext(expression("Prop. male to mat. fem."),las=3,side=2,line=2,cex=0.7,lwd=1.75)
      legend('topright', col="blue",pch=16, legend="Current F",bty='n', cex=0.8,lwd=-1)
    }

    # plot male depletion vs egg. fert rate
    if (PlotOpt==0 | PlotOpt==21) {
      ymax = 1; yint = 0.2
      xmax = 1; xint = 0.2
      plot(Res$FishMort, Res$Eq_FertRateResults,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,xmax),
           col="red",yaxt="n",xaxt="n",ylab="",xlab="")
      points(Res$FishMort[x], Res$Eq_FertRateResults[x], col="red", pch=16)
      AddAxesAndTickLabelsToPlot(xmin=NA, xmax=xmax, xint=xint, ymin=NA, ymax=ymax, yint=yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
      mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
      mtext(expression("Egg fert. rate"),las=3,side=2,line=2,cex=0.7,lwd=1.75)
      legend('topright', col="red",pch=16, legend="Current F",bty='n', cex=0.8,lwd=-1)
    }

    # plot male depletion vs eq. recruitment
    if (PlotOpt==0 | PlotOpt==22) {
      ymax = 1.0
      yint = 0.2
      plot(Res$FishMort, Res$Eq_RecResults,"l",frame.plot=F,ylim=c(0,ymax),xlim=c(0,max(Res$FishMort)),
           col="red",yaxt="n",xaxt="n",ylab="",xlab="")
      points(Current_F, Res$Eq_Rec, cex=1.2,col="red",pch=16)
      AddAxesAndTickLabelsToPlot(xmin=NA, xmax=max(Res$FishMort), xint=0.5, ymin=NA, ymax=1, yint=0.2, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
      mtext(expression(paste(plain("Equil. Recruitment"))),las=3,side=2,line=2.5,cex=0.7,lwd=1.75)
      mtext(expression(paste(italic("F") ~ (year^{-1}))),las=1,side=1,line=2,cex=0.7,lwd=1.75)
      if (!is.na(EggFertParam) &  EggFertParam < 1) {
        lines(Res$FishMort, Res$Eq_Rec_AllEggFertResults, lty="dotted")
        points(Current_F, Res$Eq_Rec_AllEggFert, cex=1.2,col="black",pch=16)
        legend('topright', col=c("black","red"),lty=c("dotted","solid"),
               legend=c("100% egg fert. rate","calc. egg fert. rate"),bty='n', cex=0.8,lwd=1.75)
      }
    }
  }



  # reset default par options
  par(.pardefault)

}

#' Get plots associated with expected sizes and weights of fish in catches
#' at different levels of mortality
#'
#' This function provides plots of expected size and weights of fish in catches,
#' for females, males and combined sexes. The expected values for a specified fishing mortality
#' are compared to other mortality values. 60 and 95 percent prediction intervals are plotted.
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
#' @param FinalSex_L50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param EstSexRatioAtLen NA  # sex ratio at length, inputted as vector
#' @param EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' @param mat_L50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_L95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtLen vector of proportion mature at length (set to NA if using age at maturity parameters)
#' @param sel_L50 logistic parameter for gear selectivity curve
#' @param sel_L95 logistic parameter for gear selectivity curve
#' @param EstGearSelAtLen gear selectivity curve inputted as vector
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_L50 logistic parameter for fish retention curve
#' @param ret_L95 logistic parameter for fish retention curve
#' @param EstRetenAtLen retention curve inputted as vector
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param Current_F level of fishing mortality of interest, e.g. current fishing mortality
#' @param FittedRes stored results of length-based per recruit analysis, can set to NA
#' @param PlotOpt 0=all plots, 1=female size, 2=male size, 3=combined sex size, 4=female weight,
#' 5=male weight, 6=combined sex weight
#' @param xmax maximum value for x axis
#' @param xint interval for x axis
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
#'                          EstMalWtAtLen=NA) # weight at length, inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at recruitment age/length
#' FinalSex_L50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_L95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtLen <- NA  # sex ratio at length inputted as vector
#' EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, ~0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_L50 <- c(250, 250) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_L95 <- c(300, 300) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
#'                           EstMalMatAtLen=NA) # maturity at length, inputted as values in data frame
#' sel_L50 <- c(250, 250) # females, males - Logistic length selectivity relationship parameters
#' sel_L95 <- c(300, 300) # females, males - Logistic length selectivity relationship parameters
#' EstGearSelAtLen <- data.frame(EstFemGearSelAtLen=NA,
#'                               EstMalGearSelAtLen=NA)
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_L50 <- c(250, 250) # females, males - Logistic fish retention at length parameters
#' ret_L95 <- c(300, 300) # females, males - Logistic fish retention at length parameters
#' EstRetenAtLen <- data.frame(EsFemtRetenAtLen=NA,
#'                             EstMalRetenAtLen=NA)
#' DiscMort <- 0.25 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort <- 4.22 / MaxModelAge # natural mortality  (year-1)
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' Current_F <- 0.2
#' PlotOpt <- 0 # 0=all plots, 1=female size, 2=male size, 3=combined sex size, 4=female weight, 5=male weight, 6=combined sex weight
#' Output_Opt = 2 # including additional length and weight outputs
#' FittedRes=GetPerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
#'                             RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
#'                             ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
#'                             EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
#'                             ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, Current_F, Output_Opt)
#' PlotPerRecruit_ExpCatchSizeDistns_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
#'                                      RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
#'                                      ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
#'                                      EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
#'                                      ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, Current_F,
#'                                      FittedRes, PlotOpt, xmax=NA, xint=NA)
#' @export
PlotPerRecruit_ExpCatchSizeDistns_LB <- function(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                                 RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                                                 ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                                                 EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                                                 ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, Current_F,
                                                 FittedRes, PlotOpt, xmax=NA, xint=NA) {

  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    Res =  FittedRes
  } else {
    Output_Opt = 2 # 1=standard output, 2=with added length and weight outputs (slower)
    Res=GetPerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                                ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                                EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                                ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, Current_F, Output_Opt)
  }
  FishMort=Res$FishMort

  if (is.na(xmax)) xmax = 2
  if (is.na(xint)) xint = 0.5

  # get relationship between fishing mortality and mean length of females in catches
  if (PlotOpt==0) {
    par(mfrow=c(2,3),mar=c(5,4,2,2))
  }
  xaxis_lab = expression(paste("Fishing mortality," ~ y^{-1}))

  # females - lengths
  if (PlotOpt==0 | PlotOpt==1) {
    yaxis_lab = "Female lengths, mm"
    ylims = Get_yaxis_scale(c(Res$MeanCatchLenResults$FemCatchLen.97.5qntl,
                              Res$MeanCatchLenResults$MalCatchLen.97.5qntl))
    ymax = ylims$ymax; yint=ylims$yint
    plot(FishMort, Res$MeanCatchLenResults[,2], "l", ylim=c(0,ymax), xlim=c(0,xmax), ylab = yaxis_lab, xlab = xaxis_lab,
         col="red", bty='n', xaxt='n', yaxt='n')
    x = c(FishMort,rev(FishMort))
    y = c(Res$MeanCatchLenResults$FemCatchLen.2.5qntl,
          rev(Res$MeanCatchLenResults$FemCatchLen.97.5qntl))
    polygon(x,y, col="pink",border=NA)
    y = c(Res$MeanCatchLenResults$FemCatchLen.20qntl,
          rev(Res$MeanCatchLenResults$FemCatchLen.80qntl))
    polygon(x,y, col="light blue",border=NA)
    lines(FishMort, Res$MeanCatchLenResults[,2],col="red")
    x=min(which(FishMort>=Current_F))
    points(Current_F, Res$MeanCatchLenResults[x,2], pch=16, cex=2, col="red")
    legend('topright', bty='n', cex=0.8, lwd=c(5,5,-1), pch = c(NA,NA,16), col=c("light blue","pink","red"),
           legend=c("20 & 80 quantiles","2.5 & 97.5 quantiles","Current F"))
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
  }

  # males - lengths
  if (PlotOpt==0 | PlotOpt==2) {
    yaxis_lab = "Male lengths, mm"
    ylims = Get_yaxis_scale(c(Res$MeanCatchLenResults$FemCatchLen.97.5qntl,
                              Res$MeanCatchLenResults$MalCatchLen.97.5qntl))
    ymax = ylims$ymax; yint=ylims$yint
    plot(FishMort, Res$MeanCatchLenResults[,3], "l", ylim=c(0,ymax), xlim=c(0,xmax), ylab = yaxis_lab,
         xlab = xaxis_lab, col="red", bty='n', xaxt='n', yaxt='n')
    x = c(FishMort,rev(FishMort))
    y = c(Res$MeanCatchLenResults$MalCatchLen.2.5qntl,
          rev(Res$MeanCatchLenResults$MalCatchLen.97.5qntl))
    polygon(x,y, col="pink",border=NA)
    y = c(Res$MeanCatchLenResults$MalCatchLen.20qntl,
          rev(Res$MeanCatchLenResults$MalCatchLen.80qntl))
    polygon(x,y, col="light blue",border=NA)
    lines(FishMort, Res$MeanCatchLenResults[,3],col="blue")
    x=min(which(FishMort>=Current_F))
    points(Current_F, Res$MeanCatchLenResults[x,3], pch=16, cex=2, col="blue")
    legend('topright', bty='n', cex=0.8, lwd=c(5,5,-1), pch = c(NA,NA,16), col=c("light blue","pink","blue"),
           legend=c("20 & 80 quantiles","2.5 & 97.5 quantiles","Current F"))
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
  }

  # combined sex - lengths
  if (PlotOpt==0 | PlotOpt==3) {
    yaxis_lab = "Comb. sex lengths, mm"
    ylims = Get_yaxis_scale(c(Res$MeanCatchLenResults$CombSexCatchLen.2.5qntl,
                              Res$MeanCatchLenResults$CombSexCatchLen.97.5qntl))
    ymax = ylims$ymax; yint=ylims$yint
    plot(FishMort, Res$MeanCatchLenResults[,4], "l", ylim=c(0,ymax), xlim=c(0,xmax), ylab = yaxis_lab,
         xlab = xaxis_lab, col="red", bty='n', xaxt='n', yaxt='n')
    x = c(FishMort,rev(FishMort))
    y = c(Res$MeanCatchLenResults$CombSexCatchLen.2.5qntl,
          rev(Res$MeanCatchLenResults$CombSexCatchLen.97.5qntl))
    polygon(x,y, col="pink",border=NA)
    y = c(Res$MeanCatchLenResults$CombSexCatchLen.20qntl,
          rev(Res$MeanCatchLenResults$CombSexCatchLen.80qntl))
    polygon(x,y, col="light blue",border=NA)
    lines(FishMort, Res$MeanCatchLenResults[,3],col="orange")
    x=min(which(FishMort>=Current_F))
    points(Current_F, Res$MeanCatchLenResults[x,3], pch=16, cex=2, col="orange")
    legend('topright', bty='n', cex=0.8, lwd=c(5,5,-1), pch = c(NA,NA,16), col=c("light blue","pink","orange"),
           legend=c("20 & 80 quantiles","2.5 & 97.5 quantiles","Current F"))
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
  }

  # females - weights
  if (PlotOpt==0 | PlotOpt==4) {
    if(max(c(Res$MeanCatchWtResults$FemCatchWt.97.5qntl,
             Res$MeanCatchWtResults$MalCatchWt.97.5qntl)<1)) {
      ScaleFact = 1000
      ylabel = "Female weights, g"
    } else {
      ScaleFact = 1
      ylabel = "Female weights, kg"
    }
    ylims = Get_yaxis_scale(c(Res$MeanCatchWtResults$FemCatchWt.97.5qntl,
                            Res$MeanCatchWtResults$MalCatchWt.97.5qntl))
    ymax = ylims$ymax*ScaleFact; yint=ylims$yint*ScaleFact
    plot(FishMort, Res$MeanCatchWtResults[,2]*ScaleFact, "l", ylim=c(0,ymax), xlim=c(0,xmax), ylab = ylabel,
         xlab = xaxis_lab, col="red", bty='n', xaxt='n', yaxt='n')
    x = c(FishMort,rev(FishMort))
    y = c(Res$MeanCatchWtResults$FemCatchWt.2.5qntl*ScaleFact,
          rev(Res$MeanCatchWtResults$FemCatchWt.97.5qntl*ScaleFact))
    polygon(x,y, col="pink",border=NA)
    y = c(Res$MeanCatchWtResults$FemCatchWt.20qntl*ScaleFact,
          rev(Res$MeanCatchWtResults$FemCatchWt.80qntl*ScaleFact))
    polygon(x,y, col="light blue",border=NA)
    lines(FishMort, Res$MeanCatchWtResults[,2]*ScaleFact,col="red")
    x=min(which(FishMort>=Current_F))
    points(Current_F, Res$MeanCatchWtResults[x,2]*ScaleFact, pch=16, cex=2, col="red")
    # abline(h=100,lwd=2, lty="dotted")
    legend('topright', bty='n', cex=0.8, lwd=c(5,5,-1), pch = c(NA,NA,16), col=c("light blue","pink","red"),
           legend=c("20 & 80 quantiles","2.5 & 97.5 quantiles","Current F"))
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
  }

  # males - weights
  if (PlotOpt==0 | PlotOpt==5) {
    if(max(c(Res$MeanCatchWtResults$FemCatchWt.97.5qntl,
             Res$MeanCatchWtResults$MalCatchWt.97.5qntl)<1)) {
      ScaleFact = 1000
      ylabel = "Male weights, g"
    } else {
      ScaleFact = 1
      ylabel = "Male weights, kg"
    }
    ylims = Get_yaxis_scale(c(Res$MeanCatchWtResults$FemCatchWt.97.5qntl,
                              Res$MeanCatchWtResults$MalCatchWt.97.5qntl))
    ymax = ylims$ymax*ScaleFact; yint=ylims$yint*ScaleFact
    plot(FishMort, Res$MeanCatchWtResults[,3]*ScaleFact, "l", ylim=c(0,ymax), xlim=c(0,xmax), ylab = ylabel,
         xlab = xaxis_lab, col="red", bty='n', xaxt='n', yaxt='n')
    x = c(FishMort,rev(FishMort))
    y = c(Res$MeanCatchWtResults$MalCatchWt.2.5qntl*ScaleFact,
          rev(Res$MeanCatchWtResults$MalCatchWt.97.5qntl*ScaleFact))
    polygon(x,y, col="pink",border=NA)
    y = c(Res$MeanCatchWtResults$MalCatchWt.20qntl*ScaleFact,
          rev(Res$MeanCatchWtResults$MalCatchWt.80qntl*ScaleFact))
    polygon(x,y, col="light blue",border=NA)
    lines(FishMort, Res$MeanCatchWtResults[,3]*ScaleFact,col="blue")
    x=min(which(FishMort>=Current_F))
    points(Current_F, Res$MeanCatchWtResults[x,3]*ScaleFact, pch=16, cex=2, col="blue")
    # abline(h=100,lwd=2, lty="dotted")
    legend('topright', bty='n', cex=0.8, lwd=c(5,5,-1), pch = c(NA,NA,16), col=c("light blue","pink","blue"),
           legend=c("20 & 80 quantiles","2.5 & 97.5 quantiles","Current F"))
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
  }

  # combined sex - weights
  if (PlotOpt==0 | PlotOpt==6) {
    if(max(c(Res$MeanCatchWtResults$FemCatchWt.97.5qntl,
             Res$MeanCatchWtResults$MalCatchWt.97.5qntl)<1)) {
      ScaleFact = 1000
      ylabel = "Comb. sex weights, g"
    } else {
      ScaleFact = 1
      ylabel = "Comb. sex weights, kg"
    }
    ylims = Get_yaxis_scale(c(Res$MeanCatchWtResults$FemCatchWt.97.5qntl,
                              Res$MeanCatchWtResults$MalCatchWt.97.5qntl))
    ymax = ylims$ymax*ScaleFact; yint=ylims$yint*ScaleFact
    plot(FishMort, Res$MeanCatchWtResults[,4]*ScaleFact, "l", ylim=c(0,ymax), xlim=c(0,xmax), ylab = "Comb. sex weight, g",
         xlab = xaxis_lab, col="red", bty='n', xaxt='n', yaxt='n')
    x = c(FishMort,rev(FishMort))
    y = c(Res$MeanCatchWtResults$CombSexCatchWt.2.5qntl*ScaleFact,
          rev(Res$MeanCatchWtResults$CombSexCatchWt.97.5qntl*ScaleFact))
    polygon(x,y, col="pink",border=NA)
    y = c(Res$MeanCatchWtResults$CombSexCatchWt.20qntl*ScaleFact,
          rev(Res$MeanCatchWtResults$CombSexCatchWt.80qntl*ScaleFact))
    polygon(x,y, col="light blue",border=NA)
    lines(FishMort, Res$MeanCatchWtResults[,4]*ScaleFact,col="orange")
    x=min(which(FishMort>=Current_F))
    points(Current_F, Res$MeanCatchWtResults[x,4]*ScaleFact, pch=16, cex=2, col="orange")
    # abline(h=100,lwd=2, lty="dotted")
    legend('topright', bty='n', cex=0.8, lwd=c(5,5,-1), pch = c(NA,NA,16), col=c("light blue","pink","orange"),
           legend=c("20 & 80 quantiles","2.5 & 97.5 quantiles","Current F"))
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
  }
}

#' Plot outputs from length-based value per recruit analysis
#'
#' This function provides several plots associated with outputs from length-based 'value per recruit' analysis.
#' The relative value of fish at different lengths is specified according to a logistic curve, with specified
#' input parameter values.
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
#' @param FinalSex_L50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param EstSexRatioAtLen NA  # sex ratio at length, inputted as vector
#' @param EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' @param mat_L50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_L95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtLen vector of proportion mature at length (set to NA if using age at maturity parameters)
#' @param sel_L50 logistic parameter for gear selectivity curve
#' @param sel_L95 logistic parameter for gear selectivity curve
#' @param EstGearSelAtLen gear selectivity curve inputted as vector
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_L50 logistic parameter for fish retention curve
#' @param ret_L95 logistic parameter for fish retention curve
#' @param EstRetenAtLen retention curve inputted as vector
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param Current_F fishing mortality
#' @param FittedRes results from value per recruit analysis, can set to NA
#' @param PlotOpt = 0 # 0=all plots, 1=rel value vs length, 2=rel value vs weight, 3=eq catch and ypr, 4=rel value per recruit
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
#' FinalSex_L50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_L95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtLen <- NA  # sex ratio at length inputted as vector
#' EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, ~0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_L50 <- c(250, 250) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_L95 <- c(300, 300) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
#'                           EstMalMatAtLen=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_L50 <- c(250, 250) # females, males - Logistic age selectivity relationship parameters
#' sel_L95 <- c(300, 300) # females, males - Logistic age selectivity relationship parameters
#' EstGearSelAtLen <- data.frame(EstFemGearSelAtLen=NA,
#'                               EstMalGearSelAtLen=NA)
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_L50 <- c(250, 250) # females, males - Logistic age fish retention at age parameters
#' ret_L95 <- c(300, 300) # females, males - Logistic age fish retention at age parameters
#' EstRetenAtLen <- data.frame(EsFemtRetenAtLen=NA,
#'                             EstMalRetenAtLen=NA)
#' DiscMort <- 0.25 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 4.22 / MaxModelAge # natural mortality  (year-1)
#' Current_F = 0.2
#' MinFishWtOfVal = 100 # size of fish below which fish has no value
#' ValScaleFact = 1.3 # amount by which value scales according to size, above minimum size of value
#' PlotOpt = 0 # 0=all plots, 1=rel value vs length, 2=rel value vs weight, 3=eq catch and ypr, 4=rel value per recruit
#' FittedRes = Get_Relative_Value_Per_Recruit_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams, RefnceAges,
#'             CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen, ReprodScale, ReprodPattern,
#'             InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen, EggFertParam, mat_L50, mat_L95,
#'             EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax, ret_L50, ret_L95, EstRetenAtLen,
#'             DiscMort, Steepness, SRrel_Type, NatMort, Current_F, MinFishWtOfVal, ValScaleFact)
#' PlotValuePerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
#'                               RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
#'                               ReprodScale, ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge,
#'                               EggFertParam,mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
#'                               ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, Current_F,
#'                               MinFishWtOfVal, ValScaleFact, FittedRes, PlotOpt)
#' @export
PlotValuePerRecruitResults_LB <- function(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                          RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                                          ReprodScale, ReprodPattern, InitRatioFem, FinalSex_A50, FinalSex_A95, EstSexRatioAtAge,
                                          EggFertParam,mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                                          ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, Current_F,
                                          MinFishWtOfVal, ValScaleFact, FittedRes, PlotOpt) {


  if (is.list(FittedRes)) {
    Res =  FittedRes
  } else {
    Res = Get_Relative_Value_Per_Recruit_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams, RefnceAges,
                                            CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen, ReprodScale, ReprodPattern,
                                            InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen, EggFertParam, mat_L50, mat_L95,
                                            EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax, ret_L50, ret_L95, EstRetenAtLen,
                                            DiscMort, Steepness, SRrel_Type, NatMort, Current_F, MinFishWtOfVal, ValScaleFact)
  }

  RelValAtLen=Res$RelValAtLen
  RelValPerRecAtFMort=Res$RelValPerRecAtFMort
  OptFishMort=Res$OptFishMort
  MaxValPerRec=Res$MaxValPerRec
  Fem_SPRResults=Res$Fem_SPRResults
  Eq_FemRelSpBiomResults=Res$Eq_FemRelSpBiomResults
  YPRResults=Res$YPRResults
  Eq_CatchResults=Res$Eq_CatchResults
  MinFishLenOfVal = Res$MinFishLenOfVal
  FMort = Current_F
  Res=CalcYPRAndSPRForFMort_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                               RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                               ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                               EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                               ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, FMort)

  WtAtLen = (Res$ModelDiag$FemWtAtLen + Res$ModelDiag$MalWtAtLen)/2



  # plot relationship between fish length and relative value of fish, considered by fishers
  if (PlotOpt==0) {
    par(mfrow=c(2,2))
  }

  if (PlotOpt==0 | PlotOpt==1) {
    xaxis_lab = "Fish length, mm"
    yaxis_lab = "Relative value"
    xlims = Get_xaxis_scale(midpt)
    xmax = max(ubnd)
    xint = xlims$xint
    ylims = Get_yaxis_scale(RelValAtLen)
    ymax = ylims$ymax
    yint = ylims$yint
    plot(midpt,RelValAtLen,"p",pch=16,xlab=xaxis_lab, ylab=yaxis_lab,yaxt="n",xaxt="n",bty="n",
         xlim=c(0,xmax),ylim=c(0,ymax),cex=0.6)
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    # abline(v=Res$MinFishLenOfVal,lty="dotted")
    legend('topright', legend=paste("Min. len. of val. =",round(MinFishLenOfVal,0),"mm"),
           lty=-1,bty='n', cex=0.8, pch=-1)
  }

  # plot relationship between fish length and relative value of fish, considered by fishers
  if (PlotOpt==0 | PlotOpt==2) {
    xaxis_lab = "Fish weight, kg"
    yaxis_lab = "Relative value"
    xlims = Get_xaxis_scale(WtAtLen)
    xmax = xlims$xmax
    xint = xlims$xint
    ylims = Get_yaxis_scale(RelValAtLen)
    ymax = ylims$ymax
    yint = ylims$yint
    plot(WtAtLen,RelValAtLen,"p",pch=16,xlab=xaxis_lab, ylab=yaxis_lab,yaxt="n",xaxt="n",bty="n",
         xlim=c(0,xmax),ylim=c(0,ymax),cex=0.6)
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    # abline(v=MinFishWtOfVal/1000,lty="dotted")
    legend('topright', legend=paste("Min. wt. of val. =",MinFishWtOfVal/1000,"kg"),
           lty=-1,bty='n', cex=0.8, pch=-1)
  }

  # plot relationship between fishing mortality and equilbrium catch
  if (PlotOpt==0 | PlotOpt==3) {
    xaxis_lab = "Fishing mortality"
    yaxis_lab = "YPR/Eq. catch"
    FishMort = seq(0,2,0.01)
    xlims = Get_xaxis_scale(FishMort)
    xmax = 2.0
    xint = xlims$xint
    Eq_Catch = Res$Eq_Catch
    ylims = Get_yaxis_scale(YPRResults)
    ymax = ylims$ymax
    yint = ylims$yint
    plot(FishMort,Eq_CatchResults,"l", xlab=xaxis_lab, ylab=yaxis_lab, yaxt="n",xaxt="n",bty="n",
         xlim=c(0,xmax),ylim=c(0,ymax))
    x=which(FishMort==OptFishMort)
    points(OptFishMort,Eq_CatchResults[x],pch=16)
    lines(FishMort,YPRResults,lty="dotted")
    points(OptFishMort,YPRResults[x],pch=16)
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    legend('topright', legend=c(bquote("Opt. F =" ~ .(OptFishMort) ~ y^-1),"Eq. Catch","YPR"),
           lty=c(NA,"solid","dotted"),bty='n', cex=0.8, pch=c(16,-1,-1))
  }

  # plot relationship between fishing mortality and relative value per recruit
  if (PlotOpt==0 | PlotOpt==4) {
    xaxis_lab = "Fishing mortality"
    yaxis_lab = "Relative value per recruit"
    FishMort = seq(0,2,0.01)
    xlims = Get_xaxis_scale(FishMort)
    xmax = xlims$xmax
    xint = xlims$xint
    ylims = Get_yaxis_scale(RelValPerRecAtFMort)
    ymax = ylims$ymax
    yint = ylims$yint
    plot(FishMort,RelValPerRecAtFMort,"l", xlab=xaxis_lab, ylab=yaxis_lab, yaxt="n",xaxt="n",bty="n",
         xlim=c(0,xmax),ylim=c(0,ymax))
    points(OptFishMort,MaxValPerRec,pch=16)
    AddAxesAndTickLabelsToPlot(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=0.8, lwdval=1.5, lineval=NA, lasval=NA)
    legend('topright', legend=bquote("Opt. F =" ~ .(OptFishMort) ~ y^-1), bty='n', cex=0.8, pch=16)
  }
}

#' Get outputs from length-based value per recruit analysis
#'
#' This function provides outputs associated with a 'value per recruit' analysis, calculated based on
#' length-based per recruit analysis (with a stock recruitment relationship to account
#' for potential impacts of fishing on recruitment). The relative value of fish at different lengths
#' is specified according to a logistic curve, with specified input parameter values.
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
#' @param FinalSex_L50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param EstSexRatioAtLen NA  # sex ratio at length, inputted as vector
#' @param EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' @param mat_L50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_L95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtLen vector of proportion mature at length (set to NA if using age at maturity parameters)
#' @param sel_L50 logistic parameter for gear selectivity curve
#' @param sel_L95 logistic parameter for gear selectivity curve
#' @param EstGearSelAtLen gear selectivity curve inputted as vector
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_L50 logistic parameter for fish retention curve
#' @param ret_L95 logistic parameter for fish retention curve
#' @param EstRetenAtLen retention curve inputted as vector
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param FMort fishing mortality
#'
#' @return relative value of fish vs length to fishers (RelValAtLen), relative value per recruit to fishers
#' of catch for specified fishing mortality (RelValPerRecAtFMort), maximum value per recruit to fishers
#' across range of fishing mortalities (MaxValPerRec), level of fishing mortality required to achieve MaxValPerRec.
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
#' FinalSex_L50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_L95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtLen <- NA  # sex ratio at length inputted as vector
#' EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, ~0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_L50 <- c(250, 250) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_L95 <- c(300, 300) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
#'                           EstMalMatAtLen=NA) # maturity at age (from age 0), inputted as values in data frame
#' sel_L50 <- c(250, 250) # females, males - Logistic age selectivity relationship parameters
#' sel_L95 <- c(300, 300) # females, males - Logistic age selectivity relationship parameters
#' EstGearSelAtLen <- data.frame(EstFemGearSelAtLen=NA,
#'                               EstMalGearSelAtLen=NA)
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_L50 <- c(250, 250) # females, males - Logistic age fish retention at age parameters
#' ret_L95 <- c(300, 300) # females, males - Logistic age fish retention at age parameters
#' EstRetenAtLen <- data.frame(EsFemtRetenAtLen=NA,
#'                             EstMalRetenAtLen=NA)
#' DiscMort <- 0.25 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 4.22 / MaxModelAge # natural mortality  (year-1)
#' Current_F <- 0.2
#' MinFishWtOfVal = 100 # size of fish below which fish has no value
#' ValScaleFact = 1.2 # amount by which value scales according to size, above mininum size of value
#' FittedRes = Get_Relative_Value_Per_Recruit_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams, RefnceAges,
#'                                               CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen, ReprodScale, ReprodPattern,
#'                                               InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen, EggFertParam, mat_L50, mat_L95,
#'                                               EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax, ret_L50, ret_L95, EstRetenAtLen,
#'                                               DiscMort, Steepness, SRrel_Type, NatMort, Current_F, MinFishWtOfVal, ValScaleFact)
#' @export
Get_Relative_Value_Per_Recruit_LB <- function(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams, RefnceAges,
                                              CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen, ReprodScale, ReprodPattern,
                                              InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen, EggFertParam, mat_L50, mat_L95,
                                              EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax, ret_L50, ret_L95, EstRetenAtLen,
                                              DiscMort, Steepness, SRrel_Type, NatMort, Current_F, MinFishWtOfVal, ValScaleFact) {

  # get weight length relationships
  Res=CalcYPRAndSPRForFMort_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                               RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                               ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                               EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                               ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, FMort=0.01)


  if (is.na(EstWtAtLen[1,1])) {
    if (WLrel_Type == 1) { # power relationship
      # already have params needed
    }
    if (WLrel_Type == 2) { # log-log relationship
      lenwt_a = exp(ln_lenwt_a)
    }
  }
  if (!is.na(EstWtAtLen[1,1])) {
    # average calculated female and male weights at lengths in case these are a different, to simplify analysis
    WtAtLen = 1000*((Res$MalWtAtLen + Res$FemWtAtLen)/2) # in g
    ln_Wt = log(WtAtLen)
    ln_len = log(midpt)
    mod1=lm(ln_Wt~ln_len)
    lenwt_a=(exp(mod1$coefficients[1]))
    lenwt_b=mod1$coefficients[2]
  }

  # calculate length at minimum weight considered to be of value, e.g. filletable size
  MinFishLenOfVal = (MinFishWtOfVal/lenwt_a)^(1/lenwt_b)

  # allocate a 'relative value' of a fish, to fishers, according to its size
  ValAtLen = MinFishLenOfVal+((midpt-MinFishLenOfVal)^ValScaleFact)
  x=which(midpt<MinFishLenOfVal)
  ValAtLen[x]=0
  # calculate value according to, relative to smallest size of fish considered to be of value
  RelValAtLen = ValAtLen / midpt
  # plot(midpt,RelValAtLen,"l")

  FishMort <- seq(0,2,0.01)
  nFVals = length(FishMort)
  RelValPerRecAtFMort = rep(NA,nFVals)
  Fem_SPRResults = rep(NA,nFVals)
  Eq_FemRelSpBiomResults = rep(NA,nFVals)
  YPRResults = rep(NA,nFVals)
  Eq_CatchResults = rep(NA,nFVals)

  for (i in 1:nFVals) {
    FMort=FishMort[i]
    Res=CalcYPRAndSPRForFMort_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                 RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                                 ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                                 EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                                 ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, FMort)

    RelValPerRecAtFMort[i] = sum(Res$ModelDiag$CombSexCatchBiomLen*RelValAtLen*Res$Eq_Rec)

    # for comparison, also save female SPR and relbiom, and ypr and eq catch
    Fem_SPRResults[i] = Res$Fem_SPR
    Eq_FemRelSpBiomResults[i] = Res$Eq_FemRelSpBiom
    YPRResults[i] = Res$YPR
    Eq_CatchResults[i] = Res$Eq_Catch

    cat("F",FishMort[i],"RelValPerRecAtFMort",RelValPerRecAtFMort[i], '\n')
  }
  x=which(RelValPerRecAtFMort==max(RelValPerRecAtFMort))
  OptFishMort = FishMort[x]
  MaxValPerRec = RelValPerRecAtFMort[x]


  Results = list(RelValAtLen = RelValAtLen,
                 RelValPerRecAtFMort = RelValPerRecAtFMort,
                 MaxValPerRec = MaxValPerRec,
                 OptFishMort = OptFishMort,
                 Fem_SPRResults = Fem_SPRResults,
                 Eq_FemRelSpBiomResults = Eq_FemRelSpBiomResults,
                 YPRResults = YPRResults,
                 Eq_CatchResults = Eq_CatchResults,
                 MinFishLenOfVal = MinFishLenOfVal)

  return(Results)

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
#' @param FinalSex_A50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_A95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param EstSexRatioAtAge NA  # sex ratio at age (from age 0) inputted as vector
#' @param EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' @param mat_A50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_A95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtAge vector of proportion mature at age (set to NA if using age at maturity parameters)
#' @param Gear_sel_A50 logistic parameter for gear selectivity curve
#' @param Gear_sel_A95 logistic parameter for gear selectivity curve
#' @param EstGearSelAtAge vector for gear selectivity at age (set to NA if using age at gear selectivity parameters)
#' @param Land_sel_A50 logistic parameter for selectivity of landings curve
#' @param Land_sel_A95 logistic parameter for selectivity of landings curve
#' @param EstGearSelAtAge vector of gear selectivity at age (set to NA if using age at gear selectivity parameters)
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_A50 logistic parameter for fish retention curve
#' @param ret_A95 logistic parameter for fish retention curve
#' @param EstRetenAtAge vector of fish retention at age (set to NA if using age at fish retention parameters)
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param PlotOpt # 1=females, 2=males, 3=combined sex
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
#' FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtAge <- NA  # sex ratio at age (from age 0) inputted as vector
#' EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' Gear_sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' Gear_sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' EstGearSelAtAge <- data.frame(FemGearSelAtAge=NA, MalGearSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' Land_sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' Land_sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' EstLandSelAtAge <- data.frame(FemSelLandAtAge=NA, MalSelLandAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
#' EstRetenAtAge <- data.frame(EstFemRetenAtAge=rep(1,MaxModelAge+1), EstMalRetenAtAge=rep(1,MaxModelAge+1)) # fish retention at age (from age 0), inputted as values in data frame
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort <- 0.2 # natural mortality  (year-1)
#' Current_F <- 0.2 # estimate of fishing mortality, e.g. from catch curve analysis
#' PlotOpt <- 1 # 1=females, 2=males, 3=combined sex
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' PlotPerRecruit_Biom_no_err_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge, lenwt_a, ln_lenwt_a,
#'                               lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern, InitRatioFem,
#'                               FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam, mat_A50, mat_A95,
#'                               EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge, Land_sel_A50, Land_sel_A95,
#'                               EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness,
#'                               SRrel_Type, NatMort, PlotOpt, RefPointPlotOpt, Current_F)
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
#' FinalSex_A50 <- 35 # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- 60 # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtAge <- NA  # sex ratio at age (from age 0) inputted as vector
#' EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_A50 <- c(20, 20) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(30, 30) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' Gear_sel_A50 <- c(20, 20) # females, males - Logistic age selectivity relationship parameters
#' Gear_sel_A95 <- c(30, 30) # females, males - Logistic age selectivity relationship parameters
#' EstGearSelAtAge <- data.frame(FemGearSelAtAge=NA, MalGearSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' Land_sel_A50 <- c(25, 25) # females, males - Logistic age selectivity relationship parameters
#' Land_sel_A95 <- c(35, 35) # females, males - Logistic age selectivity relationship parameters
#' EstLandSelAtAge <- data.frame(FemSelLandAtAge=NA, MalSelLandAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
#' EstRetenAtAge <- data.frame(EstFemRetenAtAge=rep(1,MaxModelAge+1), EstMalRetenAtAge=rep(1,MaxModelAge+1)) # fish retention at age (from age 0), inputted as values in data frame
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort <- 0.07 # natural mortality  (year-1)
#' Current_F <- 0.1 # estimate of fishing mortality, e.g. from catch curve analysis
#' PlotOpt <- 1 # 1=females, 2=males, 3=combined sex
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' PlotPerRecruit_Biom_no_err_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge, lenwt_a, ln_lenwt_a,
#'                               lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern, InitRatioFem,
#'                               FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam, mat_A50, mat_A95,
#'                               EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge, Land_sel_A50, Land_sel_A95,
#'                               EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness,
#'                               SRrel_Type, NatMort, PlotOpt, RefPointPlotOpt, Current_F)
#' @export
PlotPerRecruit_Biom_no_err_AB <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge, lenwt_a, ln_lenwt_a,
                                          lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern, InitRatioFem,
                                          FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam, mat_A50, mat_A95,
                                          EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge, Land_sel_A50, Land_sel_A95,
                                          EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness,
                                          SRrel_Type, NatMort, PlotOpt, RefPointPlotOpt, Current_F) {

  Res = GetPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge, lenwt_a, ln_lenwt_a,
                                lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern, InitRatioFem,
                                FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam, mat_A50, mat_A95,
                                EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge, Land_sel_A50, Land_sel_A95,
                                EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort,
                                Steepness, SRrel_Type, NatMort, Current_F)

  # F vs SPR and Brel
  xmax = max(Res$FishMort)
  if (NatMort < 0.15) xmax = 1.0

  if (PlotOpt == 1) { # females
    plot(Res$FishMort, Res$Eq_FemRelSpBiomResults, "l", frame.plot=F, ylim=c(0, 1), xlim=c(0, xmax),
         col="black", yaxt="n", xaxt="n", ylab="", xlab="")
    points(Current_F, Res$Eq_FemRelSpBiom, cex=1.2, col="black", pch=16)
    legend("topleft", col="black", pch = 16, lty=0, legend="Estimate - females",
           bty="n", cex=0.8, inset = 0.05)
  }
  if (PlotOpt == 2) { # males
    plot(Res$FishMort, Res$Eq_MalRelSpBiomResults, "l", frame.plot=F, ylim=c(0, 1), xlim=c(0, xmax),
         col="black", yaxt="n", xaxt="n", ylab="", xlab="")
    points(Current_F, Res$Eq_MalRelSpBiom, cex=1.2, col="black", pch=16)
    legend("topleft", col="black", pch = 16, lty=0, legend="Estimate - males",
           bty="n", cex=0.8, inset = 0.05)
  }
  if (PlotOpt == 3) { # combined sex
    plot(Res$FishMort, Res$Eq_CombSexRelSpBiomResults, "l", frame.plot=F, ylim=c(0, 1), xlim=c(0, xmax),
         col="black", yaxt="n", xaxt="n", ylab="", xlab="")
    points(Current_F, Res$Eq_CombSexRelSpBiom, cex=1.2, col="black", pch=16)
    legend("topleft", col="black", pch = 16, lty=0, legend="Estimate - comb. sex",
           bty="n", cex=0.8, inset = 0.05)
  }


  axis(1, at=seq(0, xmax, 0.5), cex.axis=1, lwd=1.75, lab=F)
  axis(2, at=seq(0, 1, 0.2), cex.axis=1, lwd=1.75, lab=F)
  axis(1, at=seq(0, xmax, 0.5), labels = seq(0, xmax, 0.5),
       cex.axis=1, line=0.5, las=1, lwd=1.5, tick=F)
  axis(2, at=seq(0, 1, 0.2), cex.axis=1, line=0.5, las=1, lwd=1.5, tick=F)
  mtext(expression(paste(plain("Relative spawning biomass"))), las=3, side=2, line=3, cex=1, lwd=1.75)
  mtext(expression(paste(italic("F") ~ (year^{-1}))), las=1, side=1, line=3, cex=1, lwd=1.75)

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
#' @param FinalSex_L50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param EstSexRatioAtLen NA  # sex ratio at length, inputted as vector
#' @param EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' @param mat_L50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at length)
#' @param mat_L95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at length)
#' @param EstMatAtLen vector of proportion mature at length (set to NA if using length at maturity parameters)
#' @param sel_L50 logistic parameter for gear selectivity curve
#' @param sel_L95 logistic parameter for gear selectivity curve
#' @param EstGearSelAtLen gear selectivity curve inputted as vector
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_L50 logistic parameter for fish retention curve
#' @param ret_L95 logistic parameter for fish retention curve
#' @param EstRetenAtLen retention curve inputted as vector
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param PlotOpt # 1=females, 2=males, 3=combined sex
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
#'                          EstMalWtAtLen=NA) # weight at length, inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at recuitment age/length
#' FinalSex_L50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_L95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtLen <- NA  # sex ratio at length inputted as vector
#' EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, ~0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_L50 <- c(250, 250) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_L95 <- c(300, 300) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
#'                           EstMalMatAtLen=NA) # maturity at length, inputted as values in data frame
#' sel_L50 <- c(250, 250) # females, males - Logistic length selectivity relationship parameters
#' sel_L95 <- c(300, 300) # females, males - Logistic length selectivity relationship parameters
#' EstGearSelAtLen <- data.frame(EstFemGearSelAtLen=NA,
#'                               EstMalGearSelAtLen=NA)
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_L50 <- c(250, 250) # females, males - Logistic fish retention at length parameters
#' ret_L95 <- c(300, 300) # females, males - Logistic fish retention at length parameters
#' EstRetenAtLen <- data.frame(EsFemtRetenAtLen=NA,
#'                             EstMalRetenAtLen=NA)
#' DiscMort <- 0.25 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort = 4.22 / MaxModelAge # natural mortality  (year-1)
#' PlotOpt <- 1 # 1=females, 2=males, 3=combined sex
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' Current_F = 0.2
#' PlotPerRecruit_Biom_no_err_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams, RefnceAges,
#'                               CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,  ReprodScale, ReprodPattern,
#'                               InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,  EggFertParam, mat_L50, mat_L95,
#'                               EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax, ret_L50, ret_L95, EstRetenAtLen,
#'                               DiscMort, Steepness, SRrel_Type, NatMort, PlotOpt, RefPointPlotOpt, Current_F)
#' @export
PlotPerRecruit_Biom_no_err_LB <- function(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams, RefnceAges,
                                          CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,  ReprodScale, ReprodPattern,
                                          InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,  EggFertParam, mat_L50, mat_L95,
                                          EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax, ret_L50, ret_L95, EstRetenAtLen,
                                          DiscMort, Steepness, SRrel_Type, NatMort, PlotOpt, RefPointPlotOpt, Current_F) {

  Output_Opt = 1 # 1=standard output, 2=with added length and weight outputs (slower)
  Res = GetPerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                                ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                                EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                                ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, Current_F, Output_Opt)

  # F vs SPR and Brel
  xmax = max(Res$FishMort)
  if (NatMort < 0.15) xmax = 1.0

  if (PlotOpt==1) { # plot females
    plot(Res$FishMort, Res$Eq_FemRelSpBiomResults, "l", frame.plot=F, ylim=c(0, 1), xlim=c(0, xmax),
         col="black", yaxt="n", xaxt="n", ylab="", xlab="")
    points(Current_F, Res$Eq_FemRelSpBiom, cex=1.2, col="black", pch=16)
    legend("topleft", col="black", pch = 16, lty=0, legend="Estimate - females",
           bty="n", cex=0.8, inset = 0.05)
  }
  if (PlotOpt==2) { # plot males
    plot(Res$FishMort, Res$Eq_MalRelSpBiomResults, "l", frame.plot=F, ylim=c(0, 1), xlim=c(0, xmax),
         col="black", yaxt="n", xaxt="n", ylab="", xlab="")
    points(Current_F, Res$Eq_MalRelSpBiom, cex=1.2, col="black", pch=16)
    legend("topleft", col="black", pch = 16, lty=0, legend="Estimate - males",
           bty="n", cex=0.8, inset = 0.05)
  }
  if (PlotOpt==3) { # plot males
    plot(Res$FishMort, Res$Eq_CombSexRelSpBiomResults, "l", frame.plot=F, ylim=c(0, 1), xlim=c(0, xmax),
         col="black", yaxt="n", xaxt="n", ylab="", xlab="")
    points(Current_F, Res$Eq_CombSexRelSpBiom, cex=1.2, col="black", pch=16)
    legend("topleft", col="black", pch = 16, lty=0, legend="Estimate - comb. sex",
           bty="n", cex=0.8, inset = 0.05)
  }

  axis(1, at=seq(0, xmax, 0.5), cex.axis=1, lwd=1.75, lab=F)
  axis(2, at=seq(0, 1, 0.2), cex.axis=1, lwd=1.75, lab=F)
  axis(1, at=seq(0, xmax, 0.5), labels = seq(0, xmax, 0.5),
       cex.axis=1, line=0.5, las=1, lwd=1.5, tick=F)
  axis(2, at=seq(0, 1, 0.2), cex.axis=1, line=0.5, las=1, lwd=1.5, tick=F)
  mtext(expression(paste(plain("Relative spawning biomass"))), las=3, side=2, line=3, cex=1, lwd=1.75)
  mtext(expression(paste(italic("F") ~ (year^{-1}))), las=1, side=1, line=3, cex=1, lwd=1.75)

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
#' @param FinalSex_A50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_A95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param EstSexRatioAtAge NA  # sex ratio at age (from age 0) inputted as vector
#' @param EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' @param mat_A50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_A95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtAge vector of proportion mature at age (set to NA if using age at maturity parameters)
#' @param Gear_sel_A50 logistic parameter for gear selectivity curve
#' @param Gear_sel_A95 logistic parameter for gear selectivity curve
#' @param EstGearSelAtAge vector of gear selectivity at age (set to NA if using age at gear selectivity parameters)
#' @param Land_sel_A50 logistic parameter for selectivity of landings curve
#' @param Land_sel_A95 logistic parameter for selectivity of landings curve
#' @param EstLandSelAtAge vector for for selectivity of landings at age (set to NA if using age at selectivity of landings parameters)
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
#' @return fishing mortality values for analysis (PerRec_FValues), female SPR estimates vs F (Fem_SPR_Vals)
#' estimates of equilbrium female spawning biomass vs F (Eq_RelFemSpBiom_Vals), estimated BMSY values (BMSY_Vals)
#' resampled female SPR values vs F (Sim_FemSPR=Sim_FemSPR), resampled equilbrium female, male and combined sex
#' spawning biomass values vs F (Sim_Eq_RelFemSpBiom, Sim_Eq_RelMalSpBiom, Sim_Eq_RelCombSexSpBiom),
#' median female, male and combined sex SPR values vs F (EstFemSPR, EstMalSPR, EstCombSexSPR), median BMSY ratio estimate
#' (EstMedBMSYratio), summary of estimates with associated 95 percent confidence limits (ResSummary_with_err).
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
#' FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtAge <- NA  # sex ratio at age (from age 0) inputted as vector
#' EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' Gear_sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' Gear_sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' EstGearSelAtAge <- data.frame(FemGearSelAtAge=NA, MalGearSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' Land_sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' Land_sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' EstLandSelAtAge <- data.frame(FemSelLandAtAge=NA, MalSelLandAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
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
#' GetPerRecruitResults_AB_with_err(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge, lenwt_a, ln_lenwt_a,
#'                                  lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern, InitRatioFem,
#'                                  FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam, mat_A50, mat_A95,
#'                                  EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge, Land_sel_A50, Land_sel_A95,
#'                                  EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness,
#'                                  Steepness_sd, SRrel_Type, NatMort, NatMort_sd, Current_F, Current_F_sd, nReps)
#' @export
GetPerRecruitResults_AB_with_err <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge, lenwt_a, ln_lenwt_a,
                                             lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern, InitRatioFem,
                                             FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam, mat_A50, mat_A95,
                                             EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge, Land_sel_A50, Land_sel_A95,
                                             EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness,
                                             Steepness_sd, SRrel_Type, NatMort, NatMort_sd, Current_F, Current_F_sd, nReps) {



  FValues = rnorm(nReps, Current_F, Current_F_sd)
  hValues = rnorm(nReps, Steepness, Steepness_sd)
  MValues = rnorm(nReps, NatMort, NatMort_sd)


  for (i in 1:nReps) {
    FMort = FValues[i]
    Steepness = hValues[i]
    NatMort = MValues[i]
    PREst = GetPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge, lenwt_a, ln_lenwt_a,
                                    lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern, InitRatioFem,
                                    FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam, mat_A50, mat_A95,
                                    EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge, Land_sel_A50, Land_sel_A95,
                                    EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort,
                                    Steepness, SRrel_Type, NatMort, Current_F)

    if (i == 1) {
      Fem_SPR_Vals = rep(0, nReps)
      Mal_SPR_Vals = rep(0, nReps)
      CombSex_SPR_Vals = rep(0, nReps)
      Eq_RelFemSpBiom_Vals = rep(0, nReps)
      Eq_RelMalSpBiom_Vals = rep(0, nReps)
      Eq_RelCombSexSpBiom_Vals = rep(0, nReps)
      BMSY_Vals = rep(0,nReps)
      Sim_FemSPR = data.frame(matrix(nrow=nReps, ncol=length(PREst$FishMort)))
      colnames(Sim_FemSPR) = PREst$FishMort
      Sim_FemSPR = as.matrix(Sim_FemSPR)
      Sim_Eq_RelFemSpBiom = as.matrix(Sim_FemSPR)
      Sim_Eq_RelMalSpBiom = Sim_Eq_RelFemSpBiom
      Sim_Eq_RelCombSexSpBiom = Sim_Eq_RelFemSpBiom
    }

    Fem_SPR_Vals[i] = PREst$Fem_SPR
    Mal_SPR_Vals[i] = PREst$Mal_SPR
    CombSex_SPR_Vals[i] = PREst$CombSex_SPR

    Eq_RelFemSpBiom_Vals[i] = PREst$Eq_FemRelSpBiom
    Eq_RelMalSpBiom_Vals[i] = PREst$Eq_MalRelSpBiom
    Eq_RelCombSexSpBiom_Vals[i] = PREst$Eq_CombSexRelSpBiom

    Sim_FemSPR[i,] = PREst$Fem_SPRResults
    Sim_Eq_RelFemSpBiom[i,] = PREst$Eq_FemRelSpBiomResults
    Sim_Eq_RelMalSpBiom[i,] = PREst$Eq_MalRelSpBiomResults
    Sim_Eq_RelCombSexSpBiom[i,] = PREst$Eq_CombSexRelSpBiomResults

    BMSY_Vals[i] = PREst$BMSY_Thresh

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

  EstMedEquilRelFemSpBiom = round(median(Eq_RelFemSpBiom_Vals),3)
  Low95EquilRelFemSpBiom = as.numeric(round(quantile(Eq_RelFemSpBiom_Vals, 0.025),3))
  Upp95EquilRelFemSpBiom = as.numeric(round(quantile(Eq_RelFemSpBiom_Vals, 0.975),3))
  EstEquilRelFemSpBiom = data.frame(EstMedEquilRelFemSpBiom=EstMedEquilRelFemSpBiom,
                                    Low95EquilRelFemSpBiom=Low95EquilRelFemSpBiom,
                                    Upp95EquilRelFemSpBiom=Upp95EquilRelFemSpBiom)

  EstMedEquilRelMalSpBiom = round(median(Eq_RelMalSpBiom_Vals),3)
  Low95EquilRelMalSpBiom = as.numeric(round(quantile(Eq_RelMalSpBiom_Vals, 0.025),3))
  Upp95EquilRelMalSpBiom = as.numeric(round(quantile(Eq_RelMalSpBiom_Vals, 0.975),3))
  EstEquilRelMalSpBiom = data.frame(EstMedEquilRelMalSpBiom=EstMedEquilRelMalSpBiom,
                                    Low95EquilRelMalSpBiom=Low95EquilRelMalSpBiom,
                                    Upp95EquilRelMalSpBiom=Upp95EquilRelMalSpBiom)

  EstMedEquilRelCombSexSpBiom = round(median(Eq_RelCombSexSpBiom_Vals),3)
  Low95EquilRelCombSexSpBiom = as.numeric(round(quantile(Eq_RelCombSexSpBiom_Vals, 0.025),3))
  Upp95EquilRelCombSexSpBiom = as.numeric(round(quantile(Eq_RelCombSexSpBiom_Vals, 0.975),3))
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
                 Eq_RelFemSpBiom_Vals=Eq_RelFemSpBiom_Vals,
                 BMSY_Vals=BMSY_Vals,
                 Sim_FemSPR=Sim_FemSPR,
                 Sim_Eq_RelFemSpBiom=Sim_Eq_RelFemSpBiom,
                 Sim_Eq_RelMalSpBiom=Sim_Eq_RelMalSpBiom,
                 Sim_Eq_RelCombSexSpBiom=Sim_Eq_RelCombSexSpBiom,
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
#' @param FinalSex_L50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param EstSexRatioAtLen NA  # sex ratio at length, inputted as vector
#' @param EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' @param mat_L50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_L95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtLen vector of proportion mature at length (set to NA if using age at maturity parameters)
#' @param sel_L50 logistic parameter for gear selectivity curve
#' @param sel_L95 logistic parameter for gear selectivity curve
#' @param EstGearSelAtLen gear selectivity curve inputted as vector
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_L50 logistic parameter for fish retention curve
#' @param ret_L95 logistic parameter for fish retention curve
#' @param EstRetenAtLen retention curve inputted as vector
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
#' estimates of equilbrium female spawning biomass vs F (Eq_RelFemSpBiom_Vals), estimated BMSY values (BMSY_Vals)
#' resampled female SPR values vs F (Sim_FemSPR=Sim_FemSPR), resampled equilbrium female, male and combined sex
#' spawning biomass values vs F (Sim_Eq_RelFemSpBiom, Sim_Eq_RelMalSpBiom, Sim_Eq_RelCombSexSpBiom),
#' median female, male and combined sex SPR values vs F (EstFemSPR, EstMalSPR, EstCombSexSPR), median BMSY ratio estimate
#' (EstMedBMSYratio), summary of estimates with associated 95 percent confidence limits (ResSummary_with_err).
#'
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
#'                          EstMalWtAtLen=NA) # weight at length, inputted as values in data frame
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' InitRatioFem <- 0.5 # Ratio of females to males at recruitment age/size
#' FinalSex_L50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_L95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtLen <- NA  # sex ratio at length inputted as vector
#' EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, ~0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_L50 <- c(250, 250) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_L95 <- c(300, 300) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
#'                           EstMalMatAtLen=NA) # maturity at length, inputted as values in data frame
#' sel_L50 <- c(250, 250) # females, males - Logistic length selectivity relationship parameters
#' sel_L95 <- c(300, 300) # females, males - Logistic length selectivity relationship parameters
#' EstGearSelAtLen <- data.frame(EstFemGearSelAtLen=NA,
#'                               EstMalGearSelAtLen=NA)
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_L50 <- c(250, 250) # females, males - Logistic fish retention at length parameters
#' ret_L95 <- c(300, 300) # females, males - Logistic fish retention at length parameters
#' EstRetenAtLen <- data.frame(EsFemtRetenAtLen=NA,
#'                             EstMalRetenAtLen=NA)
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
#'                                  RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
#'                                  ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
#'                                  EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
#'                                  ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd,
#'                                  Current_F, Current_F_sd, nReps)
#' @export
GetPerRecruitResults_LB_with_err <- function(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                             RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                                             ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                                             EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                                             ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd,
                                             Current_F, Current_F_sd, nReps) {



  FValues = rnorm(nReps, Current_F, Current_F_sd)
  hValues = rnorm(nReps, Steepness, Steepness_sd)
  MValues = rnorm(nReps, NatMort, NatMort_sd)
  Output_Opt = 1 # 1=standard output, 2=with added length and weight outputs (slower)

  for (i in 1:nReps) {
    FMort = FValues[i]
    Steepness = hValues[i]
    NatMort = MValues[i]
    PREst = GetPerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                    RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                                    ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                                    EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                                    ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, Current_F, Output_Opt)

    if (i == 1) {
      Fem_SPR_Vals = rep(0, nReps)
      Mal_SPR_Vals = rep(0, nReps)
      CombSex_SPR_Vals = rep(0, nReps)

      Eq_RelFemSpBiom_Vals = rep(0, nReps)
      Eq_RelMalSpBiom_Vals = rep(0, nReps)
      Eq_RelCombSexSpBiom_Vals = rep(0, nReps)

      BMSY_Vals = rep(0,nReps)
      Sim_FemSPR = data.frame(matrix(nrow=nReps, ncol=length(PREst$FishMort)))
      colnames(Sim_FemSPR) = PREst$FishMort
      Sim_Eq_RelFemSpBiom = as.matrix(Sim_FemSPR)
      Sim_Eq_RelMalSpBiom = Sim_Eq_RelFemSpBiom
      Sim_Eq_RelCombSexSpBiom = Sim_Eq_RelFemSpBiom
    }

    Fem_SPR_Vals[i] = PREst$Fem_SPR
    Mal_SPR_Vals[i] = PREst$Mal_SPR
    CombSex_SPR_Vals[i] = PREst$CombSex_SPR

    Eq_RelFemSpBiom_Vals[i] = PREst$Eq_FemRelSpBiom
    Eq_RelMalSpBiom_Vals[i] = PREst$Eq_MalRelSpBiom
    Eq_RelCombSexSpBiom_Vals[i] = PREst$Eq_CombSexRelSpBiom

    BMSY_Vals[i] = PREst$BMSY_Thresh
    Sim_FemSPR[i,] = PREst$Fem_SPRResults
    Sim_Eq_RelFemSpBiom[i,] = PREst$Eq_FemRelSpBiomResults
    Sim_Eq_RelMalSpBiom[i,] = PREst$Eq_MalRelSpBiomResults
    Sim_Eq_RelCombSexSpBiom[i,] = PREst$Eq_CombSexRelSpBiomResults

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

  EstMedEquilRelFemSpBiom = round(median(Eq_RelFemSpBiom_Vals),3)
  Low95EquilRelFemSpBiom = as.numeric(round(quantile(Eq_RelFemSpBiom_Vals, 0.025),3))
  Upp95EquilRelFemSpBiom = as.numeric(round(quantile(Eq_RelFemSpBiom_Vals, 0.975),3))
  EstEquilRelFemSpBiom = data.frame(EstMedEquilRelFemSpBiom=EstMedEquilRelFemSpBiom,
                                    Low95EquilRelFemSpBiom=Low95EquilRelFemSpBiom,
                                    Upp95EquilRelFemSpBiom=Upp95EquilRelFemSpBiom)

  EstMedEquilRelMalSpBiom = round(median(Eq_RelMalSpBiom_Vals),3)
  Low95EquilRelMalSpBiom = as.numeric(round(quantile(Eq_RelMalSpBiom_Vals, 0.025),3))
  Upp95EquilRelMalSpBiom = as.numeric(round(quantile(Eq_RelMalSpBiom_Vals, 0.975),3))
  EstEquilRelMalSpBiom = data.frame(EstMedEquilRelMalSpBiom=EstMedEquilRelMalSpBiom,
                                    Low95EquilRelMalSpBiom=Low95EquilRelMalSpBiom,
                                    Upp95EquilRelMalSpBiom=Upp95EquilRelMalSpBiom)

  EstMedEquilRelCombSexSpBiom = round(median(Eq_RelCombSexSpBiom_Vals),3)
  Low95EquilRelCombSexSpBiom = as.numeric(round(quantile(Eq_RelCombSexSpBiom_Vals, 0.025),3))
  Upp95EquilRelCombSexSpBiom = as.numeric(round(quantile(Eq_RelCombSexSpBiom_Vals, 0.975),3))
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
                 Eq_RelFemSpBiom_Vals=Eq_RelFemSpBiom_Vals,
                 BMSY_Vals=BMSY_Vals,
                 Sim_FemSPR=Sim_FemSPR,
                 Sim_Eq_RelFemSpBiom=Sim_Eq_RelFemSpBiom,
                 Sim_Eq_RelMalSpBiom=Sim_Eq_RelMalSpBiom,
                 Sim_Eq_RelCombSexSpBiom=Sim_Eq_RelCombSexSpBiom,
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
#' @param FinalSex_A50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_A95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param EstSexRatioAtAge NA  # sex ratio at age (from age 0) inputted as vector
#' @param EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' @param mat_A50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param mat_A95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at age)
#' @param EstMatAtAge vector of proportion mature at age (set to NA if using age at maturity parameters)
#' @param Gear_sel_A50 logistic parameter for gear selectivity curve
#' @param Gear_sel_A95 logistic parameter for gear selectivity curve
#' @param EstGearSelAtAge vector of gear selectivity at age (set to NA if using age at gear selectivity parameters
#' @param Land_sel_A50 logistic parameter for selectivity of landings curve
#' @param Land_sel_A95 logistic parameter for selectivity of landings curve
#' @param EstLandSelAtAge vector for for selectivity of landings at age (set to NA if using age at selectivity of landings parameters)
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
#' @param PlotOpt # 1=females, 2=males, 3=combined sex
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
#' FinalSex_A50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtAge <- NA  # sex ratio at age (from age 0) inputted as vector
#' EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_A50 <- c(2.5, 2.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(3.5, 3.5) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' Gear_sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' Gear_sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' EstGearSelAtAge <- data.frame(FemGearSelAtAge=NA, MalGearSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' Land_sel_A50 <- c(2.5, 2.5) # females, males - Logistic age selectivity relationship parameters
#' Land_sel_A95 <- c(3.5, 3.5) # females, males - Logistic age selectivity relationship parameters
#' EstLandSelAtAge <- data.frame(FemSelLandAtAge=NA, MalSelLandAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
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
#' PlotOpt <- 1 # 1=females, 2=males, 3=combined sex
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' nReps = 50
#' FittedRes=GetPerRecruitResults_AB_with_err(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge, lenwt_a, ln_lenwt_a,
#'                                  lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern, InitRatioFem,
#'                                  FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam, mat_A50, mat_A95,
#'                                  EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge, Land_sel_A50, Land_sel_A95,
#'                                  EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness,
#'                                  Steepness_sd, SRrel_Type, NatMort, NatMort_sd, Current_F, Current_F_sd, nReps)
#' # Plot. Note, can skip above step and set FittedRes=NA (plot function will be slower
#' PlotPerRecruit_Biom_with_err_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge, lenwt_a, ln_lenwt_a, lenwt_b,
#'                                 WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_A50,
#'                                 FinalSex_A95, EstSexRatioAtAge, EggFertParam, mat_A50, mat_A95, EstMatAtAge,
#'                                 Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge, Land_sel_A50, Land_sel_A95, EstLandSelAtAge,
#'                                 ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness, Steepness_sd, SRrel_Type,
#'                                 NatMort, NatMort_sd, Current_F, Current_F_sd, PlotOpt, RefPointPlotOpt, FittedRes, nReps,
#'                                 MainLabel=NA, xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA, ymax=NA, yint=NA)
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
#' FinalSex_A50 <- 35 # Logistic sex change relationship parameters (inflection point)
#' FinalSex_A95 <- 60 # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtAge <- NA  # sex ratio at age (from age 0) inputted as vector
#' EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_A50 <- c(20, 20) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_A95 <- c(30, 30) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtAge <- data.frame(EstFemMatAtAge=NA, EstMalMatAtAge=NA) # maturity at age (from age 0), inputted as values in data frame
#' Gear_sel_A50 <- c(20, 20) # females, males - Logistic age selectivity relationship parameters
#' Gear_sel_A95 <- c(30, 30) # females, males - Logistic age selectivity relationship parameters
#' EstGearSelAtAge <- data.frame(FemGearSelAtAge=NA, MalGearSelAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' Land_sel_A50 <- c(25, 25) # females, males - Logistic age selectivity relationship parameters
#' Land_sel_A95 <- c(35, 35) # females, males - Logistic age selectivity relationship parameters
#' EstLandSelAtAge <- data.frame(FemSelLandAtAge=NA, MalSelLandAtAge=NA) # gear selectivity at age (from age 0), inputted as values in data frame
#' ret_Pmax <- NA  # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_A50 <- NA  # females, males - Logistic age fish retention at age parameters
#' ret_A95 <- NA  # females, males - Logistic age fish retention at age parameters
#' EstRetenAtAge <- data.frame(EstFemRetenAtAge=rep(1,MaxModelAge+1), EstMalRetenAtAge=rep(1,MaxModelAge+1)) # fish retention at age (from age 0), inputted as values in data frame
#' DiscMort <- 0.0 # discard mortality (e.g. 50% released fish die = 0.5)
#' Steepness <- 0.75 # steepness parameter of the Beverton and Holt stock-recruitment relationship
#' Steepness_sd <- 0.025
#' SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
#' NatMort <- 0.07 # natural mortality  (year-1)
#' NatMort_sd <- 0.005
#' Current_F <- 0.1 # estimate of fishing mortality, e.g. from catch curve analysis
#' Current_F_sd <- 0.005
#' PlotOpt <- 1 # 1=females, 2=males, 3=combined sex
#' RefPointPlotOpt <- 1 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
#' nReps = 10
#' FittedRes=GetPerRecruitResults_AB_with_err(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge, lenwt_a, ln_lenwt_a,
#'                                  lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern, InitRatioFem,
#'                                  FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam, mat_A50, mat_A95,
#'                                  EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge, Land_sel_A50, Land_sel_A95,
#'                                  EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness,
#'                                  Steepness_sd, SRrel_Type, NatMort, NatMort_sd, Current_F, Current_F_sd, nReps)
#' # Plot. Note, can skip above step and set FittedRes=NA (plot function will be slower
#' PlotPerRecruit_Biom_with_err_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge, lenwt_a, ln_lenwt_a, lenwt_b,
#'                                 WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_A50,
#'                                 FinalSex_A95, EstSexRatioAtAge, EggFertParam, mat_A50, mat_A95, EstMatAtAge,
#'                                 Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge, Land_sel_A50, Land_sel_A95, EstLandSelAtAge,
#'                                 ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness, Steepness_sd, SRrel_Type,
#'                                 NatMort, NatMort_sd, Current_F, Current_F_sd, PlotOpt, RefPointPlotOpt, FittedRes, nReps,
#'                                 MainLabel=NA, xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA, ymax=NA, yint=NA)
#' @export
PlotPerRecruit_Biom_with_err_AB <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge, lenwt_a, ln_lenwt_a, lenwt_b,
                                            WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_A50,
                                            FinalSex_A95, EstSexRatioAtAge, EggFertParam, mat_A50, mat_A95, EstMatAtAge,
                                            Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge, Land_sel_A50, Land_sel_A95, EstLandSelAtAge,
                                            ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness, Steepness_sd, SRrel_Type,
                                            NatMort, NatMort_sd, Current_F, Current_F_sd, PlotOpt, RefPointPlotOpt, FittedRes, nReps,
                                            MainLabel, xaxis_lab, yaxis_lab, xmax, xint, ymax, yint) {

  # get BMSY reference points
  res = GetPerRecruitResults_AB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge, lenwt_a, ln_lenwt_a,
                                lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern, InitRatioFem,
                                FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam, mat_A50, mat_A95,
                                EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge, Land_sel_A50, Land_sel_A95,
                                EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort,
                                Steepness, SRrel_Type, NatMort, Current_F)

  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    Res =  FittedRes
  } else {
    Res=GetPerRecruitResults_AB_with_err(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge, lenwt_a, ln_lenwt_a,
                                         lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern, InitRatioFem,
                                         FinalSex_A50, FinalSex_A95, EstSexRatioAtAge, EggFertParam, mat_A50, mat_A95,
                                         EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge, Land_sel_A50, Land_sel_A95,
                                         EstLandSelAtAge, ret_Pmax, ret_A50, ret_A95, EstRetenAtAge, DiscMort, Steepness,
                                         Steepness_sd, SRrel_Type, NatMort, NatMort_sd, Current_F, Current_F_sd, nReps)
  }

  if (is.na(yaxis_lab)) yaxis_lab = "Relative spawning biomass"
  if (is.na(xaxis_lab)) xaxis_lab = expression(paste(italic("F") ~ (year^{-1})))
  if (is.na(xmax)) xmax = 2
  if (is.na(xint)) xint = 0.5
  if (is.na(ymax)) ymax = 1
  if (is.na(yint)) yint = 0.2


  # Plot per recruit outputs with uncertainty
  if (PlotOpt==1) { # plot females
    EqB_med = apply(Res$Sim_Eq_RelFemSpBiom,2,quantile, probs=c(0.5))
    EqB_lw = apply(Res$Sim_Eq_RelFemSpBiom,2,quantile, probs=c(0.025))
    EqB_hi = apply(Res$Sim_Eq_RelFemSpBiom,2,quantile, probs=c(0.975))
    x=which(Res$PerRec_FValues==xmax)
    plot(Res$PerRec_FValues[1:x], EqB_med[1:x], "l", frame.plot=F, ylim=c(0,ymax), xlim=c(0,xmax),
         col="red", yaxt="n", xaxt="n", ylab="", xlab="", main=MainLabel)
    polygon(c(Res$PerRec_FValues[1:x],rev(Res$PerRec_FValues[1:x])),c(EqB_lw[1:x],rev(EqB_hi[1:x])),
            col="lightpink", border="lightpink")
    lines(Res$PerRec_FValues[1:x], EqB_med[1:x],col="red")
    points(Current_F, Res$EstEquilRelFemSpBiom[1], cex=1.2, col="red", pch=16)
    lw=as.numeric(Res$EstEquilRelFemSpBiom[2]); up=as.numeric(Res$EstEquilRelFemSpBiom[3])
    arrows(Current_F, lw, Current_F, up,length=0.05, angle=90, code=3,col="red")
    legend("topleft", col="red", pch = 16, legend="Estimate - females",
           bty="n", cex=1,0, lty=0, inset = 0.05)
  }

  if (PlotOpt==2) { # plot males
    EqB_med = apply(Res$Sim_Eq_RelMalSpBiom,2,quantile, probs=c(0.5))
    EqB_lw = apply(Res$Sim_Eq_RelMalSpBiom,2,quantile, probs=c(0.025))
    EqB_hi = apply(Res$Sim_Eq_RelMalSpBiom,2,quantile, probs=c(0.975))
    x=which(Res$PerRec_FValues==xmax)
    plot(Res$PerRec_FValues[1:x], EqB_med[1:x], "l", frame.plot=F, ylim=c(0,ymax), xlim=c(0,xmax),
         col="blue", yaxt="n", xaxt="n", ylab="", xlab="", main=MainLabel)

    polygon(c(Res$PerRec_FValues[1:x],rev(Res$PerRec_FValues[1:x])),c(EqB_lw[1:x],rev(EqB_hi[1:x])),
            col="lightblue", border="lightblue")
    lines(Res$PerRec_FValues[1:x], EqB_med[1:x],col="blue")
    lw=as.numeric(Res$EstEquilRelMalSpBiom[2]); up=as.numeric(Res$EstEquilRelMalSpBiom[3])
    arrows(Current_F, lw, Current_F, up,length=0.05, angle=90, code=3,col="blue")
    points(Current_F, Res$EstEquilRelMalSpBiom[1], cex=1.2, col="blue", pch=16)
    legend("topleft", col="blue", pch = 16, legend="Estimate - males",
           bty="n", cex=1,0, lty=0, inset = 0.05)
  }

  if (PlotOpt==3) { # plot combined sex
    EqB_med = apply(Res$Sim_Eq_RelCombSexSpBiom,2,quantile, probs=c(0.5))
    EqB_lw = apply(Res$Sim_Eq_RelCombSexSpBiom,2,quantile, probs=c(0.025))
    EqB_hi = apply(Res$Sim_Eq_RelCombSexSpBiom,2,quantile, probs=c(0.975))
    x=which(Res$PerRec_FValues==xmax)
    plot(Res$PerRec_FValues[1:x], EqB_med[1:x], "l", frame.plot=F, ylim=c(0,ymax), xlim=c(0,xmax),
         col="black", yaxt="n", xaxt="n", ylab="", xlab="", main=MainLabel)

    polygon(c(Res$PerRec_FValues[1:x],rev(Res$PerRec_FValues[1:x])),c(EqB_lw[1:x],rev(EqB_hi[1:x])),
            col="lightgrey", border="lightgrey")
    lines(Res$PerRec_FValues[1:x], EqB_med[1:x],col="black")
    lw=as.numeric(Res$EstEquilRelCombSexSpBiom[2]); up=as.numeric(Res$EstEquilRelCombSexSpBiom[3])
    arrows(Current_F, lw, Current_F, up,length=0.05, angle=90, code=3,col="black")
    points(Current_F, Res$EstEquilRelCombSexSpBiom[1], cex=1.2, col="black", pch=16)
    legend("topleft", col="black", pch = 16, legend="Estimate - comb. sex",
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
#' @param FinalSex_L50 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param FinalSex_L95 logistic sex change parameter for hermaphroditic species (set to NA for gonochoristic species)
#' @param EstSexRatioAtLen NA  # sex ratio at length, inputted as vector
#' @param EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, 0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' @param mat_L50 logistic length at maturity parameter (set to NA if directly inputting proportion mature at length)
#' @param mat_L95 logistic length at maturity parameter (set to NA if directly inputting proportion mature at length)
#' @param EstMatAtLen vector of proportion mature at length (set to NA if using length at maturity parameters)
#' @param sel_L50 logistic parameter for gear selectivity curve
#' @param sel_L95 logistic parameter for gear selectivity curve
#' @param EstGearSelAtLen gear selectivity curve inputted as vector
#' @param ret_Pmax logistic parameter for fish retention curve
#' @param ret_L50 logistic parameter for fish retention curve
#' @param ret_L95 logistic parameter for fish retention curve
#' @param EstRetenAtLen retention curve inputted as vector
#' @param DiscMort proportion of fish that die following capture and release
#' @param Steepness steepness parameter of Beverton-Holt or Ricker stock-recruitment relationship
#' @param Steepness_sd standard deviation for steepness parameter
#' @param SRrel_Type 1 = Beverton-Holt, 2=Ricker stock-recruitment relationship
#' @param NatMort natural mortality
#' @param NatMort_sd standard deviation for natural mortality
#' @param Current_F estimated current fishing mortality
#' @param Current_F_sd standard deviation for estimated current fishing mortality
#' @param PlotOpt # 1=females, 2=males, 3=combined sex
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
#'                          EstMalWtAtLen=NA) # weight at length, inputted as values in data frame
#' ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
#' ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
#' InitRatioFem <- 0.5 # Ratio of females to males at recruitment age/length
#' FinalSex_L50 <- NA # Logistic sex change relationship parameters (inflection point)
#' FinalSex_L95 <- NA # Logistic sex change relationship parameters (95% of max probability)
#' EstSexRatioAtLen <- NA  # sex ratio at length inputted as vector
#' EggFertParam <- NA # (NA or from ~0.2-1) NA = no effect, ~0.2 = direct effect of popn. sex ratio changes on egg fertilisation rates, 1 = no effects
#' mat_L50 <- c(250, 250) # females, males - Logistic length (mm) at maturity relationship parameters
#' mat_L95 <- c(300, 300) # females, males - Logistic length (mm) at maturity relationship parameters
#' EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
#'                           EstMalMatAtLen=NA) # maturity at length, inputted as values in data frame
#' sel_L50 <- c(250, 250) # females, males - Logistic length selectivity relationship parameters
#' sel_L95 <- c(300, 300) # females, males - Logistic length selectivity relationship parameters
#' EstGearSelAtLen <- data.frame(EstFemGearSelAtLen=NA,
#'                               EstMalGearSelAtLen=NA)
#' ret_Pmax <- c(1.0, 1.0) # maximum retention, values lower than 1 imply discarding of fish above MLL
#' ret_L50 <- c(250, 250) # females, males - Logistic fish retention at length parameters
#' ret_L95 <- c(300, 300) # females, males - Logistic fish retention at length parameters
#' EstRetenAtLen <- data.frame(EsFemtRetenAtLen=NA,
#'                             EstMalRetenAtLen=NA)
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
#'                                            RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
#'                                            ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
#'                                            EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
#'                                            ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd,
#'                                            Current_F, Current_F_sd, nReps)
#' # Plot. Note, can skip above step and set FittedRes=NA (plot function will be slower
#' PlotPerRecruit_Biom_with_err_LB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
#'                                             lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern,
#'                                             InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen, EggFertParam, mat_A50, mat_A95,
#'                                             EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge, ret_Pmax, ret_A50, ret_A95,
#'                                             EstRetenAtAge, DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd,
#'                                             Current_F, Current_F_sd, PlotOpt, RefPointPlotOpt, FittedRes, nReps, MainLabel=NA,
#'                                             xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA, ymax=NA, yint=NA)
#' @export
PlotPerRecruit_Biom_with_err_LB <- function(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                            lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern,
                                            InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen, EggFertParam, mat_A50, mat_A95,
                                            EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge, ret_Pmax, ret_A50, ret_A95,
                                            EstRetenAtAge, DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd,
                                            Current_F, Current_F_sd, PlotOpt, RefPointPlotOpt, FittedRes, nReps, MainLabel,
                                            xaxis_lab, yaxis_lab, xmax, xint, ymax, yint) {


  # get BMSY reference points
  Output_Opt = 1 # 1=standard output, 2=with added length and weight outputs (slower)
  res=GetPerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                              RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                              ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                              EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                              ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, Current_F, Output_Opt)


  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    Res =  FittedRes
  } else {
    Res=GetPerRecruitResults_LB_with_err(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                         RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtLen,
                                         ReprodScale, ReprodPattern, InitRatioFem, FinalSex_L50, FinalSex_L95, EstSexRatioAtLen,
                                         EggFertParam, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                                         ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd,
                                         Current_F, Current_F_sd, nReps)
  }

  if (is.na(yaxis_lab)) yaxis_lab = "Relative spawning biomass"
  if (is.na(xaxis_lab)) xaxis_lab = expression(paste(italic("F") ~ (year^{-1})))
  if (is.na(xmax)) xmax = 2
  if (is.na(xint)) xint = 0.5
  if (is.na(ymax)) ymax = 1
  if (is.na(yint)) yint = 0.2

  if (PlotOpt==1) { # plot females
    EqB_med = apply(Res$Sim_Eq_RelFemSpBiom,2,quantile, probs=c(0.5))
    EqB_lw = apply(Res$Sim_Eq_RelFemSpBiom,2,quantile, probs=c(0.025))
    EqB_hi = apply(Res$Sim_Eq_RelFemSpBiom,2,quantile, probs=c(0.975))
    x=which(Res$PerRec_FValues==xmax)
    plot(Res$PerRec_FValues[1:x], EqB_med[1:x], "l", frame.plot=F, ylim=c(0,ymax), xlim=c(0,xmax),
         col="red", yaxt="n", xaxt="n", ylab="", xlab="", main=MainLabel)

    polygon(c(Res$PerRec_FValues[1:x],rev(Res$PerRec_FValues[1:x])),c(EqB_lw[1:x],rev(EqB_hi[1:x])),
            col="lightpink", border="lightpink")
    lines(Res$PerRec_FValues[1:x], EqB_med[1:x],col="red")
    points(Current_F, Res$EstEquilRelFemSpBiom[1], cex=1.2, col="red", pch=16)
    lw=as.numeric(Res$EstEquilRelFemSpBiom[2]); up=as.numeric(Res$EstEquilRelFemSpBiom[3])
    arrows(Current_F, lw, Current_F, up,length=0.05, angle=90, code=3,col="red")
    legend("topleft", col="red", pch = 16, legend="Estimate - females",
           bty="n", cex=1,0, lty=0, inset = 0.05)

  }

  if (PlotOpt==2) { # plot males
    EqB_med = apply(Res$Sim_Eq_RelMalSpBiom,2,quantile, probs=c(0.5))
    EqB_lw = apply(Res$Sim_Eq_RelMalSpBiom,2,quantile, probs=c(0.025))
    EqB_hi = apply(Res$Sim_Eq_RelMalSpBiom,2,quantile, probs=c(0.975))
    x=which(Res$PerRec_FValues==xmax)
    plot(Res$PerRec_FValues[1:x], EqB_med[1:x], "l", frame.plot=F, ylim=c(0,ymax), xlim=c(0,xmax),
         col="blue", yaxt="n", xaxt="n", ylab="", xlab="", main=MainLabel)
    polygon(c(Res$PerRec_FValues[1:x],rev(Res$PerRec_FValues[1:x])),c(EqB_lw[1:x],rev(EqB_hi[1:x])),
            col="lightblue", border="lightblue")
    lines(Res$PerRec_FValues[1:x], EqB_med[1:x],col="blue")
    lw=as.numeric(Res$EstEquilRelMalSpBiom[2]); up=as.numeric(Res$EstEquilRelMalSpBiom[3])
    arrows(Current_F, lw, Current_F, up,length=0.05, angle=90, code=3,col="blue")
    points(Current_F, Res$EstEquilRelMalSpBiom[1], cex=1.2, col="blue", pch=16)
    legend("topleft", col="blue", pch = 16, legend="Estimate - males",
           bty="n", cex=1,0, lty=0, inset = 0.05)
  }

  if (PlotOpt==3) { # plot combined sex
    EqB_med = apply(Res$Sim_Eq_RelCombSexSpBiom,2,quantile, probs=c(0.5))
    EqB_lw = apply(Res$Sim_Eq_RelCombSexSpBiom,2,quantile, probs=c(0.025))
    EqB_hi = apply(Res$Sim_Eq_RelCombSexSpBiom,2,quantile, probs=c(0.975))
    x=which(Res$PerRec_FValues==xmax)
    plot(Res$PerRec_FValues[1:x], EqB_med[1:x], "l", frame.plot=F, ylim=c(0,ymax), xlim=c(0,xmax),
         col="black", yaxt="n", xaxt="n", ylab="", xlab="", main=MainLabel)
    polygon(c(Res$PerRec_FValues[1:x],rev(Res$PerRec_FValues[1:x])),c(EqB_lw[1:x],rev(EqB_hi[1:x])),
            col="lightgrey", border="lightgrey")
    lines(Res$PerRec_FValues[1:x], EqB_med[1:x],col="black")
    lw=as.numeric(Res$EstEquilRelCombSexSpBiom[2]); up=as.numeric(Res$EstEquilRelCombSexSpBiom[3])
    arrows(Current_F, lw, Current_F, up,length=0.05, angle=90, code=3,col="black")
    points(Current_F, Res$EstEquilRelCombSexSpBiom[1], cex=1.2, col="black", pch=16)
    legend("topleft", col="black", pch = 16, legend="Estimate - comb. sex",
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

#' Get 95 percent prediction intervals for growth curves used in length-based per recruit analysis
#'
#' This function outputs 95 percent prediction intervals for growth curves used in length-based per recruit analysis
#'
#' @keywords internal
#'
#' @param nTimeSteps number of model timesteps
#' @param nLenCl number of length classes
#' @param lbnd lower bounds of length classes
#' @param ubnd upper bounds of length classes
#' @param midpt mid points bounds of length classes
#' @param Res object containing results outputted by CalcYPRAndSPRForFMort_LB function
#' @export
GetPerRecruitGrowthPredIntervals_LB <- function(nTimeSteps, nLenCl, midpt, lbnd, ubnd, Res) {

  FemLenAtAge_lw = rep(NA,nTimeSteps); FemLenAtAge_hi = rep(NA,nTimeSteps)
  MalLenAtAge_lw = rep(NA,nTimeSteps); MalLenAtAge_hi = rep(NA,nTimeSteps)

  for (i in 1:nTimeSteps) {

    # calc length class distribution for age class
    FemProbs = Res$ModelDiag$Fish_FemNPerRecAtAge[i,] / sum(Res$ModelDiag$Fish_FemNPerRecAtAge[i,])
    MalProbs = Res$ModelDiag$Fish_MalNPerRecAtAge[i,] / sum(Res$ModelDiag$Fish_MalNPerRecAtAge[i,])

    # specify number of fish in age class
    nFem = 10000; nMal = 10000
    LenInterv = (ubnd[1] - lbnd[1]) / 2

    # females
    if (!is.nan(FemProbs[1])) {
      RandFemLenCl=rep(1:nLenCl,(as.vector(rmultinom(1,nFem,FemProbs))))
      FemDat = hist(RandFemLenCl,breaks=seq(0,nLenCl,1),right=F, plot=F)
      FreqFem = as.vector(FemDat$counts)
      RandFemLenClAtAge=rep(midpt,FreqFem)
      RandFemLen = round(RandFemLenClAtAge + runif(nFem,-LenInterv, LenInterv),0)
      FemLenAtAge_lw[i] = quantile(RandFemLen,0.025); FemLenAtAge_hi[i] = quantile(RandFemLen,0.975)
    } else {
      FemLenAtAge_lw[i]=NA; FemLenAtAge_hi[i]=NA
    }
    # males
    if (!is.nan(MalProbs[1])) {
      RandMalLenCl=rep(1:nLenCl,(as.vector(rmultinom(1,nMal,MalProbs))))
      MalDat = hist(RandMalLenCl,breaks=seq(0,nLenCl,1),right=F, plot=F)
      FreqMal = as.vector(MalDat$counts)
      RandMalLenClAtAge=rep(midpt,FreqMal)
      RandMalLen = round(RandMalLenClAtAge + runif(nMal,-LenInterv, LenInterv),0)
      MalLenAtAge_lw[i] = quantile(RandMalLen,0.025); MalLenAtAge_hi[i] = quantile(RandMalLen,0.975)
    } else {
      MalLenAtAge_lw[i]=NA; MalLenAtAge_hi[i]=NA
    }
  }

  Growth_95PLs = data.frame(FemLenAtAge_lw=FemLenAtAge_lw,
                            FemLenAtAge_hi=FemLenAtAge_hi,
                            MalLenAtAge_lw=MalLenAtAge_lw,
                            MalLenAtAge_hi=MalLenAtAge_hi)

  return(Growth_95PLs)

}

#' Plot assumed error distributions for M, F and h in per recruit analysis allowing for error
#'
#' This function plots assumed error distributions for M, F and h in per recruit analysis allowing for error, based on user inputs
#' for mean and sd of these parameters. Parameters are assumed to be normally distributed.
#'
#' @param NatMort mean for natural mortality
#' @param NatMort_sd sd for natural mortality
#' @param Current_F mean for fishing mortality
#' @param Current_F_sd sd for fishing mortality
#' @param Steepness mean for steepness of stock recruitment relationship
#' @param Steepness_sd sd for steepness of stock recruitment relationship
#' @export
PlotPerRecruit_Param_Err_Distns <- function(NatMort,NatMort_sd,Current_F,Current_F_sd,Steepness,Steepness_sd) {

  .pardefault <- par(no.readonly = TRUE) # store current par settings
  set.seed(123)
  nReps = 100000

  # Assumed distribution for natural mortality
  par(mfrow=c(2,2), mar=c(5,4,1,1), oma=c(1,1,1,1))
  xaxis_lab = expression(paste(italic("M ") ~ y^{-1}))
  yaxis_lab = "Frequency"
  M_FreqDat = hist(rnorm(nReps,NatMort,NatMort_sd),plot=F)
  xlims = Get_xaxis_scale(M_FreqDat$breaks)
  xmax = xlims$xmax
  xmin = xlims$xmin
  ylims = Get_yaxis_scale(M_FreqDat$counts)
  ymax = ylims$ymax
  hist(rnorm(nReps,NatMort,NatMort_sd),xlab=xaxis_lab, ylab=yaxis_lab,
       xlim=c(xmin,xmax), ylim=c(0,ymax), main="")
  abline(v=NatMort,lwd=2)
  legend("topright", legend=c(paste("mean =",NatMort),paste("sd =",NatMort_sd)), inset=c(0.13,0),
         lty=1, cex = 0.8, bty="n", lwd=-1, col="black")


  # Assumed distribution for fishing mortality (e.g. as estimated from catch curve)
  xaxis_lab = expression(paste(italic("F ") ~ y^{-1}))
  yaxis_lab = "Frequency"
  F_FreqDat = hist(rnorm(nReps,Current_F,Current_F_sd),plot=F)
  xlims = Get_xaxis_scale(F_FreqDat$breaks)
  xmax = xlims$xmax
  xmin = xlims$xmin
  ylims = Get_yaxis_scale(F_FreqDat$counts)
  ymax = ylims$ymax
  hist(rnorm(nReps,Current_F,Current_F_sd),xlab=xaxis_lab, ylab=yaxis_lab,
       xlim=c(xmin,xmax), ylim=c(0,ymax), main="")
  abline(v=Current_F,lwd=2)
  legend("topright", legend=c(paste("mean =",Current_F),paste("sd =",Current_F_sd)), inset=c(0.13,0),
         lty=1, cex = 0.8, bty="n", lwd=-1, col="black")

  # Assumed distribution for steepness
  xaxis_lab = "h"
  yaxis_lab = "Frequency"
  h_FreqDat = hist(rnorm(nReps,Steepness,Steepness_sd),plot=F)
  xlims = Get_xaxis_scale(h_FreqDat$breaks)
  xmax = xlims$xmax
  xmin = xlims$xmin
  ylims = Get_yaxis_scale(h_FreqDat$counts)
  ymax = ylims$ymax
  hist(rnorm(nReps,Steepness,Steepness_sd),xlab=xaxis_lab, ylab=yaxis_lab,
       xlim=c(xmin,xmax), ylim=c(0,ymax), main="")
  abline(v=Steepness,lwd=2)
  legend("topright", legend=c(paste("mean =",Steepness),paste("sd =",Steepness_sd)), inset=c(0.13,0),
         lty=1, cex = 0.8, bty="n", lwd=-1, col="black")

  par(.pardefault)

}

#' Produces WA risk matrix for relative biomass estimates from per recruit analysis, and probabilities
#'
#' This function outputs a risk matrix plot Produces WA risk matrix for relative biomass estimates from
#' age or length based per recruit analysis and associated likelihoods and risk levels (for each consequence
#' level and overall)
#'
#' @param EstBrel relative biomass point estimate
#' @param EstBrel_se standard error for relative biomass estimate
#' @param Sim_BrelVals Simulated relative biomass values (set to NA, or input as vector, which will override EstBrel, EstBrel_se)
#' @param B_targ target reference value for relative biomass
#' @param B_thresh threshold reference value for relative biomass
#' @param B_lim limit reference value for relative biomass
#' @return plot of risk matrix for relative biomass, likelihoods for consequence levels, and risks for consequence
#' levels and overall
#' @examples
#' EstBrel = 0.32
#' EstBrel_se = 0.05
#' Sim_BrelVals = NA
#' B_targ = 0.4
#' B_thresh = 0.3
#' B_lim = 0.2
#' PlotRiskMatrix_PerRecruit_RelBiom(EstBrel, EstBrel_se, Sim_BrelVals,
#'                                   B_targ, B_thresh, B_lim)
#' @export
PlotRiskMatrix_PerRecruit_RelBiom <- function(EstBrel, EstBrel_se, Sim_BrelVals,
                                              B_targ, B_thresh, B_lim) {

  # get simulated values for Brel, if they are not inputted
  if (is.na(Sim_BrelVals[1])) {
    if (is.na(EstBrel_se)) {
      cat("Problem: need to specify value for EstBrel_se")
    }
    Sim_BrelVals = rnorm(500,EstBrel,EstBrel_se)
  }
  nReps = length(Sim_BrelVals)

  plot(1:6,1:6,cex=0,ylab=NA,xlab=NA,xaxt='n',yaxt='n',bty='n', xlim=c(1,6), ylim=c(1,6))
  for (i in 1:5) {
    for (j in 1:5) {
      xval=c(i,i,6,6); yval=c(j,6,6,j)
      polygon(x=xval,y=yval, col="white",border=T)
    }
  }

  text(1.5,1.7,"Major", cex=0.8)
  text(1.5,1.3,"B<Lim", cex=0.8)
  text(1.5,2.7,"High", cex=0.8)
  text(1.5,2.3,"Lim<B<Thresh", cex=0.8)
  text(1.5,3.7,"Moderate", cex=0.8)
  text(1.5,3.3,"Thresh<B<Targ", cex=0.8)
  text(1.5,4.7,"Minor", cex=0.8)
  text(1.5,4.3,"B>Targ", cex=0.8)
  text(2.5,5.7,"Remote", cex=0.8)
  text(2.5,5.3,"<5%", cex=0.8)
  text(3.5,5.7,"Unlikely", cex=0.8)
  text(3.5,5.3,"5- <20%", cex=0.8)
  text(4.5,5.7,"Possible", cex=0.8)
  text(4.5,5.3,"20- <50%", cex=0.8)
  text(5.5,5.7,"Likely", cex=0.8)
  text(5.5,5.3,">= 50%", cex=0.8)
  mtext("Consequence",las=3,side=2,line=0,cex=1.2)
  mtext("Likelihood",las=1,side=3,line=0,cex=1.2)

  # calculate likelihoods for each consequence level
  Prob_Minor = round(length(which(Sim_BrelVals > B_targ)) / nReps,2) # B>targ
  Prob_Moderate = round(length(which(Sim_BrelVals < B_targ &
                                       Sim_BrelVals > B_thresh)) / nReps,2) #thresh<B<targ
  Prob_High = round(length(which(Sim_BrelVals < B_thresh &
                                   Sim_BrelVals > B_lim)) / nReps,2) #lim<B<thresh
  Prob_Major = round(length(which(Sim_BrelVals < B_lim)) / nReps,2) #B<lim

  Brel_likelihoods <- data.frame(Prob_Minor=Prob_Minor,
                                 Prob_Moderate=Prob_Moderate,
                                 Prob_High=Prob_High,
                                 Prob_Major=Prob_Major,
                                 nReps=nReps)

  # calculate risk scores for each consequence level
  if (Prob_Minor*100 < 5) {
    Risk_Minor = "Negligble"
    xval=c(2,2,3,3); yval=c(4,5,5,4)
    polygon(x=xval,y=yval, col="lightblue",border=T)
    text(2.5,4.5,"Negligble", cex=0.8)
  }
  if (Prob_Minor*100 >= 5 & Prob_Minor*100 < 20) {
    Risk_Minor = "Negligble"
    xval=c(3,3,4,4); yval=c(4,5,5,4)
    polygon(x=xval,y=yval, col="lightblue",border=T)
    text(3.5,4.5,"Negligble", cex=0.8)
  }
  if (Prob_Minor*100 >= 20 & Prob_Minor*100 < 50) {
    Risk_Minor = "Low"
    xval=c(4,4,5,5); yval=c(4,5,5,4)
    polygon(x=xval,y=yval, col="lightgreen",border=T)
    text(4.5,4.5,"Low", cex=0.8)
  }
  if (Prob_Minor*100 >= 50) {
    Risk_Minor = "Low"
    xval=c(5,5,6,6); yval=c(4,5,5,4)
    polygon(x=xval,y=yval, col="lightgreen",border=T)
    text(5.5,4.5,"Low", cex=0.8)
  }

  if (Prob_Moderate*100 < 5)  {
    Risk_Moderate = "Negligble"
    xval=c(2,2,3,3); yval=c(3,4,4,3)
    polygon(x=xval,y=yval, col="lightblue",border=T)
    text(2.5,3.5,"Negligble", cex=0.8)
  }
  if (Prob_Moderate*100 >= 5 & Prob_Moderate*100 < 20) {
    Risk_Moderate = "Low"
    xval=c(3,3,4,4); yval=c(3,4,4,3)
    polygon(x=xval,y=yval, col="lightgreen",border=T)
    text(3.5,3.5,"Low", cex=0.8)
  }
  if (Prob_Moderate*100 >= 20 & Prob_Moderate*100 < 50) {
    Risk_Moderate = "Medium"
    xval=c(4,4,5,5); yval=c(3,4,4,3)
    polygon(x=xval,y=yval, col="yellow",border=T)
    text(4.5,3.5,"Medium", cex=0.8)
  }
  if (Prob_Moderate*100 >= 50) {
    Risk_Moderate = "Medium"
    xval=c(5,5,6,6); yval=c(3,4,4,3)
    polygon(x=xval,y=yval, col="yellow",border=T)
    text(5.5,3.5,"Medium", cex=0.8)
  }

  if (Prob_High*100 < 5) {
    Risk_High = "Low"
    xval=c(2,2,3,3); yval=c(2,3,3,2)
    polygon(x=xval,y=yval, col="lightgreen",border=T)
    text(2.5,2.5,"Low", cex=0.8)
  }
  if (Prob_High*100 >= 5 & Prob_High*100 < 20) {
    Risk_High = "Medium"
    xval=c(3,3,4,4); yval=c(2,3,3,2)
    polygon(x=xval,y=yval, col="yellow",border=T)
    text(3.5,2.5,"Medium", cex=0.8)
  }
  if (Prob_High*100 >= 20 & Prob_High*100 < 50) {
    Risk_High = "High"
    xval=c(4,4,5,5); yval=c(2,3,3,2)
    polygon(x=xval,y=yval, col="orange",border=T)
    text(4.5,2.5,"High", cex=0.8)
  }
  if (Prob_High*100 >= 50) {
    Risk_High = "High"
    xval=c(5,5,6,6); yval=c(2,3,3,2)
    polygon(x=xval,y=yval, col="orange",border=T)
    text(5.5,2.5,"High", cex=0.8)
  }

  if (Prob_Major*100 < 5) {
    Risk_Major = "Low"
    xval=c(2,2,3,3); yval=c(1,2,2,1)
    polygon(x=xval,y=yval, col="lightgreen",border=T)
    text(2.5,1.5,"Low", cex=0.8)
  }
  if (Prob_Major*100 >= 5 & Prob_Major*100 < 20) {
    Risk_Major = "Medium"
    xval=c(3,3,4,4); yval=c(1,2,2,1)
    polygon(x=xval,y=yval, col="yellow",border=T)
    text(3.5,1.5,"Medium", cex=0.8)
  }
  if (Prob_Major*100 >= 20 & Prob_Major*100 < 50) {
    Risk_Major = "Severe"
    xval=c(4,4,5,5); yval=c(1,2,2,1)
    polygon(x=xval,y=yval, col="red",border=T)
    text(4.5,1.5,"Severe", cex=0.8)
  }
  if (Prob_Major*100 >= 50) {
    Risk_Major = "Severe"
    xval=c(5,5,6,6); yval=c(1,2,2,1)
    polygon(x=xval,y=yval, col="red",border=T)
    text(5.5,1.5,"Severe", cex=0.8)
  }

  Brel_risks <- data.frame(Risk_Minor=Risk_Minor,
                           Risk_Moderate=Risk_Moderate,
                           Risk_High=Risk_High,
                           Risk_Major=Risk_Major)

  if (length(which(Brel_risks=="Severe")>0)) {
    RiskOverall = "Severe"
  } else if (length(which(Brel_risks=="High")>0)) {
    RiskOverall = "High"
  } else if (length(which(Brel_risks=="Medium")>0)) {
    RiskOverall = "Medium"
  } else if (length(which(Brel_risks=="Low")>0)) {
    RiskOverall = "Low"
  } else {
    RiskOverall = "Negligible"
  }
  Brel_risks$RiskOverall=RiskOverall

  Results = list(Brel_risks=Brel_risks,
                 Brel_likelihoods=Brel_likelihoods)

  return(Results)

}

#' Produces WA risk matrix for fishing mortality from catch curve analysis, and probabilities
#'
#' This function outputs a risk matrix plot Produces WA risk matrix for fishing mortality estimates from age-based,
#' length-based or age and length-based catch curve analysis analysis, and associated
#' likelihoods and risk levels (for each consequence level and overall)
#'
#' @param EstFMort Fishing mortality point estimate
#' @param EstFMort_se standard error for fishing mortality estimate
#' @param Sim_FVals Simulated fishing mortality values (set to NA, or input as vector, which will override EstFMort, EstFMort_se)
#' @param F_targ target reference value for fishing mortality
#' @param F_thresh threshold reference value for fishing mortality
#' @param F_lim limit reference value for fishing mortality
#' @return plot of risk matrix for fishing mortality, likelihoods for consequence levels, and risks for consequence
#' levels and overall
#' @examples
#' NatMort = 0.15
#' EstFMort = 0.15
#' EstFMort_se = 0.025
#' Sim_FVals = NA
#' F_targ = 2/3*NatMort
#' F_thresh = NatMort
#' F_lim = 3/2*NatMort
#' PlotRiskMatrix_CatchCurve_FishMort(EstFMort, EstFMort_se, Sim_FVals,
#'                                    F_targ, F_thresh, F_lim)
#' @export
PlotRiskMatrix_CatchCurve_FishMort <- function(EstFMort, EstFMort_se, Sim_FVals,
                                               F_targ, F_thresh, F_lim) {

  # get simulated values for F, if they are not inputted
  if (is.na(Sim_FVals[1])) {
    if (is.na(EstFMort_se)) {
      cat("Problem: need to specify value for EstFMort_se")
    }
    Sim_FVals = rnorm(500,EstFMort,EstFMort_se)
  }
  nReps = length(Sim_FVals)

  # calculate likelihoods for each consequence level
  Prob_Minor = round(length(which(Sim_FVals < F_targ)) / nReps,2) # B>targ
  Prob_Moderate = round(length(which(Sim_FVals > F_targ &
                                       Sim_FVals < F_thresh)) / nReps,2) #thresh<B<targ
  Prob_High = round(length(which(Sim_FVals > F_thresh &
                                   Sim_FVals < F_lim)) / nReps,2) #lim<B<thresh
  Prob_Major = round(length(which(Sim_FVals > F_lim)) / nReps,2) #B<lim

  F_likelihoods <- data.frame(Prob_Minor=Prob_Minor,
                                 Prob_Moderate=Prob_Moderate,
                                 Prob_High=Prob_High,
                                 Prob_Major=Prob_Major,
                                 nReps=nReps)

  plot(1:6,1:6,cex=0,ylab=NA,xlab=NA,xaxt='n',yaxt='n',bty='n', xlim=c(1,6), ylim=c(1,6))
  for (i in 1:5) {
    for (j in 1:5) {
      xval=c(i,i,6,6); yval=c(j,6,6,j)
      polygon(x=xval,y=yval, col="white",border=T)
    }
  }

  text(1.5,1.7,"Major", cex=0.8)
  text(1.5,1.3,"F>Lim", cex=0.8)
  text(1.5,2.7,"High", cex=0.8)
  text(1.5,2.3,"Lim<F<Thresh", cex=0.8)
  text(1.5,3.7,"Moderate", cex=0.8)
  text(1.5,3.3,"Thresh<F<Targ", cex=0.8)
  text(1.5,4.7,"Minor", cex=0.8)
  text(1.5,4.3,"F<Targ", cex=0.8)
  text(2.5,5.7,"Remote", cex=0.8)
  text(2.5,5.3,"<5%", cex=0.8)
  text(3.5,5.7,"Unlikely", cex=0.8)
  text(3.5,5.3,"5- <20%", cex=0.8)
  text(4.5,5.7,"Possible", cex=0.8)
  text(4.5,5.3,"20- <50%", cex=0.8)
  text(5.5,5.7,"Likely", cex=0.8)
  text(5.5,5.3,">= 50%", cex=0.8)
  mtext("Consequence",las=3,side=2,line=0,cex=1.2)
  mtext("Likelihood",las=1,side=3,line=0,cex=1.2)

  # calculate risk scores for each consequence level
  if (Prob_Minor*100 < 5) {
    Risk_Minor = "Negligble"
    xval=c(2,2,3,3); yval=c(4,5,5,4)
    polygon(x=xval,y=yval, col="lightblue",border=T)
    text(2.5,4.5,"Negligble", cex=0.8)
  }
  if (Prob_Minor*100 >= 5 & Prob_Minor*100 < 20) {
    Risk_Minor = "Negligble"
    xval=c(3,3,4,4); yval=c(4,5,5,4)
    polygon(x=xval,y=yval, col="lightblue",border=T)
    text(3.5,4.5,"Negligble", cex=0.8)
  }
  if (Prob_Minor*100 >= 20 & Prob_Minor*100 < 50) {
    Risk_Minor = "Low"
    xval=c(4,4,5,5); yval=c(4,5,5,4)
    polygon(x=xval,y=yval, col="lightgreen",border=T)
    text(4.5,4.5,"Low", cex=0.8)
  }
  if (Prob_Minor*100 >= 50) {
    Risk_Minor = "Low"
    xval=c(5,5,6,6); yval=c(4,5,5,4)
    polygon(x=xval,y=yval, col="lightgreen",border=T)
    text(5.5,4.5,"Low", cex=0.8)
  }

  if (Prob_Moderate*100 < 5)  {
    Risk_Moderate = "Negligble"
    xval=c(2,2,3,3); yval=c(3,4,4,3)
    polygon(x=xval,y=yval, col="lightblue",border=T)
    text(2.5,3.5,"Negligble", cex=0.8)
  }
  if (Prob_Moderate*100 >= 5 & Prob_Moderate*100 < 20) {
    Risk_Moderate = "Low"
    xval=c(3,3,4,4); yval=c(3,4,4,3)
    polygon(x=xval,y=yval, col="lightgreen",border=T)
    text(3.5,3.5,"Low", cex=0.8)
  }
  if (Prob_Moderate*100 >= 20 & Prob_Moderate*100 < 50) {
    Risk_Moderate = "Medium"
    xval=c(4,4,5,5); yval=c(3,4,4,3)
    polygon(x=xval,y=yval, col="yellow",border=T)
    text(4.5,3.5,"Medium", cex=0.8)
  }
  if (Prob_Moderate*100 >= 50) {
    Risk_Moderate = "Medium"
    xval=c(5,5,6,6); yval=c(3,4,4,3)
    polygon(x=xval,y=yval, col="yellow",border=T)
    text(5.5,3.5,"Medium", cex=0.8)
  }

  if (Prob_High*100 < 5) {
    Risk_High = "Low"
    xval=c(2,2,3,3); yval=c(2,3,3,2)
    polygon(x=xval,y=yval, col="lightgreen",border=T)
    text(2.5,2.5,"Low", cex=0.8)
  }
  if (Prob_High*100 >= 5 & Prob_High*100 < 20) {
    Risk_High = "Medium"
    xval=c(3,3,4,4); yval=c(2,3,3,2)
    polygon(x=xval,y=yval, col="yellow",border=T)
    text(3.5,2.5,"Medium", cex=0.8)
  }
  if (Prob_High*100 >= 20 & Prob_High*100 < 50) {
    Risk_High = "High"
    xval=c(4,4,5,5); yval=c(2,3,3,2)
    polygon(x=xval,y=yval, col="orange",border=T)
    text(4.5,2.5,"High", cex=0.8)
  }
  if (Prob_High*100 >= 50) {
    Risk_High = "High"
    xval=c(5,5,6,6); yval=c(2,3,3,2)
    polygon(x=xval,y=yval, col="orange",border=T)
    text(5.5,2.5,"High", cex=0.8)
  }

  if (Prob_Major*100 < 5) {
    Risk_Major = "Low"
    xval=c(2,2,3,3); yval=c(1,2,2,1)
    polygon(x=xval,y=yval, col="lightgreen",border=T)
    text(2.5,1.5,"Low", cex=0.8)
  }
  if (Prob_Major*100 >= 5 & Prob_Major*100 < 20) {
    Risk_Major = "Medium"
    xval=c(3,3,4,4); yval=c(1,2,2,1)
    polygon(x=xval,y=yval, col="yellow",border=T)
    text(3.5,1.5,"Medium", cex=0.8)
  }
  if (Prob_Major*100 >= 20 & Prob_Major*100 < 50) {
    Risk_Major = "Severe"
    xval=c(4,4,5,5); yval=c(1,2,2,1)
    polygon(x=xval,y=yval, col="red",border=T)
    text(4.5,1.5,"Severe", cex=0.8)
  }
  if (Prob_Major*100 >= 50) {
    Risk_Major = "Severe"
    xval=c(5,5,6,6); yval=c(1,2,2,1)
    polygon(x=xval,y=yval, col="red",border=T)
    text(5.5,1.5,"Severe", cex=0.8)
  }

  F_risks <- data.frame(Risk_Minor=Risk_Minor,
                           Risk_Moderate=Risk_Moderate,
                           Risk_High=Risk_High,
                           Risk_Major=Risk_Major)

  if (length(which(F_risks=="Severe")>0)) {
    RiskOverall = "Severe"
  } else if (length(which(F_risks=="High")>0)) {
    RiskOverall = "High"
  } else if (length(which(F_risks=="Medium")>0)) {
    RiskOverall = "Medium"
  } else if (length(which(F_risks=="Low")>0)) {
    RiskOverall = "Low"
  } else {
    RiskOverall = "Negligible"
  }
  F_risks$RiskOverall=RiskOverall

  Results = list(F_risks=F_risks,
                 F_likelihoods=F_likelihoods)

  return(Results)

}


