#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
double CalcNLLCondAgeAtLength_cpp(const int nLenCl, const int MaxAge,
                                  NumericMatrix ObsCatchFreqAtLengthAndAge,
                                  NumericMatrix ExpCatchPropLengthGivenAge,
                                  NumericVector ExpCatchPropAtAge) {


  //conditional age at length likelihood (based on mathematical description
  //given by Piner et al., 2015) - for age and length-based catch curve
  int i;
  int t;
  double temp;
  double NLL;
  NumericVector FractDenom(nLenCl);
  NumericMatrix ExpCatchPropAgeGivenLength(MaxAge, nLenCl);

  for (i=0; i<nLenCl; i++) {
    for (t=0; t<MaxAge; t++) {
      FractDenom(i) = FractDenom(i) + (ExpCatchPropLengthGivenAge(t,i) * ExpCatchPropAtAge(t));
    } // t

    for (t=0; t<MaxAge; t++) {
      ExpCatchPropAgeGivenLength(t,i) = (ExpCatchPropLengthGivenAge(t,i) * ExpCatchPropAtAge(t))
                                          / FractDenom(i);
    } // t
  } // i

  temp = 0.0;
  for (i=0; i<nLenCl; i++) {
    for (t=0; t<MaxAge; t++) {
      temp = temp + (ObsCatchFreqAtLengthAndAge(t,i) * log(ExpCatchPropAgeGivenLength(t,i) + 1E-4));
    } // t
  } // i

  NLL = -temp;

  return(NLL);

}

// [[Rcpp::export]]
List CalcCatches_AgeAndLengthBasedCatchCurves_cpp(NumericVector params, const double NatMort,
                                             NumericVector RecLenDist, const double InitRecNumber, const int MaxAge,
                                             const int nLenCl, NumericVector midpt,
                                             NumericVector RetAtLength, NumericVector SelAtLength, const double DiscMort,
                                             NumericMatrix LTM) {

  int i;
  int ii;
  int t;
  double FishMort;
  NumericVector SelLandAtLength(nLenCl);
  NumericVector SelDiscAtLength(nLenCl);
  NumericVector FAtLen(nLenCl); // fishing mortality at length associated with landings and discards
  NumericVector FAtLenCapt(nLenCl); // fishing mortality at length associated with all fish captures
  NumericVector FAtLenReten(nLenCl); // fishing mortality at length associated with all fish captures
  NumericVector FAtLenDisc(nLenCl); // fishing mortality at length associated with all fish captures
  NumericVector ZAtLen(nLenCl);
  NumericVector Fish_NPerRec(nLenCl);
  NumericMatrix Fish_NPerRecAtAge(MaxAge, nLenCl);
  NumericMatrix Fish_SurvPerRecAtAge(MaxAge, nLenCl);
  NumericMatrix RetCatch(MaxAge, nLenCl); // retained catches at age and length
  NumericMatrix DiscCatch(MaxAge, nLenCl); // retained catches at age and length
  NumericMatrix TotCatch(MaxAge, nLenCl); // total (released + retained) catches at age and length
  NumericVector DiscCatchAtLen(nLenCl); // total released catches at length
  NumericVector RetCatchAtLen(nLenCl); // total retained catches at length
  NumericVector TotCatchAtLen(nLenCl); // total (released + retained) catches at length

  //inverse of logit transformation
  FishMort = 1/(1+exp(-params(0)));

    // per recruit numbers surviving after natural mortality
    for (i=0; i<nLenCl; i++) {
      Fish_NPerRecAtAge(0,i) = RecLenDist(i) * InitRecNumber;

      Fish_NPerRec(i) = Fish_NPerRec(i) + Fish_NPerRecAtAge(0,i);

      // selectivity of landings
      SelLandAtLength(i) = SelAtLength(i) * RetAtLength(i);

      // selectivity of discards
      SelDiscAtLength(i) = SelAtLength(i) * (1.0 - RetAtLength(i));

      // fishing mortality at length (allowing for both retention, and discarding)
      FAtLen(i) = (SelLandAtLength(i) + (DiscMort * SelDiscAtLength(i))) * FishMort;

      // total mortality at length
      ZAtLen(i) = NatMort + FAtLen(i);

      // fishing mortality used to calculate total numbers of fish caught
      FAtLenCapt(i) = SelAtLength(i) * FishMort;

      // fishing mortality used to calculate numbers of fish caught and retained
      FAtLenReten(i) = SelLandAtLength(i) * FishMort;

      // fishing mortality used to calculated numbers of fish discarded
      FAtLenDisc(i) = SelDiscAtLength(i) * FishMort;

      // calculate retained catch at length (in numbers)
      RetCatch(0,i) = Fish_NPerRecAtAge(0,i) * (FAtLenReten(i) / ZAtLen(i)) * (1 - exp(-ZAtLen(i)));
      RetCatchAtLen(i) = RetCatchAtLen(i) + RetCatch(0,i);

      // calculate discarded catch at length  (in numbers)
      DiscCatch(0,i) = Fish_NPerRecAtAge(0,i) * (FAtLenDisc(i) / ZAtLen(i)) * (1 - exp(-ZAtLen(i)));
      DiscCatchAtLen(i) = DiscCatchAtLen(i) + DiscCatch(0,i);

      // approximate total catches (caught + released)  (in numbers)
      TotCatch(0,i) = Fish_NPerRecAtAge(0,i) * (FAtLenCapt(i) / ZAtLen(i)) * (1 - exp(-ZAtLen(i)));
      TotCatchAtLen(i) = TotCatchAtLen(i) + TotCatch(0,i);
      //std::cout << " i " << i << " RecLenDist " << RecLenDist(i) <<
      //  " SelAtLength " << SelAtLength(i) <<" FAtLen " << FAtLen(i) << std::endl;
    }

  // apply mortality to calculate survival
  for (t=1; t<MaxAge; t++) {
    for (i=0; i<nLenCl; i++) {
      if (t < MaxAge-1) {
        Fish_SurvPerRecAtAge(t,i) = Fish_NPerRecAtAge(t-1,i) * exp(-ZAtLen(i));
      } // end if
      if (t == MaxAge-1) {
        Fish_SurvPerRecAtAge(t,i) = Fish_NPerRecAtAge(t-1,i) * exp(-ZAtLen(i));
        // (1 - exp(-ZAtLen(i)));
      } // end if
      //std::cout << " t " << t << " i " << i << " Fish_SurvPerRecAtAge " << Fish_SurvPerRecAtAge(t,i) << std::endl;
    } // i


    // grow fish
    for (ii=0; ii<nLenCl; ii++) {  // starting length class
      for (i=0; i<nLenCl; i++) { // ending length class
        Fish_NPerRecAtAge(t,i) = Fish_NPerRecAtAge(t,i) + (Fish_SurvPerRecAtAge(t,ii) * LTM(i,ii));
      } // i
    } // ii

    for (i=0; i<nLenCl; i++) {
      // add numbers per recruit at length, over timesteps to get marginal length distribution
      Fish_NPerRec(i) = Fish_NPerRec(i) + Fish_NPerRecAtAge(t,i);
      // calculate catch at length for timestep (in numbers)
      RetCatch(t,i) = Fish_NPerRecAtAge(t,i) * (FAtLenReten(i) / ZAtLen(i)) * (1 - exp(-ZAtLen(i)));
      RetCatchAtLen(i) = RetCatchAtLen(i) + RetCatch(t,i);

      // calculate discarded catch at length  (in numbers)
      DiscCatch(t,i) = Fish_NPerRecAtAge(t,i) * (FAtLenDisc(i) / ZAtLen(i)) * (1 - exp(-ZAtLen(i)));
      DiscCatchAtLen(i) = DiscCatchAtLen(i) + DiscCatch(t,i);

      // approximate total catches (caught + released)
      TotCatch(t,i) = Fish_NPerRecAtAge(t,i) * (FAtLenCapt(i) / ZAtLen(i)) * (1 - exp(-ZAtLen(i)));
      TotCatchAtLen(i) = TotCatchAtLen(i) + TotCatch(t,i);
    } // i

  } // t


      return List::create(Named("SelLandAtLength") = SelLandAtLength,
                          Named("RetCatch") = RetCatch,
                          Named("RetCatchAtLen") = RetCatchAtLen,
                          Named("DiscCatch") = DiscCatch,
                          Named("DiscCatchAtLen") = DiscCatchAtLen,
                          Named("TotCatch") = TotCatch,
                          Named("TotCatchAtLen") = TotCatchAtLen,
                          Named("Fish_NPerRec") = Fish_NPerRec);

} // end function


// [[Rcpp::export]]
NumericMatrix CalcLTM_cpp(NumericVector AnnGrowthSizeInc, const double CVSizeAtAge,
                                  NumericVector lbnd, NumericVector midpt, NumericVector ubnd,
                                  const int nLenCl) {

  // set up data frame for length-transition matrix
  int i;
  int ii;
  NumericMatrix LTM(nLenCl, nLenCl);
  NumericVector MeanEndingLength(nLenCl);
  NumericVector StDev(nLenCl);
  NumericVector temp;

    MeanEndingLength = midpt + AnnGrowthSizeInc;
    StDev = MeanEndingLength * CVSizeAtAge;

      for (ii=0; ii<nLenCl; ii++) { // starting length class, upper bound - lower bound

        temp = pnorm(ubnd, MeanEndingLength(ii), StDev(ii)) -
          pnorm(lbnd, MeanEndingLength(ii), StDev(ii));

        for (i=0; i<nLenCl; i++) { // starting length class
          LTM(i,ii) = temp(i);
        }
      }

      return LTM;
}

// [[Rcpp::export]]
List UpdateGrowthAndSurvival_cpp(const int ReprodPattern, const double TimeStep, const int nTimeSteps, const int nLenCl,
                                 const double InitRatioFem, NumericMatrix RecLenDist, const double NatMort, NumericVector FemZAtLen,
                                 NumericVector MalZAtLen, NumericVector PropFemAtLen, NumericMatrix LTM_Fem,
                                 NumericMatrix LTM_Mal, NumericVector FemWtAtLen, NumericVector MalWtAtLen,
                                 NumericVector FemPropMatAtLen, NumericVector MalPropMatAtLen) {

  int t;
  int i;
  int ii;
  NumericMatrix Unfish_FemNPerRecAtAge(nTimeSteps,nLenCl);
  NumericMatrix Fish_FemNPerRecAtAge(nTimeSteps,nLenCl);
  NumericMatrix Unfish_MalNPerRecAtAge(nTimeSteps,nLenCl);
  NumericMatrix Fish_MalNPerRecAtAge(nTimeSteps,nLenCl);
  NumericMatrix Unfish_FemSurvPerRecAtAge(nTimeSteps,nLenCl);
  NumericMatrix Fish_FemSurvPerRecAtAge(nTimeSteps,nLenCl);
  NumericMatrix Unfish_MalSurvPerRecAtAge(nTimeSteps,nLenCl);
  NumericMatrix Fish_MalSurvPerRecAtAge(nTimeSteps,nLenCl);
  NumericMatrix tempUnfish_FemSurvPerRecAtAge(nTimeSteps,nLenCl);
  NumericMatrix tempFish_FemSurvPerRecAtAge(nTimeSteps,nLenCl);
  NumericMatrix tempUnfish_MalSurvPerRecAtAge(nTimeSteps,nLenCl);
  NumericMatrix tempFish_MalSurvPerRecAtAge(nTimeSteps,nLenCl);
  NumericMatrix Unfish_FemBiomPerRecAtAge(nTimeSteps,nLenCl);
  NumericMatrix Fish_FemBiomPerRecAtAge(nTimeSteps,nLenCl);
  NumericMatrix Unfish_MalBiomPerRecAtAge(nTimeSteps,nLenCl);
  NumericMatrix Fish_MalBiomPerRecAtAge(nTimeSteps,nLenCl);

  NumericVector Unfish_FemNPerRec(nLenCl);
  NumericVector Fish_FemNPerRec(nLenCl);
  NumericVector Unfish_MalNPerRec(nLenCl);
  NumericVector Fish_MalNPerRec(nLenCl);
  NumericVector Unfish_FemBiomAtAge(nTimeSteps);
  NumericVector Fish_FemBiomAtAge(nTimeSteps);
  NumericVector Unfish_MalBiomAtAge(nTimeSteps);
  NumericVector Fish_MalBiomAtAge(nTimeSteps);

  NumericVector tempVec1(nLenCl);
  NumericVector tempVec2(nLenCl);
  NumericVector tempVec3(nLenCl);
  NumericVector tempVec4(nLenCl);

  for (t=0; t<nTimeSteps; t++) {
    if(t==0) {
      for (i=0; i<nLenCl; i++) {
        Unfish_FemNPerRecAtAge(t,i) = RecLenDist(0,i) * InitRatioFem;
        Fish_FemNPerRecAtAge(t,i) = Unfish_FemNPerRecAtAge(t,i);
        Unfish_FemNPerRec(i) = Unfish_FemNPerRec(i) + Unfish_FemNPerRecAtAge(t,i);
        Fish_FemNPerRec(i) = Fish_FemNPerRec(i) + Fish_FemNPerRecAtAge(t,i);

        Unfish_MalNPerRecAtAge(t,i) = RecLenDist(1,i) * (1.0 - InitRatioFem);
        Fish_MalNPerRecAtAge(t,i) = Unfish_MalNPerRecAtAge(t,i);
        Unfish_MalNPerRec(i) = Unfish_MalNPerRec(i) + Unfish_MalNPerRecAtAge(t,i);
        Fish_MalNPerRec(i) = Fish_MalNPerRec(i) + Fish_MalNPerRecAtAge(t,i);
      } // i
    } else { // t>0

      // calculate survival
      if (t < nTimeSteps-1) {
        for (i=0; i<nLenCl; i++) {
          // females
          Unfish_FemSurvPerRecAtAge(t,i) = Unfish_FemNPerRecAtAge(t-1,i) * exp(-NatMort*TimeStep);
          Fish_FemSurvPerRecAtAge(t,i) = Fish_FemNPerRecAtAge(t-1,i) * exp(-FemZAtLen(i)*TimeStep);

          // males
          Unfish_MalSurvPerRecAtAge(t,i) = Unfish_MalNPerRecAtAge(t-1,i) * exp(-NatMort*TimeStep);
          Fish_MalSurvPerRecAtAge(t,i) = Fish_MalNPerRecAtAge(t-1,i) * exp(-MalZAtLen(i)*TimeStep);
        }
      }
      if (t == nTimeSteps-1) {
        for (i=0; i<nLenCl; i++) {
          // females
          Unfish_FemSurvPerRecAtAge(t,i) = Unfish_FemNPerRecAtAge(t-1,i) * exp(-NatMort*TimeStep) / (1 - exp(-NatMort));
          Fish_FemSurvPerRecAtAge(t,i) = Fish_FemNPerRecAtAge(t-1,i) * exp(-FemZAtLen(i)*TimeStep) / (1 - exp(-FemZAtLen(i)));
          // males
          Unfish_MalSurvPerRecAtAge(t,i) = Unfish_MalNPerRecAtAge(t-1,i) * exp(-NatMort*TimeStep) / (1 - exp(-NatMort));
          Fish_MalSurvPerRecAtAge(t,i) = Fish_MalNPerRecAtAge(t-1,i) * exp(-MalZAtLen(i)*TimeStep) / (1 - exp(-MalZAtLen(i)));
        }
      }

      if (ReprodPattern > 1) { // hermaphroditic species
        // assuming that, regardless of fishing pattern, the sex ratio at age is determined by
        // by the logistic curve describing probability of sex change at age
        for (i=0; i<nLenCl; i++) {
          tempUnfish_FemSurvPerRecAtAge(t,i) = PropFemAtLen(i) * (Unfish_FemSurvPerRecAtAge(t,i) + Unfish_MalSurvPerRecAtAge(t,i));
          tempUnfish_MalSurvPerRecAtAge(t,i) = (1 - PropFemAtLen(i)) * (Unfish_FemSurvPerRecAtAge(t,i) + Unfish_MalSurvPerRecAtAge(t,i));
          Unfish_FemSurvPerRecAtAge(t,i) = tempUnfish_FemSurvPerRecAtAge(t,i);
          Unfish_MalSurvPerRecAtAge(t,i) = tempUnfish_MalSurvPerRecAtAge(t,i);

          tempFish_FemSurvPerRecAtAge(t,i) = PropFemAtLen(i) * (Fish_FemSurvPerRecAtAge(t,i) + Fish_MalSurvPerRecAtAge(t,i));
          tempFish_MalSurvPerRecAtAge(t,i) = (1 - PropFemAtLen(i)) * (Fish_FemSurvPerRecAtAge(t,i) + Fish_MalSurvPerRecAtAge(t,i));
          Fish_FemSurvPerRecAtAge(t,i) = tempFish_FemSurvPerRecAtAge(t,i);
          Fish_MalSurvPerRecAtAge(t,i) = tempFish_MalSurvPerRecAtAge(t,i);
        } // i
      } // reproductive pattern

      // grow fish
      for (ii=0; ii<nLenCl; ii++) {  // starting length class
        for (i=0; i<nLenCl; i++) { // ending length class
          tempVec1(i) = tempVec1(i) + (Unfish_FemSurvPerRecAtAge(t,ii) * LTM_Fem(i,ii));
          tempVec2(i) = tempVec2(i) + (Unfish_MalSurvPerRecAtAge(t,ii) * LTM_Mal(i,ii));
          tempVec3(i) = tempVec3(i) + (Fish_FemSurvPerRecAtAge(t,ii) * LTM_Fem(i,ii));
          tempVec4(i) = tempVec4(i) + (Fish_MalSurvPerRecAtAge(t,ii) * LTM_Mal(i,ii));
        } // i
      } // ii

      for (i=0; i<nLenCl; i++) {
        Unfish_FemNPerRecAtAge(t,i) = tempVec1(i);
        Unfish_MalNPerRecAtAge(t,i) = tempVec2(i);
        Fish_FemNPerRecAtAge(t,i) = tempVec3(i);
        Fish_MalNPerRecAtAge(t,i) = tempVec4(i);
        tempVec1(i)=0;
        tempVec2(i)=0;
        tempVec3(i)=0;
        tempVec4(i)=0;
        Unfish_FemNPerRec(i) = Unfish_FemNPerRec(i) + Unfish_FemNPerRecAtAge(t,i);
        Unfish_MalNPerRec(i) = Unfish_MalNPerRec(i) + Unfish_MalNPerRecAtAge(t,i);
        Fish_FemNPerRec(i) = Fish_FemNPerRec(i) + Fish_FemNPerRecAtAge(t,i);
        Fish_MalNPerRec(i) = Fish_MalNPerRec(i) + Fish_MalNPerRecAtAge(t,i);
      }
    } // else
  } // t

  // calculate unfished mature biomass at age
  for (t=0; t<nTimeSteps; t++) {
    for (i=0; i<nLenCl; i++) {
      Unfish_FemBiomPerRecAtAge(t,i) = Unfish_FemNPerRecAtAge(t,i) * FemWtAtLen(i) * FemPropMatAtLen(i);
      Fish_FemBiomPerRecAtAge(t,i) = Fish_FemNPerRecAtAge(t,i) * FemWtAtLen(i) * FemPropMatAtLen(i);
      Unfish_MalBiomPerRecAtAge(t,i) = Unfish_MalNPerRecAtAge(t,i) * MalWtAtLen(i) * MalPropMatAtLen(i);
      Fish_MalBiomPerRecAtAge(t,i) = Fish_MalNPerRecAtAge(t,i) * MalWtAtLen(i) * MalPropMatAtLen(i);
      Unfish_FemBiomAtAge(t) = Unfish_FemBiomAtAge(t) + Unfish_FemBiomPerRecAtAge(t,i);
      Fish_FemBiomAtAge(t) = Fish_FemBiomAtAge(t) + Fish_FemBiomPerRecAtAge(t,i);
      Unfish_MalBiomAtAge(t) = Unfish_MalBiomAtAge(t) + Unfish_MalBiomPerRecAtAge(t,i);
      Fish_MalBiomAtAge(t) = Fish_MalBiomAtAge(t) + Fish_MalBiomPerRecAtAge(t,i);
    }
  } // t

  return List::create(Named("Unfish_FemNPerRecAtAge") = Unfish_FemNPerRecAtAge,
                      Named("Fish_FemNPerRecAtAge") = Fish_FemNPerRecAtAge,
                      Named("Unfish_MalNPerRecAtAge") = Unfish_MalNPerRecAtAge,
                      Named("Fish_MalNPerRecAtAge") = Fish_MalNPerRecAtAge,
                      Named("Unfish_FemNPerRec") = Unfish_FemNPerRec,
                      Named("Fish_FemNPerRec") = Fish_FemNPerRec,
                      Named("Unfish_MalNPerRec") = Unfish_MalNPerRec,
                      Named("Fish_MalNPerRec") = Fish_MalNPerRec,
                      Named("Unfish_FemBiomPerRecAtAge") = Unfish_FemBiomPerRecAtAge,
                      Named("Fish_FemBiomPerRecAtAge") = Fish_FemBiomPerRecAtAge,
                      Named("Unfish_MalBiomPerRecAtAge") = Unfish_MalBiomPerRecAtAge,
                      Named("Unfish_FemBiomAtAge") = Unfish_FemBiomAtAge,
                      Named("Fish_FemBiomAtAge") = Fish_FemBiomAtAge,
                      Named("Unfish_MalBiomAtAge") = Unfish_MalBiomAtAge,
                      Named("Fish_MalBiomAtAge") = Fish_MalBiomAtAge);
}
