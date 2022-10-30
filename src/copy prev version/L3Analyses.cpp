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
