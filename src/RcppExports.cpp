// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// CalcNLLCondAgeAtLength_cpp
double CalcNLLCondAgeAtLength_cpp(const int nLenCl, const int MaxAge, NumericMatrix ObsCatchFreqAtLengthAndAge, NumericMatrix ExpCatchPropLengthGivenAge, NumericVector ExpCatchPropAtAge);
RcppExport SEXP _L3Assess_CalcNLLCondAgeAtLength_cpp(SEXP nLenClSEXP, SEXP MaxAgeSEXP, SEXP ObsCatchFreqAtLengthAndAgeSEXP, SEXP ExpCatchPropLengthGivenAgeSEXP, SEXP ExpCatchPropAtAgeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type nLenCl(nLenClSEXP);
    Rcpp::traits::input_parameter< const int >::type MaxAge(MaxAgeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ObsCatchFreqAtLengthAndAge(ObsCatchFreqAtLengthAndAgeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ExpCatchPropLengthGivenAge(ExpCatchPropLengthGivenAgeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ExpCatchPropAtAge(ExpCatchPropAtAgeSEXP);
    rcpp_result_gen = Rcpp::wrap(CalcNLLCondAgeAtLength_cpp(nLenCl, MaxAge, ObsCatchFreqAtLengthAndAge, ExpCatchPropLengthGivenAge, ExpCatchPropAtAge));
    return rcpp_result_gen;
END_RCPP
}
// CalcCatches_AgeAndLengthBasedCatchCurves_cpp
List CalcCatches_AgeAndLengthBasedCatchCurves_cpp(NumericVector params, const double NatMort, NumericVector RecLenDist, const double InitRecNumber, const int MaxAge, const int nLenCl, NumericVector midpt, NumericVector RetAtLength, NumericVector SelAtLength, const double DiscMort, NumericMatrix LTM);
RcppExport SEXP _L3Assess_CalcCatches_AgeAndLengthBasedCatchCurves_cpp(SEXP paramsSEXP, SEXP NatMortSEXP, SEXP RecLenDistSEXP, SEXP InitRecNumberSEXP, SEXP MaxAgeSEXP, SEXP nLenClSEXP, SEXP midptSEXP, SEXP RetAtLengthSEXP, SEXP SelAtLengthSEXP, SEXP DiscMortSEXP, SEXP LTMSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const double >::type NatMort(NatMortSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type RecLenDist(RecLenDistSEXP);
    Rcpp::traits::input_parameter< const double >::type InitRecNumber(InitRecNumberSEXP);
    Rcpp::traits::input_parameter< const int >::type MaxAge(MaxAgeSEXP);
    Rcpp::traits::input_parameter< const int >::type nLenCl(nLenClSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type midpt(midptSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type RetAtLength(RetAtLengthSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type SelAtLength(SelAtLengthSEXP);
    Rcpp::traits::input_parameter< const double >::type DiscMort(DiscMortSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type LTM(LTMSEXP);
    rcpp_result_gen = Rcpp::wrap(CalcCatches_AgeAndLengthBasedCatchCurves_cpp(params, NatMort, RecLenDist, InitRecNumber, MaxAge, nLenCl, midpt, RetAtLength, SelAtLength, DiscMort, LTM));
    return rcpp_result_gen;
END_RCPP
}
// CalcLTM_cpp
NumericMatrix CalcLTM_cpp(NumericVector TimeStepGrowthSizeInc, const double CVSizeAtAge, NumericVector lbnd, NumericVector midpt, NumericVector ubnd, const int nLenCl);
RcppExport SEXP _L3Assess_CalcLTM_cpp(SEXP TimeStepGrowthSizeIncSEXP, SEXP CVSizeAtAgeSEXP, SEXP lbndSEXP, SEXP midptSEXP, SEXP ubndSEXP, SEXP nLenClSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type TimeStepGrowthSizeInc(TimeStepGrowthSizeIncSEXP);
    Rcpp::traits::input_parameter< const double >::type CVSizeAtAge(CVSizeAtAgeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lbnd(lbndSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type midpt(midptSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ubnd(ubndSEXP);
    Rcpp::traits::input_parameter< const int >::type nLenCl(nLenClSEXP);
    rcpp_result_gen = Rcpp::wrap(CalcLTM_cpp(TimeStepGrowthSizeInc, CVSizeAtAge, lbnd, midpt, ubnd, nLenCl));
    return rcpp_result_gen;
END_RCPP
}
// UpdateGrowthAndSurvival_cpp
List UpdateGrowthAndSurvival_cpp(const int ReprodPattern, const double TimeStep, const int nTimeSteps, const int nLenCl, const double InitRatioFem, NumericMatrix RecLenDist, const double NatMort, NumericVector FemZAtLen, NumericVector MalZAtLen, NumericVector PropFemAtLen, NumericMatrix LTM_Fem, NumericMatrix LTM_Mal, NumericVector FemWtAtLen, NumericVector MalWtAtLen, const double ReprodScale, NumericVector FemPropMatAtLen, NumericVector MalPropMatAtLen);
RcppExport SEXP _L3Assess_UpdateGrowthAndSurvival_cpp(SEXP ReprodPatternSEXP, SEXP TimeStepSEXP, SEXP nTimeStepsSEXP, SEXP nLenClSEXP, SEXP InitRatioFemSEXP, SEXP RecLenDistSEXP, SEXP NatMortSEXP, SEXP FemZAtLenSEXP, SEXP MalZAtLenSEXP, SEXP PropFemAtLenSEXP, SEXP LTM_FemSEXP, SEXP LTM_MalSEXP, SEXP FemWtAtLenSEXP, SEXP MalWtAtLenSEXP, SEXP ReprodScaleSEXP, SEXP FemPropMatAtLenSEXP, SEXP MalPropMatAtLenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type ReprodPattern(ReprodPatternSEXP);
    Rcpp::traits::input_parameter< const double >::type TimeStep(TimeStepSEXP);
    Rcpp::traits::input_parameter< const int >::type nTimeSteps(nTimeStepsSEXP);
    Rcpp::traits::input_parameter< const int >::type nLenCl(nLenClSEXP);
    Rcpp::traits::input_parameter< const double >::type InitRatioFem(InitRatioFemSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type RecLenDist(RecLenDistSEXP);
    Rcpp::traits::input_parameter< const double >::type NatMort(NatMortSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type FemZAtLen(FemZAtLenSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type MalZAtLen(MalZAtLenSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type PropFemAtLen(PropFemAtLenSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type LTM_Fem(LTM_FemSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type LTM_Mal(LTM_MalSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type FemWtAtLen(FemWtAtLenSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type MalWtAtLen(MalWtAtLenSEXP);
    Rcpp::traits::input_parameter< const double >::type ReprodScale(ReprodScaleSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type FemPropMatAtLen(FemPropMatAtLenSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type MalPropMatAtLen(MalPropMatAtLenSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateGrowthAndSurvival_cpp(ReprodPattern, TimeStep, nTimeSteps, nLenCl, InitRatioFem, RecLenDist, NatMort, FemZAtLen, MalZAtLen, PropFemAtLen, LTM_Fem, LTM_Mal, FemWtAtLen, MalWtAtLen, ReprodScale, FemPropMatAtLen, MalPropMatAtLen));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_L3Assess_CalcNLLCondAgeAtLength_cpp", (DL_FUNC) &_L3Assess_CalcNLLCondAgeAtLength_cpp, 5},
    {"_L3Assess_CalcCatches_AgeAndLengthBasedCatchCurves_cpp", (DL_FUNC) &_L3Assess_CalcCatches_AgeAndLengthBasedCatchCurves_cpp, 11},
    {"_L3Assess_CalcLTM_cpp", (DL_FUNC) &_L3Assess_CalcLTM_cpp, 6},
    {"_L3Assess_UpdateGrowthAndSurvival_cpp", (DL_FUNC) &_L3Assess_UpdateGrowthAndSurvival_cpp, 17},
    {NULL, NULL, 0}
};

RcppExport void R_init_L3Assess(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
