# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

CalcExpCatchPropIntAgeGivenLength_cpp <- function(nLenCl, nAgeCl, ExpRetCatchPropLengthGivenIntAge, ExpRetCatchPropAtIntAge) {
    .Call(`_L3Assess_CalcExpCatchPropIntAgeGivenLength_cpp`, nLenCl, nAgeCl, ExpRetCatchPropLengthGivenIntAge, ExpRetCatchPropAtIntAge)
}

CalcNLLCondAgeAtLength_cpp <- function(nLenCl, nAgeCl, ObsCatchFreqAtLengthAndIntAge, ExpRetCatchPropIntAgeGivenLength) {
    .Call(`_L3Assess_CalcNLLCondAgeAtLength_cpp`, nLenCl, nAgeCl, ObsCatchFreqAtLengthAndIntAge, ExpRetCatchPropIntAgeGivenLength)
}

CalcCatches_AgeAndLengthBasedCatchCurves_cpp <- function(params, NatMort, RecLenDist, InitRecNumber, MaxAge, TimeStep, nTimeSteps, nLenCl, midpt, RetAtLength, SelAtLength, DiscMort, LTM) {
    .Call(`_L3Assess_CalcCatches_AgeAndLengthBasedCatchCurves_cpp`, params, NatMort, RecLenDist, InitRecNumber, MaxAge, TimeStep, nTimeSteps, nLenCl, midpt, RetAtLength, SelAtLength, DiscMort, LTM)
}

CalcLTM_cpp <- function(TimeStepGrowthSizeInc, CVSizeAtAge, lbnd, midpt, ubnd, nLenCl) {
    .Call(`_L3Assess_CalcLTM_cpp`, TimeStepGrowthSizeInc, CVSizeAtAge, lbnd, midpt, ubnd, nLenCl)
}

UpdateGrowthAndSurvival_cpp <- function(ReprodPattern, TimeStep, nTimeSteps, nLenCl, InitRatioFem, RecLenDist, NatMort, FemZAtLen, MalZAtLen, PropFemAtLen, LTM_Fem, LTM_Mal, FemWtAtLen, MalWtAtLen, ReprodScale, FemPropMatAtLen, MalPropMatAtLen) {
    .Call(`_L3Assess_UpdateGrowthAndSurvival_cpp`, ReprodPattern, TimeStep, nTimeSteps, nLenCl, InitRatioFem, RecLenDist, NatMort, FemZAtLen, MalZAtLen, PropFemAtLen, LTM_Fem, LTM_Mal, FemWtAtLen, MalWtAtLen, ReprodScale, FemPropMatAtLen, MalPropMatAtLen)
}

