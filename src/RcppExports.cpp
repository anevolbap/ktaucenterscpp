// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rhoAux
NumericVector rhoAux(NumericVector x, double cc);
RcppExport SEXP _ktaucenterscpp_rhoAux(SEXP xSEXP, SEXP ccSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type cc(ccSEXP);
    rcpp_result_gen = Rcpp::wrap(rhoAux(x, cc));
    return rcpp_result_gen;
END_RCPP
}
// rhoOpt
NumericVector rhoOpt(NumericVector x, double cc);
RcppExport SEXP _ktaucenterscpp_rhoOpt(SEXP xSEXP, SEXP ccSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type cc(ccSEXP);
    rcpp_result_gen = Rcpp::wrap(rhoOpt(x, cc));
    return rcpp_result_gen;
END_RCPP
}
// psiAux
NumericVector psiAux(NumericVector x, double cc);
RcppExport SEXP _ktaucenterscpp_psiAux(SEXP xSEXP, SEXP ccSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type cc(ccSEXP);
    rcpp_result_gen = Rcpp::wrap(psiAux(x, cc));
    return rcpp_result_gen;
END_RCPP
}
// psiOpt
NumericVector psiOpt(NumericVector x, double cc);
RcppExport SEXP _ktaucenterscpp_psiOpt(SEXP xSEXP, SEXP ccSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type cc(ccSEXP);
    rcpp_result_gen = Rcpp::wrap(psiOpt(x, cc));
    return rcpp_result_gen;
END_RCPP
}
// derpsiAux
NumericVector derpsiAux(NumericVector x, double cc);
RcppExport SEXP _ktaucenterscpp_derpsiAux(SEXP xSEXP, SEXP ccSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type cc(ccSEXP);
    rcpp_result_gen = Rcpp::wrap(derpsiAux(x, cc));
    return rcpp_result_gen;
END_RCPP
}
// derpsiOpt
NumericVector derpsiOpt(NumericVector x, double cc);
RcppExport SEXP _ktaucenterscpp_derpsiOpt(SEXP xSEXP, SEXP ccSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type cc(ccSEXP);
    rcpp_result_gen = Rcpp::wrap(derpsiOpt(x, cc));
    return rcpp_result_gen;
END_RCPP
}
// normal_consistency_constants
double normal_consistency_constants(int p);
RcppExport SEXP _ktaucenterscpp_normal_consistency_constants(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(normal_consistency_constants(p));
    return rcpp_result_gen;
END_RCPP
}
// constC1
NumericVector constC1(NumericVector p);
RcppExport SEXP _ktaucenterscpp_constC1(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(constC1(p));
    return rcpp_result_gen;
END_RCPP
}
// constC2
NumericVector constC2(NumericVector p);
RcppExport SEXP _ktaucenterscpp_constC2(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(constC2(p));
    return rcpp_result_gen;
END_RCPP
}
// my_median
double my_median(NumericVector x);
RcppExport SEXP _ktaucenterscpp_my_median(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(my_median(x));
    return rcpp_result_gen;
END_RCPP
}
// Mscale
double Mscale(NumericVector u, double b, double cc);
RcppExport SEXP _ktaucenterscpp_Mscale(SEXP uSEXP, SEXP bSEXP, SEXP ccSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type cc(ccSEXP);
    rcpp_result_gen = Rcpp::wrap(Mscale(u, b, cc));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ktaucenterscpp_rhoAux", (DL_FUNC) &_ktaucenterscpp_rhoAux, 2},
    {"_ktaucenterscpp_rhoOpt", (DL_FUNC) &_ktaucenterscpp_rhoOpt, 2},
    {"_ktaucenterscpp_psiAux", (DL_FUNC) &_ktaucenterscpp_psiAux, 2},
    {"_ktaucenterscpp_psiOpt", (DL_FUNC) &_ktaucenterscpp_psiOpt, 2},
    {"_ktaucenterscpp_derpsiAux", (DL_FUNC) &_ktaucenterscpp_derpsiAux, 2},
    {"_ktaucenterscpp_derpsiOpt", (DL_FUNC) &_ktaucenterscpp_derpsiOpt, 2},
    {"_ktaucenterscpp_normal_consistency_constants", (DL_FUNC) &_ktaucenterscpp_normal_consistency_constants, 1},
    {"_ktaucenterscpp_constC1", (DL_FUNC) &_ktaucenterscpp_constC1, 1},
    {"_ktaucenterscpp_constC2", (DL_FUNC) &_ktaucenterscpp_constC2, 1},
    {"_ktaucenterscpp_my_median", (DL_FUNC) &_ktaucenterscpp_my_median, 1},
    {"_ktaucenterscpp_Mscale", (DL_FUNC) &_ktaucenterscpp_Mscale, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_ktaucenterscpp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
