// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// iouC
double iouC(arma::mat img1, arma::mat img2);
RcppExport SEXP _spat1f_iouC(SEXP img1SEXP, SEXP img2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type img1(img1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type img2(img2SEXP);
    rcpp_result_gen = Rcpp::wrap(iouC(img1, img2));
    return rcpp_result_gen;
END_RCPP
}
// mask_insersectC
double mask_insersectC(arma::mat img1, arma::mat img2);
RcppExport SEXP _spat1f_mask_insersectC(SEXP img1SEXP, SEXP img2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type img1(img1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type img2(img2SEXP);
    rcpp_result_gen = Rcpp::wrap(mask_insersectC(img1, img2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spat1f_iouC", (DL_FUNC) &_spat1f_iouC, 2},
    {"_spat1f_mask_insersectC", (DL_FUNC) &_spat1f_mask_insersectC, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_spat1f(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
