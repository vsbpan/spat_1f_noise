// #include <RcppArmadillo.h>
// #include "shared_funs.h"
// using namespace Rcpp;
// 
// // [[Rcpp::export]]
// NumericVector decay_wt(NumericVector t, double k){
//   return 1 / (1 + t * k);
// }
// 
// 
// 
// // Function to handle NA indices in x and return corresponding values
// // [[Rcpp::export]]
// NumericVector index_with_NA(NumericVector x, IntegerVector i) {
//   int n = i.length();
//   NumericVector out(n);
// 
//   for (int j = 0; j < n; ++j) {
//     int idx = i[j];
//     out[j] = (idx == NA_INTEGER) ? NA_REAL : x[idx];
//   }
// 
//   return out;
// } 
// 
// // [[Rcpp::export]]
// double delta_funC(NumericVector x, NumericVector t, int i, double k){
//   i = i - 1;
//   NumericVector tau = t[Rcpp::seq_len(i + 1) - 1];
//   tau = tau - 1;
//   NumericVector w = decay_wt(t[i] - tau, k);
//   
//   NumericVector z = index_with_NA(x, Rcpp::match(tau, t) - 1);
//   NumericVector wt_xi = z * w;
// 
//   NumericVector a = ((wt_xi[!Rcpp::is_na(wt_xi)]));
//   NumericVector b = (w[!Rcpp::is_na(wt_xi)]);
//   double delta = x[i] - (Rcpp::sum(a) / Rcpp::sum(b));
//   return delta;
// }
// 
// // [[Rcpp::export]]
// IntegerVector order_(NumericVector x) {
//   NumericVector sorted = clone(x).sort();
//   return Rcpp::match(sorted, x);
// }
// 
// 
// // [[Rcpp::export]]
// NumericVector redupC(NumericVector x, NumericVector grid, NumericVector target){
//   return x[Rcpp::match(target, grid) - 1];
// }
// 
// 
// 
// // [[Rcpp::export]]
// List memdelC(NumericVector x, NumericVector t, LogicalVector Case, double k){
//   if(x.length() != t.length()){
//     stop("Different length vectors detected. x is of length %u and t is of length %u.", x.length(), t.length());
//   }
//   
//   NumericVector x2, t2;
//   
//   x = x[order_(t) - 1];
//   x2 = x[Case];
//   t2 = t[Case];
// 
//   int n = t2.length();
//   NumericVector res = NumericVector(n, NA_REAL);
//   
// 
//   for(int i = 0; i < n; ++i){
//     res[i] = delta_funC(x2, t2, i+1, k);
//   }
//   
//   // res = redup(res, t2, t);
//   return List::create(Named("res") = res,
//                       Named("t2") = t2,
//                       Named("t") = t);
// 
// } 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
