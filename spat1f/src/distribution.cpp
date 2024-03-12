#include <RcppArmadillo.h>
using namespace Rcpp;

// Function g(omega) as in eqn 24 of Gatto 2008 Stat Comput
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector g_genvonmisesC(NumericVector omega, double kappa1, double kappa2) {
  NumericVector out = kappa1 * cos(omega) + kappa2 * cos(2 * (omega - 3.141593)); 
  return out;
}


// ratio of uniform proposal generator for generalize vonmises distribution random number generation
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector propose_genvonmises(int n, double a, double kappa1, double kappa2){
  double pi = 3.141593;
  NumericVector u = runif(n, 0, a);
  NumericVector v = runif(n, 0, 2 * pi * a);
  LogicalVector s = v > 2 * pi * u;
  u = ifelse(s, a - u, u);
  v = ifelse(s, 2 * pi * a - v, v);
  NumericVector w = g_genvonmisesC(v/u, kappa1, kappa2) / 2;
  LogicalVector cond = u <= exp(w);
  NumericVector out = v[cond]/u[cond];
  return out;
}

// runif as in R
NumericVector runif(int n, double l, double u){
  arma::vec v = arma::randu(n, arma::distr_param(l, u));
  NumericVector out = NumericVector(v.begin(), v.end());
  return out;
}



