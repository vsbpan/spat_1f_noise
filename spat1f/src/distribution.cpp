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

NumericVector conc(NumericVector x, NumericVector y){
  int nx=x.size(), n=x.size()+y.size(),i,j;
  NumericVector out=no_init(n);
  for (i=0; i<nx; ++i){ out[ i ] = x[ i ];}
  for (j=i, i=0; j<n; ++j, ++i){ out[ j ] = y[i] ;}
  return out;
}


// [[Rcpp::export]]
NumericVector rgenvonmisesC(int n, double kappa1, double kappa2, int max_try = 1000){
  double pi = 3.141593;
  IntegerVector grid = seq(0, floor(2 * pi * 100000));
  
  double a = exp(max(g_genvonmisesC(as<NumericVector>(grid) / 100000, kappa1, kappa2)));
  
  int n_out = 0;
  NumericVector x = NumericVector(0);
  NumericVector x_temp;
  
  for(int i = 1; i < (max_try + 1); ++i){
    x_temp = propose_genvonmises(n, a, kappa1, kappa2);
    x = conc(x, x_temp);
    n_out = x.length();
    if(n_out > n){
      break;
    } else {
      if(i == max_try){
        double p = (double)n_out / (double)n * 100;
        stop(
          "Function timed out after %u tries. Generated %.02f percent of requested numbers. \nConsider increasing the 'max_try' argument.", 
          max_try, p
        );
      }
    }
  }

  NumericVector out = x[seq(1, n)] - pi;

  return out;
}
