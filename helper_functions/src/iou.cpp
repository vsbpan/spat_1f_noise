#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double iouC(arma::mat img1, arma::mat img2) {
  double i = arma::accu(img1 && img2);
  double u = arma::accu(img1 || img2);
  double iou = i / u;
  return iou;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double mask_insersectC(arma::mat img1, arma::mat img2) {
  double i = arma::accu(img1 && img2);
  return i;
}