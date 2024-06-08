#ifndef SHARED_FUNS_H
#define SHARED_FUNS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

NumericVector sample_dbl(NumericVector x, int sz, bool rep = false, 
                         sugar::probs_t prob = R_NilValue);

NumericVector mod(NumericVector x, int n);

NumericVector seq_len_num(int n);

IntegerVector sample_int(int n, int sz, bool rep = false, 
                         sugar::probs_t prob = R_NilValue);

NumericVector flatten_xy(NumericVector x, NumericVector y,
                         double max_x = 1000, double max_y = 1000, 
                         double dim_x = 12, double dim_y = 12);

List pick_new_theta_xy(List sl_rand, List ta_rand, int index, int n, 
                       double direction_start, double x_start, double y_start);


#endif // SHARED_FUNCTIONS_H
