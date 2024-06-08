#include <RcppArmadillo.h>
using namespace Rcpp;

NumericVector sample_dbl(NumericVector x, int sz, bool rep = false, 
                         sugar::probs_t prob = R_NilValue){
  return sample(x, sz, rep, prob);
}


NumericVector mod(NumericVector x, int n){
  return x - floor(x/n)*n;
}

NumericVector seq_len_num(int n){
  IntegerVector vec = seq(1, n);
  NumericVector vec_num = as<NumericVector>(vec);
  return vec_num;
}

IntegerVector sample_int(int n, int sz, bool rep = false, 
                         sugar::probs_t prob = R_NilValue){ 
  NumericVector res_num = sample_dbl(seq_len_num(n), sz, rep, prob);
  IntegerVector res = as<IntegerVector>(res_num);
  return res;
}


// [[Rcpp::export]]
NumericVector flatten_xy(NumericVector x, NumericVector y,
                         double max_x = 1000, double max_y = 1000, 
                         double dim_x = 12, double dim_y = 12){
  if(x.length() != y.length()){
    stop("x and y must have the same length");
  }
  NumericVector i = (ceiling(x / max_x * dim_x) + (ceiling(y / max_y * dim_y) - 1) * dim_x);
  return i;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List pick_new_theta_xy(List sl_rand, List ta_rand, int index, int n, 
                       double direction_start, double x_start, double y_start){

  NumericVector theta_new, x_new, y_new, r_new;
  
  bool loop_success = false;
  
  for(int i = 0; i < 1001; ++i){
    NumericVector r = sample_dbl(sl_rand[index-1], n, true);
    NumericVector theta = sample_dbl(ta_rand[index-1], n, true);
    
    theta_new = mod(direction_start + theta, 6.283186);
    x_new = x_start + r * cos(theta_new);
    y_new = y_start + r * sin(theta_new);

    LogicalVector out_of_bound = ((x_new < 0) | (x_new > 1000)) | ((y_new < 0) | (y_new > 1000));
    x_new = x_new[!out_of_bound];
    y_new = y_new[!out_of_bound];
    theta_new = theta_new[!out_of_bound];
    r_new = r[!out_of_bound];
    
    if(y_new.length() > 0){ // Exit condition
      loop_success = true;
      break;
    }
  }
  
  if(!loop_success){
    stop("Timed out after 1000 tries. Try increasing n.");
  }
  

  NumericVector i_new = flatten_xy(x_new, y_new);
  List out = List::create(Named("x_new") = x_new, 
                          Named("y_new") = y_new, 
                          Named("theta_new") = theta_new,
                          Named("r_new") = r_new,
                          Named("i_new") = i_new);
  
  return out;
}


