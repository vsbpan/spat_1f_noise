#include <RcppArmadillo.h>
#include "shared_funs.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector read_valueC(NumericVector x, 
                          NumericVector y,
                          double max_x, double max_y,
                          NumericMatrix img){
  NumericVector dim_img = img.attr("dim");
  img.attr("dim") = R_NilValue;
  NumericVector img_flat = as<NumericVector>(img);
  NumericVector indices = flatten_xy(x, y, max_x, max_y, dim_img[0], dim_img[1]);
  NumericVector out = NumericVector(indices.length());
  bool is_na_here;
  for(int i = 0; i < out.length(); ++i){
    is_na_here = NumericVector::is_na(indices[i]);
    
    if(is_na_here){
      out[i] = NA_REAL;
    } else {
      out[i] = img_flat[indices[i]-1];
    }
  } 
  return out;
} 



// [[Rcpp::export]]
DataFrame add_random_stepsC(int n_draws,
                            double x_start, 
                            double y_start,
                            double direction_start,
                            List sl_rand, 
                            List ta_rand, 
                            int index,
                            double rss_coef_exp,
                            NumericVector ref_grid_flat){

  List picked_list = pick_new_theta_xy(
    sl_rand,
    ta_rand,
    index,
    n_draws,
    direction_start,
    x_start,
    y_start
  );
  
  NumericVector i = picked_list["i_new"];

  NumericVector val_new = ifelse(ref_grid_flat > 0.5, 1.0, 1.0 * rss_coef_exp);
  val_new = val_new[i - 1];
  IntegerVector ind = sample_int(val_new.length(), 1, true, val_new) - 1;
  
  NumericVector theta = picked_list["theta_new"];
  NumericVector thetai = theta[ind];
  NumericVector x = picked_list["x_new"];
  NumericVector xi = x[ind];
  NumericVector y = picked_list["y_new"];
  NumericVector yi = y[ind];
  NumericVector r = picked_list["r_new"];
  NumericVector ri = r[ind];
  
  NumericVector ii = i[ind];
  NumericVector val = ref_grid_flat[i - 1];
  double ava_toxic = mean(val);
  
  
  DataFrame out = DataFrame::create(
    Named("theta") = thetai,
    Named("x") = xi,
    Named("y") = yi,
    Named("r") = ri,
    Named("on_toxic") = ref_grid_flat[ii-1],
    Named("ava_toxic") = ava_toxic
  );

  return out;
}

// [[Rcpp::export]]
DataFrame add_random_steps_iterateC(int n,
                                    int n_draws,
                                    double x_start, 
                                    double y_start,
                                    double direction_start,
                                    double diet_start,
                                    List sl_rand, 
                                    List ta_rand, 
                                    double rss_coef,
                                    bool same_move,
                                    NumericVector ref_grid_flat){
  
  double rss_coef_exp = exp(rss_coef); 
  DataFrame tempdf = DataFrame::create(
    Named("theta") = direction_start,
    Named("x") = x_start,
    Named("y") = y_start,
    Named("r") = NA_REAL,
    Named("on_toxic") = diet_start,
    Named("ava_toxic") = NA_REAL
  );
  NumericVector theta = NumericVector(n);
  NumericVector x = NumericVector(n);
  NumericVector y = NumericVector(n);
  NumericVector r = NumericVector(n);
  NumericVector on_toxic = NumericVector(n);
  NumericVector ava_toxic = NumericVector(n);
  double temp_on_toxic;
  int index;
  int slr_len = sl_rand.length();
  int tar_len = ta_rand.length();
  
  
  for(int i = 0; i < n; ++i){
    temp_on_toxic = tempdf["on_toxic"];
    if(same_move){
      index = 1;
    } else {
      if(temp_on_toxic > 0.5){
        index = 1;
      } else {
        index = 2;
      }
    }
    
    
    if(slr_len < index){
      stop("Indexed location %u is not part  of sl_rand list of length %u.", index, slr_len);
    }
    
    if(tar_len < index){
      stop("Indexed location %u is not part of ta_rand list of length %u.", index, tar_len);
    }
    
    
    tempdf = add_random_stepsC(n_draws, tempdf["x"], tempdf["y"], tempdf["theta"],
                               sl_rand, ta_rand, 
                               index, 
                               rss_coef_exp, ref_grid_flat);
    x[i] = tempdf["x"];
    y[i] = tempdf["y"];
    r[i] = tempdf["r"];
    theta[i] = tempdf["theta"];
    on_toxic[i] = tempdf["on_toxic"];
    ava_toxic[i] = tempdf["ava_toxic"];
  }
  
  DataFrame out = DataFrame::create(
    Named("theta") = theta,
    Named("x") = x,
    Named("y") = y,
    Named("r") = r,
    Named("on_toxic") = on_toxic,
    Named("ava_toxic") = ava_toxic
  );
  
  return out;
}


