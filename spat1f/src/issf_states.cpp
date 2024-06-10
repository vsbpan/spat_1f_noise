#include <RcppArmadillo.h>
#include "shared_funs.h"
using namespace Rcpp;


DataFrame add_random_stepsC2(int n_draws,
                            double x_start, 
                            double y_start,
                            double direction_start,
                            List sl_rand, 
                            List ta_rand, 
                            int index,
                            double rss_coef_exp,
                            NumericVector ref_grid_flat,
                            double max_x = 1000, double max_y = 1000, 
                            double dim_x = 12, double dim_y = 12){
  
  List picked_list = pick_new_theta_xy(
    sl_rand,
    ta_rand,
    index,
    n_draws,
    direction_start,
    x_start,
    y_start,
    max_x,
    max_y,
    dim_x,
    dim_y
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

  DataFrame out = DataFrame::create(
    Named("theta") = thetai,
    Named("x") = xi,
    Named("y") = yi,
    Named("r") = ri,
    Named("on_toxic") = ref_grid_flat[ii-1]
  );
  return out;
}


// [[Rcpp::export]]
int pick_new_state(int state, NumericMatrix transition_mat){
  if(transition_mat.ncol() != transition_mat.nrow()){
    stop("Transition matrix is not a square matrix!");
  }
  
  NumericVector v = transition_mat(state - 1 , _ );
  IntegerVector state_new = sample_int(v.length(), 1, true, v);
  int res = (int)state_new[0];
  return res;
} 


// [[Rcpp::export]]
DataFrame add_random_steps_iterate_statesC(int n,
                                     int n_draws,
                                     double x_start, 
                                     double y_start,
                                     double direction_start,
                                     double diet_start,
                                     int state_start,
                                     List sl_rand, 
                                     List ta_rand, 
                                     double rss_coef,
                                     NumericMatrix transition_mat,
                                     bool same_move,
                                     NumericVector ref_grid_flat,
                                     double max_x = 1000, double max_y = 1000, 
                                     double dim_x = 12, double dim_y = 12){
  
  double rss_coef_exp = exp(rss_coef); 
  DataFrame tempdf = DataFrame::create(
    Named("theta") = direction_start,
    Named("x") = x_start,
    Named("y") = y_start,
    Named("r") = NA_REAL,
    Named("on_toxic") = diet_start,
    Named("state") = state_start
  ); 
  NumericVector theta = NumericVector(n);
  NumericVector x = NumericVector(n);
  NumericVector y = NumericVector(n);
  NumericVector r = NumericVector(n);
  NumericVector on_toxic = NumericVector(n);
  IntegerVector state = IntegerVector(n);
  double temp_on_toxic;
  int temp_state = state_start;
  int index;
  int slr_len = sl_rand.length();
  int tar_len = ta_rand.length();
  
  // Define nstates
  int nstates = transition_mat.nrow();
  
  for(int i = 0; i < n; ++i){
    temp_on_toxic = tempdf["on_toxic"];
    if(same_move){
      index = 1;
    } else { 
      if(temp_on_toxic > 0.5){
        index = 1; // more toxic 
      } else {  
        index = 2; // less toxic
      } 
      
      index = index + nstates * (temp_state - 1); // Repeat states 
      // expects: 
      // more toxic, state 1 
      // less toxic, state 1
      // more toxic, state 2 
      // less toxic, state 2
        
    }  
    
    
    if(slr_len < index){
      stop("Indexed location %u is not part of sl_rand list of length %u.", index, slr_len);
    }  
    
    if(tar_len < index){
      stop("Indexed location %u is not part of ta_rand list of length %u.", index, tar_len);
    }  
    
    
    tempdf = add_random_stepsC2(n_draws, tempdf["x"], tempdf["y"], tempdf["theta"],
                               sl_rand, ta_rand, 
                               index, 
                               rss_coef_exp, 
                               ref_grid_flat,
                               max_x, max_y, 
                               dim_x, dim_y);
    
    x[i] = tempdf["x"];
    y[i] = tempdf["y"];
    r[i] = tempdf["r"];
    theta[i] = tempdf["theta"];
    on_toxic[i] = tempdf["on_toxic"];
    
    state[i] = temp_state;
    temp_state = pick_new_state(temp_state, transition_mat);
    
  }  
  
  DataFrame out = DataFrame::create(
    Named("theta") = theta,
    Named("x") = x,
    Named("y") = y,
    Named("r") = r,
    Named("state") = state,
    Named("on_toxic") = on_toxic
  );
  
  return out;
}  



