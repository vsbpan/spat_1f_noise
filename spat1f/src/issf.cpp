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




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List pick_new_theta_xy(List sl_rand, List ta_rand, int index, int n, 
                       double direction_start, double x_start, double y_start){
  double pi = 3.141593;

  NumericVector theta_new;
  NumericVector x_new;
  NumericVector y_new;
  NumericVector r_new;
  
  if((sl_rand.length() < index) | (index < 1) | (ta_rand.length() < index)){
    stop("Indexed location %u does not exist in list of length %u and %u.", index, sl_rand.length(), ta_rand.length());
  }

  for(int i = 0; i < 1001; ++i){
    NumericVector r = sample_dbl(sl_rand[index-1], n, true);
    NumericVector theta = sample_dbl(ta_rand[index-1], n, true);

    theta_new = mod(direction_start + theta, 2 * pi);
    x_new = x_start + r * cos(theta_new);
    y_new = y_start + r * sin(theta_new);

    LogicalVector out_of_bound = ((x_new < 0) | (x_new > 1000)) | ((y_new < 0) | (y_new > 1000));
    x_new = x_new[!out_of_bound];
    y_new = y_new[!out_of_bound];
    theta_new = theta_new[!out_of_bound];
    r_new = r[!out_of_bound];
    if(y_new.length() > 0){ // Exit condition
      break;
    } else {
      if(i == 1000){
        stop("Timed out after 1000 tries. Try increasing n.");
      }
    }
  }

  NumericVector i_new = flatten_xy(x_new, y_new);
  List out = List::create(Named("x_new") = x_new, 
                          Named("y_new") = y_new, 
                          Named("theta_new") = theta_new,
                          Named("r_new") = r_new,
                          Named("i_new") = i_new);
  
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
  NumericVector temp_on_toxic = NumericVector(n);
  IntegerVector ind = IntegerVector(1);
  int index;
  
  for(int i = 0; i < n; ++i){
    temp_on_toxic = tempdf["on_toxic"];
    if(same_move){
      ind = 1;
    } else {
      ind = ifelse(temp_on_toxic > 0.5, 1, 2);
    }
    
    index = ind[0];
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


