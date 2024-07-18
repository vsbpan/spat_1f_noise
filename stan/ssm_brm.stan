data {
  int<lower=1> N;  // total number of observations
  int<lower=1> N_y;  // number of observations
  vector[N_y] Y_y;  // response variable
  int<lower=0> Nmi_y;  // number of missings
  int<lower=1> Jmi_y[Nmi_y];  // positions of missings
  int<lower=1> N_ylatent;  // number of observations
  vector[N_ylatent] Y_ylatent;  // response variable
  int<lower=0> Nmi_ylatent;  // number of missings
  int<lower=1> Jmi_ylatent[Nmi_ylatent];  // positions of missings
  // data needed for ARMA correlations 
  int<lower=0> Kar_ylatent;  // AR order
  int<lower=0> J_lag_ylatent[N_ylatent]; // number of lags per observation
  int<lower=0> N_predictors; // number of predictors
  matrix[N_y, N_predictors] X; // design matrix without the intercept
  int has_lag; // has lag
  real<lower=0> sigma_y_prior;
  real<lower=0> sigma_ylatent_prior;
}
transformed data {
  int max_lag_ylatent = Kar_ylatent; // Compute max lag
  int n_start = 1 + max_lag_ylatent;
}
parameters {
  vector[Nmi_y] Ymi_y;  // estimated missings y
  real<lower=0> sigma_y;  // dispersion parameter for observation process 
  vector[Nmi_ylatent] Ymi_ylatent;  // estimated missings y_latent (all missing)
  vector[Kar_ylatent] ar_ylatent;  // autoregressive coefficients 
  real<lower=0> sigma_ylatent;  // dispersion parameter for random walk process
  real Intercept_ylatent;
  vector[N_predictors] beta_x_ylatent; // beta coef for predictors
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += normal_lpdf(sigma_y | 0,sigma_y_prior) - 1 * normal_lccdf(0 | 0,sigma_y_prior); // half-normal prior for sigma
  lprior += normal_lpdf(ar_ylatent | 0, 0.1); // normal prior for lag coefs
  lprior += normal_lpdf(sigma_ylatent | 0,sigma_ylatent_prior) - 1 * normal_lccdf(0 | 0,sigma_ylatent_prior); // half-normal prior for sigma
  lprior += normal_lpdf(Intercept_ylatent | 0,1); // normal prior for intrinsic growth rate
  lprior += normal_lpdf(beta_x_ylatent | 0,1); // normal prior for marginal effect of x on intrinsic growth rate
}
model {
    // Initialize
    
    // vector combining observed and missing responses
    vector[N_y] Yl_y = Y_y;
    // vector combining observed and missing responses
    vector[N_ylatent] Yl_ylatent = Y_ylatent;
    // initialize linear predictor term
    vector[N_y] mu_y = rep_vector(0.0, N_y);
    
    
    // matrix storing lagged values
    matrix[N_ylatent, max_lag_ylatent] lag_ylatent = rep_matrix(0, N_ylatent, max_lag_ylatent);
    
    // initialize linear predictor term
    vector[N_ylatent] mu_ylatent = rep_vector(0.0, N_ylatent); // start from zero
    Yl_y[Jmi_y] = Ymi_y;
    Yl_ylatent[Jmi_ylatent] = Ymi_ylatent;
    
    // Loop through data
    // Latent y relation to mu_y
    for (n in 1:N_y) {
      // add more terms to the linear predictor
      mu_y[n] += Yl_ylatent[n];
    }
    
    
    // Loop through data, dropping the first J_lag_ylatent observation
    // Autoregressive struture for ylantent
    for (n in n_start:N_ylatent) {
      
      // If no lag, y_latent reduces to a pure white noise process, not random walk
      if(has_lag){
        // Compute ylatent for lag i
        for (i in 1:J_lag_ylatent[n]) {
          lag_ylatent[n, i] = Yl_ylatent[n - i];
        }
        
        // ar_ylatent can be modified as: 
        // ar_ylatent = X %*% beta
        // add lag effects to latent y
        mu_ylatent[n] += lag_ylatent[n, 1] + lag_ylatent[n, 1:Kar_ylatent] * ar_ylatent; 
      }
      
      // drift for latent y
      mu_ylatent[n] += Intercept_ylatent;
      
      if(N_predictors > 0){
        mu_ylatent[n] += X[n, 1:N_predictors] * beta_x_ylatent; 
      }
      
    }
    
    // Observation process
    target += normal_lpdf(Yl_y | mu_y, sigma_y);
    
    // Latent y random walk
    target += normal_lpdf(Yl_ylatent | mu_ylatent, sigma_ylatent); // mu_ylatent is set to zero
    
  // priors including constants
  target += lprior;
}



