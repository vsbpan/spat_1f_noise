# Append generalized von mises distribution param estimators. na_as_zero should generally be set to TRUE as that they don't get dropped when used in tandem with :moved
append_genvonmises_estimators <- function(x, na_as_zero = FALSE){
  if(na_as_zero){
    x %>% 
      mutate(
        cos_theta_pi = ifelse(is.na(theta_rel), 0, cos(theta_rel + pi)),
        cos_2theta = ifelse(is.na(theta_rel), 0, cos(2 * (theta_rel))),
      )
  } else {
    x %>% 
      mutate(
        cos_theta_pi = cos(theta_rel + pi),
        cos_2theta = cos(2 * (theta_rel)),
      )
  }
}


# Append zigamma distribution param estimators
append_zigamma_estimators <- function(x){
  x %>% 
    mutate(
      sl = r_threshed,
      logsl = ifelse(r_threshed > 0, r_threshed, 0)
    )
}

# Append gamma distribution param estimators
append_gamma_estimators <- function(x){
  x %>% 
    mutate(
      sl = r,
      logsl = log(r)
    )
}

# Append lognormal distribution param estimators
append_lnorm_estimators <- function(x){
  x %>% 
    mutate(
      logsl = log(r),
      logslsq = log(r)^2
    )
}

# Append invgamma distribution param estimators
append_invgamma_estimators <- function(x){
  x %>% 
    mutate(
      logsl = log(r),
      invsl = 1/r
    )
}


# Append zigamma zero inflation estimator
append_moved <- function(x, thresh = NULL){
  if(is.null(thresh)){
    thresh <- attr(x, "r_thresh")
  }
  x %>% 
    mutate(
      moved = ifelse(r > thresh, 1, 0)
    )
}

# Append von mises param estimator. na_as_zero should generally be set to TRUE
append_vonmises_estimators <- function(x, na_as_zero = FALSE){
  if(na_as_zero){
    x %>% 
      mutate(
        cos_theta = ifelse(is.na(threa_rel), 0, cos(theta_rel))
      )
  } else {
    x %>% 
      mutate(
        cos_theta = cos(theta_rel)
      )
  }
}

# General append estimator wrapper for various distributions
append_estimators <- function(x, na_as_zero = FALSE){
  
  sl_est_method <- switch(attr(x, "sl")$name, 
                          "gamma" = append_gamma_estimators,
                          "lnorm" = append_lnorm_estimators,
                          "invgamma" = append_invgamma_estimators,
                          "zigamma" = append_zigamma_estimators)
  
  ta_est_method <- switch(attr(x, "ta")$name, 
                          "vonmises" = append_vonmises_estimators,
                          "genvonmises" = append_genvonmises_estimators)
  
  x %>% 
    sl_est_method() %>% 
    ta_est_method(na_as_zero = na_as_zero) %>% 
    append_moved()
  
}

interaction_grid <- function(estimators, int_list){
  init_flat <- apply(expand.grid(int_list), 1, function(x){
    paste0(x, collapse = ":")
  })
  
  out <- lapply(estimators, function(x){
    sprintf("%s:%s",
            init_flat, 
            x)
  })
  out
}


# default distribution estimator. called in pick_default_estimators()
.gamma_default_estimator <- function(){
  out <- list(
    "shape" = "logsl",
    "scale" = "sl"
  )
}

.lnorm_default_estimator <- function(){
  list(
    "mulog" = "logsl",
    "sdlog" = "logslsq"
  )
}

.invgamma_default_estimator <- function(){
  list(
    "shape" = "logsl",
    "scale" = "invsl"
  )
}

.zigamma_default_estimator <- function(){
  list(
    "shape" = "logsl",
    "scale" = "sl",
    "p" = "moved"
  )
}

.vonmises_default_estimator <- function(){
  list(
    "kappa" = "cos_theta"
  )
}

.genvonmises_default_estimator <- function(){
  list(
    "kappa1" = "cos_theta_pi",
    "kappa2" = "cos_2theta"
  )
}


pick_default_estimators <- function(dist, int_list = NULL){
  fcall <- sprintf(".%s_default_estimator", dist)
  if(isFALSE(is.function(get0(fcall)))){
    stop(
      sprintf("Could not find function '%s'. Please define in the global environment.", fcall)
    )
  }
  
  res <- do.call(fcall, list())
  
  if(!is.null(int_list)){
    res <- interaction_grid(res, int_list)
  } 
  
  return(res)
}



# General update distribution wrapper using param estimators  
update_distr <- function(x, 
                         sl_estimators = NULL,
                         ta_estimators = NULL){
  
  sl_dist_name <- x$sl$name 
  ta_dist_name <- x$ta$name
  
  if(is.null(sl_estimators)){
    sl_estimators <- pick_default_estimators(sl_dist_name)
  }
  
  if(is.null(ta_estimators)){
    ta_estimators <- pick_default_estimators(ta_dist_name)
  }
  
  map_estimators <- function(estimators, coefs){
    names(estimators) <- paste0(names(estimators), "_estimator")
    estimators <- lapply(estimators, function(estimator_name){
      coefs[estimator_name]
    })
    return(estimators)
  }
  
  x$sl_updated <- do.call(paste0("update_",sl_dist_name), 
                          c(
                            list(x$sl), 
                            map_estimators(sl_estimators, x$model$coefficients)
                          ))
  
  x$ta_updated <- do.call(paste0("update_",ta_dist_name), 
                          c(
                            list(x$ta), 
                            map_estimators(ta_estimators, x$model$coefficients)
                          ))
  
  return(x)
}
