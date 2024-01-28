# Add n random steps to the data frame, drawing from fitted sl and ta distirbutions
add_random_steps <- function(
    data,
    n = 20L,
    id = data$step_id,
    x_start = data$x1,
    y_start = data$y1,
    direction_start = data$theta_abs,
    sl_distr = fit_gamma(data$r),
    ta_distr = fit_genvonmises(data$theta_rel),
    sl_rand = NULL,
    ta_rand = NULL,
    keep_obs = TRUE){
  
  
  stopifnot(!any(duplicated(id)))
  nr <- length(id)
  stopifnot(
    nr == length(x_start),
    nr == length(y_start),
    nr == length(direction_start)
  )
  
  if(!is.null(ta_rand)){
    ta_rand <- ta_rand
  } else {
    ta_rand <- amt::random_numbers(ta_distr, n = 10e5)
  }
  
  if(!is.null(sl_rand)){
    sl_rand <- sl_rand
  } else {
    sl_rand <- amt::random_numbers(sl_distr, n = 10e5)
  }
  
  
  out <- lapply(seq_along(id), function(i){
    r <- sl_rand[sample.int(n = length(sl_rand), size = n)]
    theta <- ta_rand[sample.int(n = length(ta_rand), size = n)]
    
    x1 <- x_start[i]
    y1 <- y_start[i]
    theta_abs <- direction_start[i]
    
    out <- data.frame(
      "step_id" = id[i],
      "r" = r,
      "r_threshed" = r,
      "theta_abs" = theta_abs, 
      "theta_rel" = theta,
      "x1" = x1,
      "x2" = x1 + r * cos(theta + theta_abs),
      "y1" = y1,
      "y2" = y1 + r * sin(theta + theta_abs),
      "case" = FALSE
    )
    return(out)
  }) %>% 
    do.call("rbind", .)
  
  if(keep_obs){
    out <- plyr::rbind.fill(cbind(data, "case" = TRUE), out) %>%
      arrange(step_id)
  } else {
    out <- arrange(out, step_id)
  }
  
  attr(out, "r_thresh") <- attr(data, "r_thresh")
  attr(out, "sl") <- sl_distr
  attr(out, "ta") <- ta_distr
  attr(out, "call") <- match.call()
  return(out)
}


# Flag destinations that are out of bound
flag_invalid_steps <- function(data, x = data$x2, y = data$y2, 
                               dim_xy = c(1000, 1000), remove_invalid = FALSE){
  
  data$valid_step <- is_between(x, c(0,dim_xy[1])) & is_between(y, c(0,dim_xy[2]))
  if(remove_invalid){
    data <- data %>% dplyr::filter(valid_step)
  }
  return(data)
}

.issf_fit_internal <- function(formula, data, ...){
  clogit(
    formula,
    data = data,
    #method = "efron",
    ...
  )
}

refit_issf <- function(object, newdata = NULL){
  if(is.null(newdata)){
    newdata <- object$data
  }
  object$model <- .issf_fit_internal(formula = object$formula, data = newdata)
  object$data <- newdata
  return(object)
}

issf <- function(formula, data, 
                 zero_inflation_estimator = "moved",
                 scale_estimator = "sl", 
                 shape_estimator = "logsl", 
                 kappa_estimator = "cos_thta" ,
                 kappa1_estimator = "cos_theta_pi",
                 kappa2_estimator = "cos_2theta",
                 remove_invalid = FALSE,
                 ...){
  m <- .issf_fit_internal(formula, data, ...)
  out <- list(
    "model" = m,
    "sl" = attr(data, "sl"),
    "ta" = attr(data, "ta"),
    "data" = data,
    "call" = match.call(), 
    "formula" = formula
  )
  
  out <- update_distr(out, 
                      zero_inflation_estimator = zero_inflation_estimator,
                      scale_estimator = scale_estimator, 
                      shape_estimator = shape_estimator,
                      kappa_estimator = kappa_estimator, 
                      kappa1_estimator = kappa1_estimator,
                      kappa2_estimator = kappa2_estimator,
                      remove_invalid = remove_invalid)
  class(out) <- c("issf_fit","list")
  
  return(out)
}


print.issf_fit <- function(x, ...){
  print(summary(x$model))
}

registerS3method("print", "issf_fit", print.issf_fit)


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

append_gamma_estimators <- function(x){
  x %>% 
    mutate(
      sl = r,
      logsl = log(r)
    )
}

append_moved <- function(x, thresh = NULL){
  if(is.null(thresh)){
    thresh <- attr(x, "r_thresh")
  }
  x %>% 
    mutate(
      moved = ifelse(r > thresh, 1, 0)
    )
}

append_vonmises_estimators <- function(x, na_as_zero){
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


append_estimators <- function(x, na_as_zero = FALSE){
  sl_est_method <- switch(attr(x, "sl")$name, 
                          "gamma" = append_gamma_estimators)
  
  ta_est_method <- switch(attr(x, "ta")$name, 
                          "vonmises" = append_vonmises_estimators,
                          "genvonmises" = append_genvonmises_estimators)
  
  x %>% 
    sl_est_method() %>% 
    ta_est_method(na_as_zero = na_as_zero) %>% 
    append_moved()
  
}





update_distr <- function(x, 
                         zero_inflation_estimator = "moved",
                         scale_estimator = "sl", 
                         shape_estimator = "logsl",
                         kappa_estimator = "cos_theta",
                         kappa1_estimator = "cos_theta_pi", 
                         kappa2_estimator = "cos_2theta",
                         remove_invalid = FALSE){
  
  sl_dist_name <- x$sl$name 
  ta_dist_name <- x$ta$name
  
  
  
  if(ta_dist_name == "vonmises"){
    x$ta_updated <- amt::update_vonmises(
      x$ta, 
      beta_cos_ta = x$model$coefficients[kappa_estimator])
  } else {
    if(ta_dist_name == "genvonmises"){
      x$ta_updated <- update_genvonmises(
        x$ta, 
        beta_cos_theta_pi = x$model$coefficients[kappa1_estimator],
        beta_cos_2theta = x$model$coefficients[kappa2_estimator])
    }
  }
  
  
  if(sl_dist_name == "gamma"){
    x$sl_updated <- update_gamma(x$sl, 
                                beta_sl = x$model$coefficients[scale_estimator], 
                                beta_log_sl = x$model$coefficients[shape_estimator])
  } else {
    if(sl_dist_name == "zigamma"){
      x$sl_updated <- update_zigamma_ss(x$sl, 
                                        beta_sl = x$model$coefficients[scale_estimator], 
                                        beta_log_sl = x$model$coefficients[shape_estimator])
      x$data_updated <- resample_data(x, remove_invalid = remove_invalid)
      x <- refit_issf(x, newdata = x$data_updated)
      x$sl_updated <- update_zigamma_p(x$sl_updated,
                                     x$model$coefficients[zero_inflation_estimator])
    } else {
      # if(sl_dist_name == "cengamma"){
      #   x$sl_updated <- update_cengamma(x$sl, 
      #                                beta_sl = x$model$coefficients[scale_estimator], 
      #                                beta_log_sl = x$model$coefficients[shape_estimator])
      # }
    }
  }
  
  return(x)
}



resample_data <- function(object, n = NULL, remove_invalid = FALSE){
  if(is.null(n)){
    n <- as.list(attr(object$data, "call"))$n
  }
  # Kind of a hack right now. Need to clean up later to generalize
  object$data %>% 
    filter(case) %>% 
    add_random_steps(n = n, sl_distr = object$sl_updated, ta_distr = object$ta_updated) %>% 
    flag_invalid_steps(remove_invalid = remove_invalid) %>% 
    append_estimators()
}




# .test_p_adjust <- function(true_p, p){
#   test_d <- data.frame(
#     "step_id" = 1:3000, 
#     "r" = rzigamma(3000, p = true_p, shape = 1, scale = 10), 
#     "theta_rel" = random_numbers(make_vonmises_distr(5), 3000),
#     "p" = p
#   ) %>% 
#     add_random_steps(n = 30L, 
#                      id = step_id,
#                      x_start =  rep(NA, 3000),
#                      y_start = rep(NA, 3000),
#                      direction_start = rep(NA, 3000), 
#                      sl_distr = make_zigamma(unique(p), 1, 10), 
#                      ta_distr = make_vonmises_distr(2)) %>% 
#     mutate(
#       cos_theta = cos(theta_rel),
#       sl = r, 
#       logsl = ifelse(r == 0, 0, log(r)),
#       moved = ifelse(r == 0, 0, 1)
#     )
#   # issf(
#   #   case ~ 
#   #     moved+ 
#   #     (sl + logsl):moved + cos_theta : moved + 
#   #     strata(step_id),
#   #   data = test_d
#   # )$model
#     # coef() %>% 
#     # .[c(1)]
#   test_d
# }
# # When the tentative distributions match in all params except for p, use update zigamma() to get true updated distribution. 
# # Back transform directly from the betas very imprecise otherwise
# 
# 



