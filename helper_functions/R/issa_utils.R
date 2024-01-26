# Add n random steps to the data frame, drawing from fitted sl and ta distirbutions
add_random_steps <- function(
    data,
    n = 20L,
    id = step_id,
    x_start = x1,
    y_start = y1,
    direction_start = theta_abs,
    sl_distr = amt::fit_distr(r_threshed, "gamma"),
    ta_distr = fit_genvonmises(theta_rel)){
  
  .expose_columns_interal()
  
  
  stopifnot(!any(duplicated(id)))
  nr <- length(id)
  stopifnot(
    nr == length(x_start),
    nr == length(y_start),
    nr == length(direction_start)
  )
  
  r_rand <- amt::random_numbers(sl_distr, n = 10e5)
  theta_rand <- amt::random_numbers(ta_distr, n = 10e5)
  
  out <- lapply(seq_along(id), function(i){
    r <- r_rand[sample.int(n = 10e5, size = n)]
    theta <- theta_rand[sample.int(n = 10e5, size = n)]
    
    x1 <- x_start[i]
    y1 <- y_start[i]
    theta_abs <- direction_start[i]
    
    out <- data.frame(
      "step_id" = id[i],
      "r" = r,
      "r_threshed" = r,
      "theta_abs" = theta_abs, 
      "theta_rel" = theta,
      "cos_theta" = cos(theta + theta_abs),
      "x1" = x1,
      "x2" = x1 + r * cos(theta + theta_abs),
      "y1" = y1,
      "y2" = y1 + r * sin(theta + theta_abs),
      "case" = FALSE
    )
    return(out)
  }) %>% 
    do.call("rbind", .)
  
  
  out <- plyr::rbind.fill(cbind(data, "case" = TRUE), out) %>%
    arrange(step_id)
  
  attr(out, "sl") <- sl_distr
  attr(out, "ta") <- ta_distr
  return(out)
}


# Flag destinations that are out of bound
flag_invalid_steps <- function(data, x = x2, y = y2, 
                               dim_xy = c(1000, 1000), remove = FALSE){
  .expose_columns_interal()
  
  data$valid_step <- is_between(x, c(0,dim_xy[1])) & is_between(y, c(0,dim_xy[2]))
  if(remove){
    data <- data %>% dplyr::filter(valid_step)
  }
  return(data)
}


issf <- function(formula, data, 
                 scale_estimator = "sl", 
                 shape_estimator = "logsl", 
                 kappa1_estimator = "cos_theta_pi",
                 kappa2_estimator = "cos_2theta",
                 ...){
  m <- clogit(
    formula,
    data = data,
    ...
  )
  out <- list(
    "model" = m,
    "sl" = attr(data, "sl"),
    "ta" = attr(data, "ta"),
    "data" = data
  )
  
  out <- update_distr(out, 
                      scale_estimator = scale_estimator, 
                      shape_estimator = shape_estimator,
                      kappa1_estimator = kappa1_estimator,
                      kappa2_estimator = kappa2_estimator)
  class(out) <- c("issf_fit","list")
  
  return(out)
}


print.issf_fit <- function(x, ...){
  print(summary(x$model))
}

registerS3method("print", "issf_fit", print.issf_fit)


append_genvonmises_estimators <- function(x){
  x %>% 
    mutate(
      cos_theta_pi = cos(theta_rel + pi),
      cos_2theta = cos(2 * (theta_rel)),
    )
}

update_distr <- function(x, scale_estimator = "sl", shape_estimator = "logsl", 
                         kappa1_estimator = "cos_theta_pi", 
                         kappa2_estimator = "cos_2theta"){
  x$sl_updated <- amt::update_gamma(x$sl, 
                                    beta_sl = x$model$coefficients[scale_estimator], 
                                    beta_log_sl = x$model$coefficients[shape_estimator])
  x$ta_updated <- update_genvonmises(x$ta, 
                                     beta_cos_theta_pi = x$model$coefficients[kappa1_estimator],
                                     beta_cos_2theta = x$model$coefficients[kappa2_estimator])
  return(x)
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



