# Add n random steps to the data frame, drawing from fitted sl and ta distributions. If sl_rand or ta_rand is supplied, then the vectors are randomly sampled for number number generation. Paired if set to TRUE uses the same random number index for sl_rand and ta_rand. Keep_obs if FALSE, drops the observed data. thresh_r sets r to zero if r is below thresh_r. turn angle theta is set to zero if r == 0
add_random_steps <- function(
    data,
    n = 20L,
    id = data$step_id,
    x_start = data$x1,
    y_start = data$y1,
    direction_start = data$theta_abs,
    sl_distr = fit_lnorm(data$r),
    ta_distr = fit_genvonmises(data$theta_rel),
    sl_rand = NULL,
    ta_rand = NULL,
    keep_obs = TRUE, 
    thresh_r = NULL,
    paired = FALSE){
  
  
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
    ta_rand <- rdist(ta_distr, n = 10e5)
  }
  
  if(!is.null(sl_rand)){
    sl_rand <- sl_rand
  } else {
    sl_rand <- rdist(sl_distr, n = 10e5)
  }
  
  if(!is.null(thresh_r)){
    sl_rand[sl_rand < thresh_r] <- 0
  }
  
  
  out <- lapply(seq_along(id), function(i){
    if(paired){
      index_theta <- index_r <- sample.int(n = length(sl_rand), size = n)
    } else {
      index_r <- sample.int(n = length(sl_rand), size = n)
      index_theta <- sample.int(n = length(ta_rand), size = n)
    }
    r <- sl_rand[index_r]
    theta <- ta_rand[index_theta]
    
    
    x1 <- x_start[i]
    y1 <- y_start[i]
    theta_abs <- direction_start[i]
    theta <- ifelse(r > 0, theta, 0)
    
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
    
    out <- cbind(out, data[rep(i, n), !names(data) %in% names(out), drop = FALSE])
    
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

# wrapper for survival::clogit
.issf_fit_internal <- function(formula, data, ...){
  clogit(
    formula,
    data = data,
    ...
  )
}

# Refit issf with newdata. If is null, use the same old data
refit_issf <- function(object, newdata = NULL){
  if(is.null(newdata)){
    newdata <- object$data
  }
  object$model <- .issf_fit_internal(formula = object$formula, data = newdata)
  object$data <- newdata
  return(object)
}

# Wrapper for clogit that update the distributions
issf <- function(formula, data, 
                 mulog_estimator = "logsl",
                 sdlog_estimator = "logslsq",
                 zero_inflation_estimator = "moved",
                 scale_estimator = "sl", 
                 shape_estimator = "logsl", 
                 kappa_estimator = "cos_thta" ,
                 kappa1_estimator = "cos_theta_pi",
                 kappa2_estimator = "cos_2theta",
                 remove_invalid = TRUE,
                 adjust = TRUE,
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
  if(adjust){
    out <- update_distr(out, 
                        mulog_estimator = mulog_estimator,
                        sdlog_estimator = sdlog_estimator,
                        zero_inflation_estimator = zero_inflation_estimator,
                        scale_estimator = scale_estimator, 
                        shape_estimator = shape_estimator,
                        kappa_estimator = kappa_estimator, 
                        kappa1_estimator = kappa1_estimator,
                        kappa2_estimator = kappa2_estimator,
                        remove_invalid = remove_invalid)
  }
  
  class(out) <- c("issf_fit","list")
  
  return(out)
}

# S3 method for print and summary of issf_fit objects
print.issf_fit <- function(x, ...){
  print(summary(x$model))
  invisible(summary(x$model))
}

# S3 method registration
registerS3method("print", "issf_fit", print.issf_fit)
registerS3method("summary", "issf_fit", print.issf_fit)


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
                          "zigamma" = append_zigamma_estimators)
  
  ta_est_method <- switch(attr(x, "ta")$name, 
                          "vonmises" = append_vonmises_estimators,
                          "genvonmises" = append_genvonmises_estimators)
  
  x %>% 
    sl_est_method() %>% 
    ta_est_method(na_as_zero = na_as_zero) %>% 
    append_moved()
  
}


# General update distribution wrapper using param estimators  
update_distr <- function(x, 
                         mulog_estimator = "logsl",
                         sdlog_estimator = "logslsq",
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
    x$ta_updated <- update_vonmises(x$ta, beta_cos_ta = x$model$coefficients[kappa_estimator])
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
      warning("Zigamma not implemented")
      # x$sl_updated <- update_zigamma_ss(x$sl, 
      #                                   beta_sl = x$model$coefficients[scale_estimator], 
      #                                   beta_log_sl = x$model$coefficients[shape_estimator])
      # x$data_updated <- resample_data(x, remove_invalid = remove_invalid)
      # x <- refit_issf(x, newdata = x$data_updated)
      # x$sl_updated <- update_zigamma_p(x$sl_updated,
      #                                x$model$coefficients[zero_inflation_estimator])
    } else {
      if(sl_dist_name == "lnorm"){
        x$sl_updated <- update_lnorm(x$sl, 
                                     beta_log_sl = x$model$coefficients[mulog_estimator], 
                                     beta_log_sl_sq = x$model$coefficients[sdlog_estimator])
      }
    }
  }
  
  return(x)
}

# Simulate available steps from updated sl and ta distributions
resample_ava_steps <- function(model, n = 1L){
  model$data %>%
    filter(case) %>% 
    dplyr::select(step_id, x1, x2, y1, y2, theta_abs, r, theta_rel) %>% 
    add_random_steps(n = n,
                     sl_distr = model$sl_updated[[1]],
                     ta_distr = model$ta_updated[[1]]
    ) %>%
    flag_invalid_steps(remove = TRUE)
}



# Internal function for resampling data in refitting issf in two shot estimation of zigamma. Currently defunct bc zigamma support is mostly removed.  
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


# Compute the utilization distribution area (level_ud * 100 % isopleth area). conf is the confidence intervals. Uses kernel density estimate method from ctmm. Fits a continuous time movement model (iid bivariate gaussian)
ud_area <- function(x, y, level_ud = 0.95,conf = 0.95){
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]
  y <- y[keep]
  
  d <- data.frame(
    lon = x, 
    lat = y, 
    timestamp = Sys.time() + seq_along(x), 
    ID = 1) %>% 
    ctmm::as.telemetry() %>% 
    suppressMessages()
  d$x <- d$longitude
  d$y <- d$latitude
  
  res <- ctmm::akde(d, ctmm::ctmm.fit(d)) %>% 
    ctmm::SpatialPolygonsDataFrame.UD(
      level.UD = level_ud, 
      conf.level = conf) %>% 
    sf::st_as_sf() %>% 
    sf::st_area() %>% 
    as.vector()
  
  names(res) <- c("lower", "estimate", "upper")
  res[c(2,1,3)]
  # Result in pixels
}

# Uses fitted issf_fit object to do simulation. If use_observed = TRUE, then the raw data used to fit the model is used instead for simulation. 
iterate_random_steps <- function(start, issf_fit, 
                                 n = 100, 
                                 index = 1,
                                 paired = FALSE, 
                                 use_observed = FALSE,
                                 r_thresh = attr(issf_fit$data, "r_thresh")){
  out.list <- vector(mode = "list", length = n+1)
  out.list[[1]] <- start
  stopifnot(nrow(start) == 1)
  
  if(use_observed){
    ra <- filter(issf_fit$data, case & !is.na(r) & !is.na(theta_rel))$theta_rel
    rr <- filter(issf_fit$data, case & !is.na(r) & !is.na(theta_rel))$r
    r_thresh <- NULL
  } else {
    ra <- rdist(issf_fit$ta_updated, 10^5)[[index]]
    rr <- rdist(issf_fit$sl_updated, 10^5)[[index]]
  }
  
  new_x2 <- start$x2
  new_y2 <- start$y2
  for(i in seq_len(n)){
    new <- TRUE
    while(new | !is_between(new_x2, c(0, 1000)) | 
          !is_between(new_y2, c(0, 1000))
    ){
      out.list[[i+1]] <- out.list[[i]] %>% 
        add_random_steps(
          data = .,
          n = 1, 
          keep_obs = FALSE, 
          id = start$step_id + i, 
          x_start = .$x2,
          y_start = .$y2,
          direction_start = (.$theta_abs + .$theta_rel) %% (2 * pi),
          sl_distr = NULL, 
          ta_distr = NULL, 
          ta_rand = ra,
          sl_rand = rr, 
          paired = paired, 
          thresh_r = NULL)
      new <- FALSE
      new_x2 <- out.list[[i+1]]$x2
      new_y2 <- out.list[[i+1]]$y2
    }
  }
  do.call("rbind.fill", out.list)
}

# Create starting object for iterate_random_steps()
creat_start <- function(x, y, theta){
  data.frame("step_id" = 1, "r" = NA, "r_threshed" = NA, "theta_abs" = theta, "theta_rel" = 0, "x1" = NA, "x2" = x, "y1" = NA, "y2" = y)
}






