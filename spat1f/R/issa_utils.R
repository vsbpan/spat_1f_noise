# Add n random steps to the data frame, drawing from fitted sl and ta distributions. If sl_rand or ta_rand is supplied, then the vectors are randomly sampled for number number generation. Paired if set to TRUE uses the same random number index for sl_rand and ta_rand. Keep_obs if FALSE, drops the observed data. thresh_r sets r to zero if r is below thresh_r. turn angle theta is set to zero if r == 0
add_random_steps <- function(
    data,
    n = 20L,
    id = data$step_id,
    x_start = data$x1,
    y_start = data$y1,
    direction_start = data$theta_abs,
    state = data$state,
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
    ta_rand <- rdist(ta_distr, n = 1e5)
  }
  
  if(!is.null(sl_rand)){
    sl_rand <- sl_rand
  } else {
    sl_rand <- rdist(sl_distr, n = 1e5)
  }
  
  if(!is.null(thresh_r)){
    sl_rand[sl_rand < thresh_r] <- 0
  }
  
  if(!is.list(sl_rand)){
    sl_rand <- list(sl_rand)
  }
  
  if(!is.list(ta_rand)){
    ta_rand <- list(ta_rand)
  }
  
  stopifnot(length(sl_rand) == length(ta_rand))
  
  if(is.null(state) || length(sl_rand) == 1){
    index <- rep(1, nr)
  } else {
    index <- state
  }
  
  out <- lapply(seq_along(id), function(i){
    
    
    if(paired){
      index_theta <- index_r <- sample.int(n = length(sl_rand[[index[i]]]), size = n)
    } else {
      index_r <- sample.int(n = length(sl_rand[[index[i]]]), size = n)
      index_theta <- sample.int(n = length(ta_rand[[index[i]]]), size = n)
    }
    r <- sl_rand[[index[i]]][index_r]
    theta <- ta_rand[[index[i]]][index_theta]
    
    
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
                 sl_estimators = NULL,
                 ta_estimators = NULL,
                 update = TRUE,
                 keep_data = TRUE,
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
  if(update){
    out <- update_distr(out, 
                        sl_estimators = sl_estimators,
                        ta_estimators = ta_estimators)
  }
  
  if(!keep_data){
    out$data <- NULL
  }
  
  class(out) <- c("issf_fit","list")
  
  return(out)
}

# S3 method for coef of issf_fit objects
coef.issf_fit <- function(x, se = FALSE, ...){
  if(se){
    return(
      list(
        "estimate" = coef(x$model),
        "se" = sqrt(diag(x$model$var))
      )
    )
  } else {
    return(
      coef(x$model)
    )
  }
}

# S3 method for print and summary of issf_fit objects
print.issf_fit <- function(x, ...){
  print(summary(x$model))
  invisible(summary(x$model))
}

# S3 method registration
registerS3method("print", "issf_fit", print.issf_fit)
registerS3method("summary", "issf_fit", print.issf_fit)
registerS3method("coef", "issf_fit", coef.issf_fit)



# Simulate available steps from updated sl and ta distributions
resample_ava_steps <- function(model, n = 1L){
  model$data %>%
    filter(case) %>% 
    dplyr::select(step_id, x1, x2, y1, y2, theta_abs, r, theta_rel, state) %>% 
    add_random_steps(n = n,
                     sl_distr = model$sl_updated,
                     ta_distr = model$ta_updated
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






