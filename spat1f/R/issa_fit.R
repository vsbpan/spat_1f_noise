mledist <- function(x, dist, init, amt_format = FALSE, param_names = NULL, ...){
  if(!is.null(param_names)){
    names(init) <- param_names
    if(!all(param_names %in% names(formals(sprintf("d%s", dist))))){
      extra <- param_names[!param_names %in% names(formals(sprintf("d%s", dist)))]
      stop(sprintf("Cannot find the argument '%s' in the function d%s", extra, dist))
      
    }
  }
  
  
  
  fit <- fitdistrplus::mledist(x, dist, start = as.list(init), 
                               fix.arg = NULL, checkstartfix = TRUE, 
                               ...)
  names(fit)[1] <- "par"
  
  if(is.na(fit$value) || !is.finite(is.na(fit$value))){
    fit$convergence <- 1
  } else {
    fit$se <- tryCatch(sqrt(diag(solve(fit$hessian))), 
                       error = function(e){
                         e
                       })
    if(inherits(fit$se, "error") || any(is.na(fit$se))){
      if(inherits(fit$se, "error")){
        fit$message <- paste0(as.character(fit$se), collapse = "\n")
      } else {
        fit$message <- "Cannot estimate standard error."
      }
      
      fit$convergence <- 2
      fit$se <- rep(NA_real_, length(fit$par))
    } else {
      fit$message <- ""
    }
  }
  
  if(amt_format){
    if(is.null(param_names)){
      stop("Must supply 'param_names' if amt_format = TRUE.")
    }
    
    if(fit$convergence != 0){
      message(sprintf("Did not converge for %s. \nCode: %s\nMessage: %s\n",
                      dist,
                      fit$convergence,
                      fit$message))
      vcov <- matrix(rep(NA, length(param_names)^2),
                     nrow = length(param_names),
                     ncol = length(param_names))
    } else {
      vcov <- solve(fit$hessian) 
    }
    
    params <- as.list(fit$par)
    
    names(params) <- param_names
    rownames(vcov) <- colnames(vcov) <- param_names
    
    out <- list("name" = dist,
                "params" = params,
                "vcov" = vcov) # VCV matrix on transformed scale
    
    return(out)
  } else {
    return(fit)
  }
}



# Nice wrapper for fitting generalized von Mises distribution using MLE. 
fit_genvonmises <- function(x, 
                            init = "auto",
                            na.rm = TRUE, 
                            ...){
  
  if(na.rm){
    x <- x[!is.na(x)]
  }
  
  if(any(init == "auto")){
    init[1] <- 1
    init[2] <- 0.5
    init <- as.numeric(init)
  }
  
  param_names <- c("kappa1", "kappa2")
  dist <- "genvonmises"
  out <- mledist(x, dist, init = init, amt_format = TRUE, param_names = param_names, 
                 lower = c(0, 0), optim.method = "L-BFGS-B", ...)
  
  class(out) <- c(sprintf("%s_distr", dist), "ta_distr", "amt_distr", "list")
  return(out)
  
}


# Fit zigamma returning the same format as `amt::fit_distr()`
fit_zigamma <- function(x,
                        init = "auto",
                        na.rm = TRUE, 
                        ...){
  
  if(na.rm){
    x <- x[!is.na(x)]
  }
  
  if(any(init == "auto")){
    m1 <- mean(x[x>0])
    m2 <- mean((x[x>0])^2)
    init[1] <- mean(c(x == 0, TRUE))
    init[2] <- (2 * m2 - m1^2)/(m2 - m1^2)
    init[3] <- m1 * m2/(m2 - m1^2)
    init <- as.numeric(init)
  }
  
  param_names <- c("p", "shape", "scale")
  
  dist <- "zigamma"
  out <- mledist(x, dist, init = init, amt_format = TRUE, param_names = param_names, ...)
  
  class(out) <- c(sprintf("%s_distr", dist), "sl_distr", "amt_distr", "list")
  return(out)
}


# Fit inverse gamma distribution with MLE
fit_invgamma <- function(x, 
                         init = "auto",
                         na.rm = TRUE, 
                         ...){
  
  if(na.rm){
    x <- x[!is.na(x)]
  }
  
  if(any(init == "auto")){
    m1 <- mean(x)
    m2 <- mean(x^2)
    init[1] <- (2 * m2 - m1^2)/(m2 - m1^2)
    init[2] <- m1 * m2/(m2 - m1^2)
    init <- as.numeric(init)
  }
  
  param_names <- c("shape", "scale")
  
  dist <- "invgamma"
  out <- mledist(x, dist, init = init, amt_format = TRUE, param_names = param_names, ...)
  
  class(out) <- c(sprintf("%s_distr", dist), "sl_distr", "amt_distr", "list")
  return(out)
}






# Standardized wrapper for fitting gamma distribution
fit_gamma <- function(x, 
                      init = "auto",
                      na.rm = TRUE, 
                      ...){
  if(na.rm){
    x <- x[!is.na(x)]
  }
  
  if(any(init == "auto")){
    m <- mean(x)
    v <- var(x)
    init[1] <- m^2/v
    init[2] <- m/v
    init <- as.numeric(init)
  }
  
  param_names <- c("shape", "scale")

    dist <- "gamma"
  out <- mledist(x, dist, init = init, amt_format = TRUE, param_names = param_names, ...)
  
  class(out) <- c(sprintf("%s_distr", dist), "sl_distr", "amt_distr", "list")
  return(out)
  
}


# Standardized wrapper for fitting lognormal distribution
fit_lnorm <- function(x, 
                      init = "auto",
                      na.rm = TRUE, 
                      ...){
  if(na.rm){
    x <- x[!is.na(x)]
  }
  
  if(any(init == "auto")){
    
    n <- length(x)
    lx <- log(x)
    sd0 <- sqrt((n - 1)/n) * sd(lx)
    ml <- mean(lx)

    init[1] <- ml
    init[2] <- sd0
    init <- as.numeric(init)
  }
  
  param_names <- c("meanlog", "sdlog")
  
  dist <- "lnorm"
  out <- mledist(x, dist, init = init, amt_format = TRUE, param_names = param_names, ...)
  
  class(out) <- c(sprintf("%s_distr", dist), "sl_distr", "amt_distr", "list")
  return(out)
}

# Standardized wrapper for fitting vonmises distribution
fit_vonmises <- function(x, 
                         init = "auto",
                         na.rm = TRUE, 
                         ...){
  if(na.rm){
    x <- x[!is.na(x)]
  }
  
  if(any(init == "auto")){
    V <- mean(cos(x - atan2(sum(sin(x)), sum(cos(x)))))
    if (V > 0) {
      init <- circular:::A1inv(V)
    }
    else {
      init <- 0
    }
    init <- as.numeric(init)
  }
  
  param_names <- c("kappa")
  
  dist <- "vonmises"
  out <- mledist(x, dist, init = init, amt_format = TRUE, param_names = param_names, 
                 lower = 0, optim.method = "L-BFGS-B", ...)
  
  class(out) <- c(sprintf("%s_distr", dist), "ta_distr", "amt_distr", "list")
  return(out)
}

# function for making any fitting function with paste0("d",dist) defined
# auto_init takes a function that takes x and infer the initial values for optimization
make_fit <- function(dist, params, auto_init = NULL){
  den_fun <- sprintf("d%s", dist)
  if(is.null(get0(den_fun))){
    stop(sprintf("Cannot find '%s()' anywhere. Please define the density function. ", den_fun))
  }
  
  if(!all(params %in% names(formals(den_fun)))){
    extra <- params[!params %in% names(formals(den_fun))]
    stop(sprintf("Cannot find the argument '%s' in the function %s", extra, den_fun))
    
  }
  
  
  f <- function(x, 
                init = "auto",
                na.rm = TRUE, 
                ...){
    if(na.rm){
      x <- x[!is.na(x)]
    }
    
    if(any(init == "auto")){
      if(!is.null(auto_init)){
        init <- auto_init(x)
      } else {
        init <- rep(1, length(params))
      }
      init <- as.numeric(init)
    }
    
    param_names <- params
    out <- mledist(x, dist, init = init, amt_format = TRUE, param_names = param_names, ...)
    
    class(out) <- c(sprintf("%s_distr", dist), "sl_distr", "ta_distr", "amt_distr", "list")
    return(out)
  }
  assign(sprintf("fit_%s", dist), value = f, envir = globalenv())
  message(sprintf("Successfully created the function '%s()' in the global environment.", sprintf("fit_%s", dist)))
  return(invisible(NULL))
}


# Update gamma distribution params using fitted betas from ISSA
update_gamma <- function(dist, scale_estimator, shape_estimator){
  stopifnot(length(scale_estimator) == length(shape_estimator))
  lapply(seq_along(scale_estimator), function(i){
    amt::update_gamma(dist, scale_estimator[i], shape_estimator[i])
  })
  
}


# Update generalized von Mises params using fitted betas from ISSA
update_genvonmises <- function(dist, kappa1_estimator, kappa2_estimator){
  stopifnot(length(kappa1_estimator) == length(kappa2_estimator))
  lapply(seq_along(kappa1_estimator), function(i){
    new_kappa1 <- unname(dist$params$kappa1 + kappa1_estimator[i])
    new_kappa2 <- unname(dist$params$kappa2 + kappa2_estimator[i])
    make_genvonmises(new_kappa1, new_kappa2)
  })
}

# Update lognormal distribution params using fitted betas from ISSA
update_lnorm <- function(dist, mulog_estimator, sdlog_estimator ){
  stopifnot(length(mulog_estimator) == length(sdlog_estimator ))
  lapply(seq_along(mulog_estimator), function(i){
    amt::update_lnorm(dist, beta_log_sl = mulog_estimator[i], beta_log_sl_sq = sdlog_estimator [i])
  })
}

# Update vonmises distribution params using fitted betas from ISSA
update_vonmises <- function(dist, kappa_estimator){
  lapply(seq_along(kappa_estimator), function(i){
    amt::update_vonmises(
      dist, beta_cos_ta = kappa_estimator[i]
    )
  })
}


# Update generalized von Mises params using fitted betas from ISSA
update_invgamma <- function(dist, scale_estimator, shape_estimator){
  stopifnot(length(scale_estimator) == length(shape_estimator))
  lapply(seq_along(scale_estimator), function(i){
    new_shape <- unname(dist$params$shape - shape_estimator[i])
    new_scale <- unname(dist$params$scale - scale_estimator[i])
    make_invgamma(new_shape, new_scale)
  })
}


update_zigamma <- function(dist, p_estimator, scale_estimator, shape_estimator){
  stopifnot(
    length(p_estimator) == length(scale_estimator),
    length(scale_estimator) == length(shape_estimator)
  )
  warning("update_zigamma() not implemented.")
  lapply(seq_along(p_estimator), function(i){
    make_zigamma(p = NA, shape = NA, scale = NA)
  })
}


comp_dist <- function(x, dist = c("gamma", "lnorm", "invgamma"), criterion = "AICc"){
  x <- x[!is.na(x)]
  
  l <- lapply(dist, function(d){
    fit <- do.call(
      paste0("fit_", d),
      list(x)
    )
    
    res <- do.call(
      paste0("d", d),
      c(
        list(x),
        fit$params
      )
    )
    res <- sum(log(res))
    n <- length(x)
    
    class(res) <- "logLik"
    attr(res, "nall") <- n
    attr(res, "nobs") <- n
    attr(res, "df") <- length(fit$params)
    
    return(res)
  })
  
  out <- lapply(
    seq_along(l), function(i){
      cri <- do.call(criterion, l[i])
      df <- attributes(l[[i]])$df
      data.frame(model = dist[i], cri = cri, df = df)
    }
  ) %>% 
    do.call("rbind", .)
  
  out$dcri <- out$cri - min(out$cri)
  out <- out[order(out$cri), ]
  names(out)[2] <- criterion
  names(out)[4] <- paste0("d",criterion)
  
  return(out)
}



