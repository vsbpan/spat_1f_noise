# Nice wrapper for fitting generalized von Mises distribution using MLE. 
fit_genvonmises <- function(x, 
                            method = c("Nelder-Mead", "BFGS", "nlminb", "nlm"),
                            init = c(1, 0.5),
                            lower = c(0, 0),
                            upper = c(Inf, Inf),
                            parscale = c(1000, 1000),
                            na.rm = TRUE){
  
  method <- match.arg(method)
  
  if(na.rm){
    x <- x[!is.na(x)]
  }
  
  param_names <- c("kappa1", "kappa2")
  
  fit <- optim2(
    init = init * parscale,
    fn = function(theta){
      -sum(dgenvonmises(x,
                        kappa1 = (theta[1]) / parscale[1], # Transform to deal with numerical issue
                        kappa2 = (theta[2]) / parscale[2], 
                        log = TRUE))
    },
    lower = lower * parscale,
    upper = upper * parscale,
    method = method
  )
  
  if(fit$convergence != 0){
    message(sprintf("Did not converge.\nCode: %s\nMessage: %s", fit$convergence, fit$message))
    vcov <- matrix(rep(NA, length(param_names)^2),
                   nrow = length(param_names),
                   ncol = length(param_names))
  } else {
    vcov <- solve(fit$hessian) / crossprod(x = t(matrix(parscale)))
  }
  
  params <- as.list(fit$par)
  params[[1]] <- (params[[1]]) / parscale[1] # Back transform
  params[[2]] <- (params[[2]]) / parscale[2]
  
  names(params) <- param_names
  rownames(vcov) <- colnames(vcov) <- param_names
  
  out <- list("name" = "genvonmises",
              "params" = params,
              "vcov" = vcov) # VCV matrix on transformed scale
  
  class(out) <- c("genvonmises_distr", "ta_distr", "amt_distr", "list")
  return(out)
  
}

# Is a hack that fits zigamma in two parts. 
fit_zigamma2 <- function(x, na.rm = TRUE){
  if(na.rm){
    x <- x[!is.na(x)]
  }
  
  zero <- ifelse(x > 0, 0, 1)
  x2 <- x[!as.logical(zero)]
  
  pfit <- optim2(c(0.1), fn = function(theta){
    -sum(dbinom(zero, size = 1, prob = theta/1000, log = TRUE))
  },
  method = "Brent",
  lower = 0,
  upper = 1000,
  silent = TRUE)
  
  
  param_names <- c("p", "shape", "scale")
  
  params_gamma <- amt::fit_distr(x2, "gamma")$params
  
  params <- list(
    pfit$par / 1000,
    params_gamma[[1]],
    params_gamma[[2]]
  )
  names(params) <- param_names
  
  out <- list("name" = "zigamma",
              "params" = params,
              "vcov" = NA) # VCV matrix on transformed scale
  
  class(out) <- c("zigamma_distr", "sl_distr", "amt_distr", "list")
  return(out)
}

# Fit zigamma returning the same format as `amt::fit_distr()`
fit_zigamma <- function(x,
                        method = c("Nelder-Mead", "BFGS", "nlminb", "nlm"),
                        init = c(0.1, 1, 30),
                        lower = c(0, 0, 0),
                        upper = c(1, Inf, Inf),
                        parscale = c(1000, 10, 100),
                        na.rm = TRUE){
  
  method <- match.arg(method)
  
  if(na.rm){
    x <- x[!is.na(x)]
  }
  
  param_names <- c("p", "shape", "scale")
  
  fit <- optim2(
    init = init * parscale,
    fn = function(theta){
      -sum(dzigamma(x,
                    p = (theta[1]) / parscale[1], # Transform to deal with numerical issue
                    shape = (theta[2]) / parscale[2],
                    scale = (theta[3]) / parscale[3], log = TRUE))
    },
    lower = lower * parscale,
    upper = upper * parscale,
    method = method
  )
  
  if(fit$convergence != 0){
    message(sprintf("Did not converge.\nCode: %s\nMessage: %s", fit$convergence, fit$message))
    vcov <- matrix(rep(NA, length(param_names)^2),
                   nrow = length(param_names),
                   ncol = length(param_names))
  } else {
    vcov <- solve(fit$hessian) / crossprod(x = t(matrix(parscale)))
  }
  
  params <- as.list(fit$par)
  params[[1]] <- (params[[1]]) / parscale[1] # Back transform
  params[[2]] <- (params[[2]]) / parscale[2]
  params[[3]] <- (params[[3]]) / parscale[3]
  
  names(params) <- param_names
  rownames(vcov) <- colnames(vcov) <- param_names
  
  out <- list("name" = "zigamma",
              "params" = params,
              "vcov" = vcov) # VCV matrix on transformed scale
  
  class(out) <- c("zigamma_distr", "sl_distr", "amt_distr", "list")
  return(out)
}


# Fit inverse gamma distribution with MLE
fit_invgamma <- function(x, 
                         method = c("Nelder-Mead", "BFGS", "nlminb", "nlm"),
                         init = c(1, 0.5),
                         lower = c(0, 0),
                         upper = c(Inf, Inf),
                         parscale = c(1000, 1000),
                         na.rm = TRUE){
  
  method <- match.arg(method)
  
  if(na.rm){
    x <- x[!is.na(x)]
  }
  
  param_names <- c("shape", "scale")
  
  fit <- optim2(
    init = init * parscale,
    fn = function(theta){
      -sum(dinvgamma(x,
                     shape = (theta[1]) / parscale[1], # Transform to deal with numerical issue
                     scale = (theta[2]) / parscale[2], 
                     log = TRUE))
    },
    lower = lower * parscale,
    upper = upper * parscale,
    method = method
  )
  
  if(fit$convergence != 0){
    message(sprintf("Did not converge.\nCode: %s\nMessage: %s", fit$convergence, fit$message))
    vcov <- matrix(rep(NA, length(param_names)^2),
                   nrow = length(param_names),
                   ncol = length(param_names))
  } else {
    vcov <- solve(fit$hessian) / crossprod(x = t(matrix(parscale)))
  }
  
  params <- as.list(fit$par)
  params[[1]] <- (params[[1]]) / parscale[1] # Back transform
  params[[2]] <- (params[[2]]) / parscale[2]
  
  names(params) <- param_names
  rownames(vcov) <- colnames(vcov) <- param_names
  
  out <- list("name" = "invgamma",
              "params" = params,
              "vcov" = vcov) # VCV matrix on transformed scale
  
  class(out) <- c("invgamma_distr", "ta_distr", "amt_distr", "list")
  return(out)
}


# Fit inverse gamma distribution with MLE
fit_frechet <- function(x, 
                         method = c("Nelder-Mead", "BFGS", "nlminb", "nlm"),
                         init = c(1, 0.5),
                         lower = c(0, 0),
                         upper = c(Inf, Inf),
                         parscale = c(1000, 1000),
                         na.rm = TRUE){
  
  method <- match.arg(method)
  
  if(na.rm){
    x <- x[!is.na(x)]
  }
  
  param_names <- c("shape", "scale")
  
  fit <- optim2(
    init = init * parscale,
    fn = function(theta){
      -sum(dfrechet(x,
                     shape = (theta[1]) / parscale[1], # Transform to deal with numerical issue
                     scale = (theta[2]) / parscale[2], 
                     log = TRUE))
    },
    lower = lower * parscale,
    upper = upper * parscale,
    method = method
  )
  
  if(fit$convergence != 0){
    message(sprintf("Did not converge.\nCode: %s\nMessage: %s", fit$convergence, fit$message))
    vcov <- matrix(rep(NA, length(param_names)^2),
                   nrow = length(param_names),
                   ncol = length(param_names))
  } else {
    vcov <- solve(fit$hessian) / crossprod(x = t(matrix(parscale)))
  }
  
  params <- as.list(fit$par)
  params[[1]] <- (params[[1]]) / parscale[1] # Back transform
  params[[2]] <- (params[[2]]) / parscale[2]
  
  names(params) <- param_names
  rownames(vcov) <- colnames(vcov) <- param_names
  
  out <- list("name" = "frechet",
              "params" = params,
              "vcov" = vcov) # VCV matrix on transformed scale
  
  class(out) <- c("frechet_distr", "ta_distr", "amt_distr", "list")
  return(out)
}



# update_zigamma_p <- function(dist, beta_move){
#   new_p <- unname(plogis(qlogis(dist$params$p) - beta_move))
#   make_zigamma(new_p, dist$params$shape, dist$params$scale)
# }
# 
# update_zigamma_ss <- function(dist, beta_sl, beta_log_sl){
#   new_shape <- unname(dist$params$shape + beta_log_sl)
#   new_scale <- unname(1/((1/dist$params$scale) - beta_sl))
#   make_zigamma(dist$params$p, new_shape, new_scale)
# }

# Standardized wrapper for fitting gamma distribution
fit_gamma <- function(x){
  amt::fit_distr(x, "gamma")
}


# Standardized wrapper for fitting lognormal distribution
fit_lnorm <- function(x){
  amt::fit_distr(x, "lnorm")
}

# Standardized wrapper for fitting vonmises distribution
fit_vonmises <- function(x){
  amt::fit_distr(x, "vonmises")
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
update_vonmises <- function(dist, kappa_esimator){
  lapply(seq_along(kappa_esimator), function(i){
    amt::update_vonmises(
      dist, beta_cos_ta = kappa_esimator[i]
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


comp_dist <- function(x, dist = c("gamma", "lnorm", "invgamma", "frechet"), criterion = "AICc"){
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



