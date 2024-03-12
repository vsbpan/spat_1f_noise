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


update_genvonmises <- function(dist, beta_cos_theta_pi, beta_cos_2theta){
  stopifnot(length(beta_cos_theta_pi) == length(beta_cos_2theta))
  lapply(seq_along(beta_cos_theta_pi), function(i){
    new_kappa1 <- unname(dist$params$kappa1 + beta_cos_theta_pi[i])
    new_kappa2 <- unname(dist$params$kappa2 + beta_cos_2theta[i])
    make_genvonmises(new_kappa1, new_kappa2)
  })
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

fit_gamma <- function(x){
  amt::fit_distr(x, "gamma")
}

update_gamma <- function(dist, beta_sl, beta_log_sl){
  stopifnot(length(beta_sl) == length(beta_log_sl))
  lapply(seq_along(beta_sl), function(i){
    amt::update_gamma(dist, beta_sl[i], beta_log_sl[i])
  })
  
}

fit_lnorm <- function(x){
  amt::fit_distr(x, "lnorm")
}

update_lnorm <- function(dist, beta_log_sl, beta_log_sl_sq){
  stopifnot(length(beta_log_sl) == length(beta_log_sl_sq))
  lapply(seq_along(beta_log_sl), function(i){
    amt::update_lnorm(dist, beta_log_sl = beta_log_sl[i], beta_log_sl_sq = beta_log_sl_sq[i])
  })
}

update_vonmises <- function(dist, beta_cos_ta){
  lapply(seq_along(beta_cos_ta), function(i){
    amt::update_vonmises(
      dist, beta_cos_ta = beta_cos_ta[i]
    )
  })

}


comp_dist <- function(x){
  ln_fit <- amt:::fit_distr(x, "lnorm")
  gamma_fit <- amt:::fit_distr(x, "gamma")
  
  d_ln <- do.call(
    paste0("d",ln_fit$name),
    c(
      list(x), 
      ln_fit$params
    )
  )
  d_gamma <- do.call(
    paste0("d", gamma_fit$name),
    c(
      list(x), 
      gamma_fit$params
    )
  )
  out <- c("gamma" = sum(log(d_gamma)), "lnorm" = sum(log(d_ln)))
  out <- c(out, "delta_AIC" = unname(diff(out)))
  c(out, "lnorm_wins" = out["delta_AIC"] > 0)
}



