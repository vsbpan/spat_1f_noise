
update_cengamma <- function(dist, beta_sl, beta_log_sl){
  new_shape <- unname(dist$params$shape + beta_log_sl)
  new_scale <- unname(1/((1/dist$params$scale) - beta_sl))
  make_cengamma(new_shape, new_scale, dist$params$thresh)
}

fit_cengamma <- function(x, 
                         thresh = 10,
                         method = c("Nelder-Mead", "BFGS", "nlminb", "nlm"),
                         init = c(1, 30),
                         lower = c(0, 0),
                         upper = c(Inf, Inf),
                         parscale = c(10, 100),
                         na.rm = TRUE){
  
  method <- match.arg(method)
  
  if(na.rm){
    x <- x[!is.na(x)]
  }
  
  param_names <- c("shape", "scale", "thresh")
  
  fit <- optim2(
    init = init * parscale,
    fn = function(theta){
      -sum(dcengamma(x,
                     shape = (theta[1]) / parscale[1], # Transform to deal with numerical issue
                     scale = (theta[2]) / parscale[2], 
                     thresh = thresh,
                     log = TRUE))
    },
    lower = lower * parscale,
    upper = upper * parscale,
    method = method
  )
  
  if(fit$convergence != 0){
    message(sprintf("Did not converge.\nCode: %s\nMessage: %s", fit$convergence, fit$message))
    vcov <- matrix(rep(NA, (length(param_names)-1)^2),
                   nrow = length(param_names)-1,
                   ncol = length(param_names)-1)
  } else {
    vcov <- solve(fit$hessian) / crossprod(x = t(matrix(parscale)))
  }
  
  params <- as.list(fit$par)
  params[[1]] <- (params[[1]]) / parscale[1] # Back transform
  params[[2]] <- (params[[2]]) / parscale[2]
  params[[3]] <- thresh
  
  names(params) <- param_names
  rownames(vcov) <- colnames(vcov) <- param_names[1:2]
  
  out <- list("name" = "cengamma",
              "params" = params,
              "vcov" = vcov) # VCV matrix on transformed scale
  
  class(out) <- c("cengamma_distr", "ta_distr", "amt_distr", "list")
  return(out)
  
}








fit_cengamma2 <- function(x, 
                          thresh = 10,
                          method = c("Nelder-Mead", "BFGS", "nlminb", "nlm"),
                          init = c(1, 30),
                          lower = c(0, 0),
                          upper = c(Inf, Inf),
                          parscale = c(10, 100),
                          na.rm = TRUE){
  
  method <- match.arg(method)
  
  if(na.rm){
    x <- x[!is.na(x)]
  }
  
  param_names <- c("shape", "scale", "thresh")
  
  fit <- optim2(
    init = init * parscale,
    fn = function(theta){
      -sum(dcengamma2(x,
                      shape = (theta[1]) / parscale[1], # Transform to deal with numerical issue
                      scale = (theta[2]) / parscale[2], 
                      thresh = thresh,
                      log = TRUE))
    },
    lower = lower * parscale,
    upper = upper * parscale,
    method = method
  )
  
  if(fit$convergence != 0){
    message(sprintf("Did not converge.\nCode: %s\nMessage: %s", fit$convergence, fit$message))
    vcov <- matrix(rep(NA, (length(param_names)-1)^2),
                   nrow = length(param_names)-1,
                   ncol = length(param_names)-1)
  } else {
    vcov <- solve(fit$hessian) / crossprod(x = t(matrix(parscale)))
  }
  
  params <- as.list(fit$par)
  params[[1]] <- (params[[1]]) / parscale[1] # Back transform
  params[[2]] <- (params[[2]]) / parscale[2]
  params[[3]] <- thresh
  
  names(params) <- param_names
  rownames(vcov) <- colnames(vcov) <- param_names[1:2]
  
  out <- list("name" = "cengamma",
              "params" = params,
              "vcov" = vcov) # VCV matrix on transformed scale
  
  class(out) <- c("cengamma_distr", "ta_distr", "amt_distr", "list")
  return(out)
  
}


append_cengamma_estimators <- function(x, thresh){
  x %>% 
    mutate(
      sl = r,
      logsl = ifelse(r > thresh, log(r), 0)
    )
}



rcengamma <- function(n, shape, scale, thresh){
  x <- rgamma(n, shape = shape, scale = scale)
  ifelse(x > thresh, x, 0)
}

dcengamma <- function(x, shape, scale, thresh, log = FALSE){
  s <- x > thresh
  p0 <- pgamma(thresh, shape = shape, scale = scale)
  d <- ifelse(!s, 
              p0, 
              NA)
  d[s] <- (1-p0) * dgamma(x = x[s], shape = shape, scale = scale)
  
  if(log){
    d <- log(d)
  }
  return(d)
}

dcengamma2 <- function(x, shape, scale, thresh, log = FALSE){
  x <- x[x > thresh]
  
  p0 <- pgamma(thresh, shape = shape, scale = scale)
  d <- dgamma(x = x, shape = shape, scale = scale) / (1-p0)
  
  if(log){
    d <- log(d)
  }
  return(d)
}



rbivonmises <- function(n, kappa1, kappa2){
  x <- circular::rmixedvonmises(n = n, 
                                mu1 = circular::circular(pi), 
                                mu2 = circular::circular(0), 
                                kappa1 = kappa1, 
                                kappa2 = kappa2, 
                                prop = 0.5) %>% 
    as.numeric()
  x - pi
}

dbivonmises <- function(x, kappa1, kappa2, log = FALSE){
  d <- circular::dmixedvonmises(circular::circular(x + pi),
                                mu1 = circular::circular(pi), 
                                mu2 = circular::circular(0), 
                                kappa1 = kappa1, 
                                kappa2 = kappa2, 
                                prop = 0.5) %>% 
    as.numeric()
  if(log){
    d <- log(d)
  }
  return(d)
}

make_bivonmises <- function(kappa1, kappa2){
  param_names <- c("kappa1", "kappa2")
  
  params <- as.list(c(kappa1, kappa2))
  names(params) <- param_names
  
  out <- list("name" = "bivonmises",
              "params" = params,
              "vcov" = NULL)
  
  class(out) <- c("bivonmises_distr", "ta_distr", "amt_distr", "list")
  return(out)
}


make_cengamma <- function(shape, scale, thresh){
  param_names <- c("shape", "scale", "thresh")
  
  params <- as.list(c(shape, scale, thresh))
  names(params) <- param_names
  
  out <- list("name" = "cengamma",
              "params" = params,
              "vcov" = NULL)
  
  class(out) <- c("cengamma_distr", "sl_distr", "amt_distr", "list")
  return(out)
}



fit_bivonmises <- function(x, 
                           method = c("Nelder-Mead", "BFGS", "nlminb", "nlm"),
                           init = c(1, 1),
                           lower = c(0, 0),
                           upper = c(Inf, Inf),
                           parscale = c(100, 100),
                           na.rm = TRUE){
  
  method <- match.arg(method)
  
  if(na.rm){
    x <- x[!is.na(x)]
  }
  
  param_names <- c("kappa1", "kappa2")
  
  fit <- optim2(
    init = init * parscale,
    fn = function(theta){
      -sum(dbivonmises(x,
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
  
  out <- list("name" = "bivonmises",
              "params" = params,
              "vcov" = vcov) # VCV matrix on transformed scale
  
  class(out) <- c("bivonmises_distr", "ta_distr", "amt_distr", "list")
  return(out)
  
}
