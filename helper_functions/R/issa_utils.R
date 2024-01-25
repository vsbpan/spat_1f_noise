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

print.issf_fit <- function(x, ...){
  print(summary(x$model))
}

registerS3method("print", "issf_fit", print.issf_fit)


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


rgenvonmises <- function(n, kappa1, kappa2, max_try = 500){
  g <- function(omega){
    kappa1 * cos(omega) + kappa2 * cos(2 * (omega - pi))
  }
  # m <- max(g(seq(0, 2 * pi, by = 0.0001)))
  # n_out <- 0
  # n_first_shot <- ceiling(n * 2L)
  # x <- c()
  # 
  # while(n_out < n){
  #   u <- runif(n_first_shot, 0, 2 * pi)
  #   v <- runif(n_first_shot, 0, m)
  #   x_temp <- u[v <= g(u)]
  #   x <- c(x, x_temp)
  #   n_out <- length(x)
  # }
  # return(x[seq_len(n)] - pi)
  
  # Standard ratio-of-uniforms algorithm from 10.1016/j.stamet.2006.11.003
  a <- exp(max(g(seq(0, 2 * pi, by = 0.00001))))
  n_out <- 0
  n_first_shot <- ceiling(n * 5L)
  x <- c()

  for (i in seq_len(max_try)){
    u <- runif(n_first_shot, 0, a)
    v <- runif(n_first_shot, 0, 2 * pi * a)
    s <- v > 2*pi*u
    u[s] <- a-u[s]
    v[s] <- 2 * pi * a - v[s]
    w <- g(v/u)/2
    cond <- u <= exp(w)
    
    x_temp <- v[cond]/u[cond]
    x <- c(x, x_temp)
    n_out <- length(x)
    if(n_out > n){
      break
    } else {
      if(i == max_try){
        stop(sprintf("Function timed out after %s tries. Generated %s percent of requested numbers.", max_try, n_out / n * 100))
      }
    }
  }
  return(x[seq_len(n)] - pi)
}


dgenvonmises <- function(x, kappa1, kappa2, log = FALSE){
  num <- exp(kappa1 * cos(x + pi) + kappa2 * cos(2 * x))
  den <- integrate(function(x) {
    exp(kappa1 * cos(x + pi) + kappa2 * cos(2 * x))
  },lower =  0,upper =  2 * pi)$value
  
  d <- num / den
  
  if(log){
    d <- log(d)
  }
  return(d)
}


make_genvonmises <- function(kappa1, kappa2){
  param_names <- c("kappa1", "kappa2")
  
  params <- as.list(c(kappa1, kappa2))
  names(params) <- param_names
  
  out <- list("name" = "genvonmises",
              "params" = params,
              "vcov" = NULL)
  
  class(out) <- c("genvonmises_distr", "ta_distr", "amt_distr", "list")
  return(out)
}


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
  new_kappa1 <- unname(dist$params$kappa1 + beta_cos_theta_pi)
  new_kappa2 <- unname(dist$params$kappa2 + beta_cos_2theta)
  make_genvonmises(new_kappa1, new_kappa2)
}

append_genvonmises_estimators <- function(x){
  x %>% 
    mutate(
      cos_theta_pi = cos(theta_rel + pi),
      cos_2theta = cos(2 * (theta_rel)),
    )
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
# update_zigamma <- function(dist, beta_move){
#   new_p <- unname(plogis(qlogis(dist$params$p) - beta_move))
#   make_zigamma(new_p, dist$params$shape, dist$params$scale)
# }
# 
# make_zigamma <- function(p, shape, scale){
#   param_names <- c("p", "shape", "scale")
#   
#   params <- as.list(c(p, shape, scale))
#   names(params) <- param_names
#   
#   out <- list("name" = "zigamma", 
#               "params" = params, 
#               "vcov" = NULL) 
#   
#   class(out) <- c("zigamma_distr", "sl_distr", "amt_distr", "list")
#   return(out)
# }


# Is a hack that fits zigamma in two parts. 
# fit_zigamma2 <- function(x, na.rm = TRUE){
#   if(na.rm){
#     x <- x[!is.na(x)]
#   }
#   
#   zero <- ifelse(x > 0, 0, 1)
#   x2 <- x[!as.logical(zero)]
#   
#   pfit <- optim2(c(0.1), fn = function(theta){
#     -sum(dbinom(zero, size = 1, prob = theta/1000, log = TRUE))
#   }, 
#   method = "Brent", 
#   lower = 0, 
#   upper = 1000, 
#   silent = TRUE)
#   
#   
#   param_names <- c("p", "shape", "scale")
#   
#   params_gamma <- amt::fit_distr(x2, "gamma")$params
#   
#   params <- list(
#     pfit$par / 1000,
#     params_gamma[[1]],
#     params_gamma[[2]]
#   )
#   names(params) <- param_names
#   
#   out <- list("name" = "zigamma", 
#               "params" = params, 
#               "vcov" = NA) # VCV matrix on transformed scale
#   
#   class(out) <- c("zigamma_distr", "sl_distr", "amt_distr", "list")
#   return(out)
# }

# Zero hurdle (inflated) gamma distribution random number generator
rzigamma <- function(n, p, shape, scale = 1/rate, rate){
  x <- rbinom(n, size = 1, prob = 1-p)
  x[x == 1] <- rgamma(sum(x), shape = shape, scale = scale)
  return(x)
}

# Zero hurdle (inflated) gamma distribution density function
dzigamma <- function (x, p, shape, scale = 1/rate, rate, log = FALSE) {
  x_zero <- x == 0
  d <- ifelse(x_zero, p, NA)
  d[!x_zero] <- dgamma(x[!x_zero], shape = shape, scale = scale) * (1 - p)
  if(log){
    d <- log(d)
  }
  return(d)
}

# Fit zigamma returning the same format as `amt::fit_distr()`
# fit_zigamma <- function(x, 
#                         method = c("Nelder-Mead", "BFGS", "nlminb", "nlm"),
#                         init = c(0.1, 1, 30),
#                         lower = c(0, 0, 0),
#                         upper = c(1, Inf, Inf),
#                         parscale = c(1000, 10, 100),
#                         na.rm = TRUE){
#   
#   method <- match.arg(method)
#   
#   if(na.rm){
#     x <- x[!is.na(x)]
#   }
#   
#   param_names <- c("p", "shape", "scale")
#   
#   fit <- optim2(
#     init = init * parscale, 
#     fn = function(theta){
#       -sum(dzigamma(x, 
#                     p = (theta[1]) / parscale[1], # Transform to deal with numerical issue
#                     shape = (theta[2]) / parscale[2], 
#                     scale = (theta[3]) / parscale[3], log = TRUE))
#     }, 
#     lower = lower * parscale,
#     upper = upper * parscale,
#     method = method
#   )
#   
#   if(fit$convergence != 0){
#     message(sprintf("Did not converge.\nCode: %s\nMessage: %s", fit$convergence, fit$message))
#     vcov <- matrix(rep(NA, length(param_names)^2), 
#                    nrow = length(param_names), 
#                    ncol = length(param_names))
#   } else {
#     vcov <- solve(fit$hessian) / crossprod(x = t(matrix(parscale)))
#   }
#   
#   params <- as.list(fit$par)
#   params[[1]] <- (params[[1]]) / parscale[1] # Back transform
#   params[[2]] <- (params[[2]]) / parscale[2]
#   params[[3]] <- (params[[3]]) / parscale[3]
#   
#   names(params) <- param_names
#   rownames(vcov) <- colnames(vcov) <- param_names
#   
#   out <- list("name" = "zigamma", 
#               "params" = params, 
#               "vcov" = vcov) # VCV matrix on transformed scale
#   
#   class(out) <- c("zigamma_distr", "sl_distr", "amt_distr", "list")
#   return(out)
# }
