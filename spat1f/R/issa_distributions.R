# General wrapper for amt_dist type object for getting density function
ddist <- function(dist, x_max = 1500, len = 100, return_x = TRUE){
  
  f <- function(dist){
    dist_name <- dist$name
    
    if(inherits(dist, "ta_distr")){
      x <- seq(-pi, pi, len = len)
    } else {
      x <- seq(0, x_max, len = len)
    }
    
    d <- do.call(
      paste0("d",dist_name),
      c(
        list(x), 
        dist$params
      )
    )
    
    d[!is.finite(d)] <- NA_real_
    
    if(return_x){
      data.frame("x" = x, "density" = d)
    }
  }
  
  if(!any(class(dist) %in% c("ta_distr", "sl_distr","amt_distr"))){
    lapply(dist, function(x){
      f(x)
    })
  } else {
    f(dist)
  }
  
}

# General wrapper for amt_dist type object for getting random number generation
rdist <- function(dist, n, ...){
  dots <- as.list(match.call(expand.dots = TRUE))[-c(1:3)]
  if(!any(class(dist) %in% c("ta_distr", "sl_distr","amt_distr"))){
    lapply(
      dist,
      function(x){
        do.call(paste0("r", x$name), c(list(n = n), x$params, dots)) 
      }
    ) 
  } else {
    do.call(paste0("r", dist$name), c(list(n = n), dist$params, dots)) 
  }
}




# Random number generation for generalized vonmises distribution
rgenvonmises <- function(n, kappa1, kappa2, max_try = 1000){
  rgenvonmisesC(n, kappa1, kappa2, max_try)
}
# .g_genvonmises <- function(omega, kappa1, kappa2){
#   kappa1 * cos(omega) + kappa2 * cos(2 * (omega - pi))
# }
# Standard ratio-of-uniforms algorithm from 10.1016/j.stamet.2006.11.003
# if(is.null(a)){
#   a <- exp(max(g_genvonmisesC(seq(0, 2 * pi, by = 0.00001), kappa1, kappa2))) 
# }
# n_out <- 0
# x <- c()
# 
# for (i in seq_len(max_try)){
#   # u <- runif(n_first_shot, 0, a)
#   # v <- runif(n_first_shot, 0, 2 * pi * a)
#   # s <- v > 2*pi*u
#   # u[s] <- a-u[s]
#   # v[s] <- 2 * pi * a - v[s]
#   # w <- g_genvonmisesC(v/u, kappa1, kappa2)/2
#   # cond <- u <= exp(w)
#   # x_temp <- v[cond]/u[cond]
#   x_temp <- propose_genvonmises(n, a, kappa1, kappa2)
#   x <- c(x, x_temp)
#   n_out <- length(x)
#   if(n_out > n){
#     break
#   } else {
#     if(i == max_try){
#       stop(sprintf("Function timed out after %s tries. Generated %s percent of requested numbers. \nConsider increasing the 'max_try' argument.", max_try, n_out / n * 100))
#     }
#   }
# }
# return(x[seq_len(n)] - pi)




# Density function for generalized von mises distribution
dgenvonmises <- function(x, kappa1, kappa2, log = FALSE){
  num <- exp(g_genvonmisesC(x + pi, kappa1, kappa2))
  den <- integrate(function(x) {
    exp(g_genvonmisesC(x + pi, kappa1, kappa2))
  },lower =  0,upper =  2 * pi)$value
  
  d <- num / den
  
  if(log){
    d <- log(d)
  }
  return(d)
}

# create amt dist for inverse gamma
make_invgamma <- function(shape, scale){
  param_names <- c("shape", "scale")
  
  params <- as.list(c(shape, scale))
  names(params) <- param_names
  
  out <- list("name" = "invgamma",
              "params" = params,
              "vcov" = NULL)
  
  class(out) <- c("invgamma_distr", "sl_distr", "amt_distr", "list")
  return(out)
}


# create amt dist for generalized von mises 
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

# create amt dist for generalized von mises 
make_vonmises <- function(kappa){
  param_names <- c("kappa")
  
  params <- as.list(c(kappa))
  names(params) <- param_names
  
  out <- list("name" = "vonmises",
              "params" = params,
              "vcov" = NULL)
  
  class(out) <- c("vonmises_distr", "ta_distr", "amt_distr", "list")
  return(out)
}

# wrapper for amt dist gamma
make_gamma <- function(shape, scale){
  amt::make_gamma_distr(shape = shape, scale = scale)
}

# Wrapper for amt dist lnorm
make_lnorm <- function(meanlog, sdlog){
  amt::make_lnorm_distr(meanlog = meanlog, sdlog = sdlog)
}


# make amt dist 
make_zigamma <- function(p, shape, scale){
  param_names <- c("p", "shape", "scale")
  
  params <- as.list(c(p, shape, scale))
  names(params) <- param_names
  
  out <- list("name" = "zigamma",
              "params" = params,
              "vcov" = NULL)
  
  class(out) <- c("zigamma_distr", "sl_distr", "amt_distr", "list")
  return(out)
}

# make amt dist 
make_unif <- function(min = -pi, max = pi){
  param_names <- c("min", "max")
  
  params <- as.list(c(min, max))
  names(params) <- param_names
  
  out <- list("name" = "unif",
              "params" = params,
              "vcov" = NULL)
  
  class(out) <- c("unif_distr", "ta_distr", "amt_distr", "list")
  return(out)
}

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


# kth moment of absolute value of genvonmises dist
genvonmises_abs_moments <- function(kappa1, kappa2, k = 1){
  if(is.na(kappa1) | is.na(kappa2) | is.na(k)){
    return(NA)
  } else {
    integrate(
      function(x){
        dgenvonmises(x, kappa1, kappa2) * abs(x^k)
      }, lower = -pi, upper = pi
    )$value 
  }
}

rinvgamma <- function(n, shape, scale){
  1 / rgamma(n, shape = shape, scale = scale)
}

dinvgamma <- function(x, shape, scale, log = FALSE){
  extraDistr::dinvgamma(x, alpha = shape, beta = scale, log = log)
}



rvonmises <- function(n, kappa){
  x <- as.numeric(circular::rvonmises(n = n, mu = circular::circular(pi), kappa = kappa))
  x <- x %%(2 * pi)
  ifelse(x > base::pi, x - (2 * base::pi), x)
}

dvonmises <- function(x, kappa, log = FALSE){
  # x <- (x + pi)
  circular::dvonmises(circular::circular(x), 
                      mu = circular::circular(-pi), 
                      kappa = kappa, 
                      log = log)
}

#### Add beta prime distribution


