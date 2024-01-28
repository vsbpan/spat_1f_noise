rgenvonmises <- function(n, kappa1, kappa2, a = NULL, max_try = 500){
  .g_genvonmises <- function(omega, kappa1, kappa2){
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
  if(is.null(a)){
    a <- exp(max(.g_genvonmises(seq(0, 2 * pi, by = 0.00001), kappa1, kappa2))) 
  }
  n_out <- 0
  n_first_shot <- ceiling(n * 5L)
  x <- c()
  
  for (i in seq_len(max_try)){
    u <- runif(n_first_shot, 0, a)
    v <- runif(n_first_shot, 0, 2 * pi * a)
    s <- v > 2*pi*u
    u[s] <- a-u[s]
    v[s] <- 2 * pi * a - v[s]
    w <- .g_genvonmises(v/u, kappa1, kappa2)/2
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


make_gamma <- function(shape, scale){
  amt::make_gamma_distr(shape = shape, scale = scale)
}

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



