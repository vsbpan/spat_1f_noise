
# Density function of double pareto lognormal distribution
ddpln <- function(x, alpha, beta, meanlog, sdlog, log = FALSE){
  distributionsrd::ddoubleparetolognormal(x = x, shape1 = alpha, shape2 = beta, meanlog = meanlog, sdlog = sdlog, log = log)
}

# random number generation of double pareto lognormal distribution
rdpln <- function(n, alpha, beta, meanlog, sdlog){
  distributionsrd::rdoubleparetolognormal(n, shape1 = alpha, shape2 = beta, meanlog = meanlog, sdlog = sdlog)
}


dfrechet <- function(x, shape, scale, log = FALSE){
  extraDistr::dfrechet(x = x, lambda = shape, mu = 0, sigma = scale, log = log)
}

rfrechet <- function(n, shape, scale){
  extraDistr::dfrechet(n, lambda = shape, mu = 0, sigma = scale)
}


fit_dpln <- function(x, 
                     init = "auto",
                     na.rm = TRUE, 
                     ...){
  if(na.rm){
    x <- x[!is.na(x)]
  }
  
  if(any(init == "auto")){
    init[1] <- 2
    init[2] <- 2
    init[4] <- sd(log(x))
    init <- as.numeric(init)
    init[3] <- mean(log(x)) - 1/init[1] - 1/init[2]
  }
  
  param_names <- c("alpha", "beta", "meanlog", "sdlog")
  dist <- "dpln"
  out <- mledist(x, dist, init = init, amt_format = TRUE, param_names = param_names, ...)
  
  class(out) <- c(sprintf("%s_distr", dist), "sl_distr", "amt_distr", "list")
  return(out)
  
}

# Fit inverse gamma distribution with MLE
fit_frechet <- function(x, 
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
  
  dist <- "frechet"
  out <- mledist(x, dist, init = init, amt_format = TRUE, param_names = param_names, ...)
  
  class(out) <- c(sprintf("%s_distr", dist), "sl_distr", "amt_distr", "list")
  return(out)
}
