dbigamma <- function(x , a, b, d, log = FALSE){
  # res <- (1 + (1- d * x)^2)  * b^a / gamma(a) * x^(a-1) * exp(-b*x)
  # / (2 + a * d / b * ((1+a) * d/b - 2) )
  
  res <- (d * x)^2  * b^a / gamma(a) * x^(a-1) * exp(-b*x)
  
  if(log){
    res <- log(res)
  }
  return(res)
}

x <- seq(0, 20, by = 0.01)

data.frame(
  log(x),
  "y" = dbigamma(x, a = 5, b = 0.7, 0.09)
) %>% 
  plot()
  





