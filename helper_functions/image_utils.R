plot_image_guide <- function(img, col = "red", cex = 0.8, main = NULL){
  ny <- nrow(img)
  nx <- ncol(img)
  
  d <- expand.grid(
    "y" = seq_len(ny),
    "x" = seq_len(nx)
  ) %>% 
    mutate(
      label = paste0(LETTERS[y],x)
    ) %>% 
    mutate(
      y = y/ny*(ny-1) + 0.5,
      x = x/nx*(nx-1) + 0.5
    )
  
  plot(img, main = main)
  text(x = d$x, y = d$y, col = col, label = d$label, cex = cex)
}


image_flatten <- function(img){
  as.numeric(c(img))
}

image_unflatten <- function(x){
  n <- length(x)
  ns <- sqrt(n)
  stopifnot(ns%%1 == 0)
  mat <- matrix(x, nrow = ns, ncol = ns)
  
  as.cimg(mat)
}



