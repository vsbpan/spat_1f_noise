# Takes a matrix and threshold by the amplitude in the frequency domain
fourier_thresh <- function(mat, quantile){
  r.mat <- c(Mod(mat))
  thresh <- quantile(r.mat, probs = quantile)
  mat[r.mat < thresh] <- 0
  return(mat)
} 

# Takes a matrix and run a low or high pass filter through the spectrum
fourier_pass <- function(mat, r, side = c("lower","higher", "between", "outside")){
  side <- match.arg(side)
  mat <- fftshift(mat)
  
  if(side %in% c("outside","between")){
    stopifnot(length(r) == 2)
  } else {
    stopifnot(length(r) == 1)
  }
  
  if(!all(r < round(dim(mat)/2))){
    stop("r too high. Try r < ", floor(min(dim(mat))/2))
  }
  
  R <- expand.grid( "x" = seq_len(nrow(mat)), 
                    "y" = seq_len(ncol(mat))) %>% 
    cbind("xm"= round((nrow(mat))/2), 
          "ym"= round((ncol(mat))/2)) %>% 
    with(., sqrt((x-xm)^2 + (y-ym)^2))
  mat[switch(side, 
             "lower" = R > r, 
             "higher" = R < r, 
             "outside" = R >= r[1] & R <= r[2], 
             "between" = R < r[1] | R > r[2])] <- 0
  mat <- fftshift(mat, inverse = TRUE)
  return(mat)
}


# Rearrange the matrix into the highest frequency starts in the middle of the matrix, not the top left of the matrix, as is the output of mvfft(). Setting inverse = TRUE turns it back. 
fftshift <- function(input_matrix, dim = -1, inverse = FALSE) {
  rows <- dim(input_matrix)[1]    
  cols <- dim(input_matrix)[2]
  
  if(inverse){
    foo <- floor
  } else {
    foo <- ceiling
  }
  
  swap_up_down <- function(input_matrix) {
    rows_half <- foo(rows/2)
    return(rbind(input_matrix[((rows_half+1):rows), (1:cols)], input_matrix[(1:rows_half), (1:cols)]))
  }
  
  swap_left_right <- function(input_matrix) {
    cols_half <- foo(cols/2)
    return(cbind(input_matrix[1:rows, ((cols_half+1):cols)], input_matrix[1:rows, 1:cols_half]))
  }
  
  
  if (dim == -1) {
    if(inverse){
      input_matrix <- swap_left_right(input_matrix)
      input_matrix <- swap_up_down(input_matrix)
    } else {
      input_matrix <- swap_up_down(input_matrix)
      input_matrix <- swap_left_right(input_matrix) 
    }
  }
  else if (dim == 1) {
    input_matrix <- swap_up_down(input_matrix)
  }
  else if (dim == 2) {
    input_matrix <- swap_left_right(input_matrix)
  }
  else {
    stop("Invalid dimension parameter")
  }
  return(input_matrix)
}

# Retrun fft transformed image back as a cimg object
fft_img_inv <- function(x, as.cimg = TRUE){
  inverted.mat <- Re(fft(x,inverse = TRUE)) / length(x)
  inverted.mat[inverted.mat < 0] <- 0
  if(as.cimg){
    inverted.mat <- as.cimg(inverted.mat)
  }
  return(inverted.mat)
}


# Plot the amplitude or phase of an image on the frequency domain
fft_img_plot <- function(x, type = c("magnitude", "phase"), 
                         trans, 
                         shift = TRUE, plot = TRUE, as.cimg = FALSE){
  
  type <- match.arg(type)
  
  if(missing(trans)){
    trans <- switch(type, 
                    "magnitude" = function(x){log(x)}, 
                    "phase" = function(x){x})
  }
  
  if(shift){
    x <- fftshift(x)
  }
  
  if(is.null(trans) || (!is.function(trans) && is.na(trans))){
    trans <- function(x){
      x
    }
  } else {
    trans <- match.fun(trans)
  }
  
  foo <- switch(type, 
                "magnitude" = function(x){
                  Mod(x)
                }, 
                "phase" = function(x){
                  Arg(x)
                })
  
  img.mat <- apply(x, 1,function(z){
    trans(foo(z))
  })
  
  if(plot){
    plot(as.cimg(img.mat))
  }
  
  if(as.cimg){
    img.mat <- as.cimg(img.mat)
  }
  
  invisible(img.mat)
}



# Simulate image by drawing nfreq frequencies 
sim_img <- function(dim = c(150,150), 
                    nfreq = 10, 
                    zr = 1 / (rgamma(nfreq, 2.7699452,0.5074743) / max(dim)),
                    ztheta = runif(nfreq, 0, 2 * pi), 
                    alpha = 1.3, 
                    normalize = TRUE){
  x <- seq(0, 1, length.out = dim[1])
  y <- seq(0, 1, length.out = dim[2])
  zx <- zr * cos(ztheta)
  zy <- zr * sin(ztheta)
  za <- vapply(1/zr, function(x){(runif(1,0,x)^(alpha/2))}, FUN.VALUE = numeric(1))
  za <- za/sum(za)
  A <- tcrossprod(zx,x)
  B <- tcrossprod(zy,y)
  m <- crossprod((za * sin(A)), (za * cos(B))) + crossprod((za * cos(A)), (za * sin(B)))
  
  m <- m - min(m) 
  if(normalize){
    m <- m / max(m)
  }
  return(as.cimg(m))
}




# Simulate image using spectral synthesis method
syn_spec <- function(n = 100, beta = 1, plot = TRUE, threshold = TRUE, invert = TRUE,...){
  stopifnot(n%%2==0)
  epsilon <- rnorm((n^2))
  expand.grid("x" = seq.int(n)-floor(n/2), "y" = seq.int(n)-floor(n/2)) %>% 
    mutate("r" = sqrt((x)^2 + (y)^2), "epsilon" = epsilon) %>% 
    mutate(A = 1/(r)^(beta/2) * epsilon, phi = runif(n^2, 0, 2*pi)) %>% 
    mutate(A = ifelse(r == 0, 0, A)) %>% 
    mutate(z = A * exp((0+1i) * phi)) %>% 
    with(., matrix(data = z, ncol = n, nrow = n)) %>% 
    fftshift(inverse = TRUE) -> m
  
  if(invert){
    m <- fft_img_inv(m)
    if(threshold){
      
      m <- as.cimg(m[] + rnorm(prod(dim(m)), 0, 0.00001))
      
      m <- m %>% imeval(~. > median(c(.)))
    }
  }
  
  if(plot){
    if(invert){
      m %>% as.cimg() %>% plot(interpolate = FALSE)
    } else {
      fft_img_plot(m,...)
    }
    
  }
  invisible(m)
}


# Computes center of mass (mean) under periodic boundary conditions
com_periodic <- function(x, y){
  foo <- function(z,n){
    if(missing(n)){
      n <- length(z) 
    }
    theta <- z / n * 2 * pi
    thetamu<- atan2(-mean(sin(theta)), -mean(cos(theta))) + pi
    mu <- n * thetamu / (2 * pi) 
    return(mu)
  }
  if(missing(y)){
    return(foo(x))
  } else {
    return(c(foo(x, n = sqrt(length(x))) , foo(y, n = sqrt(length(y)))))
  }
}

# Computes the variance under periodic boundary conditions
var_periodic <- function(x,y){
  mu <- com_ref(x, y)
  n<-sqrt(length(x))
  var_x <- mean(((mu[1] - x) %% n)^2)
  var_y <- mean(((mu[2] - y) %% n)^2)
  return(c(sqrt(var_x^2 + var_y^2)))
}


