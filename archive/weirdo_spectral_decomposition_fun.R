ma <- function(x, n = 5){
  stats::filter(x, rep(1 / n, n), sides = 2)
}

thin_rows <- function(x, n = 10L){
  n <- round(n)
  x[seq.int(nrow(x)/n) * n,]
}


min_peaks <- function(x, thresh, p.check = FALSE){
  if(is.null(x) || any(is.na(x))){
    return(NA_real_)
  }
  x <- x[order(x[,1]),]
  mmax <- amax2(x[,2],thresh = 0)
  mmin <- amax2(-x[,2],thresh = 0)
  x <- x[sort(c(mmin,mmax,nrow(x),1)),]
  z <- amax2(-x[,2], thresh = thresh)
  out <- x[z,]
  
  if(!isFALSE(p.check)){
    out <- out[out[,2] < p.check,]
  }
  return(out[,1])
}

amax2<- function(x, thresh = 0){
  a1 <- c(0, x - thresh, 0)
  a2 <- c(x, 0, 0)
  a3 <- c(0, 0, x)
  e <- which((a1 >= a2 & a1 > a3)[2:(length(x))])
  if (!is.na(e[1] == 1)) 
    if (e[1] == 1) 
      e <- e[-1]
  if (length(e) == 0) 
    e <- NaN
  return(e)
}

get_peaks <- function(img, FAP = 0.001, clean = TRUE){
  
  if(is.cimg(img)){
    mat <-img[,,1,1] 
  } else {
    mat <- img
  }
  
  min_freq <- 1 / min(dim(mat)/2)
  FT <- spec2D(mat, FAP = !clean)
  a <- FT %>% 
    mutate(fr = sqrt(fx^2 + fy^2),
           ftheta = atan(fy/fx) * 180 / pi)
  
  s <- attr(img, "px.size")$size
  if(clean){
    a <- a %>% 
      filter(fr > min_freq) %>% 
      filter(fr > 0 & ftheta > 0) %>% 
      mutate(period = 1/fr * s) %>%
      mutate(fr = round(fr/min(fr), 2)) %>% 
      group_by(fr) %>% 
      summarise(period = mean(period),
                PSD = sum(PSD), # total power
                nr = n()) %>% 
      ungroup() %>% 
      mutate(p = get_FAP(PSD, FAP = FAP)) %>% 
      dplyr::select(-fr)
  }
  return(as.data.frame(a))
}


smooth_spec <- function(x, sp = 0.1, k = 100, m = NA, lower_cut_off = 0, 
                        plot = FALSE, thresh = 1, p.check = 0.001, 
                        method = c("ma", "gam","none"), window = 1, by = FALSE, 
                        return.PSD = FALSE){
  #warning("Not a good idea. Don't use me!")
  method <- match.arg(method)
  l <- nrow(x)
  x <- x %>% 
    mutate(log.p = log10(p),
           log.period = log10(period)) %>% 
    filter(log.period > lower_cut_off) %>% 
    arrange(log.period)
  x <- x %>%
    mutate(log.p = ifelse(is.infinite(log.p) & sign(log.p) == -1,
                          min(x$log.p[is.finite(x$log.p)]) - 10,
                          ifelse(is.infinite(log.p) & sign(log.p) == 1,
                                 1,
                                 log.p)))
  if(is.na(k) | k >= nrow(x)){
    k <- nrow(x) - 1
  }
  
  if(method == "ma"){
    if(isFALSE(by)){
      smoothed.p <- data.frame(
        "log10period" = x$log.period, 
        "PSD" = ma(x$PSD, n = window)) %>% 
        filter(!is.na(PSD)) %>% 
        mutate(log10p.pred = log10(get_FAP(PSD, N = l)))
    } else {
      if(is.na(by)){
        by <- min(abs(diff(x$period)), na.rm = TRUE)
      }
      zx <- seq(min(x$period), max(x$period),
                by = by)
      zy <- approxfun(x = x$period,
                      y= x$log.p,
                      method = "linear")(zx)
      dz <- data.frame("log.p" = zy,
                       "log.period" = log10(zx))
      smoothed.p <- data.frame(
        "log10period" = dz$log.period, 
        "log10p.pred" = ma(dz$log.p, window)) %>% 
        filter(!is.na(log10p.pred))
      message(paste0("Smoothing with 'by' = ",by))
    }
  } else if(method == "gam"){
    m <- mgcv::gam(log.p ~ 
                     s(log.period, 
                       bs = "tp", 
                       k = k, 
                       fx = FALSE, 
                       m = m,
                       sp = sp), 
                   data = x)
    
    new.x <- seq(min(x$log.period), max(x$log.period), length.out = 10000)
    smoothed.p <- data.frame(
      "log10period" = new.x, 
      "log10p.pred" = mgcv::predict.gam(m, newdata = data.frame("log.period" = new.x))) 
  } else if(method == "none"){
    smoothed.p <- data.frame(
      "log10period" = x$log.period, 
      "log10p.pred" = log10(x$p))
  }
  
  if(plot){
    peaks <- na.omit(min_peaks(smoothed.p %>% 
                                 dplyr::select(log10period, log10p.pred), 
                               thresh = thresh, p.check = p.check))
    g <- x %>% 
      ggplot(aes(x = log10(period), y = log10(p))) + 
      geom_line() + 
      geom_point() + 
      geom_line(data = smoothed.p, 
                aes(x = log10period, y = (log10p.pred)), 
                color = "red") + 
      geom_vline(data = data.frame("peaks" =  peaks), 
                 aes(xintercept = peaks), 
                 color = "blue", 
                 linetype = 3, 
                 size = 1) + 
      theme_bw(base_size = 15) + 
      labs(subtitle = paste0(length(peaks), " peaks detected"))
    print(g)
  }
  if(!return.PSD){
    smoothed.p <- smoothed.p %>% dplyr::select(-PSD)
  }
  invisible(smoothed.p)
}



spec2D <- function(x, FAP = 0.001){
  if(!is.matrix(x)){
    x <- as.matrix(x)
  }
  if(var(c(x)) == 0){
    stop("matrix contrains zero variance")
  }
  
  x.fft <- fft(x)
  x.fft <- fftshift(x.fft)
  nx <- nrow(x)
  ny <- ncol(x)
  n <- nx * ny
  center <- ceiling(c(nx, ny) / 2)
  fy <- (seq.int(ny) - center[2]) / ny
  fx <- (seq.int(nx) - center[1]) / nx
  amplitude <- c(Mod(x.fft)) / n
  out <- cbind(expand.grid("fx"=fx, "fy" = fy) ,amplitude)
  out$PSD <- out$amplitude^2
  if(!isFALSE(FAP)){
    out$p <- get_FAP(out$PSD, FAP = FAP) 
  }
  #out <- out[!(out$fx == 0 & out$fy == 0), ]
  return(out)
}


get_FAP <- function(x, N = NA, method = c("Lomb-Scargle","Chisquare"), FAP = 0.001){
  method <- match.arg(method)
  x <- as.vector(x)
  if(sum(x) <=0){
    return(rep(NA_real_, length(x)))
  }
  
  if(is.na(N)){
    N <- length(x)
  }
  
  x <- x / sum(x)
  
  if(method == "Lomb-Scargle"){
    M <- -6.362 + 1.193 * N + 0.00098 * N^2
    tmp <- (1 - x)^((N - 3)/2)
    p <- M * tmp
    if (any(p > FAP)){
      p[p > FAP] <- 1 - (1 - tmp[p > FAP])^M
    }
  } else if(method == "Chisquare"){
    # For spatial point process
    p <- pchisq(x * 2, df = 2, lower.tail = FALSE)
  }
  return(p)
}

get_pseudo_peaks <- function(x1, x2 = NULL, smoothed = FALSE, window = 5){
  if(is.null(x2) && is.list(x1)){
    x2 <- x1[[2]]
    x1 <- x1[[1]]
  }
  
  z <- get_peaks(x1) %>% 
    rename_all(.funs = function(x) {
      ifelse(x != "period",  paste0(x,"1"), x) 
    }) %>% 
    left_join(get_peaks(x2) %>% 
                rename_all(.funs = function(x){
                  ifelse(x != "period",  paste0(x,"2"), x) 
                }), by = "period") 
  
  if(smoothed){
    z <- z %>% 
      mutate(PSD1 = ma(PSD1, n = window), 
             PSD2 = ma(PSD2, n = window)) %>% 
      filter(!is.na(PSD1) & !is.na(PSD2))
  }
  
  z <- z %>% 
    mutate(diff = PSD1 - PSD2) %>% 
    mutate(epsilon = min(diff[diff>0])/2) %>%
    mutate(epsilon = ifelse(epsilon < 0 | !is.finite(epsilon), 0,epsilon)) %>% 
    mutate(diff = ifelse(diff < 0, epsilon, diff))
  z$p_diff <- get_FAP(z$diff)
  return(z)
}

get_null_img <- function(x, px_size, thr = "99.9"){
  img1 <- x %>% crop_leaf(invalid = 0, thr = thr) %>% 
    add_px_size(px_size)
  img2 <- x %>% 
    crop_leaf(invalid = NA, thr = thr) %>% 
    imeval(~!is.na(.)) %>% 
    as.cimg() %>% 
    add_px_size(px_size)
  return(list(img1, img2))
}

estimate_alpha <- function(data, x = "period", y= "PSD", method = c("OLS", "MA","SMA"), 
                           plot = FALSE, thresh = 0.1){
  data <- as.data.frame(data)
  data$x <- data[,x]
  data$y <- data[,y]
  method <- match.arg(method)
  data <- data %>% filter(x>thresh)
  z <- suppressMessages({lmodel2::lmodel2(log(y) ~ log(x), 
                                          data = data)})
  if(plot){
    plot(z, method = method)
  }
  
  r <- switch(method, 
              "MA" = 2, 
              "OLS" = 1, 
              "SMA" = 3)
  
  out <- c("estimate" = z$regression.results[r,3], 
           "lower" = z$confidence.intervals[r, 4], 
           "upper" = z$confidence.intervals[r, 5])
  return(out)
}


plot_spectra <- function(data, x = "period", y= "PSD"){
  data <- as.data.frame(data)
  data$x <- data[,x]
  data$y <- data[,y]
  
  g <- data %>% 
    ggplot(aes(x = x, y = y )) + 
    geom_line() + 
    labs(y = "PSD", x = "Period") + 
    theme_bw(base_size = 15) + 
    scale_x_continuous(trans = "log10") + 
    scale_y_continuous(trans = "log10")
  
  return(g)
}


acfFFT <- function(data, y= "PSD"){
  data <- as.data.frame(data)
  data$y <- data[,y]
  z <- fft(data$y, inverse = TRUE)
  return(z)
}
