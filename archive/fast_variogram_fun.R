variogram <- function(X, summarise = TRUE, na.rm = TRUE, max_eval = 5000, res = 1, diag = FALSE, ...){
  # Algorithm n^2 complexity
  Z <- melt(X, drop = TRUE)
  if(na.rm){
    Z <- Z[!is.na(Z$val), ]
  }
  
  nr <- nrow(Z)
  if(nr > max_eval){
    Z <- Z[sample(seq_len(nr), size = max_eval, replace = FALSE),]
  }
  
  variogram_calc(Z[,(names(Z)!= "val")], Z$val, ...)
}


variogram_calc <- function(X, val, summarise = TRUE, res = 1, diag = FALSE){ 

  dist_vec <- Rfast::upper_tri(Rfast::Dist(x = X, method = "euclidean", square = FALSE), 
                               suma = FALSE, diag = diag)
  group_vec <- seq(min(dist_vec),max(dist_vec), res)
  
  dist_d <- cbind(
    "distance" = dist_vec, 
    "gamma" = Rfast::upper_tri(Rfast::vecdist(as.numeric(val))^2, suma = FALSE, diag = diag), 
    "group" = findInterval(dist_vec,group_vec)
  ) 
  #cat(length(dist_vec))
  if(summarise){
    gmean <- Rfast::group(dist_d[,"gamma"], dist_d[,"group"],method = "mean") 
    out <- data.frame("distance" = group_vec[Rfast::sort_unique(dist_d[,"group"])], 
                      "gamma_mean" = gmean)
    return(out)
  } else {
    return(dist_d)
  }
}


Moran_I <- function(X, type = c("global", "local"), q = 1, z_score = FALSE, max_eval = 5000){
  Z <- melt(X, drop = TRUE)
  
  nr <- nrow(Z)
  if(nr > max_eval){
    Z <- Z[sample(seq_len(nr), size = max_eval, replace = FALSE),]
  }
  
  n <- length(Z$val)
  w_mat <- Rfast::Dist(Z[,(names(Z)!= "val")])^-q
  diag(w_mat) <- NA
  type <- match.arg(type)
  
  
  if(type == "local"){
    b1 <- sum(w_mat^2, na.rm = TRUE)
    
    Z$val <- lapply(seq_along(Z$val), function(i){
      val <- Z$val[-i]
      w <- w_mat[,i][-i]
      
      I <- sum(w * (val - mean(val, na.rm = TRUE)), na.rm = TRUE) * 
        Rfast::Mad(val, method = "mean", na.rm = TRUE) / 
        Rfast::Var(val, na.rm = TRUE)
      
      if(z_score){
        E_I <- -sum(w, na.rm = TRUE) / (n-1)
        b2 <- Rfast::kurt(val[!is.na(val)])
        V_I <- (n - b2) / (n-1) * sum(w^2, na.rm = TRUE) - (2*b2 - n) / (n-1) / (n-2) * b1 - E_I^2
        z <- (I - E_I) / sqrt(V_I)
        return(z)
      } else {
        return(I)
      }
    }) %>% do.call("c",.)
    
    return(unmelt(Z)) 
    
  } else {
    if(type == "global"){
      W <- sum(w_mat, na.rm = TRUE)
      
      I <- sum((tcrossprod(Z$val - mean(Z$val, na.rm = TRUE)) * 
             w_mat), na.rm = TRUE) / W / (Rfast::Var(Z$val, na.rm = TRUE))
      if(z_score){
        E_I <- -1/(n - 1)
        
        # w_{ij} = w_{ji} for distance weights
        s1 <- sum(w_mat^2, na.rm = TRUE) * 2
        s2 <- sum(Rfast::colsums(w_mat, na.rm = TRUE)^2)
        s3 <- Rfast::kurt(Z$val[!is.na(Z$val)])
        s4 <- (n^2 + 3*n + 3) * s1 - n*s2 + 3*W^2
        s5 <- (n^2 - n) * s1 - 2*n*s2 + 6*W^2
        V_I <- (n * s4 - s3 * s5) / (n-1) / (n-2) / (n-3) / W^2 - E_I^2
        
        z <- (I - E_I) / sqrt(V_I)
        return(z)
      } else {
        return(I)
      }
    }
  }
}
