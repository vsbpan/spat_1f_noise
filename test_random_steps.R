iterate_random_steps2 <- function(issf_fit, 
                                 start = make_start2(),
                                 n = 100,
                                 ref_grid = dummy_spec(),
                                 rss_coef = 0){
  
  out.list <- vector(mode = "list", length = n+1)
  out.list[[1]] <- start
  stopifnot(nrow(start) == 1)
  
  
  ref_grid_flat <- c(ref_grid)
  ref_grid_flat <- ifelse(ref_grid_flat == 1, 1, 1 * exp(rss_coef))
  
  
  ra <- rdist(issf_fit$ta_updated, 10^5)
  rr <- rdist(issf_fit$sl_updated, 10^5)
  
  for(i in seq_len(n)){
    out.list[[i+1]] <- add_random_steps2(
      n = 100L, 
      x_start = out.list[[i]]$x,
      y_start = out.list[[i]]$y,
      direction_start = out.list[[i]]$theta,
      ta_rand = ra,
      sl_rand = rr,
      ref_grid_flat = ref_grid_flat)
  }
  do.call("rbind.fill", out.list)
}

add_random_steps2 <- function(
    n = 20L,
    x_start = NULL,
    y_start = NULL,
    direction_start = NULL,
    sl_rand = NULL,
    ta_rand = NULL,
    ref_grid_flat = NULL){
  
  
  index <- 1
  
  loop <- TRUE
  
  while(loop){
    index_r <- sample.int(n = length(sl_rand[[index]]), size = n)
    index_theta <- sample.int(n = length(ta_rand[[index]]), size = n)
    
    r <- sl_rand[[index]][index_r]
    theta <- ta_rand[[index]][index_theta]
    
    theta_new <- (direction_start + theta) %% (2 * pi)
    
    x_new <- x_start + r * cos(theta_new)
    y_new <- y_start + r * sin(theta_new)
    
    out_of_bound <- (x_new < 0 | x_new > 1000) | (y_new < 0 | y_new > 1000)
    
    x_new <- x_new[!out_of_bound]
    y_new <- y_new[!out_of_bound]
    theta_new <- theta_new[!out_of_bound]
    
    if(length(y_new) > 0){ # Exit condition
      loop <- FALSE
    }
  }
  
  i_new <- (ceiling(y_new / 1000 * 12) + floor(x_new / 1000 * 12) * 12)
  val_new <- ref_grid_flat[i_new]
  i <- sample.int(length(val_new), size = 1, replace = TRUE, prob = val_new)
  
    out <- data.frame(
      "theta" = theta_new[i],
      "x" = x_new[i],
      "y" = y_new[i]
    )
  
  return(out)
}

make_start2 <- function(theta = 0, x = 500, y = 500){
  data.frame("theta" = theta, "x" = x, "y" = y)
}

spec <- dummy_spec()
spec <- c(1:12 %% 2,
          (1:12+1) %% 2) %>% 
  rep(6) %>% 
  matrix(nrow = 12, ncol = 12) %>% 
  as.cimg()



iterate_random_steps2(issf_fit, 
                      start = make_start2(0,500,500), 
                      n = 5000, 
                      ref_grid = spec, 
                      rss_coef = 0.5) %>% 
  mutate(head_x = x, head_y = y) %>% 
  plot_track_overlay(
    repID = "Foo", 
    colored_track = FALSE, 
    trt_spec = spec, 
    plot_elements = "track"
  )

