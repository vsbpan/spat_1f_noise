iterate_random_steps2 <- function(issf_fit, 
                                  start = make_start2(),
                                  n = 100,
                                  ref_grid = dummy_spec(),
                                  rss_coef = 0){
  stopifnot(nrow(start) == 1)
  
  
  ref_grid_flat <- c(ref_grid)
  
  ra <- rdist(issf_fit$ta_updated, 10^5)
  rr <- rdist(issf_fit$sl_updated, 10^5)
  
  sim <- add_random_steps_iterateC(
    n = n, 
    n_draws = 100L, 
    x_start = start$x,
    y_start = start$y, 
    direction_start = start$theta, 
    diet_start = start$on_toxic, 
    sl_rand = ra, 
    ta_rand = rr, 
    rss_coef = rss_coef, 
    ref_grid_flat = ref_grid_flat
  )
  
  return(rbind.fill(start, sim))
  
  # for(i in seq_len(n)){
  #   out.list[[i+1]] <- add_random_steps2(
  #     n = 100L, 
  #     x_start = out.list[[i]]$x,
  #     y_start = out.list[[i]]$y,
  #     direction_start = out.list[[i]]$theta,
  #     diet_start = out.list[[i]]$on_toxic,
  #     ta_rand = ra,
  #     sl_rand = rr,
  #     rss_coef = rss_coef,
  #     ref_grid_flat = ref_grid_flat)
  # }
  # do.call("rbind", out.list[-1])
}

add_random_steps2 <- function(
    n = 20L,
    x_start = NULL,
    y_start = NULL,
    direction_start = NULL,
    diet_start = 1,
    sl_rand = NULL,
    ta_rand = NULL,
    rss_coef = NULL,
    ref_grid_flat = NULL){
  
  
  index <- ifelse(diet_start == 1, 1, 2)
  
  # picked_list <- pick_new_theta_xy(
  #   sl_rand, 
  #   ta_rand, 
  #   index = index, 
  #   n = n, 
  #   direction_start = direction_start, 
  #   x_start = x_start, 
  #   y_start = y_start)
  # 
  # 
  # val_new <- ifelse(ref_grid_flat == 1, 1, 1 * exp(rss_coef))[picked_list$i_new]
  # i <- sample.int(length(val_new), size = 1, replace = TRUE, prob = val_new)
  # 
  # out <- data.frame(
  #   "theta" = picked_list$theta_new[i],
  #   "x" = picked_list$x_new[i],
  #   "y" = picked_list$y_new[i],
  #   "r" = picked_list$r_new[i],
  #   "on_toxic" = ref_grid_flat[picked_list$i_new[i]],
  #   "ava_toxic" = mean(ref_grid_flat[picked_list$i_new])
  # )
  # 
  out <- add_random_stepsC(n = n, 
                           x_start = x_start, 
                           y_start = y_start, 
                           direction_start = direction_start, 
                           sl_rand = sl_rand, 
                           ta_rand = ta_rand, 
                           index = index, 
                           rss_coef_exp = exp(rss_coef), 
                           ref_grid_flat = ref_grid_flat)
  return(out)
}

make_start2 <- function(theta = 0, x = 500, y = 500, on_toxic = 1){
  data.frame("theta" = theta, "x" = x, "y" = y, "on_toxic" = on_toxic)
}
