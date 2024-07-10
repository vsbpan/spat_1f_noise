source("spat1f/init.R")

# Define simulation function
sim_fun <- function(sample_period = 1, rss = 0, k = 1){
  # Generate random spectrum under white noise
  spec <- syn_spec(n = 12, beta = 0, plot = FALSE)
  
  # Simulate fake movement
  d <- iterate_random_steps_states(
    ta_sl_list = list(
      "sl" = list(
        make_gamma(1, 90),    # More toxic
        make_gamma(1, 90 / k) # Less toxic
      ),
      "ta" = list(
        make_unif(), # More toxic
        make_unif() # Less toxic
      )
    ), 
    transition_mat = diag(1),
    start = make_start2(0,500,500, 1), 
    n = 2000 * sample_period, 
    ref_grid = spec, 
    same_move = FALSE, 
    rss_coef = rss 
  )
  
  # Drop out
  d <- d[seq_len(nrow(d)) %% sample_period == 0,]
  
  # Add random steps and append estimators
  d <- d %$% 
    move_seq(x, y) %>% 
    add_random_steps(n = 100, 
                     sl_distr = make_gamma(1, 90), 
                     ta_distr = make_unif()) %>% 
    flag_invalid_steps(remove_invalid = TRUE) %>% 
    mutate(
      less_toxic_start = 1 - (round(read_value(x1, y1, 
                                               ref_img = spec, transform = FALSE))),
      less_toxic_end = 1 - (round(read_value(x2, y2, 
                                             ref_img = spec, transform = FALSE)))
    ) %>% 
    append_estimators() 
  
  # Fit issf
  m <- issf(
    case ~ 
      less_toxic_end + 
      as.factor(less_toxic_start):sl +
      as.factor(less_toxic_start):logsl +
      strata(step_id), 
    data = d, 
    keep_data = FALSE, 
    update = FALSE
  ) 
  
  # Compute estimated params
  z <- unname(coef(m))
  k <- exp(-diff(log(unname(1/((1/m$sl$params$scale) - z[2:3])))))
  
  # Return result
  return(c("rss" = z[1], "k" = k))
}



# Define simulation grid
sim_d <- data.frame(
  "rss" = c(0, log(3), 0, log(3)),
  "k" = c(3, 1, 3, 1),
  "sample_period" = c(1, 1, 3, 3)
) %>% 
  rep_data.frame(200)

# Simulate
out <- pb_par_lapply(seq_len(nrow(sim_d)), 
                     FUN = function(i, sim_fun, sim_d){
                       sample_period <- sim_d[i, "sample_period"]
                       rss <- sim_d[i, "rss"]
                       k <- sim_d[i, "k"]
                       
                       res <- sim_fun(sample_period = sample_period,
                                      rss = rss, 
                                      k = k)
                       
                       return(c(
                         res, 
                         "period" = sample_period,
                         "rss_true" = rss, 
                         "k_true" = k
                       ))
                     }, 
                     sim_fun = sim_fun,
                     sim_d = sim_d,
                     cores = 8, 
                     inorder = FALSE)

res <- do.call("rbind", out) %>% 
  as.data.frame() %>% 
  group_by(
    rss_true, k_true, period
  ) %>% 
  mutate(
    scenario_id = cur_group_id()
  )

# write_csv(res, "simulation/validate_issf_estimation_sim.csv")




  

