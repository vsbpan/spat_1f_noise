source("spat1f/init.R")

for (j in 1:10){
  sim_d <- expand.grid(
    "b" = c(5, 0, -5),
    "rss" = c(log(0.75), 0, log(1.333333)), 
    "scale" = c(8.3333, 25, 83.333),
    "k" = c(0.2, 1, 5)
  ) %>% 
    rep_data.frame(40)
  
  out <- pb_par_lapply(
    seq_len(nrow(sim_d)), function(i, sim_d){
      issf_fit2 <- list(
        "ta_updated" = list(
          make_genvonmises(kappa1 = 0.4, kappa2 = 0.3),
          make_genvonmises(kappa1 = 0.4, kappa2 = 0.3)
        ),
        "sl_updated" = list(
          make_gamma(shape = 0.8, scale = sim_d[i, "scale"]),
          make_gamma(shape = 0.8, scale = sim_d[i, "scale"] / sim_d[i, "k"]) 
        )
      )
      
      spec <- as.cimg(syn_spec(12, sim_d[i,"b"], plot = FALSE))
      x <- iterate_random_steps2(issf_fit2, 
                                 start = make_start2(0,500,500, 1), 
                                 n = 1000, 
                                 ref_grid = spec, 
                                 rss_coef = sim_d[i,"rss"])
      x <- x[, c("on_toxic", "r", "ava_toxic")]
      
      data.frame(
        sim_d[i, ,drop = FALSE],
        as.data.frame(t(colMeans(x, na.rm = TRUE)))
      )
    },
    sim_d = sim_d,
    cores = 8, 
    inorder = FALSE
  )
  
  write_csv(do.call("rbind", out), file = sprintf("simulation/check_points/check_point_%s.csv",j))
}



# out <- list.files("simulation/check_points", full.names = TRUE) %>% 
#   lapply(function(x){
#     suppressMessages(read_csv(x, progress = FALSE))
#   }) %>% 
#   do.call("rbind", .)

# write_csv(out, "simulation/move_rules_sim.csv")

# file.remove(list.files("simulation/check_points", full.names = TRUE))



out