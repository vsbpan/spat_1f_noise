source("spat1f/init.R")

rm(list = ls()[grepl("subm",ls())])


#### Test mediator hypotheses ####

og_set <- find_significant_pairs(sem_summary$coefficients, p = 0.05)
by_data <- expand.grid("var" = c("beta_numeric_scale", "var_high"),
                       "target" = c("RGR_scale"),
                       "only" = c("sl_mean_obs_log_scale", "area_herb_log_scale", 
                                  "mean_toxic_conc_scale", "var_toxic_12_scale", 
                                  "select_scale"), 
                       "cat_size" = c(-2, 2), 
                       "exclude" = c(NA)) %>% 
  rbind(
    expand.grid("var" = c("beta_numeric_scale", "var_high"),
                "target" = c("sl_mean_obs_log_scale", "area_herb_log_scale", 
                             "mean_toxic_conc_scale", "var_toxic_12_scale", 
                             "select_scale"),
                "only" = NA, 
                "cat_size" = c(-2, 2), 
                "exclude" = c(NA)),
    expand.grid("var" = c("sl_mean_obs_log_scale", "area_herb_log_scale", 
                          "mean_toxic_conc_scale", "var_toxic_12_scale", 
                          "select_scale"),
                "target" = "RGR_scale",
                "only" = NA, 
                "cat_size" = c(-2, 2), 
                "exclude" = c(NA))
  )

out_list <- pb_par_lapply(
  1:10,
  function(j, sem_fit, by_data, og_set){
    sem_fit_booti <- bootSEM(sem_fit, nboot = 1, cores = 1, silent = TRUE)[[1]]
    v <- lapply(seq_len(nrow(by_data)), function(i){
      SEM_pred_coef(sem_fit_booti, 
                    var = by_data[i,"var"],
                    target = by_data[i,"target"], 
                    cat_size = by_data[i,"cat_size"],
                    exclude = by_data[i,"exclude"],
                    only = by_data[i,"only"],
                    og_set = og_set)
    }) %>% 
      do.call("c", .)
    data.frame(by_data, "val" = v)
  },
  sem_fit = sem_fit,
  by_data = by_data,
  og_set = og_set,
  cores = 6
)

out_list_d <- out_list %>% 
  do.call("rbind", .)

# write_csv(out_list_d, "cleaned_data/SEM_sim_hypotheses.csv")


#### Figure S? ####
og_set <- find_significant_pairs(sem_summary$coefficients, p = 0.05)
by_data <- expand.grid("var" = c("beta_numeric_scale", "var_high"),
                       "only" = c("sl_mean_obs_log_scale", "area_herb_log_scale",
                                  "mean_toxic_conc_scale", "var_toxic_12_scale", NA),
                       "cat_size" = c(-2, 2),
                       "exclude" = c(NA, "ava_qual"))

out_list <- pb_par_lapply(
  1:500,
  function(j, sem_fit, by_data, og_set){
    sem_fit_booti <- bootSEM(sem_fit, nboot = 1, cores = 1, silent = TRUE)[[1]]
    v <- lapply(seq_len(nrow(by_data)), function(i){
      SEM_pred_coef(sem_fit_booti,
                    var = by_data[i,"var"],
                    target = "RGR_scale",
                    cat_size = by_data[i,"cat_size"],
                    exclude = by_data[i,"exclude"],
                    only = by_data[i,"only"],
                    og_set = og_set)
    }) %>%
      do.call("c", .)
    data.frame(by_data, "val" = v)
  },
  sem_fit = sem_fit,
  by_data = by_data,
  og_set = og_set,
  cores = 6
)

out_list_d <- out_list %>%
  do.call("rbind", .)
# write_csv(out_list_d, "cleaned_data/SEM_sim_node_removal.csv")


