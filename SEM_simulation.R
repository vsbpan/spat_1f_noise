source("helper_functions/init.R")
# sem_fit_boot <- bootSEM(sem_fit, nboot = 200, cores = 8)
# saveRDS(sem_fit_boot, "invisible/sem_fit_boot.rds")
sem_fit_boot <- readRDS("invisible/sem_fit_boot.rds")

og_set <- find_significant_pairs(sem_summary$coefficients, p = 0.05)
by_data <- expand.grid("var" = c("beta_red", "beta_white", "var_high"), 
                       "only" = c("sl_mean_obs_log_scale", "area_herb_log_scale", 
                                  "mean_toxic_conc_scale", "var_toxic_12_scale", NA), 
                       "cat_size" = c(-2, 2), 
                       "exclude" = c(NA, "ava_mean_toxin_scale"))

names(sem_fit$data)

out_list <- list()
for (j in seq_along(sem_fit_boot)){
  cat(sprintf("\nProcessing boot %s\n", j))
  v <- pb_par_lapply(seq_len(nrow(by_data)), 
                     function(i, by_data, og_set, sem_fit_booti){
                       SEM_pred_coef(sem_fit_booti, 
                                     var = by_data[i,"var"],
                                     target = "RGR_scale", 
                                     cat_size = by_data[i,"cat_size"],
                                     exclude = by_data[i,"exclude"],
                                     only = by_data[i,"only"],
                                     og_set = og_set) 
                     }, cores = 1, 
                     by_data = by_data, 
                     og_set = og_set, 
                     sem_fit_booti = sem_fit_boot[[j]]) %>% 
    do.call("c", .)
  out_list[[j]] <- data.frame(by_data, "val" = v)
}

out_list_d <- out_list %>% 
  do.call("rbind", .)
#write_csv(out_list_d, "cleaned_data/SEM_sim.csv")



