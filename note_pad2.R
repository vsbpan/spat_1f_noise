source("helper_functions/init.R")
source("helper_functions/init_analysis.R")

IDs <- fetch_repID()

# event_list <- IDs %>%
#   pb_par_lapply(
#     function(x){
#       out <- fetch_events(x)
#       l <- fetch_data_dict(x)
#       trt_spec <- fetch_trt_spec(x, .ref_data = get("ref_data", parent.frame()))
#       out$head_high_trt <- read_value(
#         x = out$head_x,
#         y = out$head_y,
#         dim_xy = get_dim(l),
#         ref_img = trt_spec)
#       out$head_high_trt <- read_value(
#         x = out$centroid_x,
#         y = out$centroid_y,
#         dim_xy = get_dim(l),
#         ref_img = trt_spec)
#       out <- spat1f::insert_gaps(out)
#       return(out)
#     }, cores = 8,
#     export_fun_only = FALSE
#   ) %>%
#   append_name(paste0("rep",IDs))
# saveRDS(event_list, "cleaned_data/events_list.rds")
event_list <- readRDS("cleaned_data/events_list.rds") %>% 
  lapply(function(x) left_join(x, ref_data, by = "repID"))

move_list <- event_list %>% 
  lapply(function(d){
    cbind(
      move_seq(x = d$head_x, y = d$head_y, 1000/12 * 3), # t to t + 1 r and theta
      d[-nrow(d),], # Use the meta data at t
      append_name(d[-1, "head_high_trt"], "head_high_trt_t1") # head_high_trt at t+1
      ) 
  })


fetch_events(50) %$% 
  move_seq(head_x, head_y)


move_list %>% 
  do.call("rbind",.) %>% 
  filter(r / 1000 * 12 > 1) %>% 
  ggplot(aes(x = log(r), y = abs(rad2degree(theta_rel)))) + 
  geom_point() + 
  geom_smooth()




move_list %>% 
  do.call("rbind", .) %>% 
  filter(r / 1000 * 12 > 1) %>% 
  mutate(
    beta = ifelse(is.na(beta), "constant", as.character(beta))
  ) %>% 
  filter(!is.na(theta_rel)) %>% 
  mutate(
       theta_rel = nearest_bin(rad2degree(theta_rel), seq(-180,180, by = 30))
  ) %>% 
  group_by(var_trt, theta_rel) %>% 
  summarise(
    count = n()
  ) %>% 
  group_by(var_trt) %>% 
  mutate(
    tot = sum(count),
    p = count / tot
  ) %>% 
  ggplot(aes(x = (theta_rel), y = p, fill = var_trt)) +
  geom_col(position = "dodge") + 
  #geom_errorbar(position = "dodge", stat = "summary") + 
  labs(x = "Turn angle (degrees)", y = "Probability", fill = "Var trt") + 
  theme_bw(base_size = 15) -> g2


move_list %>% 
  do.call("rbind", .) %>% 
  filter(r / 1000 * 12 > 1) %>% 
  mutate(
    beta = ifelse(is.na(beta), "constant", as.character(beta))
  ) %>% 
  filter(!is.na(theta_rel)) %>% 
  mutate(
    theta_rel = nearest_bin(rad2degree(theta_rel), seq(-180,180, by = 30))
  ) %>% 
  group_by(beta, theta_rel) %>% 
  summarise(
    count = n()
  ) %>% 
  group_by(beta) %>% 
  mutate(
    tot = sum(count),
    p = count / tot
  ) %>% 
  ggplot(aes(x = (theta_rel), y = p, fill = beta)) +
  geom_col(position = "dodge") + 
  #geom_errorbar(position = "dodge", stat = "summary") + 
  labs(x = "Turn angle (degrees)", y = "Probability", fill = expression(beta ("Autocor."))) + 
  theme_bw(base_size = 15) ->g1


ggarrange(g1, g2, nrow = 2)

move_list %>% 
  do.call("rbind", .)  %>% 
  mutate(
    beta = ifelse(is.na(beta), "constant", as.character(beta))
  ) %>% 
  mutate(lr = log(r)) %>% 
  filter(!is.na(lr)) %>% 
  mutate(
    lr = nearest_bin(lr, seq(0, 8, by = 1))
  ) %>% 
  group_by(beta, lr) %>% 
  summarise(
    count = n()
  ) %>% 
  group_by(beta) %>% 
  mutate(
    tot = sum(count),
    p = count / tot
  ) %>% 
  ggplot(aes(x = exp(lr) / 1000 * 120, y = p, fill = beta)) +
  geom_col(position = "dodge") + 
  scale_x_continuous(trans = "log10") + 
  labs(x = "Step length (mm)", y = "Probability", fill = expression(beta)) + 
  theme_bw(base_size = 15) ->g2

ggarrange(g1, g2, common.legend = TRUE, legend = "top")



d_move2 <- move_list %>% 
  do.call("rbind", .) %>% 
  filter(
    !is_gap
  ) %>% 
  mutate(
    tz = ifelse(step_id < 24 * 6, "0h - 24h", "24h - 120h"),
    r = r 
  )  %>% 
  group_by(
    rep_id, tz
  ) %>% 
  summarise(
    mean_lr = mean(log(r), na.rm = TRUE),
    var_lr = var((r), na.rm = TRUE),
    q10_lr = quantile(na.omit(log(r)), prob = 0.10),
    q90_lr = quantile(na.omit(log(r)), prob = 0.90),
    mean_theta_abs = mean(abs(theta_rel), na.rm = TRUE), 
    var_theta_abs  = var(abs(theta_rel), na.rm = TRUE)
  ) %>% 
  left_join(
    ref_data %>% mutate(rep_id = as.character(rep_id)), by = "rep_id"
  ) %>% 
  filter(
    !is.na(rep_id)
  )

d_move2 %>% 
  ggplot(aes(x = factor(var_trt), y = mean_lr)) + 
  geom_point(position = position_jitter(width = 0.2)) + 
  geom_pointrange(stat = "summary")


d_move2 %>% 
  filter(mean_trt == "1 mg/g") %>% 
  ggplot(aes(x = factor(beta), y = exp(q90_lr) / 1000 * 120, color = var_trt)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6)) + 
  geom_pointrange(stat = "summary", position = position_dodge(0.6)) + 
  facet_wrap(~tz) + 
  scale_y_continuous(trans = "log10") + 
  labs(x = expression(beta), y = "Q90 step length (mm)", color = "Var trt") + 
  theme_bw(base_size = 15) + 
  theme(legend.position = "top")

d_move2 %>% 
  filter(mean_trt == "1 mg/g") %>% 
  ggplot(aes(x = factor(beta), y = mean_theta_abs, color = var_trt)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6)) + 
  geom_pointrange(stat = "summary", position = position_dodge(0.6))+ 
  facet_wrap(~tz) + 
  labs(x = expression(beta), y = "Mean turn angle bias (radians)", color = "Var trt") + 
  theme_bw(base_size = 15) + 
  theme(legend.position = "top")


d_move2 %>% 
  filter(mean_trt == "1 mg/g") %>% 
  ggplot(aes(x = (cat_pre_wt), y = mean_theta_abs, color = factor(var_trt))) + 
  geom_point() + 
  geom_smooth(method = "lm")+ 
  facet_wrap(~tz) + 
  scale_x_continuous(trans = "log10") + 
  labs(x = "Cat pre-weight (g)", y = "Mean turn angle bias (radians)", color = "Var trt") + 
  theme_bw(base_size = 15) + 
  theme(legend.position = "top") -> g1


d_move2 %>% 
  filter(mean_trt == "1 mg/g") %>% 
  ggplot(aes(x = (cat_pre_wt), y = q10_lr, color = factor(beta))) + 
  geom_point() + 
  geom_smooth(method = "lm")+ 
  facet_wrap(~tz) + 
  scale_x_continuous(trans = "log10") + 
  labs(x = "Cat pre-weight (g)", y = "Q10 step length (mm)", color = "Var trt") + 
  theme_bw(base_size = 15) + 
  theme(legend.position = "top") -> g2;g2


ggpubr::ggarrange(g1, g2, common.legend = TRUE, nrow = 2, legend = "top")


d2 <- move_list %>% 
  do.call("rbind", .) %>% 
  filter(
    !is_gap
  ) %>% 
  mutate(
    tz = ifelse(step_id < 24 * 6, "0h - 24h", "24h - 120h")
  ) %>% 
  filter(
    var_trt != "constant"
  )


names(d2)



m <- glmmTMB(
  r ~ 
    tz + var_trt * log(cat_pre_wt) + 
    factor(head_high_trt) * beta + 
    (1|session_id) + (1|rep_id),  
  family = Gamma(link = "log"),
  data = d2 %>% 
    mutate(
      turn_bias = abs(theta_rel)
    ) %>% 
    filter(
      cat_dead_cam_end == 0 & pupated_cam_end == 0
    )
); summary(m)

m <- glmmTMB(
  turn_bias ~ 
    tz + beta + var_trt + log(cat_pre_wt) + 
    factor(head_high_trt) + 
    (1|session_id) + (1|rep_id),  
  family = gaussian(),
  data = d2 %>% 
    mutate(
      turn_bias = abs(theta_rel)
    ) %>% 
    filter(
      cat_dead_cam_end == 0 & pupated_cam_end == 0
    )
); summary(m)

DHARMa::simulateResiduals(m) %>% plot()

plot_model(m, 
           type = "pred", 
           terms = c("cat_pre_wt[all]", "tz")) +
  theme_bw(base_size = 15) + 
  labs(title = NULL, 
       subtitle = NULL, 
       y = "Turn bias", x = "Cat pre weight", color = "Time") + 
  theme(legend.position = "top")















event_list %>% 
  do.call("rbind",.) %>% 
  .$rank




moving_var <- function(x, n){
  if(length(x) < n){
    return(NA)
  }
  append(rep(NA, n),
         vapply(seq_len(length(x) - n), FUN = function(i){
           var(x[i:(i+n)])
         }, FUN.VALUE = numeric(1)))
}

d_event2 <- event_list %>% 
  do.call("rbind",.) %>% 
  mutate(
    tz = ifelse(rank < (24 * 10), "0h - 24h", "24h - 120h")
  ) %>% 
  group_by(rep_id, tz) %>% 
  summarise(
    tot_sd = sqrt(var(head_x, na.rm = TRUE) + var(head_y, na.rm = TRUE)),
    mean_quality = mean(head_high_trt, na.rm = TRUE), 
    var_quality = mean(moving_var(head_high_trt, n = 120), na.rm = TRUE)
  )%>% 
  left_join(
    ref_data %>% mutate(rep_id = as.character(rep_id)), by = "rep_id"
  ) %>% 
    filter(!is.na(var_trt))

parse_conc <- function(x) {
  as.numeric(gsub(" mg/g","", x))
}


d_event2 <- d_event2 %>% 
  filter(mean_trt == "1 mg/g") %>% 
  mutate(
    mean_quality = ifelse(var_trt == "constant", 
                          1 , 
                          (1 - mean_quality) * parse_conc(low_diet) + mean_quality * parse_conc(high_diet))
  )

d_event2 %>% 
  filter(
    #pupated_cam_end == 0
  ) %>% 
  filter(!is.na(RGR)) %>% 
  ggplot(aes(x = var_quality, y = RGR)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~tz) + 
  labs(x = "12 hour window variance in diet toxin (mg / g)^2") + 
  theme_bw(base_size = 15)

event_d <- event_list %>% 
  do.call("rbind",.) %>% 
  mutate(
    tz = ifelse(rank < (24 * 10), "0h - 24h", "24h - 120h")
  ) %>% 
  filter(
   tz == "24h - 120h" 
  )


rr <- lapply(
  seq(5, 1200, len = 1200/5),
  function(w){
    tryCatch(event_d %>% 
      group_by(rep_id) %>% 
      summarise(
        var_quality = mean(moving_var(head_high_trt, n = w), na.rm = TRUE)
      ) %>% 
      left_join(
        ref_data %>% 
          mutate(rep_id = as.character(rep_id)) %>% 
          select(RGR, rep_id), by = "rep_id"
      ) %>% 
      select(
        var_quality, RGR
      ) %>% 
      cor(use = "co") %>% 
      .[1,2] %>% 
      .^2, error = function(e){
        NA
      }) ->o
    cat(w,"\r")
    return(o)
  }
) %>% do.call("c",.)


data.frame("r2" = rr, "w" = seq(5, 1200, len = 1200 / 5)) %>% 
  ggplot(aes(x = w, y = rr)) + 
  geom_point() 



ref_data %>% 
  filter(var_trt == "constant") %>% 
  ggplot(aes(x = parse_conc(mean_trt), y= RGR)) + 
  geom_point() + 
  geom_smooth()
















d_event2 %>% 
  filter(mean_trt == "1 mg/g") %>% 
  ggplot(aes(x = factor(beta), y = tot_sd, color = var_trt)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6)) + 
  geom_pointrange(stat = "summary", position = position_dodge(0.6)) + 
  facet_wrap(~tz) + 
  labs(x = expression(beta), y = "Spatial SD", color = "Var trt") + 
  theme_bw(base_size = 15) + 
  scale_y_continuous(trans = "log10") +
  theme(legend.position = "top") -> g1


d_event2 %>% 
  filter(mean_trt == "1 mg/g") %>% 
  ggplot(aes(x = log(cat_pre_wt), y = tot_sd, color = factor(var_trt))) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  # ggpubr::stat_cor(
  #   aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~"))
  # )+ 
  facet_wrap(~tz)



m <- glmmTMB(
  log(tot_sd) ~ 
    tz + beta + var_trt + log(cat_pre_wt) + 
    (1|session_id),  
  family = gaussian(),
  data = d_event2 %>% 
    filter(var_trt != "constant") %>% 
    filter(
      
    )
); summary(m)


plot_model(m, 
           type = "pred", 
           terms = c("cat_pre_wt[all]", "var_trt")) +
  theme_bw(base_size = 15) + 
  labs(title = NULL, 
       subtitle = NULL, 
       y = "Spatial SD", x = "Cat pre weight", color = "var_trt") + 
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +
  theme(legend.position = "top")






d_event <- event_list %>% 
  do.call("rbind",.) %>% 
  mutate(
    tz = ifelse(rank < (24 * 6), "0h - 24h", "24h - 120h")
  )






d_event2 %>% 
  ggplot(aes(x = factor(beta), y = mean_quality, color = var_trt)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6)) + 
  geom_pointrange(stat = "summary", position = position_dodge(0.6)) + 
  facet_wrap(~tz) + 
  labs(x = expression(beta), y = "Mean diet toxin (mg / g)", color = "Var trt") + 
  theme_bw(base_size = 15) + 
  theme(legend.position = "top") 



d_event2 %>% 
  filter(mean_trt == "1 mg/g") %>% 
  ggplot(aes(x = (cat_pre_wt), y = mean_quality, color = factor(beta))) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_wrap(~tz) + 
  scale_x_continuous(trans = "log10") +
  labs(x = "Cat pre weight", y = "Mean diet toxin (mg / g)", color = expression(beta)) + 
  theme_bw(base_size = 15) + 
  theme(legend.position = "top") 




m <- glmmTMB(
  mean_quality ~ 
    var_trt + tz + beta * log(cat_pre_wt) + 
    (1|session_id),  
  family = gaussian(),
  data = d_event2 %>% 
    filter(var_trt != "constant") 
); summary(m)

plot_model(m, 
           type = "pred", 
           terms = c("cat_pre_wt[all]", "beta")) +
  theme_bw(base_size = 15) + 
  labs(title = NULL, 
       subtitle = NULL, 
       y = "Mean diet toxin (mg /g)", x = "Cat pre weight", color = "beta") + 
  scale_x_continuous(trans = "log10") +
  theme(legend.position = "top")




































ref_data %>% 
  filter(mean_trt == "1 mg/g") %>% 
  ggplot(aes(x = beta, y = RGR)) + 
  geom_point(position = position_jitter(height = 0, width = 0.2)) + 
  geom_pointrange(stat = "summary")

ref_data %>% 
  #filter(var_trt = "constant") %>% 
  filter(mean_trt == "1 mg/g") %>% 
  ggplot(aes(x = var_trt, y = RGR, color = as.factor(beta))) + 
  geom_point(position = position_jitterdodge(jitter.height = 0, jitter.width = 0.2, dodge.width = 0.5)) + 
  geom_pointrange(stat = "summary", position = position_dodge(width = 0.5)) + 
  theme_bw(base_size = 15)  + 
  labs(x = "var_trt", color = expression(beta))


ref_data %>% 
  filter(mean_trt == "1 mg/g") %>% 
  ggplot(aes(x = beta, y = cat_dead_cam_end)) + 
  geom_point(position = position_jitter(height = 0)) + 
  geom_pointrange(stat = "summary") + 
  theme_bw(base_size = 15)  + 
  labs(y = "Dead within 5 days", x = expression(beta))


rgr_m <- glmmTMB(
  RGR ~ 
    log(cat_pre_wt) * (as.factor(beta) * var_trt) + pupated_cam_end + (1|session_id), 
  family = gaussian(), 
  data = ref_data %>% 
    filter(!is.na(cat_pre_wt))
); summary(rgr_m)

plot_model(rgr_m, type = "pred", terms = c("cat_pre_wt[all]","var_trt", "beta")) + 
  labs(title = NULL, 
       subtitle = NULL, 
       y = "RGR (h^-1)", x = "Cat pre weight", color = "var_trt") + 
  scale_x_continuous(trans = "log10") +
  theme(legend.position = "top") +
  theme_bw(base_size = 15)


plot_model(rgr_m, type = "pred", terms = c("cat_pre_wt[all]", "var_trt"))


death_m <- glmmTMB(
  cat_dead_cam_end ~ 
    log(cat_pre_wt) + var_trt + as.factor(beta) + (1|session_id), 
  family = binomial(), 
  data = ref_data %>% 
    filter(!is.na(cat_pre_wt))
); summary(death_m)

plot_model(death_m, type = "pred", terms = c("beta")) + 
  theme_bw(base_size = 15) + 
  labs(title = NULL, 
       subtitle = NULL, 
       y = "Prob. dead in 5 days", x = "Beta")




pup_m <- glmmTMB(
  pupated ~ 
    log(cat_pre_wt) + var_trt + as.factor(beta) + (1|session_id), 
  family = binomial(), 
  data = ref_data %>% 
    filter(!is.na(cat_pre_wt))
); summary(pup_m)

plot_model(pup_m, type = "pred", terms = c("beta")) + 
  theme_bw(base_size = 15) + 
  labs(title = NULL, 
       subtitle = NULL, 
       y = "Prop. Pupated", x = "Beta")


move_list %>% 
  do.call("rbind",.) %>% 
  select(r) %>% 
  summary()





d_event2 %>% 
  filter(
    #pupated_cam_end == 0
  ) %>% 
  filter(!is.na(RGR)) %>% 
  ggplot(aes(x = var_quality, y = RGR)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~tz) + 
  labs(x = "Hourly variance in diet toxin (mg / g)^2") + 
  theme_bw(base_size = 15)







