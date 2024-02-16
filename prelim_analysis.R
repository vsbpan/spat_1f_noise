source("helper_functions/init_analysis.R")

d <- ref_data %>% 
  filter(error == 0) %>% 
  mutate(
    cat_pre_wt_log = log(cat_pre_wt),
    cat_pre_wt_log_scale = as.numeric(scale(log(cat_pre_wt))),
    cat_size = ifelse(cat_pre_wt < median(cat_pre_wt), "small", "big"),
    cat_big = ifelse(cat_size == "big", 1, 0),
    has_var = ifelse(var_trt != "constant", 1, 0),
    mean_toxic_conc_scale = as.numeric(scale(mean_toxic_conc)),
    mean_toxic_scale = as.numeric(scale(mean_toxic)),
    area_herb_log_scale = as.numeric(scale(log(area_herb+1))),
    var_toxic_12_scale = as.numeric(scale(var_toxic_12)),
    on_toxic_logit_scale = as.numeric(scale(adjust_prop(on_toxic, 
                                                        trans = "emp", 
                                                        nudge.method = "none", 
                                                        na.action = "ignore"))),
    toxic_select_scale = as.numeric(scale(toxic)),
    sl_mean_obs_log_scale = as.numeric(scale(log(sl_mean_obs))),
    ava_mean_toxin_scale = as.numeric(scale(ava_mean_toxin)),
    cat_pre_wt_log_scale_sq = as.numeric(scale(log(cat_pre_wt)^2))
  )
  # mutate(
  #   `toxic:start_toxicyes` = ifelse(`toxic:start_toxicyes` < -10, NA, `toxic:start_toxicyes`),
  #   switch_pref = (`toxic:start_toxicno` - `toxic:start_toxicyes` )/2,
  #   mean_pref = (`toxic:start_toxicno` + `toxic:start_toxicyes`) / 2
  # ) %>% 
  # rename(
  #   toxic_start_no = `toxic:start_toxicno`,
  #   toxic_start_yes = `toxic:start_toxicyes`
  # )



plot_track_overlay(repID=problem_ids[10])

fetch_events(84)$score %>% hist()

fetch_events(84) %>% 
  filter(false_cluster)

plot_track_overlay(
  events = fetch_events(84) %>% 
    clean_events() %>% 
    filter(score > 0.9) %>% 
    insert_gaps()
)


plot_track_overlay(repID=problem_ids[11])

fetch_events(problem_ids_2[5]) %>% 
  clean_events() %>% 
  filter(score > 0.9) %>% 
  flag_false_cluster(bin_size = 60, r_thresh = 200, r = 60) %>% 
  filter(!false_cluster) %>%
  flag_false_cluster(bin_size = 60, r_thresh = 200, r = 60) %>% 
  filter(!false_cluster) %>%
  plot_track_overlay()



fetch_data_dict(84) %>% plot(frame = 606)

z2 %>% 
  mutate(
    sus = sqrt((x1 - z3$cx1)^2 + (y1 - z3$cy1)^2)
  ) %>% 
  View()

d %>% 
  filter(mean_trt_numeric == 1) %>% 
  ggplot(aes(x = sex, y = sl_mean_obs, color = sex, fill = sex)) + 
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6)) + 
  geom_pointrange(stat = "summary", position = position_dodge(0.6), color = "black") + 
  theme_bw(base_size = 15) + 
  scale_y_continuous(trans = "log2") + 
  labs(y = "Mean step length (cm)")


d %>% 
  ggplot(aes(x = sex, y = toxic, color = sex, fill = sex)) + 
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6)) + 
  geom_pointrange(stat = "summary", position = position_dodge(0.6), color = "black") + 
  theme_bw(base_size = 15) + 
  labs(y = "Preference for toxic diet")

d %>% 
  filter(var_trt != "constant") %>% 
  filter(`toxic:start_toxicyes` > -10) %>% 
  select(-today) %>% 
  gather(key = start, value = val, `toxic:start_toxicyes`:`toxic:start_toxicno`) %>% 
  ggplot(aes(x = beta, 
             y = val, 
             color = start, 
             fill = start)) + 
  geom_point(
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6), color = "black") + 
  geom_pointrange(stat = "summary", position = position_dodge(0.6)) + 
  theme_bw(base_size = 15) 

d %>% 
  filter(mean_trt == "1 mg/g") %>% 
  #filter(toxic_start_yes > -10) %>% 
  ggplot(aes(x = var_trt, 
             y = toxic,  # (no - yes)/2 = switch preference
             # (no + yes)/2 = toxic preference
             color = beta, 
             fill = beta)) + 
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6)) + 
  geom_pointrange(stat = "summary", position = position_dodge(0.6), color = "black") + 
  theme_bw(base_size = 15) + 
  facet_wrap(~cat_size)

d %>% 
  filter(var_trt != "constant") %>% 
  filter(abs(toxic) < 10) %>% 
  ggplot(aes(x = var_trt, 
             y =  toxic,  # (no - yes)/2 = switch preference
             # (no + yes)/2 = toxic preference
             color = beta, 
             fill = beta)) + 
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6)) + 
  geom_pointrange(stat = "summary", position = position_dodge(0.6), color = "black") + 
  theme_bw(base_size = 15)



d %>% 
  filter(var_trt != "constant") %>% 
  #filter(`toxic:start_toxicyes` > -10) %>% 
  ggplot(aes(y = on_toxic, 
             x = log(sl_mean_obs),  # (no - yes)/2 = switch preference
             # (no + yes)/2 = toxic preference
             color = beta, 
             fill = beta)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_bw(base_size = 15)

d %>% 
  filter(mean_trt == "1 mg/g") %>% 
  mutate(cat_size = ifelse(cat_pre_wt > median(cat_pre_wt), "big", "small")) %>% 
  ggplot(aes(x = beta, y = RGR, color = beta, fill = beta)) + 
  geom_pointrange(stat = "summary", position = position_dodge(0.6), color = "black") +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6)
    ) + 
  facet_wrap(~ cat_size)


names(d)



d %>% 
  #left_join(w, by = "rep_id") %>% 
  filter(var_trt != "constant") %>% 
  #filter(mean_trt == "1 mg/g") %>% 
  ggplot(aes(x = log(var_toxic_12), y = RGR, color = beta)) + 
  geom_point() + 
  #geom_pointrange(aes(ymax = toxic_est + toxic_se, ymin = toxic_est - toxic_se)) + 
  geom_smooth(method = "lm") 
  #scale_y_continuous(trans = "log")
  geom_abline(aes(intercept = 0, slope = 1))


  
  
glmmTMB(
    RGR ~ 
      scale(mean_toxic_conc) + 
      scale(log(area_herb+1)) + 
      cat_pre_wt_log_scale + 
      scale(sl_mean_obs) + 
      scale(var_toxic_12) + 
      (1|session_id),
    data = d %>% 
      filter(
        mean_trt_numeric == 1 & var_trt != "constant"
      )
  ) %>% summary()

m<- glmmTMB(
  on_toxic ~ 
    var_trt + beta + cat_pre_wt_log_scale + 
    (1|session_id),
  data = d %>% 
    filter(
      var_trt != "constant"
    )
);summary(m)
plot_model(m,type = "eff", terms = c("cat_pre_wt_log_scale", "beta")) 

m <- glmmTMB(
  area_herb ~ 
    var_trt * beta + cat_pre_wt_log_scale + 
    (1|session_id),
  family = Gamma(link = "log"),
  data = d %>% 
    filter(
      var_trt != "constant"
    )
);summary(m)
plot_model(m, type = "eff", terms = c("var_trt", "beta")) 





glmmTMB(
  prob_move_obs ~ 
    var_trt + beta * cat_pre_wt_log_scale + 
    (1|session_id),
  family = beta_family(),
  data = d %>% 
    filter(
      var_trt != "constant"
    )
) %>% plot_model(type = "eff", terms = c("cat_pre_wt_log_scale", "beta")) 

glmmTMB(
  log(sl_mean_obs) ~ 
    var_trt + beta * cat_pre_wt_log_scale + 
    (1|session_id),
  family = gaussian(),
  data = d %>% 
    filter(
      var_trt != "constant"
    )
) %>% plot_model(type = "eff", terms = c("cat_pre_wt_log_scale", "beta")) 


glmmTMB(
  RGR ~ 
    (var_trt + beta) * cat_pre_wt_log_scale + 
    #mean_toxic_conc + cat_pre_wt_log_scale + 
    #log(sl_mean_obs) + log(area_herb + 1) +
    #var_toxic_12 + 
    (1|session_id),
  data = d %>% 
    filter(
      mean_trt_numeric == 1
    ) 
) %>% summary()
  plot_model(type = "eff", terms = c("cat_pre_wt[all]", "beta"), show.data = TRUE) + 
  scale_x_continuous(trans = "log10")



d$pupation_time %>% hist()

m <- glmmTMB(
  log(cat_pre_wt) ~ 
    (var_trt + beta) + sex +
    (1|session_id),
  data = d, 
  family = gaussian(),
); summary(m)

m <- glmmTMB(
  pupation_time ~ 
    (var_trt * beta) + cat_pre_wt_log_scale + 
    (1|session_id),
  data = d %>% 
    filter(var_trt != "constant"), 
  family = Gamma(link = "log"),
); summary(m)

m <- glmmTMB(
  on_toxic_conc ~ 
    var_trt + beta * cat_pre_wt_log_scale +  
    (1 | session_id),
  data = d,
  family = gaussian(),
); summary(m)

plot_model(m, type = "eff", terms = c("beta", "var_trt"), show.data = TRUE, dodge = TRUE)


m <- glmmTMB(
  var_toxic_12 ~ 
    cat_pre_wt_log_scale + var_trt + cat_pre_wt_log_scale + beta + 
    #(var_trt + beta) * cat_pre_wt_log_scale +
    (1 | session_id),
  data = d,
  family = gaussian(),
); summary(m)

m <- glmmTMB(
  mean_toxic ~ 
    ava_mean_toxin * mean_pref + cat_pre_wt_log_scale + var_trt + 
    #(var_trt + beta) * cat_pre_wt_log_scale +
    (1 | session_id),
  data = d %>% 
    filter(
      `toxic:start_toxicyes` > -10
    ) %>% 
    mutate(
      switch_pref = `toxic:start_toxicno` - `toxic:start_toxicyes`,
      mean_pref = `toxic:start_toxicno` + `toxic:start_toxicyes`
    ), 
  family = gaussian(),
); summary(m)


check_model(m)

DHARMa::simulateResiduals(m) %>% plot()

plot_model(m, type = "eff", terms = c("ava_switch_toxin", "switch_pref"), show.data = TRUE)
plot_model(m, type = "eff", terms = c("ava_mean_toxin", "mean_pref"), show.data = TRUE)


m <- glmmTMB(
  adult_time ~ 
    var_trt + beta * cat_pre_wt_log_scale +
    (1 | session_id),
  data = d, 
  family = Gamma(link = "log"),
); summary(m)


marginal_effects(m, terms = c("cat_pre_wt_log_scale", "beta")) %>% 
  ggplot(aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), 
             y = yhat)) + 
  geom_ribbon(aes(fill = beta, ymin = lower, ymax = upper), alpha = 0.2) + 
  geom_line(aes(color = beta), size = 2)  + 
  theme_bw(base_size = 15) + 
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") + 
  geom_point(
    data = m$frame, aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), 
                        y = exp(m$frame[,1]), 
                        color = beta),
    alpha = 0.5,
    size = 3
  ) + 
  scale_color_brewer(type = "qual", aesthetics = c("fill", "color"))




DHARMa::simulateResiduals(m) %>% plot()

coxph(
  Surv(adult_time, observed_dead) ~ 
    (var_trt + beta) * log(cat_pre_wt) +
    strata(session_id), 
  data = d
) %>% summary()

m <- coxph(
  Surv(surv_time, observed_dead) ~ 
    var_trt + beta + cat_pre_wt_log_scale + I(cat_pre_wt_log_scale^2) + 
    strata(session_id), 
  data = d %>% 
    filter(var_trt != "constant") %>% 
    mutate(
      beta = as.factor(as.character(beta))
    )
);summary(m)

cox.zph(m)
drop1(m, test = "Chisq")


m <- glmmTMB(
  time_to_eclosure ~ 
    var_trt + beta * cat_pre_wt_log_scale +
    (1|session_id),
  data = d %>% 
    filter(var_trt != "constant"), 
  family = Gamma(link = "log"),
); summary(m)

m <- glmmTMB(
  eclosed ~ 
    var_trt + beta + cat_pre_wt_log_scale +
    (1|session_id),
  data = d%>% 
    filter(var_trt != "constant") %>% 
    mutate(
      eclosed = ifelse(deformed_adult == 1, 0, eclosed)
    ), 
  family = binomial(),
); summary(m)


m <- glmmTMB(
  time_to_pupation ~ 
    var_trt + beta * cat_pre_wt_log_scale + I(cat_pre_wt_log_scale^2) + 
    #var_toxic_12 + mean_toxic_conc + cat_pre_wt_log_scale + 
    (1|session_id),
  data = d %>% 
    filter(var_trt != "constant"), 
  family = Gamma(link = "log"),
); summary(m)

plot_model(m , type = "eff", terms = c("cat_pre_wt[all]", "beta")) + 
  scale_y_log10()

car::Anova(m, type = "III")


m <- glmmTMB(
  log(pupal_weight) ~ 
    (var_trt + beta) + cat_pre_wt_log_scale +
    (1|session_id),
  data = d, 
  family = gaussian(),
); summary(m)


problem_ids





# 45, 46, 49

ggarrange(
  plot_track_overlay(repID = 45),
  plot_track_overlay(repID = 46),
  plot_track_overlay(repID = 49),
  plot_track_overlay(repID = 61),
  common.legend = FALSE, ncol = 2, nrow = 2, legend = "right"
)

