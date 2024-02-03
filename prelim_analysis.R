source("helper_functions/init_analysis.R")

d <- ref_data %>% 
  filter(error == 0) %>% 
  mutate(
    cat_pre_wt_log = log(cat_pre_wt),
    cat_pre_wt_log_scale = as.numeric(scale(log(cat_pre_wt)))
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
             y = log(sl_mean_obs),  # (no - yes)/2 = switch preference
             # (no + yes)/2 = toxic preference
             color = beta, 
             fill = beta)) + 
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6)) + 
  geom_pointrange(stat = "summary", position = position_dodge(0.6), color = "black") + 
  theme_bw(base_size = 15)

d %>% 
  filter(var_trt != "constant") %>% 
  filter(toxic_start_yes > -10) %>% 
  ggplot(aes(x = var_trt, 
             y = switch_pref,  # (no - yes)/2 = switch preference
             # (no + yes)/2 = toxic preference
             color = beta, 
             fill = beta)) + 
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6)) + 
  geom_pointrange(stat = "summary", position = position_dodge(0.6), color = "black") + 
  theme_bw(base_size = 15)



d %>% 
  filter(var_trt != "constant") %>% 
  filter(`toxic:start_toxicyes` > -10) %>% 
  ggplot(aes(y = RGR, 
             x = `toxic:start_toxicno` - `toxic:start_toxicyes`,  # (no - yes)/2 = switch preference
             # (no + yes)/2 = toxic preference
             color = beta, 
             fill = beta)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_bw(base_size = 15)

d %>% 
  filter(mean_trt == "1 mg/g") %>% 
  mutate(cat_size = ifelse(cat_pre_wt > 0.03, "big", "small")) %>% 
  ggplot(aes(x = var_trt, y = RGR, color = beta, fill = beta)) + 
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
  ggplot(aes(x = log(cat_pre_wt), y = log(sl_mean_obs), color = beta)) + 
  geom_point() + 
  #geom_pointrange(aes(ymax = toxic_est + toxic_se, ymin = toxic_est - toxic_se)) + 
  geom_smooth(method = "lm") 
  #scale_y_continuous(trans = "log")
  geom_abline(aes(intercept = 0, slope = 1))


  
  
glmmTMB(
    RGR ~ 
      scale(mean_toxic_conc) + scale(log(area_herb + 1)) + cat_pre_wt_log_scale + 
      scale(sl_mean_obs) + 
      (1|session_id),
    data = d %>% 
      filter(
        var_trt != "constant"
      )
  ) %>% summary()

glmmTMB(
  mean_toxic_conc ~ 
    var_trt + beta * cat_pre_wt_log_scale + 
    (1|session_id),
  data = d %>% 
    filter(
      var_trt != "constant"
    )
) %>% plot_model(type = "eff", terms = c("cat_pre_wt_log_scale", "beta")) 

glmmTMB(
  log(sl_mean_obs) ~ 
    var_trt + beta * cat_pre_wt_log_scale + 
    (1|session_id),
  data = d %>% 
    filter(
      var_trt != "constant"
    )
) %>% plot_model(type = "eff", terms = c("cat_pre_wt_log_scale", "beta")) 

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
  kappa1 ~ 
    var_trt + beta * cat_pre_wt_log_scale + 
    (1|session_id),
  family = gaussian(),
  data = d %>% 
    filter(
      var_trt != "constant"
    )
) %>% summary() #%>% plot_model(type = "eff", terms = c("cat_pre_wt_log_scale", "beta")) 


glmmTMB(
  RGR ~ 
    (var_trt + beta) * cat_pre_wt_log_scale + 
    (1|session_id),
  data = d
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
  data = d, 
  family = Gamma(link = "log"),
); summary(m)

m <- glmmTMB(
  toxic ~ 
    var_trt + beta + cat_pre_wt_log_scale +  
    (1 | session_id),
  data = d,
  family = gaussian(),
); summary(m)

plot_model(m, type = "eff", terms = c("beta", "var_trt"), show.data = TRUE, dodge = TRUE)


m <- glmmTMB(
  var_toxic_12 ~ 
    ava_switch_toxin * switch_pref + cat_pre_wt_log_scale + var_trt + 
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

coxph(
  Surv(surv_time, observed_dead) ~ 
    var_trt + (beta) * cat_pre_wt_log_scale +
    strata(session_id), 
  data = d %>% 
    filter(var_trt != "constant")
) %>% summary()

m <- glmmTMB(
  time_to_eclosure ~ 
    var_trt + beta * log(cat_pre_wt) +
    (1|session_id),
  data = d, 
  family = Gamma(link = "log"),
); summary(m)

m <- glmmTMB(
  eclosed ~ 
    var_trt + beta + cat_pre_wt_log_scale +
    (1|session_id),
  data = d, 
  family = binomial(),
); summary(m)


m <- glmmTMB(
  time_to_pupation ~ 
    var_trt + beta * log(cat_pre_wt) + 
    (1|session_id),
  data = d, 
  family = Gamma(link = "log"),
); summary(m)
plot_model(m , type = "eff", terms = c("cat_pre_wt[all]", "beta")) + 
  scale_y_log10()


m <- glmmTMB(
  log(pupal_weight) ~ 
    (var_trt + beta) + cat_pre_wt_log_scale +
    (1|session_id),
  data = d, 
  family = gaussian(),
); summary(m)






