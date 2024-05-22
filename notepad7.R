source("spat1f/init_analysis.R")



plot_track_overlay(repID = 104, plot_elements = "track", colored_track = "states")
subm_rgr <- glmmTMB(
  RGR_scale ~ 
    cat_pre_wt_log_scale + 
    cat_pre_wt_log_scale_sq +
    var_trt + beta_numeric_scale * cat_pre_wt_log_scale +
    (1|session_id),
  data = d2
); summary(subm_rgr)

subm_rgr <- glmmTMB(
  RGR_scale ~ 
    cat_pre_wt_log_scale + 
    cat_pre_wt_log_scale_sq +
    var_toxic_12_scale + 
    sl_mean_obs_log_scale + 
    mean_toxic_conc_scale + 
    area_herb_log_scale +
    beta_numeric_scale * cat_pre_wt_log_scale + 
    var_trt * (sl_mean_obs_log_scale) + 
    prop_explore_logit_scale * (beta_numeric_scale + var_trt) + 
    (1|session_id),
  data = d %>% 
    mutate(
      beta_numeric_scale = as.numeric(scale(as.numeric(beta)))
    )
); summary(subm_rgr)

check_collinearity(subm_rgr)

plot_model(subm_rgr, type = "eff", terms = c("prop_explore_logit_scale", "beta_numeric_scale "))


subm_var_toxic <- glmmTMB(
  var_toxic_12_scale ~ 
    beta_numeric_scale + 
    var_high + 
    sl_mean_obs_log_scale + 
    (1|session_id),
  data = d2
); summary(subm_var_toxic)


subm_sl <- glmmTMB(
  sl_mean_obs_log_scale ~ 
    beta_numeric_scale + 
    # on_toxic_logit_scale + 
    var_high +  
    cat_pre_wt_log_scale + 
    cat_pre_wt_log_scale_sq + 
    (1|session_id),
  data = d2,
); summary(subm_sl)

subm_toxin_ingested <- glmmTMB(
  mean_toxic_conc_scale ~ 
    beta_numeric_scale * cat_pre_wt_log_scale + 
    var_high * cat_pre_wt_log_scale + 
    cat_pre_wt_log_scale_sq +  
    (1|session_id),
  data = d2 
); summary(subm_toxin_ingested)


subm_on_toxic <- glmmTMB(
  on_toxic_logit_scale ~ 
    `state1:less_toxic` + `state2:less_toxic` + 
    var_high + beta_numeric_scale + 
    (1|session_id),
  family = gaussian(),
  data = d2
); summary(subm_on_toxic)


subm_herb <- glmmTMB(
  area_herb_log_scale ~ 
    beta_numeric_scale + var_high * cat_pre_wt_log_scale + 
    (1|session_id),
  data = d2,
); summary(subm_herb)

subm_select <- glmmTMB(
  `state1:less_toxic` ~ 
    beta_numeric_scale + var_high + cat_pre_wt_log_scale + 
    (1|session_id),
  data = d2,
); summary(subm_select)

subm_select <- glmmTMB(
  `state2:less_toxic` ~ 
    beta_numeric_scale + var_high + cat_pre_wt_log_scale + 
    (1|session_id),
  data = d2,
); summary(subm_select)

subm_sl <- glmmTMB(
  log(sl_mean_pred1) ~ 
    beta_numeric_scale + var_high + cat_pre_wt_log_scale + 
    (1|session_id),
  data = d2 %>% 
    filter(sl_mean_pred1 < 1000),
); summary(subm_sl)

# subm_ava <- glmmTMB(
#   ava_qual_logit_scale ~ 
#     var_high + beta_numeric_scale + cat_pre_wt_log_scale + 
#     (1|session_id),
#   family = gaussian(),
#   data = d2,
# ); summary(subm_ava)

subm_exp <- glmmTMB(
  prop_explore_logit_scale ~ 
    var_high + beta_numeric_scale + cat_pre_wt_log_scale + 
    (1|session_id),
  family = gaussian(),
  data = d2,
); summary(subm_exp)


plot_model(subm_exp, type = "eff", terms = c("cat_pre_wt_log_scale", "beta_numeric_scale", "var_high"))


d2 %>% 
  filter(`state2:less_toxic` < 5) %>% 
  ggplot(aes(x = cat_pre_wt_log_scale, y = `state2:less_toxic`)) + 
  geom_point() + 
  geom_smooth(method = "lm")


View(d2)

d2 %>% 
  # filter(sl_mean_pred1 < 1000) %>% 
  ggplot(aes(x = cat_pre_wt_log_scale, y = log(sl_mean_pred1))) + 
  geom_point() + 
  geom_smooth(method = "lm")


z <- fetch_events(121) %>% 
  clean_events() %$% 
  move_seq(head_x, head_y)

plot_track_overlay(repID = 150, colored_track = "states", plot_elements = "track")


fit <- z %>% 
  as.moveData() %>% 
  fit_HMM()

x <- ((z$r)[viterbi(fit) == 2])

z$theta_rel %>% 
  loghist(log.p = FALSE,log.x = FALSE, 
              draw_dist = list(
                fit_vonmises(x), 
                fit_genvonmises(x)
              ), nclass = 25) + 
  scale_x_ta()




x %>% loghist(log.p = FALSE, 
              draw_dist = list(
                fit_gamma(x), 
                fit_lnorm(x),
                fit_invgamma(x),
                fit_dpln(x)
              ), nclass = 50)

comp_dist(x, dist = c("lnorm", "gamma", "dpln","invgamma"))



make_fit("wald", params = c("mu", "lambda"), auto_init = function(x){
  c(mean(x), sd(x))
})















