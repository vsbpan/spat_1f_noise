source("spat1f/init_analysis.R")



plot_track_overlay(repID = 104, plot_elements = "track", colored_track = "states")

subm_rgr <- glmmTMB(
  RGR_scale ~ 
    cat_pre_wt_log_scale + 
    cat_pre_wt_log_scale_sq + 
    sl_mean_obs_log_scale +
    prop_explore_logit_scale *  
    var_toxic_12_scale + 
    mean_toxic_conc_scale + 
    area_herb_log_scale + 
    beta_numeric_scale + 
    (1|session_id),
  data = d2
); summary(subm_rgr)

check_collinearity(subm_rgr)

plot_model(subm_rgr, type = "eff", terms = c("var_toxic_12_scale", "prop_explore_logit_scale"))


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
    var_high +  
    cat_pre_wt_log_scale + 
    I(cat_pre_wt_log_scale^2) + 
    (1|session_id),
  data = d2,
); summary(subm_sl)

plot_model(subm_sl, type = "eff", terms = c("cat_pre_wt_log_scale"))


subm_toxin_ingested <- glmmTMB(
  mean_toxic_conc_scale ~ 
    beta_numeric_scale + cat_pre_wt_log_scale + 
    var_high * cat_pre_wt_log_scale + 
    cat_pre_wt_log_scale_sq +  
    (1|session_id),
  data = d2 
); summary(subm_toxin_ingested)


subm_on_toxic <- glmmTMB(
  on_toxic_logit_scale ~ 
    scale(`state1:less_toxic`) + scale(`state2:less_toxic`) + 
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
  `state1:less_toxic` ~ 
    beta_numeric_scale + var_high + cat_pre_wt_log_scale + 
    (1|session_id),
  data = d2 %>% 
    filter(abs(`state1:less_toxic`) < 5),
); summary(subm_select)

subm_sl <- glmmTMB(
  log(sl_mean_pred4) ~ 
    beta_numeric_scale + var_high + cat_pre_wt_log_scale + 
    (1|session_id),
  data = d2 %>% 
    filter(),
); summary(subm_sl)

# subm_ava <- glmmTMB(
#   ava_qual_logit_scale ~ 
#     var_high + beta_numeric_scale + cat_pre_wt_log_scale + 
#     (1|session_id),
#   family = gaussian(),
#   data = d2,
# ); summary(subm_ava)

expand.grid(c("state1","state2"),c("toxic", "less_toxic"))

subm_exp <- glmmTMB(
  prop_explore_logit_scale ~ 
    var_high + beta_numeric_scale * cat_pre_wt_log_scale + 
    (1|session_id),
  family = gaussian(),
  data = d2,
); summary(subm_exp)


plot_model(subm_exp, type = "eff", terms = c("cat_pre_wt_log_scale", "beta_numeric_scale"))


d2 %>% 
  filter(sl_mean_pred1 > 0) %>% 
  filter(abs(`state2:less_toxic`) < 5) %>% 
  #filter(`state2:less_toxic` < 5) %>% 
  ggplot(aes(x = cat_pre_wt_log_scale, y = `state2:less_toxic`)) + 
  geom_point() + 
  geom_smooth(method = "lm")


View(d2)
d2 %>% filter(
  sl_mean_pred1 < 0
) %>% View()
plot_track_overlay(repID = 95, plot_elements = "track", colored_track = "states")


d2 %>% 
  filter(sl_mean_pred1 > 0) %>% 
  # filter(sl_mean_pred1 < 1000) %>% 
  ggplot(aes(x = log(sl_mean_pred1), y = log(sl_mean_pred3))) +
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() + 
  geom_smooth(method = "lm")


d2 %>% 
  filter(sl_mean_pred1 > 0) %>% 
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









.get_prop_state1_engine(55)

fetch_events(55) %>% 
  clean_events() %$%
  move_seq(head_x, head_y) %>% 
  as.moveData() %>% 
  fit_HMM() %>% 
  viterbi() %>% 
  vapply(function(x){
    x == 1
  }, numeric(1)) %>% 
  mean()




