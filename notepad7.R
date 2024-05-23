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
  ava_qual ~ 
    var_high + beta_numeric_scale * cat_pre_wt_log_scale +
    (1|session_id),
  family = beta_family(),
  data = d2
); summary(subm_on_toxic)

plot_model(subm_on_toxic, type = "eff", terms = c("cat_pre_wt_log_scale", "beta_numeric_scale"))


subm_herb <- glmmTMB(
  area_herb_log_scale ~ 
    beta_numeric_scale + var_high * cat_pre_wt_log_scale + 
    (1|session_id),
  data = d2,
); summary(subm_herb)

subm_select <- glmmTMB(
  s2.less_toxic.est ~ 
    beta_numeric_scale + var_high + cat_pre_wt_log_scale + prop_explore_logit_scale + 
    (1|session_id),
  data = d3,
); summary(subm_select)

subm_select <- glmmTMB(
  s1.less_toxic.est ~ 
    beta_numeric_scale + var_high + cat_pre_wt_log_scale + prop_explore_logit_scale + 
    (1|session_id),
  data = d3,
); summary(subm_select)

subm_sl <- glmmTMB(
  log(scale_1) ~ 
    cat_pre_wt_log_scale + 
    (1|session_id),
  data = d3 %>% 
    filter(),
); summary(subm_sl)

subm_sl <- glmmTMB(
  log(k1) ~ 
    cat_pre_wt_log_scale * beta + var_trt + 
    (1|session_id),
  data = d3 %>% 
    filter(),
); summary(subm_sl)
plot_model(subm_sl, type = "eff", terms = c("cat_pre_wt_log_scale", "beta"))

contrast_by_pre_wt(subm_sl, "beta")

# subm_ava <- glmmTMB(
#   ava_qual_logit_scale ~ 
#     var_high + beta_numeric_scale + cat_pre_wt_log_scale + 
#     (1|session_id),
#   family = gaussian(),
#   data = d2,
# ); summary(subm_ava)

expand.grid(c("toxic", "less_toxic"), c("state1","state2"))

subm_exp <- glmmTMB(
  prop_explore_logit_scale ~ 
    var_high + beta_numeric_scale * cat_pre_wt_log_scale + 
    (1|session_id),
  family = gaussian(),
  data = d2,
); summary(subm_exp)


plot_model(subm_exp, type = "eff", terms = c("cat_pre_wt_log_scale", "beta_numeric_scale"))


