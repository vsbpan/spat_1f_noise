source("spat1f/init_analysis.R")



subm_sl <- glmmTMB(
  log(scale_4) ~ 
    cat_pre_wt_log_scale + 
    (1|session_id),
  data = d3 %>% 
    filter(),
); summary(subm_sl)

subm_sl <- glmmTMB(
  log(k2) ~ 
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





















d <- ref_data %>% 
  filter(error == 0) %>% 
  mutate(
    cat_pre_wt_log = log(cat_pre_wt),
    cat_pre_wt_log_scale = as.numeric(scale(log(cat_pre_wt))),
    cat_size = ifelse(cat_pre_wt < median(cat_pre_wt), "small", "big"),
    has_var = ifelse(var_trt != "constant", 1, 0),
    mean_toxic_conc_scale = as.numeric(scale(mean_toxic_conc)),
    mean_toxic_scale = as.numeric(scale(mean_toxic)),
    area_herb_log_scale = as.numeric(scale(log(area_herb+1))),
    var_toxic_12_scale = as.numeric(scale(var_toxic_12)),
    on_toxic_logit_scale = as.numeric(scale(adjust_prop(on_toxic, 
                                                        trans = "emp", 
                                                        nudge.method = "none", 
                                                        na.action = "ignore"))),
    ava_qual_scale = as.numeric(scale(ava_qual)),
    ava_qual_logit_scale = as.numeric(scale(adjust_prop(ava_qual, 
                                                        trans = "emp", 
                                                        nudge.method = "none", 
                                                        na.action = "ignore"))),
    sl_mean_obs_log_scale = as.numeric(scale(log(sl_mean_obs))),
    cat_pre_wt_log_scale_sq = as.numeric(scale(log(cat_pre_wt)^2)),
    RGR_scale = as.numeric(scale(RGR)), 
    sl_mean_pred1 = shape_1 * scale_1,
    sl_mean_pred2 = shape_2 * scale_2,
    prop_explore_logit = qlogis(prop_explore),
    prop_explore_logit_scale = as.numeric(scale(prop_explore_logit))
  )

# Subset of data.frame for more detailed analyses
d3 <- d %>% 
  filter(
    var_trt != "constant"
  ) %>% 
  filter(
    sl_mean_pred2 < 2000) %>% 
  filter_at(vars(contains("shape_"), contains("scale_[0-9]")), 
            all_vars(. > 0)) %>% 
  filter(
    s1.less_toxic.se < 10 &
      s2.less_toxic.se < 10
  ) %>% 
  mutate(
    var_high = ifelse(var_trt == "high_var", 1, 0),
    beta_red = ifelse(as.character(beta) == 5, 1, 0),
    beta_white = ifelse(as.character(beta) == 0, 1, 0),
    beta_red_cat_pre_wt_log_scale = beta_red * cat_pre_wt_log_scale,
    beta_white_cat_pre_wt_log_scale = beta_white * cat_pre_wt_log_scale,
    var_high_cat_pre_wt_log_scale = var_high * cat_pre_wt_log_scale,
    beta_numeric_scale = as.numeric(scale(as.numeric(beta)))
  )





o <- read_csv("cleaned_data/event_derivative_arrestment.csv", progress = FALSE) %>% 
  rename_all(.funs = function(x){
    o <- gsub("state","s",gsub("estimate", "est", gsub(":","\\.", x)))
    vapply(o, function(z){
      z <- str_split_1(z, pattern = "__")
      z2 <- z
      z2[1] <- z[length(z)]
      z2[length(z)] <- z[1]
      paste0(z2, collapse = ".")
    }, FUN.VALUE = character(1))
  }) %>% 
  left_join(
    read_csv("cleaned_data/event_derivative_select.csv", progress = FALSE) %>% 
      rename_all(.funs = function(x){
        o <- gsub("state","s",gsub("estimate", "est", gsub(":","\\.", x)))
        vapply(o, function(z){
          z <- str_split_1(z, pattern = "__")
          z2 <- z
          z2[1] <- z[length(z)]
          z2[length(z)] <- z[1]
          paste0(z2, collapse = ".")
        }, FUN.VALUE = character(1))
      }) %>% 
      mutate(
        sl_mean_pred1 = shape_1 * scale_1,
        sl_mean_pred2 = shape_2 * scale_2
      ) %>% 
      filter(
        sl_mean_pred2 < 2000) %>% 
      filter_at(vars(contains("shape_"), contains("scale_[0-9]")), 
                all_vars(. > 0)) %>% 
      filter(
        s1.less_toxic.se < 10 &
          s2.less_toxic.se < 10
      ) %>% 
      dplyr::select(rep_id, contains("less_toxic")),
    by = "rep_id"
  )

write_csv(o, "cleaned_data/event_derivative_master.csv")