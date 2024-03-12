eff_exp <- function(est,se){
  signif(100 *(exp(c(est, est - se * 1.96, est + se*1.96)) - 1 ), 2)
}

eff <- function(est,se){
  signif((c(est, est - se * 1.96, est + se*1.96)), 2)
}

eff(-0.21,0.10)

eff_exp(0.21, 0.0915)


glmmTMB(RGR_scale ~ var_toxic_12_scale + cat_pre_wt_log_scale * (var_trt + beta) + I(cat_pre_wt_log_scale^2) + (1|session_id), data = d2) %>% summary()#car::Anova(type = "II")# plot_model(type = "eff", terms = c("cat_pre_wt_log_scale", "var_toxic_12"))


m <- glmmTMB(mean_toxic_conc_scale  ~ 
          toxic_select_scale + scale(on_toxic) * cat_pre_wt_log_scale + 
            cat_pre_wt_log_scale * (var_trt) + beta + 
          (1|session_id), data = d2);summary(m) #car::Anova(type = "II")# 

plot_model(m, type = "eff", terms = c("on_toxic"), show.data = TRUE)

check_collinearity(m)





s <- 0.005786767



0.0027 / s



d3 <- d %>% 
  filter(var_trt == "constant")


m <- glmmTMB(
  sl_mean_obs ~ 
    mean_trt_numeric * cat_pre_wt_log_scale + (1|session_id), 
  data = d3,
  family = Gamma(link = "log")
);summary(m) 
plot_model(m, type = "eff", terms = c("cat_pre_wt_log_scale", "mean_trt_numeric"))


contrast_by_pre_wt(m, "mean_trt_numeric")

emmeans::emmeans(m, 
                 pairwise ~ mean_trt_numeric | cat_pre_wt_log_scale, 
                 var = "mean_trt_numeric", 
                 at = list(cat_pre_wt_log_scale = c(-2, 2), 
                           mean_trt_numeric = c(0, 0.5)))






















ID <- 93

d <- fetch_events(ID) %>%
  clean_events(ref_data = ref_data) %>%
  insert_gaps() %>% 
  mutate(
    head_x = ifelse(score >=0.9, head_x, NA),
    head_y = ifelse(score >=0.9, head_y, NA)
  ) %$%
  move_seq(head_x, head_y, r_thresh = 0, inherit.theta = FALSE) %>%
  filter(!is.na(r)) %>%
    add_random_steps(n = 100L,
                     sl_distr = fit_gamma(.$r),
                     ta_distr = fit_genvonmises(.$theta_rel)
    ) %>%
    flag_invalid_steps(remove = TRUE) %>%
    mutate(
      # start_toxic = ifelse(
      #   read_value(x1, y1, c(1000, 1000),
      #              ref_img = fetch_trt_spec(ID, .ref_data = ref_data, quiet = TRUE)) == 1,
      #   "yes",
      #   "no"
      # ),
      toxic = read_value(x2, y2, c(1000, 1000),
                         ref_img = fetch_trt_spec(ID, .ref_data = ref_data, quiet = TRUE))
    ) %>%
    append_estimators(na_as_zero = TRUE)

out <- issf(
  case ~
    toxic + 
    (cos_theta_pi + cos_2theta) +
    (sl + logsl) +
    #logsl + logslsq + 
    strata(step_id),
  data = d,
  shape_estimator = c("logsl"),
  scale_estimator = c("sl"),
  kappa1_estimator = "cos_theta_pi",
  kappa2_estimator = "cos_2theta"
)
out

reload()


d %>% filter(case) %>% .$r %>% loghist(nclass = 20)
d %>% filter(case) %>% .$r %>% log() %>% hist(nclass = 20)


d %>% filter(case) %>% .$r %>% log %>% hist(nclass = 50)
rdist(out$sl_updated, 10000)[[1]] %>% log() %>% hist(nclass = 50)


creat_start(500, 500, 0) %>% 
  iterate_random_steps(issf_fit = out, n = 1000) %>% 
  plot_track(x2, y2)

creat_start(500, 500, 0) %>% 
  iterate_random_steps(issf_fit = out, n = 3000, paired = TRUE, use_observed = TRUE) %>% 
  plot_track(x2, y2)
