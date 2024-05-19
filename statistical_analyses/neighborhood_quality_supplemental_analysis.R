source("spat1f/init_analysis.R")

# Some quick cleaning
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
    select_scale = as.numeric(scale(less_toxic)),
    ava_qual_scale = as.numeric(scale(ava_qual)),
    ava_qual_logit_scale = as.numeric(scale(adjust_prop(ava_qual, 
                                                        trans = "emp", 
                                                        nudge.method = "none", 
                                                        na.action = "ignore"))),
    sl_mean_obs_log_scale = as.numeric(scale(log(sl_mean_obs))),
    cat_pre_wt_log_scale_sq = as.numeric(scale(log(cat_pre_wt)^2)),
    RGR_scale = as.numeric(scale(RGR)), 
    sl_mean_pred = shape * scale
  )


library(mgcv)

ava_beta_m <- glmmTMB(
  log(shape) ~ 
    var_trt + cat_pre_wt_log_scale + beta_numeric_scale + 
    (1|session_id),
  family = gaussian(),
  data = d2,
); summary(ava_beta_m)


ID <- 81

fetch_events(ID) %>% 
  clean_events() %$%
  move_seq(head_x, head_y) %>% 
  filter(!is.na(r)) %>% 
  add_random_steps(n = 100L, # Simulate random available steps
                   sl_distr = fit_gamma(.$r), # Fit gamma step length dist
                   ta_distr = fit_genvonmises(.$theta_rel) # Generalized von Mises turn angle dist
  ) %>%
  flag_invalid_steps(remove = TRUE) %>%
  mutate(
    toxic = read_value(x2, y2, c(1000, 1000),
                       ref_img = fetch_trt_spec(ID, 
                                                .ref_data = ref_data, 
                                                quiet = TRUE)), 
    less_toxic = 1 - toxic
  ) %>%
  append_estimators(na_as_zero = TRUE) %>% 
  issf(
    case ~
      less_toxic +
      (cos_theta_pi + cos_2theta) +  
      (sl + logsl) + 
      strata(step_id),
    data = .) -> issf_fit



spat1f::iterate_random_steps(
  spat1f::creat_start(500, 500, 0),
  issf_fit, n = 3000
) %>% 
  plot_track(x1, y1, type = "track")




select_m <- glmmTMB(
  less_toxic ~ 
    var_trt + as.numeric(beta) * cat_pre_wt_log_scale + 
    (1|session_id),
  family = gaussian(),
  data = d2,
); summary(select_m)

ava_beta_m <- glmmTMB(
  ava_qual ~ 
    var_trt + (less_toxic + cat_pre_wt_log_scale) * as.numeric(beta) + 
    (1|session_id),
  data = d2,
  family = beta_family()
); summary(ava_beta_m)



# Resource Selection Strength (RSS) coefficient model
select_m <- glmmTMB(
  less_toxic ~ 
    var_trt + beta * cat_pre_wt_log_scale + 
    (1|session_id),
  family = gaussian(),
  data = d,
); summary(select_m)

# Average neighborhood quality model. 
ava_beta_m <- glmmTMB(
  ava_qual ~ 
    var_trt + (cat_pre_wt_log_scale + less_toxic) * beta_numeric_scale + 
    (1|session_id),
  family = beta_family(link = "logit"),
  data = d2,
  
); summary(ava_beta_m)


# Post hoc contrast and means 
emmeans::emmeans(ava_beta_m, 
                 pairwise ~ beta | cat_pre_wt_log_scale, 
                 var = "beta", 
                 at = list(cat_pre_wt_log_scale = c(-2, 2)))

emmeans::emtrends(ava_beta_m, var = "cat_pre_wt_log_scale", specs = "beta")

# Backtransform to proportion
plogis(0.6432)
plogis(0.2066)
plogis(-0.0141)


# Plots
g1 <- marginal_effects(ava_beta_m, terms = c("cat_pre_wt_log_scale","beta")) %>% 
  ggplot(aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = yhat)) +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed", color = "black", size = 1.5) + 
  geom_ribbon(aes(ymax = upper, ymin = lower, fill = beta), alpha = 0.2) + 
  geom_line(aes(color = beta), linewidth = 2) + 
  geom_point(
    data = ava_beta_m$frame,
    aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = ava_qual, color = beta),
    size = 3, alpha = 0.8
  ) + 
  theme_bw(base_size = 15) + 
  scale_x_continuous(trans = "log10", labels = fancy_scientific) +
  theme(legend.position = "top") +
  scale_color_brewer(type = "qual", aesthetics = c("fill", "color")) + 
  labs(x = "Pre-weight (g)", y = "Neighborhood diet quality (Prop.)",
       color = expression(beta), fill = expression(beta))


g2 <- marginal_effects(select_m, terms = c("cat_pre_wt_log_scale","beta")) %>% 
  ggplot(aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = exp(yhat))) +
  geom_hline(aes(yintercept = 1), linetype = "dashed", color = "black", size = 1.5) + 
  geom_ribbon(aes(ymax = exp(upper), ymin = exp(lower), fill = beta), alpha = 0.2) + 
  geom_line(aes(color = beta), linewidth = 2) + 
  geom_point(
    data = select_m$frame,
    aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = exp(less_toxic), color = beta),
    size = 3, alpha = 0.8
  ) + 
  theme_bw(base_size = 15) + 
  scale_x_continuous(trans = "log10", labels = fancy_scientific) +
  theme(legend.position = "top") +
  scale_color_brewer(type = "qual", aesthetics = c("fill", "color")) + 
  labs(x = "Pre-weight (g)", y = "Relative selection towards \nless toxic diet (Odds)",
       color = expression(beta), fill = expression(beta)) + 
  scale_y_continuous(trans = "log10")

g_bind <- ggarrange(g1, g2, common.legend = TRUE, legend = "top",labels = "AUTO")

ggsave(filename = "graphs/manuscript1_figures/selection_nighborhood_means.png", dpi = 600, height = 4.5, width = 8)






