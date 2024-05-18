source("spat1f/init_analysis.R")

d <- ref_data %>% 
  filter(mean_trt == "1 mg/g") %>%  # Filter only trials with the correct concentration
  mutate(cat_pre_wt_log_scale = as.numeric(scale(log(cat_pre_wt))),
         cat_pre_wt_log = log(cat_pre_wt))

# Simple RGR model as before, but dropping beta 
m <- glmmTMB(
  RGR ~ 
    var_trt * cat_pre_wt_log_scale + I(cat_pre_wt_log_scale^2) + (1|session_id),
  data = d
); summary(m)

car::Anova(m, type = "III") 

emmeans::emmeans(m, 
                 trt.vs.ctrl ~ var_trt | cat_pre_wt_log_scale, 
                 var = "var_trt", 
                 at = list(cat_pre_wt_log_scale = c(-2, 2)))

# Same thing here but with time to pupation
m2.1 <- glmmTMB(
  time_to_pupation ~ 
    var_trt * cat_pre_wt_log_scale + I(cat_pre_wt_log_scale^2) + (1|session_id),
  data = d, 
  family = Gamma(link = "log"),
); summary(m2.1)

car::Anova(m2.1, type = "III")

# Drop interaction
m2 <- glmmTMB(
  time_to_pupation ~ 
    var_trt + cat_pre_wt_log_scale + I(cat_pre_wt_log_scale^2) + (1|session_id),
  data = d, 
  family = Gamma(link = "log"),
); summary(m2)

car::Anova(m2, type = "II")

emmeans::emmeans(m2, 
                 trt.vs.ctrl ~ var_trt | cat_pre_wt_log_scale, 
                 var = "var_trt", 
                 at = list(cat_pre_wt_log_scale = c(-2, 2)))

plot_model(type = "eff", terms = c("cat_pre_wt[all]", "var_trt"))







# Manuscript plot
g1 <- marginal_effects(m, terms = c("cat_pre_wt_log_scale","var_trt")) %>% 
  ggplot(aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = yhat)) + 
  geom_ribbon(aes(ymax = upper, ymin = lower, fill = var_trt), alpha = 0.2) + 
  geom_line(aes(color = var_trt), linewidth = 2) + 
  geom_point(
    data = m$frame,
    aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = RGR, color = var_trt),
    size = 3, alpha = 0.8
  ) + 
  theme_bw(base_size = 15) + 
  scale_x_continuous(trans = "log10", labels = fancy_scientific) +
  theme(legend.position = "top") +
  scale_color_discrete(type = c("navy","steelblue","skyblue"), 
                       aesthetics = c("color","fill"), 
                       label = c("High", "Low","Constant")) + 
  labs(x = "Pre-weight (g)", y = expression(RGR~(hour^-1)),
       color = "Variation", fill = "Variation");g1

g2 <- marginal_effects(m2, terms = c("var_trt")) %>% 
  mutate(
    var_trt = factor(as.character(var_trt), levels = c("constant", "low_var","high_var"))
  ) %>% 
  ggplot(aes(x = (var_trt), y = yhat)) + 
  geom_point(
    data = m2$frame,
    aes(x = factor(as.character(var_trt), levels = c("constant", "low_var","high_var")), 
        y = time_to_pupation, color = var_trt),
    position = position_jitter(width = 0.2),
    size = 3, alpha = 0.8
  )  + 
  geom_pointrange(
    aes(ymax = upper, ymin = lower), 
    color = "black",
    size = 1, linewidth = 2, shape = 1,
  ) +
  theme_bw(base_size = 15) + 
  scale_y_continuous(trans = "log2") +
  theme(legend.position = "top") +
  scale_x_discrete(labels = c("Constant","Low","High")) + 
  scale_color_discrete(type = rev(c("steelblue","navy","skyblue")), 
                       aesthetics = c("color"), 
                       label = c("Constant","High","Low")) + 
  labs(x = "Variation", y = "Time to pupation (days)",
       color = "Variation"); g2


g_bind <- ggarrange(g1, g2, common.legend = TRUE, labels = "AUTO", label.x = 0.05)
ggsave("graphs/constant_performance.png",g_bind, dpi = 600, width = 8, height = 5)



