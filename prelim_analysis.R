source("helper_functions/init_analysis.R")

ref_data <- ref_data %>% filter(error == 0)

names(ref_data)
ref_data %>% 
  filter(cat_post_wt > 0.6) %>% 
  select(rep_id, pupated, eclosed) %>% View()


ref_data %>% 
  filter(mean_trt_numeric == 1) %>% 
  ggplot(aes(x = sex, y = shape * scale / 1000 * 12, color = sex, fill = sex)) + 
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6)) + 
  geom_pointrange(stat = "summary", position = position_dodge(0.6), color = "black") + 
  theme_bw(base_size = 15) + 
  scale_y_continuous(trans = "log2") + 
  labs(y = "Mean step length (cm)")


ref_data %>% 
  ggplot(aes(x = sex, y = toxic, color = sex, fill = sex)) + 
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6)) + 
  geom_pointrange(stat = "summary", position = position_dodge(0.6), color = "black") + 
  theme_bw(base_size = 15) + 
  labs(y = "Preference for toxic diet")



ref_data %>% 
  ggplot(aes(x = var_trt, y = RGR, color = beta, fill = beta)) + 
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6)) + 
  geom_pointrange(stat = "summary", position = position_dodge(0.6), color = "black")



glmmTMB(
  RGR ~ 
    (var_trt + beta) * scale(log(cat_pre_wt)) +
    (1|session_id),
  data = ref_data
) %>% 
  plot_model(type = "eff", terms = c("cat_pre_wt[all]", "beta"), show.data = TRUE) + 
  scale_x_continuous(trans = "log10")











