source("spat1f/init.R")
compute_delta2 <- function(hypothesis, res){
  temp <- hypothesis %>% 
    left_join(res, by = c("k", "rss","b"), relationship = "many-to-many") %>% 
    group_by(cat_size, id) %>% 
    summarise(
      delta = diff(on_toxic)
    )
  out1 <- temp %>% 
    mutate(cat_size = ifelse(cat_size == "s", "Small", "Large"))
  return(out1)
}


grid <- list(
  "No arrestment\nNo immigration" = data.frame("k" = c(1, 1, 1, 1), 
                  "rss" = c(1, 1, 1, 1),
                  "cat_size" = c("s", "s", "l", "l"),
                  "b" = c(-5, 5, -5, 5)),
  "Size depend. arrestment\nNo immigration" = data.frame("k" = c(0.25, 0.25, 4, 4), 
             "rss" = c(1,1, 1, 1),
             "cat_size" = c("s", "s", "l", "l"),
             "b" = c(-5, 5, -5, 5)),
  "No arrestment\nSize depend. immigration" = data.frame("k" = c(1, 1, 1, 1), 
             "rss" = c(0.8,0.8, 1.25, 1.25),
             "cat_size" = c("s", "s", "l", "l"),
             "b" = c(-5, 5, -5, 5)),
  "Size depend. arrestment\nSize depend. immigration" = data.frame("k" = c(0.25, 0.25, 4, 4), 
             "rss" = c(0.8,0.8, 1.25, 1.25),
             "cat_size" = c("s", "s", "l", "l"),
             "b" = c(-5, 5, -5, 5)),
  "Size & β depend. arrestment\nNo immigration" = data.frame("k" = c(1, 0.25, 1, 4), 
             "rss" = c(1,1, 1, 1),
             "cat_size" = c("s", "s", "l", "l"),
             "b" = c(-5, 5, -5, 5)),
  "Size & β depend. arrestment\nSize depend. immigration" = data.frame("k" = c(1, 0.25, 1, 4), 
             "rss" = c(0.8,0.8, 1.25, 1.25),
             "cat_size" = c("s", "s", "l", "l"),
             "b" = c(-5, 5, -5, 5))
)


grid_d <- lapply(seq_along(grid), function(i){
  compute_delta2(grid[[i]], z) %>% 
    cbind("model" = names(grid)[i])
}) %>% 
  do.call("rbind", .)



g <- epred_draws(m, expand.grid("cat_pre_wt_log_scale" = c(-2, 2), 
                           "beta" = c(-5, 5),
                           "var_trt" = c("low_var","high_var"),
                           "session_id" = "foo")) %>% 
  group_by(cat_pre_wt_log_scale, beta, draw) %>% 
  summarise(val = mean(val)) %>% 
  group_by(
    cat_pre_wt_log_scale, draw
  ) %>% 
  summarise(
    delta = diff(val)
  ) %>% 
  mutate(
    cat_size = factor(ifelse(cat_pre_wt_log_scale < 0, "Small", "Large"), levels = c("Small", "Large"))
  ) %>% 
  ggplot(aes(y = delta)) + 
  tidybayes::stat_interval(.width = c(0.99, .95, 0.66, 0.5, 0.33, 0.05, 0.01), size = 1000, alpha = 0) + 
  geom_hline(yintercept = 0, color = "violetred", size = 1.5, linetype = "dashed") +
  stat_halfeye(
    data = grid_d %>% 
      mutate(
        model = factor(model, level = rev(names(grid))),
        cat_size = factor(cat_size, levels = c("Small", "Large"))
      ),
    aes(x = model, y = delta),
    alpha = 0
  ) + 
  scale_color_brewer() + 
  facet_grid(~cat_size, scales = "free") +
  coord_flip() + 
  theme_bw(base_size = 15) + 
  theme(legend.position = "none") + 
  guides(color = guide_legend(override.aes = list(size = 1)))  + 
  labs(x = "Model", y = "Marginal effect of clusteredness\n(Proportion time on toxic diet)");g


ggsave("graphs/table_1_model_comparison_as_figure_empty2.png", dpi = 600, width = 7, height = 5)






ggarrange(g1 + 
            scale_color_brewer(label = c("Dispersed", "Random", "Clustered"), 
                               aesthetics = c("fill", "color"), 
                               type = "qual") + 
            labs(color = expression(Clusteredness~(beta)),
                 fill = expression(Clusteredness~(beta))), 
          g3, common.legend = TRUE)



syn_spec(1000, beta = 3, threshold = FALSE, plot = FALSE) %>% 
  plot(axes = FALSE)


imlist(
  syn_spec(100, beta = -5, threshold = FALSE, plot = FALSE),
  syn_spec(100, beta = 0, threshold = FALSE, plot = FALSE),
  syn_spec(100, beta = 1, threshold = FALSE, plot = FALSE),
  syn_spec(100, beta = 3, threshold = FALSE, plot = FALSE),
  syn_spec(100, beta = 5, threshold = FALSE, plot = FALSE),
  syn_spec(100, beta = 10, threshold = FALSE, plot = FALSE)
) %>% 
  plot.imlist(interpolate = FALSE, main.panel = c(""), axes = FALSE)


syn_spec(100, beta = 1, threshold = FALSE, plot = FALSE, invert = FALSE) %>% 
  fft_img_plot(type = "mag", trans = function(x){
    x + 0.001
  })



x <- syn_spec(1000, beta = 3, threshold = FALSE, plot = FALSE, invert = TRUE) %>%
  .[1, ] %>% 
  fft() %>% 
  Mod() %>% 
  .[seq_len(floor(length(.)/2))]


data.frame("f" = seq_along(x), "mag" = x) %>%
  .[-1, ] %>% 
  ggplot(aes(x = f, y = mag)) + 
  geom_line(size = 0.6, color = "steelblue") + 
  scale_xy_log(axis = "y") + 
  scale_xy_log(axis = "x") + 
  labs(x = "Frequency", y = "Variance") + 
  theme_bw(base_size = 25)

imlist(
  syn_spec(12, beta = -50, threshold = TRUE, plot = FALSE),
  syn_spec(12, beta = 0, threshold = TRUE, plot = FALSE),
  syn_spec(12, beta = 50, threshold = TRUE, plot = FALSE)
) %>% 
  plot.imlist(interpolate = FALSE, main.panel = c(""), axes = FALSE)






rgr_m2 <- glmmTMB(
  RGR ~ 
    cat_pre_wt_log_scale +
    I(cat_pre_wt_log_scale^2) +  # residual plot show significant non-linearity
    (1|session_id),
  data = d %>% 
    filter(
      var_trt != "constant"  & !is.na(RGR)
    )
); summary(rgr_m2)



rgr_m3 <- glmmTMB(
   resid ~ 
    cat_pre_wt_log_scale * (var_trt + beta) + (1|session_id),
  data = d %>% 
    filter(
      var_trt != "constant" & !is.na(RGR)
    ) %>% 
    mutate(
      resid = residuals(rgr_m2) / sigma(rgr_m2)
    )
); summary(rgr_m3)



marginal_effects(rgr_m3, terms = c("cat_pre_wt_log_scale","beta")) %>% 
  ggplot(aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = yhat)) + 
  geom_ribbon(aes(ymax = upper, ymin = lower, fill = beta), alpha = 0.2) + 
  geom_line(aes(color = beta), linewidth = 2) + 
  geom_point(
    data = rgr_m3$frame,
    aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = resid, color = beta),
    size = 3, alpha = 0.8
  ) + 
  theme_bw(base_size = 15) + 
  theme(legend.position = "none") +
  scale_y_continuous(breaks = seq(-3, 3, by = 1)) + 
  scale_x_continuous(trans = "log10", labels = fancy_scientific) +
  scale_color_brewer(type = "qual", aesthetics = c("fill","color")) + 
  labs(x = "Pre-weight (g)", y = expression(Scaled~RGR~residual~(SD)),
       color = expression(beta), fill = expression(beta))







marginal_effects(m_select_s1, terms = c("cat_pre_wt_log_scale", "beta")) %>% 
  cbind(state = "1") %>% 
  ggplot(aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = exp(yhat), 
             color = state, linetype = beta)) + 
  geom_ribbon(aes(ymax = exp(upper), ymin = exp(lower), fill = state, color = NULL), alpha = 0.2) + 
  geom_line(linewidth = 2) + 
  geom_point(
    data = m_select_s1$frame %>% 
      cbind(state = "1"),
    aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = exp(s1.less_toxic_end.est), shape = beta),
    size = 3, alpha = 0.5
  ) +  
  theme_bw(base_size = 15) + 
  geom_hline(aes(yintercept = 1), color = "black", linetype = "dashed", linewidth = 1) + 
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0.3, 3)) +
  scale_y_continuous(trans = "log10", breaks = c(1, 0.3, 3)) + 
  scale_x_continuous(trans = "log10", labels = fancy_scientific) +
  scale_color_manual(values = rev(.getPalette(2)), 
                     aesthetics = c("color", "fill"), 
                     labels = c("Exploration", "Resting/feeding")) + 
  scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
  labs(x = "Pre-weight (g)", y = "Less toxic diet relative \nimmigration strength (odds ratio)")


marginal_effects(m_select_s2, terms = c("cat_pre_wt_log_scale", "beta")) %>%
  cbind(state = "2") %>%
  ggplot(aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = exp(yhat), 
             color = state, linetype = beta)) + 
  geom_ribbon(aes(ymax = exp(upper), ymin = exp(lower), fill = state, color = NULL), alpha = 0.2) + 
  geom_line(linewidth = 2) + 
  geom_point(
    data = m_select_s2$frame %>%
      cbind(state = "2"),
    aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = exp(s2.less_toxic_end.est), shape = beta),
    size = 3, alpha = 0.5
  ) +
  theme_bw(base_size = 15) + 
  geom_hline(aes(yintercept = 1), color = "black", linetype = "dashed", linewidth = 1) + 
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0.3, 3)) +
  scale_y_continuous(trans = "log10", breaks = c(1, 0.3, 3)) + 
  scale_x_continuous(trans = "log10", labels = fancy_scientific) +
  scale_color_manual(values = rev(.getPalette(2)[1]), 
                     aesthetics = c("color", "fill"), 
                     labels = c("Exploration", "Resting/feeding")) + 
  scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
  labs(x = "Pre-weight (g)", y = "Less toxic diet relative \nimmigration strength (odds ratio)")






marginal_effects(m_arrest_s1, terms = c("cat_pre_wt_log_scale", "beta")) %>% 
  cbind(state = "1") %>% 
  # rbind(
  #   marginal_effects(m_arrest_s2, terms = c("cat_pre_wt_log_scale", "beta")) %>% 
  #     cbind(state = "2")
  # ) %>% 
  ggplot(aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = (yhat), color = state)) + 
  geom_ribbon(aes(ymax = (upper), ymin = lower, fill = state, linetype = beta), alpha = 0.2) + 
  geom_line(aes(linetype = beta), linewidth = 2) + 
  geom_point(
    data = m_arrest_s1$frame %>% 
      cbind(state = "1"),
    aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = exp(`log(k1)`), shape = beta),
    size = 3, alpha = 0.3
  ) +  
  # geom_point(
  #   data = m_arrest_s2$frame %>% 
  #     cbind(state = "2"),
  #   aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = exp(`log(k2)`), shape = beta),
  #   size = 3, alpha = 0.3
  # ) + 
  theme_bw(base_size = 15) + 
  geom_hline(aes(yintercept = 1), color = "black", linetype = "dashed", linewidth = 1) + 
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0.1, 10)) +
  scale_y_continuous(trans = "log10", breaks = c(1, 0.3, 3, 0.1, 10)) + 
  scale_x_continuous(trans = "log10", labels = fancy_scientific) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
  scale_color_manual(values = rev(.getPalette(2)), 
                     aesthetics = c("color", "fill"), 
                     labels = c("Exploration", "Resting/feeding")) + 
  labs(x = "Pre-weight (g)", y = "Less toxic diet arresetment \nstrength (ratio)", 
       linetype = expression(beta), shape = expression(beta)) + 
  guides(color = guide_legend())



marginal_effects(m_arrest_s2, terms = c("cat_pre_wt_log_scale", "beta")) %>% 
       cbind(state = "2") %>% 
  ggplot(aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = (yhat), color = state)) + 
  geom_ribbon(aes(ymax = (upper), ymin = lower, fill = state, linetype = beta), alpha = 0.2) + 
  geom_line(aes(linetype = beta), linewidth = 2) + 
  geom_point(
    data = m_arrest_s2$frame %>%
      cbind(state = "2"),
    aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = exp(`log(k2)`), shape = beta),
    size = 3, alpha = 0.3
  ) +
  theme_bw(base_size = 15) + 
  geom_hline(aes(yintercept = 1), color = "black", linetype = "dashed", linewidth = 1) + 
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0.1, 10)) +
  scale_y_continuous(trans = "log10", breaks = c(1, 0.3, 3, 0.1, 10)) + 
  scale_x_continuous(trans = "log10", labels = fancy_scientific) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
  scale_color_manual(values = rev(.getPalette(2)[1]), 
                     aesthetics = c("color", "fill"), 
                     labels = c("Exploration", "Resting/feeding")) + 
  labs(x = "Pre-weight (g)", y = "Less toxic diet arresetment \nstrength (ratio)", 
       linetype = expression(beta), shape = expression(beta)) + 
  guides(color = guide_legend())













x <- iterate_random_steps_states(
  list(
    "ta_updated" = list(
      make_unif(),
      make_unif(),
      make_genvonmises(kappa1 = 0.4, kappa2 = 0.3),
      make_genvonmises(kappa1 = 0.4, kappa2 = 0.3)
    ),
    "sl_updated" = list(
      make_gamma(shape = 0.8, 
                 scale = 100), # More toxic diet
      make_gamma(shape = 0.8, 
                 scale = 100),  # Less toxic diet
      make_gamma(shape = 0.8, 
                 scale = 5), # More toxic diet
      make_gamma(shape = 0.8, 
                 scale = 5)  # Less toxic diet
    )
  ), 
  start = make_start2(0,500,500, 1), 
  n = 2000, 
  ref_grid = syn_spec(n = 12, beta = 5, plot = FALSE), 
  max_xy = c(1000, 1000), 
  transition_mat = matrix(
    c(0.9, 0.1, 
      0.04, 0.96), 
    byrow = TRUE,
    ncol = 2,
    nrow = 2
  ),
  rss_coef = 0.4)


plot_track_overlay(x, plot_elements = "track")

plot_d <- attr(x, "ref_grid") %>% 
  as.matrix() %>% 
  melt()

g <- plot_d %>% 
  ggplot(aes(x = dim1, y = dim2)) + 
  geom_tile(
    fill = ifelse(plot_d$val == 1, "white", "grey"),
    alpha = 0.5)

g <- g +
  theme_void() + 
  geom_point(data = data.frame(), aes(x = 0.5, y = 0.5), alpha = 0) + 
  geom_point(data = data.frame(), aes(x = 12 + 0.5, y = 12 + 0.5), alpha = 0) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) 


d <- x %>% 
  mutate(
    x = x / 1000 * 12 + 0.5,
    y = y / 1000 * 12 + 0.5,
    t = seq_along(x)
  )



for(i in c(2 * 1:1000)){
  cat(sprintf("Processing %s of %s\r", i, 2000))
  
  xy <- d %>% 
    filter(t < (i + 1)) %>% 
    .[nrow(.), c("x", "y")] %>% 
    unlist()
  
  g_temp <- g + 
    geom_path(
      data = d %>% 
        filter(t < (i + 1)),
      aes(x = x, y = y), 
      color = "black"
    ) + 
    ggimage::geom_image(
      data = data.frame(
        "img" = "graphs/methods_figure/cat.PNG",
        "x" = xy[1],
        "y" = xy[2]
      ), 
      aes(image = img, x = x, y = y), size = 0.1
      )
  
  
  jpeg(sprintf("invisible/test_video_overlay/image-%s.jpg",i/2), quality = 100, width = 500, height = 500)
  print(g_temp)
  dev.off()
}

imager::make.video(fps = 20, 
                   pattern = "image-%d.jpg", 
                   dname = abs_path("invisible/test_video_overlay"), 
                   fname = abs_path(paste0("simulated_track_select",".mp4")), 
                   verbose = TRUE)
