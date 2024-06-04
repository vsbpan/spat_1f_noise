library(tidybayes)
source("spat1f/init.R")



#### Performance plot #####


g_bind <- ggarrange(
  g2 + 
    labs(fill = "Variation", color = "Variation", x = "Pre-weight (g)"), 
  NULL,
  g1 + 
    labs(y = "", x = "Pre-weight (g)") + 
    theme(plot.background = element_rect(fill='transparent', colour = "transparent")),  
  g4 + 
    theme(legend.position = "none") + 
    scale_x_discrete(breaks = c("high_var","low_var"), label = c("high", "low")), 
  NULL,
  g3 + labs(y = "", x = "Pre-weight (g)")  + 
    theme(legend.position = "none") + 
    theme(plot.background = element_rect(fill='transparent', colour = "transparent")), 
  align = "v", 
  heights = c(1,0.8), 
  widths = c(1, -0.1, 1, 1,-0.1, 1),
  ncol = 3, nrow = 2, 
  labels = c("a", "", "b", "c", "", "d"), 
  label.x = 0.15, 
  label.y = c(0.85,0,0.85, 1,0, 1)
); g_bind
#ggsave("graphs/perfromance_plot.png",g_bind, dpi = 600, width = 6.5, height = 6.5, bg = "white")




#### Methods figure insets ####

source("spat1f/init_analysis.R")

# Example specs
fetch_trt_meta() %>% 
  mutate_if(is.numeric, .funs = function(x){
    ifelse(x == 1, 1, 0.5)
  }) %>% 
  mutate(beta = gsub("beta|_.*","",syn_id)) %>% 
  arrange(beta) %>% 
  group_by(beta) %>% 
  slice_sample(n = 2) %>% 
  trt_meta_as_list() %>% 
  lapply(
    seq_along(.),
    function(i,l){
      fn <- names(l)[i]
      jpeg(paste0("graphs/methods_figure/", fn, ".jpg"), width = 5, height = 5, units = "cm", res = 600)
      plot_image_guide(l[[i]], mar = c(0,0,0,0), col = rgb(1,1,1,0))
      dev.off()
    },
    l = .
  )


# Example raw image
img <- fetch_image(55, 100, type = "raw")
jpeg(paste0("graphs/", "processed_rep55_rank100_raw", ".jpg"), 
     width = 15, height = 15, units = "cm", res = 600)
plot(img)
dev.off()

# Exampled cropped image
img <- fetch_image(55, 100)
jpeg(paste0("graphs/", "processed_rep55_rank100_cropped", ".jpg"), 
     width = 15, height = 15, units = "cm", res = 600)
plot(img, axes = FALSE)
dev.off()

# Example predicted image
img <- fetch_image(55, 100)
jpeg(paste0("graphs/", "processed_rep55_rank100_prediction", ".jpg"), 
     width = 15, height = 15, units = "cm", res = 600)


fetch_trt_spec(55) %>% 
  flip_xy() %>% 
  resize(1000, 1000) %>% 
  imdraw(
    img, ., opacity = 0.1
  ) %>% 
  plot.cimg(axes = FALSE)
fetch_data_dict(55)[100] %>% 
  plot(add = TRUE)
dev.off()

# Example track
g <- fetch_trt_spec(55) %>% 
  flip_xy() %>% 
  as.matrix() %>% 
  melt() %>% 
  ggplot(aes(x = dim1, y = 12 - dim2)) + 
  geom_tile(
    aes(fill = as.factor(val)),
    alpha = 0.5) +
  geom_path(
    data = fetch_events(55) %>%
      clean_events(keep_sus = TRUE) %>%
      insert_gaps() %>% 
      mutate(
        head_x = head_x / 1000 * 12 + 0.5,
        head_y = 12-(head_y / 1000 * 12 + 0.5)
      ),
    aes(x = head_x, y = head_y)
  ) +
  theme_void() + 
  scale_fill_discrete(type = c("grey","white")) +
  theme(legend.position = "none") +
  geom_point(aes(x = 0.5, y = 0.5), alpha = 0) + # Force the edges to align
  geom_point(aes(x = 12.5, y = 12.5), alpha = 0) # Force the edges to align

# ggsave("graphs/rep55_tracks.jpg", g, width = 7.5, height = 8, dpi = 600)

# Example track segmented
g <- fetch_trt_spec(55) %>% 
  flip_xy() %>% 
  as.matrix() %>% 
  melt() %>% 
  ggplot(aes(x = dim1, y = 12 - dim2)) + 
  geom_tile(
    aes(fill = as.factor(val)),
    alpha = 0.5) +
  geom_path(
    data = fetch_events(55) %>%
      clean_events(keep_sus = TRUE) %>%
      mutate(
        head_x = head_x / 1000 * 12 + 0.5,
        head_y = 12-(head_y / 1000 * 12 + 0.5)
      ) %>% 
      .[-nrow(.),],
    aes(x = head_x, y = head_y, group = 1, 
        color = as.character(
          viterbi(
            fit_HMM(
              as.moveData(fetch_events(55) %>%
                            clean_events(keep_sus = TRUE) %$% 
                            move_seq(head_x, head_y)))))
    )
  ) +
  theme_void() + 
  scale_fill_discrete(type = c("grey","white")) +
  scale_color_discrete(type = rev(.getPalette(2))) +
  theme(legend.position = "none") +
  geom_point(aes(x = 0.5, y = 0.5), alpha = 0) + # Force the edges to align
  geom_point(aes(x = 12.5, y = 12.5), alpha = 0);g # Force the edges to align


# ggsave("graphs/rep55_tracks_segmented.jpg", g, width = 7.5, height = 8, dpi = 600)


#### Herbivory mask processing figure #####

img1 <- fast_load_image("graphs/methods_figure/processed_rep81__cam20_s0_rank1.jpg")
img2 <- fast_load_image("graphs/methods_figure/processed_rep81__cam20_s441162_rank1225.jpg")
mask <- (fetch_data_dict(81)[1225] %>% get_mask())[[1]]

hue_diff <- (HUE(img1) - HUE(img2)) 
hue_diff2 <- hue_diff
hue_diff2[imager::grow(mask, x = 20)] <- quantile(hue_diff2, probs = 0.1)
spe_img <- hue_diff2 %>% imagerExtra::SPE(lamda = 0.001,range = c(0,1))
thr_img <- spe_img %>% threshold2(thr = "otsu")
img_final <- thr_img %>% shrink(10)


jpeg(paste0("graphs/methods_figure/", "herbivory_detection", ".jpg"), width = 15, height = 10, units = "cm", res = 600)
plot.imlist(imlist(img2, HUE(img2), hue_diff,
                   spe_img, as.cimg(thr_img), as.cimg(img_final)), 
            main.panel = c("Raw image", "HUE", "Delta HUE censored",
                           "SPE enhanced", "Otsu thresholded", "Eroded final"), 
            axes = FALSE, interpolate = FALSE)
dev.off()



#### Generalized von Mises Distribution plot #####


x <- fetch_events(81) %$%
  move_seq(head_x, head_y) %>% 
  .$theta_rel %>% 
  na.omit() %>% 
  as.vector()

grid <- seq(-pi, pi, len = 30)

data.frame(
  "x" = nearest_bin(x, grid)
) %>% 
  group_by(x) %>% 
  tally() %>% 
  mutate(
    p = n / sum(n)
  ) %>% 
  ggplot(
    aes(x = x, y = p / mean(diff(grid)))
  ) + 
  geom_col(
   fill = "steelblue", color = "navy" 
  ) + 
  geom_line(
    data = ddist(fit_genvonmises(x)),
    aes(x = x, y = density),
    size = 2, 
    color = "tomato"
  ) +
  geom_line(
    data = data.frame(
      "x" = grid, 
      "density" = as.vector(circular::dvonmises(x = grid, 
                                                mu = -pi, 
                                                kappa = amt::fit_distr(x-pi, "vonmises")$params$kappa))
    ),
    size = 2,
    aes(x = x, y = density),
    color = "black",
    linetype = "dashed"
  ) + 
  theme_bw(base_size = 15) + 
  scale_x_continuous(breaks = c(-pi, -pi/2, 0, pi/2, pi), 
                     labels = c(expression(-pi), expression(-pi/2), 0, expression(pi/2), expression(pi)) 
  ) + 
  labs(x = "Turn angles (radians)", y = "Density") -> g

ggsave("graphs/generalized_von_mises_fit.png",g, dpi = 600)






#### ISSF Params ####


g1 <- marginal_effects(m_select_s1, terms = c("cat_pre_wt_log_scale")) %>% 
  cbind(state = "1") %>% 
  rbind(marginal_effects(m_select_s2, terms = c("cat_pre_wt_log_scale")) %>% 
          cbind(state = "2")) %>% 
  ggplot(aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = exp(yhat), color = state)) + 
  geom_ribbon(aes(ymax = exp(upper), ymin = exp(lower), fill = state, color = NULL), alpha = 0.2) + 
  geom_line(linewidth = 2) + 
  geom_point(
    data = m_select_s1$frame %>% 
      cbind(state = "1"),
    aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = exp(s1.less_toxic.est)),
    size = 3, alpha = 0.5
  ) +  
  geom_point(
    data = m_select_s2$frame %>% 
      cbind(state = "2"),
    aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = exp(s2.less_toxic.est)),
    size = 3, alpha = 0.5
  ) + 
  theme_bw(base_size = 15) + 
  geom_hline(aes(yintercept = 1), color = "black", linetype = "dashed", linewidth = 1) + 
  theme(legend.position = "top") +
  scale_y_continuous(trans = "log10") + 
  scale_x_continuous(trans = "log10", labels = fancy_scientific) +
  scale_color_manual(values = rev(.getPalette(2)), 
                     aesthetics = c("color", "fill"), 
                     labels = c("Exploration", "Resting/feeding")) + 
  labs(x = "Pre-weight (g)", y = "Less toxic diet selection \nstrength (odds ratio)")


g2 <- marginal_effects(m_arrest_s1, terms = c("cat_pre_wt_log_scale", "beta")) %>% 
  cbind(state = "1") %>% 
  rbind(
    marginal_effects(m_arrest_s2, terms = c("cat_pre_wt_log_scale", "beta")) %>% 
      cbind(state = "2")
  ) %>% 
ggplot(aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = (yhat), color = state)) + 
  geom_ribbon(aes(ymax = (upper), ymin = lower, fill = state, linetype = beta), alpha = 0.2) + 
  geom_line(aes(linetype = beta), linewidth = 2) + 
  geom_point(
    data = m_arrest_s1$frame %>% 
      cbind(state = "1"),
    aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = exp(`log(k1)`), shape = beta),
    size = 3, alpha = 0.3
  ) +  
  geom_point(
    data = m_arrest_s2$frame %>% 
      cbind(state = "2"),
    aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = exp(`log(k2)`), shape = beta),
    size = 3, alpha = 0.3
  ) + 
  theme_bw(base_size = 15) + 
  geom_hline(aes(yintercept = 1), color = "black", linetype = "dashed", linewidth = 1) + 
  theme(legend.position = "top", legend.box = "verticle", legend.margin=margin(-1,-1,-1,-1)) +
  scale_y_continuous(trans = "log10") + 
  scale_x_continuous(trans = "log10", labels = fancy_scientific) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
  scale_color_manual(values = rev(.getPalette(2)), 
                     aesthetics = c("color", "fill"), 
                     labels = c("Exploration", "Resting/feeding")) + 
  labs(x = "Pre-weight (g)", y = "Less toxic diet arresetment \nstrength (ratio)", 
       linetype = expression(beta), shape = expression(beta)) + 
  guides(color = guide_legend())



g3 <- marginal_effects(m_sl_pred_s1, terms = c("cat_pre_wt_log_scale")) %>% 
  cbind(state = "1") %>% 
  rbind(
    marginal_effects(m_sl_pred_s2, terms = c("cat_pre_wt_log_scale")) %>% 
      cbind(state = "2")
  ) %>% 
  ggplot(aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = (yhat) / 1000 * 12, color = state)) + 
  geom_ribbon(aes(ymax = (upper)/ 1000 * 12, ymin = lower/ 1000 * 12, color = NULL, fill = state), 
              alpha = 0.2) + 
  geom_line(linewidth = 2) + 
  geom_point(
    data = m_sl_pred_s1$frame %>% 
      cbind(state = "1"),
    aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), 
        y = exp(`log(sl_mean_pred1)`) / 1000 * 12),
    size = 3, alpha = 0.5
  ) +  
  geom_point(
    data = m_sl_pred_s2$frame %>% 
      cbind(state = "2"),
    aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), 
        y = exp(`log(sl_mean_pred3)`) / 1000 * 12),
    size = 3, alpha = 0.5
  ) +  
  theme_bw(base_size = 15) + 
  theme(legend.position = "top") +
  scale_y_continuous(trans = "log10") + 
  scale_x_continuous(trans = "log10", labels = fancy_scientific) +
  scale_color_brewer(type = "qual", aesthetics = c("color","fill")) + 
  scale_color_manual(values = rev(.getPalette(2)), 
                     aesthetics = c("color", "fill"), 
                     labels = c("Exploration", "Resting/feeding")) + 
  labs(x = "Pre-weight (g)", y = "Selection free mean step \nlength on toxic diet (cm)")






foo <- function(d, draw_dist, label){
  n <- length(draw_dist)
  den_data <- vector(mode = "list", length = n)
  
  for (i in seq_len(n)){
    dist_name <- draw_dist[[i]]$name
    
    den <- do.call(
      paste0("d",dist_name),
      c(
        list(d$x), 
        draw_dist[[i]]$params
      )
    )
    
    den[!is.finite(den)] <- NA_real_
    
    p <- den * d$x 
    
    den_data[[i]] <- data.frame(
      "x" = d$x, 
      "p" = p,
      "state" = label[i]
    )
  }
  do.call("rbind", den_data)
}


# Arbitrary caterpillar -- same one as in the demo

ID <- 55
temp <- fetch_events(ID) %>% 
  clean_events(ref_data = ref_data) %$%
  move_seq(head_x, head_y) %>% 
  mutate(toxic = read_value(x2, y2, ref_img = fetch_trt_spec(ID, .ref_data = ref_data)))
temp$state <- viterbi(fit_HMM(as.moveData(temp)))
temp <- temp %>% filter(!is.na(r))
x1 <- as.numeric(na.omit(temp$r[temp$state == 1 & temp$toxic == 0]))
x2 <- as.numeric(na.omit(temp$r[temp$state == 2 & temp$toxic == 0]))
x3 <- as.numeric(na.omit(temp$r[temp$state == 1 & temp$toxic == 1]))
x4 <- as.numeric(na.omit(temp$r[temp$state == 2 & temp$toxic == 1]))
hist_d <- c(x1, x2, x3, x4) %>% 
  log() %>% 
  hist(plot = FALSE, nclass = 50)
hist_d$x <- exp(hist_d$mids)
hist_d$pred_dens <- foo(hist_d,
                        list(
                          fit_gamma(x1),
                          fit_gamma(x2),
                          fit_gamma(x3),
                          fit_gamma(x4)
                        ), 
                        c("1_less", "2_less","1_more","2_more")) %>% 
  mutate(
    p = case_when(
      state == "1_less" ~ p * mean(temp$state == 1 & temp$toxic == 0, na.rm = TRUE),
      state == "2_less" ~ p * mean(temp$state == 2 & temp$toxic == 0, na.rm = TRUE),
      state == "1_more" ~ p * mean(temp$state == 1 & temp$toxic == 1, na.rm = TRUE),
      state == "2_more" ~ p * mean(temp$state == 2 & temp$toxic == 1, na.rm = TRUE)
    )
  )
hist_d$pred_dens <- hist_d$pred_dens %>% 
  rbind(
    hist_d$pred_dens %>% 
      group_by(x) %>% 
      summarise(p = sum(p), state = "total")
  )
hist_d$pred_dens <- hist_d$pred_dens %>% 
  mutate(
    is_total = state == "total",
    lt = ifelse(state == "total", "dotted", 
                ifelse(
                  grepl("less", state), "solid", "dashed"
                )),
    state = ifelse(grepl("1", state), "1", "2")
  )



g4 <- data.frame(
  "x" = hist_d$x,
  "p" = hist_d$density
) %>% 
  ggplot(aes(x = x / 1000 * 12, y = p)) + 
  geom_col(fill = "lightgrey") + 
  scale_xy_log(axis = "x") + 
  theme_bw(base_size = 15) + 
  geom_line(
    data = hist_d$pred_dens %>% 
      filter(!is_total),
    aes(color = state, linetype = lt), 
    size = 1.5,
    alpha = 0.8
  ) +
  geom_line(
    data = hist_d$pred_dens %>% 
      filter(is_total),
    aes(linetype = lt), 
    size = 1.5,
    color = "black"
  ) + 
  scale_linetype_identity() + 
  scale_color_manual(values = c(rev(.getPalette(2)), "black")) +
  labs(x = "Observed step lengths (cm)", y = "Density", color = "State");g4





g_final <- ggarrange(g2, g1, g3,
                     ncol = 1,
          common.legend = TRUE, 
          align = "hv", 
          labels = "AUTO",
          legend = "top")

ggsave("graphs/manuscript1_figures/issf_params.png", g_final, dpi = 600, width = 4, 
       height = 11, bg = "white")


#### Main text simplified ISSF simulation result ####
# out <- read_csv("simulation/move_rules_sim.csv")
# 
# out <- out %>% 
#   mutate(
#     scale = round(scale / 1000 * 12, 2),
#     rss = as.factor(round(exp(rss), 2))
#   )
# 
# out %>% 
#   filter(b != 0) %>% 
#   filter(scale == 1) %>% 
#   group_by(b, rss, k) %>%
#   dplyr::select(on_toxic) %>% 
#   mutate(
#     id = seq_along(on_toxic)
#   ) %>% 
#   arrange(b, rss, k, id) %>% 
#   group_by(rss, k, id) %>% 
#   summarise(
#     delta = diff(on_toxic)
#   ) %>% 
#   ggplot(aes(x = as.factor(k), y = delta, color = factor(rss))) +
#   geom_hline(aes(yintercept = 0), color = "grey", 
#              size = 1, linetype = "dashed") + 
#   geom_point(position = position_jitterdodge(jitter.width = 0.2, 
#                                              jitter.height = 0, 
#                                              dodge.width = 0.5), 
#              alpha = 0.2, 
#              size = 2) + 
#   geom_point(stat = "summary", 
#              position = position_dodge(width = 0.5),
#              shape = "—",
#              size = 4,
#              color = "black",
#              aes(group = factor(rss))) + 
#   theme_bw(base_size = 15) +
#   theme(legend.position = "top") + 
#   guides(colour = guide_legend(
#     override.aes = list(alpha = 1, size = 4), 
#     title.position = "top"
#   )) + 
#   labs(x = "Less toxic diet arresetment strength (ratio)", 
#        y = expression(atop(Change~"in"~proprtion~time~on, paste("more toxic diet",~(beta[5]-beta[-5])))),
#        color = "Less toxic diet selection strength (odds)") + 
#   scale_color_brewer(type = "seq") -> g;g
# 
# 
# 
# 
# ggsave("graphs/manuscript1_figures/issf_simplified_sim_result.png", g, dpi = 600, width = 5, height = 4.5)
# 


#### Supplement full ISSF simulation result ####
out <- read_csv("simulation/move_rules_sim_reduce.csv")

out <- out %>% 
  mutate(
    scale = round(scale / 1000 * 12, 2),
    rss = as.factor(round(exp(rss), 2))
  )

k_lab <- unique(out$k)
names(k_lab) = sprintf("k = %s", k_lab)
scale_lab <- unique(out$scale)
names(scale_lab) = sprintf("scale = %s", scale_lab)
rss_lab <- unique(out$rss)
names(rss_lab) = sprintf("rss = %s", rss_lab)


out %>% 
  ggplot(aes(x = as.factor(k), y = on_toxic, color = factor(b))) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, 
                                             jitter.height = 0, 
                                             dodge.width = 0.5), 
             alpha = 0.2,
             size = 1) + 
  geom_point(stat = "summary", 
             position = position_dodge(width = 0.5),
             shape = "—",
             size = 3,
             color = "black",
             aes(group = factor(b))) + 
  theme_bw(base_size = 15) +
  guides(colour = guide_legend(
    override.aes = list(alpha = 1, size = 3), 
  )) + 
  facet_grid(scale ~ rss, labeller = labeller(
    rss = reverse_names(rss_lab), scale = reverse_names(scale_lab)
  )) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  labs(x = "Less toxic diet arresetment strength (ratio)", 
       y = "Proportion time on toxic diet", 
       color = expression(beta)) + 
  scale_color_brewer(type = "qual")->g;g


ggsave("graphs/manuscript1_figures/issf_on_toxic_sim_result.png", g, dpi = 600, width = 8, height = 6)




#### On toxic observed plot ####

m <- readRDS("invisible/fitted_models/on_toxic_brm.rds")


g <- brms::conditional_effects(m,effects = c("cat_pre_wt_log_scale:beta"), re_formula = NA)[[1]] %>% 
  ggplot(aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = estimate__, color = beta)) + 
  geom_ribbon(aes(ymax = upper__, ymin = lower__, color = NULL, fill = beta), alpha = 0.2) + 
  geom_line(size = 2) + 
  geom_point(
    data = m$data,
    aes(x = unscalelog(d$cat_pre_wt_log)(cat_pre_wt_log_scale), y = on_toxic),
    size = 5, 
    alpha = 0.5
  ) + 
  theme_bw(base_size = 15) + 
  scale_x_continuous(trans = "log10", labels = fancy_scientific) +
  scale_color_brewer(type = "qual", aesthetics = c("fill", "color")) + 
  labs(x = "Pre-weight (g)", y = "Proprtion of time on toxic diet", color = expression(beta), fill = expression(beta)) +
  theme(legend.position = "top")


ggsave("graphs/manuscript1_figures/observed_on_toxic.png", g, dpi = 600, width = 5, height = 5)
