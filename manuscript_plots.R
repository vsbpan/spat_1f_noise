library(tidybayes)

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







sem_sim_d <- read_csv("cleaned_data/SEM_sim.csv")




sem_sim_d <- sem_sim_d %>% 
  mutate(
    exclude = ifelse(is.na(exclude), "total", exclude),
    only = ifelse(is.na(only), "total", only),
    cat_size = factor(ifelse(cat_size < 0, "small", "large"), 
                      level = rev(c("small", "large"))),
    only = factor(only, 
                  level = rev(c("total", "var_toxic_12_scale", "mean_toxic_conc_scale", "sl_mean_obs_log_scale","area_herb_log_scale"))),
    var = case_when(
      var == "beta_numeric_scale" ~ "Clusteredness~(beta)",
      var == "var_high" ~ "Var.~high~vs.~Var.~low"
    )
  )

sem_sim_d %>% 
  ggplot(
    aes(
      x = only, y = val, color = cat_size, shape = exclude
    )
  ) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") + 
  stat_pointinterval(
    position = position_dodge(width = 0.7), stroke = 2
  ) + 
  scale_x_discrete(label = rev(c("Total", "Temporal var. toxin", "Mean toxin ingested","Mean step length","Diet consumption"))) +
  coord_flip() + 
  theme_bw(base_size = 15) +
  theme(legend.position = "top", legend.box="vertical", legend.margin=margin(b = -5, t = -5)) +
  scale_shape_discrete(label = c("yes", "no")) +
  facet_wrap(~ var, labeller = label_parsed, scales = "free_x") +
  scale_color_discrete(type = c("steelblue","limegreen")) + 
  labs(y = "Standardized total effect on RGR", 
       x = "Proximal mediator", color = "Pre-weight", 
       shape = "Neighborhood diet quality constant") + 
  geom_text(
    data = sem_sim_d %>% 
      group_by(var,cat_size,exclude,only) %>% 
      summarise(
        lower = quantile(val, probs = 0.025),
        upper = quantile(val, probs = 0.97)
      ) %>% 
      ungroup() %>% 
      mutate(
        sig = (lower < 0 & upper < 0) | (lower > 0 & upper > 0),
        star = ifelse(sig, "*", "ns")
      ) %>% 
      group_by(var) %>% 
      mutate(
        max_val = max(upper)
      ),
    aes(label = star, y = max_val * 1.3, fill = cat_size, color = NULL),
    size = 4,
    position = position_dodge(width = 0.7)
  ) + 
  geom_text(
    data = sem_sim_d %>% 
      group_by(var,cat_size,exclude,only) %>% 
      summarise(
        lower = quantile(val, probs = 0.025),
        upper = quantile(val, probs = 0.975)
      ) %>% 
      ungroup() %>% 
      mutate(
        sig = (lower < 0 & upper < 0) | (lower > 0 & upper > 0),
        star = ifelse(sig, "*", "ns")
      ) %>% 
      group_by(var) %>% 
      mutate(
        max_val = max(upper)
      ),
    aes(label = star, y = max_val * 1.4, group = cat_size, color = NULL),
    alpha = 0,
    size = 4,
    position = position_dodge(width = 0.7)
  ) + 
  scale_y_continuous(labels = fancy_linear) -> g3;g3


ggsave("graphs/SEM_forest.png", g3, width = 6, height = 5.5, dpi = 600)






















sem_sim_d <- read_csv("cleaned_data/SEM_sim2.csv")


sem_sim_d <- sem_sim_d %>% 
  mutate(
    exclude = ifelse(is.na(exclude), "total", exclude),
    only = ifelse(is.na(only), "total", only),
    cat_size = factor(ifelse(cat_size < 0, "small", "large"), 
                      level = rev(c("small", "large"))),
    only = factor(only, 
                  level = rev(c("var_toxic_12_scale", "mean_toxic_conc_scale", "sl_mean_obs_log_scale","area_herb_log_scale","ava_mean_toxin_scale","select_scale"))),
    var = case_when(
      var == "beta_numeric_scale" ~ "Clusteredness~(beta)",
      var == "var_high" ~ "Var.~high~vs.~Var.~low"
    )
  )

sem_sim_d %>% 
  ggplot(
    aes(
      x = only, y = val, color = cat_size
    )
  ) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") + 
  stat_pointinterval(
    position = position_dodge(width = 0.7), stroke = 2
  ) + 
  scale_x_discrete(label = rev(c("Temporal var. toxin", "Mean toxin ingested","Mean step length","Diet consumption", "Neighborhood diet quality", "Diet selectivity"))) +
  coord_flip() + 
  theme_bw(base_size = 15) +
  theme(legend.position = "top", legend.box="vertical", legend.margin=margin(b = -5, t = -5)) +
  scale_shape_discrete(label = c("yes", "no")) +
  facet_wrap(~ var, labeller = label_parsed, scales = "free_x") +
  scale_color_discrete(type = c("steelblue","limegreen")) + 
  labs(y = "Standardized total effect on RGR", 
       x = "Proximal mediator", color = "Pre-weight") + 
  geom_text(
    data = sem_sim_d %>% 
      group_by(var,cat_size,only) %>% 
      summarise(
        lower = quantile(val, probs = 0.025),
        upper = quantile(val, probs = 0.975)
      ) %>% 
      ungroup() %>% 
      mutate(
        sig = (lower < 0 & upper < 0) | (lower > 0 & upper > 0),
        star = ifelse(sig, "*   ", "ns")
      ) %>% 
      group_by(var) %>% 
      mutate(
        max_val = max(upper)
      ) %>% 
      left_join(
        expand.grid(
          "only" = c('select_scale', 'ava_mean_toxin_scale', 'area_herb_log_scale', 
                     'sl_mean_obs_log_scale', 'mean_toxic_conc_scale', 'var_toxic_12_scale'),
          "var" = c("Clusteredness~(beta)", "Var.~high~vs.~Var.~low")
        ) %>% 
          cbind(
            "expected" = c("-", "+", "+","-","-","+",
                           "+", "+", "+", "-", "+", "-")
          )
      ),
    aes(label = paste0("(",expected, ")   ", star), y = max_val * 1.4, fill = cat_size, color = NULL),
    size = 4,
    position = position_dodge(width = 0.7)
  ) +
  geom_text(
    data = sem_sim_d %>% 
      group_by(var,cat_size,only) %>% 
      summarise(
        lower = quantile(val, probs = 0.025),
        upper = quantile(val, probs = 0.975)
      ) %>% 
      ungroup() %>% 
      mutate(
        sig = (lower < 0 & upper < 0) | (lower > 0 & upper > 0),
        star = ifelse(sig, "*", "ns")
      ) %>% 
      group_by(var) %>% 
      mutate(
        max_val = max(upper)
      ),
    aes(label = "2122223", y = max_val * 1.7, group = cat_size, color = NULL),
    alpha = 0,
    size = 4,
    position = position_dodge(width = 0.7)
  ) + 
  scale_y_continuous(labels = fancy_linear)



#### Methods figure insets ####

source("helper_functions/init_analysis.R")


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


img <- spat1f::fast_load_image("graphs/methods_figure/processed_rep81__cam20_s80016_rank223.jpg")
jpeg(paste0("graphs/methods_figure/", "processed_rep81__cam20_s80016_rank223_prediction", ".jpg"), 
     width = 15, height = 15, units = "cm", res = 600)
fetch_trt_spec(81) %>% 
  flip_xy() %>% 
  resize(1000, 1000) %>% 
  imdraw(
    img, ., opacity = 0.1
  ) %>% 
  plot.cimg(axes = FALSE)
fetch_data_dict(81)[223] %>% 
  plot(add = TRUE)
dev.off()

g <- fetch_trt_spec(81) %>% 
  flip_xy() %>% 
  as.matrix() %>% 
  melt() %>% 
  ggplot(aes(x = dim1, y = 12 - dim2)) + 
  geom_tile(
    aes(fill = as.factor(val)),
    alpha = 0.5) +
  geom_path(
    data = fetch_events(81) %>%
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

# ggsave("graphs/methods_figure/rep81_tracks.jpg", g, width = 7.5, height = 8, dpi = 600)


