
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
  syn_spec(300, beta = 1, threshold = FALSE, plot = FALSE),
  syn_spec(300, beta = 5, threshold = FALSE, plot = FALSE)
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




