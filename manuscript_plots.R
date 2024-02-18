

g_bind <- ggarrange(g1, g2 + labs(y = "", fill = "Variation", color = "Variation"), 
          g3 + theme(legend.position = "none"), 
          g4 + labs(y = "")  + theme(legend.position = "none"), 
          align = "v", heights = c(1,0.8), labels = "AUTO", label.x = 0.1, label.y = c(0.85,0.85, 1, 1))
ggsave("graphs/figure2.png",g_bind, dpi = 600, width = 6.5, height = 6.5)







sem_sim_d <- read_csv("cleaned_data/SEM_sim.csv")


library(tidybayes)

sem_sim_d <- sem_sim_d %>% 
  mutate(
    exclude = ifelse(is.na(exclude), "total", exclude),
    only = ifelse(is.na(only), "total", only),
    cat_size = factor(ifelse(cat_size < 0, "small", "large"), 
                      level = rev(c("small", "large"))),
    only = factor(only, 
                  level = rev(c("total", "var_toxic_12_scale", "mean_toxic_conc_scale", "sl_mean_obs_log_scale","area_herb_log_scale"))),
    var = case_when(
      var == "beta_red" ~ "beta==5",
      var == "beta_white" ~ "beta==0",
      var == "var_high" ~ "Var.~high"
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
    position = "dodge", stroke = 2
  ) + 
  scale_x_discrete(label = rev(c("Total", "Var. toxin", "Mean toxin","Step length","Consumption"))) +
  coord_flip() + 
  theme_bw(base_size = 15) +
  theme(legend.position = "top") +
  scale_shape_discrete(label = c("yes", "no")) +
  facet_wrap(~ var, labeller = label_parsed, scales = "free_x") +
  scale_color_discrete(type = c("steelblue","limegreen")) + 
  labs(y = "Standardized total effect on RGR", 
       x = "Proximal mediator", color = "Pre-weight", 
       shape = "Availability constant") + 
  geom_text(
    data = sem_sim_d %>% 
      group_by(var,cat_size,exclude,only) %>% 
      summarise(
        lower = quantile(val, probs = 0.025),
        upper = quantile(val, probs = 0.975)
      ) %>% 
      ungroup() %>% 
      mutate(
        sig = lower & upper < 0 | lower & upper > 0,
        star = ifelse(sig, "*", "ns")
      ) %>% 
      group_by(var) %>% 
      mutate(
        max_val = max(upper)
      ),
    aes(label = star, y = max_val * 1.3, fill = cat_size, color = NULL),
    size = 4,
    position = position_dodge(width = 1)
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
        sig = lower & upper < 0 | lower & upper > 0,
        star = ifelse(sig, "*", "ns")
      ) %>% 
      group_by(var) %>% 
      mutate(
        max_val = max(upper)
      ),
    aes(label = star, y = max_val * 1.4, group = cat_size, color = NULL),
    alpha = 0,
    size = 4,
    position = position_dodge(width = 1)
  ) + 
  scale_y_continuous(labels = fancy_linear) -> g3


#ggsave("graphs/SEM_forest.png", g3, width = 8, height = 6, dpi = 600)






