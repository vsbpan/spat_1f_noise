source("spat1f/init.R")

for (j in 1:10){
  sim_d <- expand.grid(
    "b" = c(5, 0, -5),
    "rss" = c(log(0.75), 0, log(1.333333)), 
    "scale" = c(8.3333, 25, 83.333),
    "k" = c(0.2, 1, 5)
  ) %>% 
    rep_data.frame(40)
  
  out <- pb_par_lapply(
    seq_len(nrow(sim_d)), function(i, sim_d){
      issf_fit2 <- list(
        "ta_updated" = list(
          make_genvonmises(kappa1 = 0.4, kappa2 = 0.3),
          make_genvonmises(kappa1 = 0.4, kappa2 = 0.3)
        ),
        "sl_updated" = list(
          make_gamma(shape = 0.8, 
                     scale = sim_d[i, "scale"] * sim_d[i, "k"]), # More toxic diet
          make_gamma(shape = 0.8, 
                     scale = sim_d[i, "scale"])  # Less toxic diet
        )
      )
      
      spec <- as.cimg(syn_spec(12, sim_d[i,"b"], plot = FALSE))
      x <- iterate_random_steps2(issf_fit2, 
                                 start = make_start2(0,500,500, 1), 
                                 n = 1000, 
                                 ref_grid = spec, 
                                 rss_coef = sim_d[i,"rss"])
      x <- x[, c("on_toxic", "r", "ava_toxic")]
      
      data.frame(
        sim_d[i, ,drop = FALSE],
        as.data.frame(t(colMeans(x, na.rm = TRUE)))
      )
    },
    sim_d = sim_d,
    cores = 8, 
    inorder = FALSE
  )
  
  write_csv(do.call("rbind", out), file = sprintf("simulation/check_points/check_point_%s.csv",j))
}



# out <- list.files("simulation/check_points", full.names = TRUE) %>%
#   lapply(function(x){
#     suppressMessages(read_csv(x, progress = FALSE))
#   }) %>%
#   do.call("rbind", .)
# 
# write_csv(out, "simulation/move_rules_sim.csv")
# 
# file.remove(list.files("simulation/check_points", full.names = TRUE))



out <- read_csv("simulation/move_rules_sim.csv")


out <- out %>% 
  mutate(
    scale = round(scale / 1000 * 12, 2),
    rss = as.factor(round(exp(rss), 2))
  )

reverse_names <- function(x){
  val <- x
  nms <- names(x)
  names(nms) <- val
  return(nms)
}

k_lab <- unique(out$k)
names(k_lab) = sprintf("k = %s", k_lab)
scale_lab <- unique(out$scale)
names(scale_lab) = sprintf("scale = %s", scale_lab)


out %>% 
  ggplot(aes(x = (rss), y = 1 - ava_toxic, color = factor(b))) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, 
                                             jitter.height = 0, 
                                             dodge.width = 0.5), 
             alpha = 0.2) + 
  geom_pointrange(stat = "summary", 
                  position = position_dodge(width = 0.5), 
                  color = "black",
                  aes(group = factor(b)), linewidth = 1) + 
  facet_grid(k ~ scale, labeller = labeller(
    k = reverse_names(k_lab), scale = reverse_names(scale_lab)
  )) + 
  theme_bw() + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  labs(x = "Less toxic diet selection strength (Odds)", 
       y = "Neighborhood quality (Proportion)", 
       color = expression(beta)) + 
  scale_color_brewer(type = "qual")->g;g



out %>% 
  ggplot(aes(x = as.factor(rss), y = 1 - on_toxic, color = factor(b))) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, 
                                             jitter.height = 0, 
                                             dodge.width = 0.5), 
             alpha = 0.2) + 
  geom_pointrange(stat = "summary", 
                  position = position_dodge(width = 0.5), 
                  color = "black",
                  aes(group = factor(b)), linewidth = 1) + 
  facet_grid(k ~ scale, labeller = labeller(
    k = reverse_names(k_lab), scale = reverse_names(scale_lab)
  )) + 
  theme_bw() + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  labs(x = "Less toxic diet selection strength (Odds)", 
       y = "Time on less toxic diet (Proportion)", 
       color = expression(beta)) + 
  scale_color_brewer(type = "qual")->g;g




out %>% 
  filter(scale == 1) %>% 
  ggplot(aes(x = as.factor(rss), y = r / 1000 * 12, color = factor(b))) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, 
                                             jitter.height = 0, 
                                             dodge.width = 0.5), 
             alpha = 0.2) + 
  geom_pointrange(stat = "summary", 
                  position = position_dodge(width = 0.5), 
                  color = "black",
                  aes(group = factor(b)), linewidth = 1) + 
  facet_grid(k ~ ., labeller = labeller(
    k = reverse_names(k_lab), scale = reverse_names(scale_lab), 
  ), scale = "free") + 
  theme_bw() + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  labs(x = "Less toxic diet selection strength (LogOdds)", 
       y = "Average step length (cm)", 
       subtitle = "Fixed at scale = 1",
       color = expression(beta)) + 
  scale_color_brewer(type = "qual") -> g;g



out %>% 
  filter(k == 1 & rss == 1) %>% 
  group_by(scale, b) %>% 
  do(as.data.frame(t(probe_dist(.$on_toxic, probe = c("mean","se")))))




k_lab <- unique(out$k)
names(k_lab) = sprintf("k = %s", k_lab)
scale_lab <- unique(out$scale)
names(scale_lab) = sprintf("scale = %s", scale_lab)
rss_lab <- unique(out$rss)
names(rss_lab) = sprintf("rss = %s", rss_lab)

out %>% 
  filter(scale == 1) %>% 
  ggplot(aes(x = as.factor(k), y = on_toxic, color = factor(b))) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, 
                                             jitter.height = 0, 
                                             dodge.width = 0.5), 
             alpha = 0.2) + 
  geom_pointrange(stat = "summary", 
                  position = position_dodge(width = 0.5), 
                  color = "black",
                  aes(group = factor(b)), linewidth = 1) + 
  facet_grid(rss ~ scale, labeller = labeller(
    rss = reverse_names(rss_lab), scale = reverse_names(scale_lab)
  )) + 
  theme_bw() + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  labs(x = "Less toxic diet arresetment strength (ratio)", 
       y = "Proportion time on toxic diet", 
       color = expression(beta)) + 
  scale_color_brewer(type = "qual")



out %>% 
  filter(b != 0) %>% 
  filter(scale == 1) %>% 
  group_by(b, rss, k) %>%
  dplyr::select(on_toxic) %>% 
  mutate(
    id = seq_along(on_toxic)
  ) %>% 
  arrange(b, rss, k, id) %>% 
  group_by(rss, k, id) %>% 
  summarise(
    delta = -diff(on_toxic)
  ) %>% 
  ggplot(aes(x = as.factor(k), y = delta, color = factor(rss))) +
  geom_hline(aes(yintercept = 0), color = "grey", 
             size = 1, linetype = "dashed") + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, 
                                             jitter.height = 0, 
                                             dodge.width = 0.5), 
             alpha = 0.2, 
             size = 2) + 
  geom_point(stat = "summary", 
             position = position_dodge(width = 0.5),
             shape = "—",
             size = 4,
             color = "black",
             aes(group = factor(rss))) + 
  theme_bw(base_size = 15) +
  theme(legend.position = "top") + 
  guides(colour = guide_legend(
    override.aes = list(alpha = 1, size = 4), 
    title.position = "top"
    )) + 
  labs(x = "Less toxic diet arresetment strength (ratio)", 
       y = expression(atop(Change~"in"~proprtion~time~on, paste("less toxic diet",~(beta[5]-beta[-5])))),
       color = "Less toxic diet selection strength (odds)") + 
  scale_color_brewer(type = "seq") -> g;g













