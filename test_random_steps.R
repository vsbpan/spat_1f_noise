source("spat1f/init_analysis.R")

ID <- 81

fetch_events(ID) %>% 
  clean_events() %$%
  move_seq(head_x, head_y) %>% 
  filter(!is.na(r)) %>% 
  add_random_steps(n = 100L,
                   sl_distr = fit_invgamma(.$r),
                   ta_distr = fit_genvonmises(.$theta_rel)
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
      (invsl + logsl) + 
      strata(step_id),
    data = .) -> issf_fit




spec <- dummy_spec()
spec <- c(1:12 %% 2,
          (1:12+1) %% 2) %>% 
  rep(6) %>% 
  matrix(nrow = 12, ncol = 12) %>% 
  as.cimg()

fetch_trt_spec(5) %>% plot()

spec <- fetch_trt_spec(ID)

iterate_random_steps2(
  ta_sl_list = list(
    "sl" = list(
      make_gamma(1, 30),
      make_gamma(1, 10)
    ),
    "ta" = list(
      make_genvonmises(0.48, 0.267),
      make_genvonmises(0.48, 0.267)
    )
  ), 
  start = make_start2(0,500,500, 1), 
  n = 5000, 
  ref_grid = spec, 
  same_move = FALSE, 
  rss_coef = 0.5
) %>% 
  mutate(head_x = x, head_y = y) %>% 
  plot_track_overlay(
    repID = "Foo", 
    colored_track = "none", 
    trt_spec = flip_xy(spec), 
    plot_elements = "track"
  )

plot_track_overlay(
  repID = ID, 
  colored_track = FALSE, 
  trt_spec = spec, 
  plot_elements = "track"
)

reload()




ID <- sample(setdiff(fetch_repID(), problem_ids), size = 1)
x <- fetch_events(ID) %>% 
  clean_events() %$% 
  move_seq(head_x, head_y) %>% 
  .$r %>% 
  na.omit()

plot_track_overlay(repID = ID)

comp_dist(x)
x %>% 
  loghist(nclass = 50, log.p = FALSE, geom = "col", 
          draw_dist = list(fit_frechet(x), fit_invgamma(x), fit_lnorm(x), fit_gamma(x))
          )

extraDistr::r(10000, 5) %>% loghist(log.p = FALSE)


comp_dist(x, dist = c("frechet","invgamma","gamma","lnorm"))






survival_plot(x)







issf_fit <- list(
  "ta_updated" = list(
    make_genvonmises(kappa1 = 0.4867337, kappa2 = 0.2466265),
    make_genvonmises(kappa1 = 0.4867337, kappa2 = 0.2466265)
  ),
  "sl_updated" = list(
    make_gamma(shape = 0.5717441, scale = 50),
    make_gamma(shape = 0.5717441, scale = 30)
  )
)



x



for (j in 1:5){
  sim_d <- expand.grid(
    "b" = c(5, 0, -5),
    "rss" = c(0, 0.5, 1), 
    "scale" = c(10, 30, 100),
    "k" = c(0.33, 1, 3)
  ) %>% 
    rep_data.frame(40)
  
  out <- pb_par_lapply(
    seq_len(nrow(sim_d)), function(i, sim_d){
      issf_fit2 <- list(
        "ta_updated" = list(
          make_genvonmises(kappa1 = 0.4867337, kappa2 = 0.2466265),
          make_genvonmises(kappa1 = 0.4867337, kappa2 = 0.2466265)
        ),
        "sl_updated" = list(
          make_gamma(shape = 0.5717441, scale = sim_d[i, "scale"]),
          make_gamma(shape = 0.5717441, scale = sim_d[i, "scale"] / sim_d[i, "k"]) 
        )
      )
      
      spec <- as.cimg(syn_spec(12, sim_d[i,"b"], plot = FALSE))
      x <- iterate_random_steps2(issf_fit2, 
                                 start = make_start2(0,500,500, 1), 
                                 n = 10000, 
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
  
  write_csv(do.call("rbind", out), file = sprintf("check_point_%s.csv",j))
}



d <- read_csv("cleaned_data/issf_sim_experiment_(k_rss_sl_beta).csv")

reverse_names <- function(x){
  val <- x
  nms <- names(x)
  names(nms) <- val
  return(nms)
}

k_lab <- unique(d$k)
names(k_lab) = sprintf("k = %s", k_lab)
scale_lab <- unique(d$scale)
names(scale_lab) = sprintf("scale = %s", scale_lab)


d %>% 
  ggplot(aes(x = as.factor(rss), y = 1 - ava_toxic, color = factor(b))) + 
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
  labs(x = "Less toxic diet selection strength (LogOdds)", 
       y = "Neighborhood quality (Proportion)", 
       color = expression(beta)) + 
  scale_color_brewer(type = "qual")->g;g

ggsave("graphs/issf_sim_neighborhood_quality.png", plot = g, dpi = 400, height = 7, width = 8)


d %>% 
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
  labs(x = "Less toxic diet selection strength (LogOdds)", 
       y = "Time on less toxic diet (Proportion)", 
       color = expression(beta)) + 
  scale_color_brewer(type = "qual")->g;g

ggsave("graphs/issf_sim_on_less_toxic.png", plot = g, dpi = 400, height = 7, width = 8)



d %>% 
  filter(k == 1) %>% 
  ggplot(aes(x = as.factor(rss), y = r / 1000 * 12, color = factor(b))) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, 
                                             jitter.height = 0, 
                                             dodge.width = 0.5), 
             alpha = 0.2) + 
  geom_pointrange(stat = "summary", 
                  position = position_dodge(width = 0.5), 
                  color = "black",
                  aes(group = factor(b)), linewidth = 1) + 
  facet_grid(scale ~ ., labeller = labeller(
    k = reverse_names(k_lab), scale = reverse_names(scale_lab), 
  ), scale = "free") + 
  theme_bw() + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  labs(x = "Less toxic diet selection strength (LogOdds)", 
       y = "Average step length (cm)", 
       subtitle = "Fixed at k = 1",
       color = expression(beta)) + 
  scale_color_brewer(type = "qual") -> g;g

ggsave("graphs/issf_sim_step_length2.png", plot = g, dpi = 400, height = 7, width = 8)






