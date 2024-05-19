source("spat1f/init_analysis.R")

ID <- 81

fetch_events(ID) %>% 
  clean_events() %$%
  move_seq(head_x, head_y) %>% 
  filter(!is.na(r)) %>% 
  add_random_steps(n = 100L,
                   sl_distr = fit_gamma(.$r),
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
      (sl + logsl) + 
      strata(step_id),
    data = .) -> issf_fit




spec <- dummy_spec()
spec <- c(1:12 %% 2,
          (1:12+1) %% 2) %>% 
  rep(6) %>% 
  matrix(nrow = 12, ncol = 12) %>% 
  as.cimg()

fetch_trt_spec(5) %>% plot()

iterate_random_steps2(issf_fit, 
                      start = make_start2(0,500,500, 1), 
                      n = 5000, 
                      ref_grid = spec, 
                      rss_coef = 0.9) %>% 
  mutate(head_x = x, head_y = y) %>% 
  plot_track_overlay(
    repID = "Foo", 
    colored_track = FALSE, 
    trt_spec = spec, 
    plot_elements = "track"
  )



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

microbenchmark::microbenchmark(
  iterate_random_steps2(issf_fit, 
                        start = make_start2(0,500,500, on_toxic = 1), 
                        n = 5000, 
                        ref_grid = as.cimg(syn_spec(12, 5, plot = FALSE)), 
                        rss_coef = 0.9), times = 10
)
x <- iterate_random_steps2(issf_fit, 
                           start = make_start2(0,500,500, on_toxic = 1), 
                           n = 5000, 
                           ref_grid = as.cimg(syn_spec(12, 5, plot = FALSE)), 
                           rss_coef = 0.9)


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



d %>% 
  mutate(
    scale = sprintf("scale = %s", scale),
    k = sprintf("k = %s", k)
  ) %>% 
  ggplot(aes(x = as.factor(rss), y = 1 - ava_toxic, color = factor(b))) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, 
                                             jitter.height = 0, 
                                             dodge.width = 0.5)) + 
  geom_pointrange(stat = "summary", 
                  position = position_dodge(width = 0.5), 
                  color = "black",
                  aes(group = factor(b)), linewidth = 1) + 
  facet_wrap(~scale + k) + 
  theme_bw() + 
  labs(x = "Less toxic diet selection strength (LogOdds)", 
       y = "Neighborhood quality (Proportion)", 
       color = expression(beta)) + 
  scale_color_brewer(type = "qual")


d %>% 
  ggplot(aes(x = as.factor(rss), y = 1 - on_toxic, color = factor(b))) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, 
                                             jitter.height = 0, 
                                             dodge.width = 0.5)) + 
  geom_pointrange(stat = "summary", 
                  position = position_dodge(width = 0.5), 
                  color = "black",
                  aes(group = factor(b)), linewidth = 1) + 
  facet_wrap(~scale + k)


d %>% 
  ggplot(aes(x = as.factor(rss), y = r, color = factor(b))) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, 
                                             jitter.height = 0, 
                                             dodge.width = 0.5)) + 
  geom_pointrange(stat = "summary", 
                  position = position_dodge(width = 0.5), 
                  color = "black",
                  aes(group = factor(b)), linewidth = 1) + 
  facet_wrap(~scale + k, scales = "free")



