reload()
spec <- dummy_spec()
spec <- c(1:12 %% 2,
          (1:12+1) %% 2) %>% 
  rep(6) %>% 
  matrix(nrow = 12, ncol = 12) %>% 
  as.cimg()

fetch_trt_spec(5) %>% plot()

iterate_random_steps2(issf_fit, 
                      start = make_start2(0,500,500), 
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





mean(x$ava_qual, na.rm = TRUE)
mean(x$on_toxic, na.rm = TRUE)


sim_d <- expand.grid(
  "b" = c(5, 0, -5),
  "rss" = c(0, 0.5, 1, 3)
) %>% 
  rep_data.frame(30)

out <- pb_par_lapply(
  seq_len(nrow(sim_d)), function(i, issf_fit, sim_d){
    spec <- as.cimg(syn_spec(12, sim_d[i,"b"], plot = FALSE))
    x <- iterate_random_steps2(issf_fit, 
                               start = make_start2(0,500,500), 
                               n = 10000, 
                               ref_grid = spec, 
                               rss_coef = sim_d[i,"rss"])
    data.frame(
      "b" = sim_d[i,"b"],
      "rss" = sim_d[i,"rss"],
      "ava_toxic" = mean(x$ava_toxic, na.rm = TRUE),
      "on_toxic" = mean(x$on_toxic, na.rm = TRUE)
    )
  },
  issf_fit = issf_fit,
  sim_d = sim_d,
  cores = 6, 
  inorder = FALSE
)


do.call("rbind", out) %>% 
  ggplot(aes(x = as.factor(b), y = 1 - ava_qual)) + 
  geom_pointrange(stat = "summary") + 
  geom_point(position = position_jitter(width = 0.2, height = 0))


do.call("rbind", out) %>% 
  ggplot(aes(x = as.factor(b), y = on_toxic)) + 
  geom_pointrange(stat = "summary") + 
  geom_point(position = position_jitter(width = 0.2, height = 0))





