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



out <- pb_par_lapply(
  rep(c(5, 0, -5), 30), function(b, issf_fit){
    spec <- as.cimg(syn_spec(12, b, plot = FALSE))
    x <- iterate_random_steps2(issf_fit, 
                               start = make_start2(0,500,500), 
                               n = 10000, 
                               ref_grid = spec, 
                               rss_coef = 0.9)
    data.frame(
      "b" = b,
      "ava_qual" = mean(x$ava_qual, na.rm = TRUE),
      "on_toxic" = mean(x$on_toxic, na.rm = TRUE)
    )
  },
  issf_fit = issf_fit,
  cores = 6, 
  inorder = FALSE
)


do.call("rbind", out)







