source("spat1f/init_analysis.R")



plot_track_overlay(repID = 81)

ID <- 81

d0 <- fetch_events(ID) %>%
  clean_events(ref_data = ref_data) %>%
  insert_gaps() %>% 
  mutate(
    head_x = ifelse(score >=0.9, head_x, NA),
    head_y = ifelse(score >=0.9, head_y, NA)
  )


d <- d0 %$%
  move_seq(head_x, head_y, r_thresh = 0, inherit.theta = FALSE) %>% 
  bind_cols(
    d0[-nrow(d0), c("size_px")]
  ) %>%
  filter(!is.na(r))


d <- d %>%
  add_random_steps(n = 100L,
                   sl_distr = fit_gamma(.$r),
                   ta_distr = fit_genvonmises(.$theta_rel)
  ) %>%
  flag_invalid_steps(remove = TRUE) %>%
  mutate(
    start_toxic = ifelse(
      read_value(x1, y1, c(1000, 1000),
                 ref_img = fetch_trt_spec(ID, .ref_data = ref_data, quiet = TRUE)) == 1,
      "yes",
      "no"
    ),
    toxic = read_value(x2, y2, c(1000, 1000),
                       ref_img = fetch_trt_spec(ID, .ref_data = ref_data, quiet = TRUE))
  ) %>%
  append_estimators(na_as_zero = TRUE) %>% 
  mutate(
    time = step_id / 1000
  )


scale(unique(d$size_px))



m <- issf(
  case ~
    toxic + 
    (cos_theta_pi + cos_2theta) +
    (sl + logsl) + 
    sl:start_toxic + logsl:start_toxic + 
    strata(step_id),
  data = d,
  shape_estimator = c("logsl", "logsl:start_toxicyes"),
  scale_estimator = c("sl","sl:start_toxicyes"),
  kappa1_estimator = "cos_theta_pi",
  kappa2_estimator = "cos_2theta"
); summary(m)

m$sl_updated
plot_track_overlay(repID = 83)




fetch_events(90) %>%
  clean_events(ref_data = ref_data) %>%
  insert_gaps() %>% 
  mutate(
    head_x = ifelse(score >=0.9, head_x, NA),
    head_y = ifelse(score >=0.9, head_y, NA)
  ) %$%
  move_seq(head_x, head_y, r_thresh = 0, inherit.theta = FALSE) %>% 
  select(r, theta_rel) %>%
  mutate(r = log(r)) %>% 
  as.matrix() %>% 
  scale() %>% 
  na.omit() %>% 
  princomp() %>% 
  vegan::scores(choice = 1) %>% 
  unlist() %>% 
  na.omit() %>% 
  as.vector() %>%
  ts() %>% 
  acf()




































