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


# more toxic state 1
# less toxic state 1
# more toxic state 2
# less toxic state 2


spec <- rbinom(50^2, 1, 0.9) %>% image_unflatten()
spec <- syn_spec(beta = 3, plot = FALSE)
mat <- matrix(c(0.9, 0.1, 0.03, 0.97), nrow = 2, ncol = 2, byrow = TRUE)
mat <- matrix(1)

iterate_random_steps_states(ta_sl_list = list(
  "sl" = list(
    make_gamma(1, 150),
    make_gamma(1, 50),
    make_gamma(1, 6),
    make_gamma(1, 6)
  ),
  "ta" = list(
    make_unif(),
    make_unif(),
    make_genvonmises(0.48, 0.267),
    make_genvonmises(0.48, 0.267)
  )
), 
  n = 1000, 
  ref_grid = syn_spec(beta = 3, plot = FALSE), 
  rss_coef = 0.5,
  transition_mat =  diag(2), 
  max_xy = c(5000,5000)
) %>% 
  plot_track_overlay(
    colored_track = "none", 
    plot_elements = "track"
  )




reload()


