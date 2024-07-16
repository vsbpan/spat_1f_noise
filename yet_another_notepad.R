source("spat1f/init_analysis.R")
ID <- 55
d <- fetch_events(ID) %>% # Get mask-R-CNN instances
  clean_events(ref_data = ref_data) %$% # Clean those instance
  move_seq(head_x, head_y, r_thresh = 0, inherit.theta = FALSE) # compute the step length and turn angl

d$state <- as.factor(viterbi(fit_HMM(as.moveData(d))))

d$r_lag1 <- lag(d$r, n = 1)
d$r_lag2 <- lag(d$r, n = 2)
d$r_lag3 <- lag(d$r, n = 3)
d$r_lag4 <- lag(d$r, n = 4)
d$r_lag5 <- lag(d$r, n = 5)
d$r_lag6 <- lag(d$r, n = 6)
d$r_lag10 <- lag(d$r, n = 10)
d$theta_lag1 <- lag(d$theta_rel, n = 1)
d$size_px <- fetch_events(ID) %>%
  clean_events(ref_data = ref_data) %>% 
  .$size_px %>% 
  .[-1]

d <- d %>% 
  mutate(
    toxic = read_value(x1, y1, c(1000, 1000),
                       ref_img = fetch_trt_spec(ID, 
                                                .ref_data = ref_data, 
                                                quiet = TRUE)), 
    less_toxic = 1 - toxic,
    less_toxic_f = as.factor(less_toxic)
  )
d$less_toxic_f_lag1 <- lag(d$less_toxic_f, n = 1)
d$less_toxic_f_lag2 <- lag(d$less_toxic_f, n = 2)
d$less_toxic_f_lag3 <- lag(d$less_toxic_f, n = 3)
d$less_toxic_f_lag4 <- lag(d$less_toxic_f, n = 4)
d$less_toxic_f_lag5 <- lag(d$less_toxic_f, n = 5)
d$less_toxic_f_lag6 <- lag(d$less_toxic_f, n = 6)

d$less_toxic_h1 <- d$less_toxic %>% roll_vapply(11, function(x) mean(x, na.rm = TRUE))
d$less_toxic_h1_lag1 <- lag(d$less_toxic_h1, n = 10)
d$less_toxic_h1_lag2 <- lag(d$less_toxic_h1, n = 20)
d$less_toxic_h1_lag3 <- lag(d$less_toxic_h1, n = 30)
d$less_toxic_h1_lag4 <- lag(d$less_toxic_h1, n = 40)
d$less_toxic_h1_lag5 <- lag(d$less_toxic_h1, n = 50)
d$less_toxic_h1_lag6 <- lag(d$less_toxic_h1, n = 60)
d$less_toxic_h1_lag7 <- lag(d$less_toxic_h1, n = 70)
d$less_toxic_h1_lag8 <- lag(d$less_toxic_h1, n = 80)



d <- d %>%
  filter(!is.na(r)) %>%  # throw out steps where r is NA
  add_random_steps(n = 100L, # Simulate random available steps
                   sl_distr = fit_gamma(.$r), # Fit gamma step length dist
                   ta_distr = fit_genvonmises(.$theta_rel) # Generalized von Mises turn angle dist
  ) %>%
  flag_invalid_steps(remove = TRUE) %>% # Throw out any random step that is outside of the arena
  mutate(
    # Find the diet value at the simulated destination location
    toxic = read_value(x2, y2, c(1000, 1000),
                       ref_img = fetch_trt_spec(ID, 
                                                .ref_data = ref_data, 
                                                quiet = TRUE)), 
    less_toxic = 1 - toxic,
    less_toxic_f = as.factor(less_toxic)
  ) %>%
  append_estimators(na_as_zero = TRUE) # Append correct step length and turn angle estimators



out <- issf( # Fit the issf
  case ~ 
    (cos_theta_pi + cos_2theta) +
    state:sl:less_toxic_f + 
    #(state:less_toxic_f:sl + state:less_toxic_f:logsl) + 
    #state:less_toxic +
    state:less_toxic_f:sl:less_toxic_f_lag1 + 
    state:less_toxic_f:sl:less_toxic_f_lag2 + 
    #state:less_toxic_h1_lag1:sl + 
    #state:less_toxic_h1_lag6 + 
    strata(step_id), # Stratify be step ID
  data = d ,
  sl_estimators = pick_default_estimators(
    "gamma", 
    list(
      c("less_toxic_f0", "less_toxic_f1"),
      c("state1", "state2")
    )
  ), 
  ta_estimators = pick_default_estimators(
    "genvonmises", 
    list(
      c("less_toxic_f0", "less_toxic_f1"),
      c("state1", "state2")
    )
  ), 
  keep_data = TRUE
); summary(out); AIC(out$model)




d %>% 
  filter(case) %>%
  ggplot(aes(x = sl, y = lag(sl))) + 
  geom_path() + 
  scale_x_log10() + 
  scale_y_log10() + 
  facet_wrap(~state)


d


z <- d %>% filter(case)


library(biwavelet)

fetch_events(81) %>%
  clean_events(ref_data = ref_data, keep_sus = TRUE, score_thresh = 0, insert_gaps = TRUE) %>% 
  mutate(
    h = floor(rank / 5)
  ) %>% 
  group_by(h) %>% 
  summarise(
    head_x = mean(head_x, na.rm = TRUE),
    head_y = mean(head_y, na.rm = TRUE)
  ) %>% 
  as.data.frame() %>% 
  .[-nrow(.),] %>% 
  as.matrix() %>% 
  wt() %>% 
  biwavelet::plot.biwavelet()













