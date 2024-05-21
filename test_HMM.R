
d2 <- fetch_events(55) %>% 
  clean_events() %$% 
  move_seq(head_x, head_y)

hmm_fit <- fit_HMM(as.moveData(d2))
plot(hmm_fit, type = "sl")
plot(hmm_fit, type = "ta")
plot(hmm_fit, type = "resid")
plot(hmm_fit, animals = 1)

plot_track_overlay(repID = 81, colored_track = "states", plot_elements = "track")




ID <- 50
d2 <- fetch_events(ID) %>% 
  clean_events() %$% 
  move_seq(head_x, head_y)

hmm_fit <- fit_HMM(as.moveData(d2))

d2$state <- as.factor(viterbi(hmm_fit))

d2 <- d2 %>%
  filter(!is.na(r)) %>% 
  add_random_steps(n = 100L, # Simulate random available steps
                   sl_distr = fit_lnorm(.$r), # Fit gamma step length dist
                   ta_distr = fit_genvonmises(.$theta_rel) # Generalized von Mises turn angle dist
  ) %>%
  flag_invalid_steps(remove = TRUE) %>% # Throw out any random step that is outside of the arena
  mutate(
    # Find the diet value at the simulated destination location
    toxic = read_value(x2, y2, c(1000, 1000),
                       ref_img = fetch_trt_spec(ID, 
                                                .ref_data = ref_data, 
                                                quiet = TRUE)), 
    less_toxic = 1 - toxic
  ) %>%
  append_estimators(na_as_zero = TRUE) # Append correct step length and turn angle estimators


fit <- issf(
  case ~
    less_toxic + 
    (state:cos_theta_pi + state:cos_2theta) +
    (state:logslsq + state:logsl) +
    strata(step_id),
  data = d2, 
  sl_estimators = lapply(.lnorm_default_estimator(), function(x){
    sprintf("%s:%s", c("state1", "state2"), x)
  }), 
  ta_estimators = lapply(.genvonmises_default_estimator(), function(x){
    sprintf("%s:%s", c("state1", "state2"), x)
  })
)

plot_track_overlay(repID = 55, colored_track = "states", plot_elements = "track")






