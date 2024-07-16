source("spat1f/init_dev.R")

data_prep <- function(ID, ref_data = get("ref_data", envir = globalenv())){
  d4 <- fetch_events(ID) %>% # Get mask-R-CNN instances
    clean_events(ref_data = ref_data) %$% # Clean those instance
    move_seq(head_x, head_y, r_thresh = 0, inherit.theta = FALSE) # compute the step length and turn angl
  d4$state <- as.factor(viterbi(fit_HMM(as.moveData(d4))))
  d4 <- d4 %>%
    filter(!is.na(r)) %>%  # throw out steps where r is NA
    add_random_steps(n = 100L, # Simulate random available steps
                     sl_distr = fit_gamma(.$r), # Fit gamma step length dist
                     ta_distr = fit_genvonmises(.$theta_rel) # Generalized von Mises turn angle dist
    ) %>%
    flag_invalid_steps(remove = TRUE) %>% # Throw out any random step that is outside of the arena
    mutate(
      # Find the diet value at the simulated destination location
      toxic = read_value(x1, y1, c(1000, 1000),
                         ref_img = fetch_trt_spec(ID, 
                                                  .ref_data = ref_data, 
                                                  quiet = TRUE)), 
      less_toxic = 1 - toxic,
      less_toxic_f = as.factor(less_toxic),
      is_state1 = ifelse(state == 1, 1, 0)
    ) %>%
    append_estimators(na_as_zero = TRUE)  # Append correct step length and turn angle estimators
}


m <- issf( 
  case ~ 
    state:cos_theta_pi:less_toxic_f + 
    state:cos_2theta:less_toxic_f + 
    state:sl:less_toxic_f +
    state:logsl:less_toxic_f +
    is_state1:sl:memdel(less_toxic, 
                        step_id, 
                        case,
                        weight_FUN = FUN,
                        k = 0.01) + 
    strata(step_id),
  data = d4, 
  update = FALSE, 
  keep_data = TRUE
);summary(m);AIC(m) 


m <- issf( 
  case ~ 
    state:cos_theta_pi:less_toxic_f + 
    state:cos_2theta:less_toxic_f + 
    state:sl:less_toxic_f +
    state:logsl:less_toxic_f +
    is_state1:sl:memdel(less_toxic, 
                        step_id, 
                        case,
                        weight_FUN = FUN) + 
    strata(step_id),
  data = d4, 
  update = FALSE, 
  keep_data = TRUE
);summary(m);AIC(m) 


FUN <- function(t, k = 0.1){
  1 / (1 + t * k)
}

debug(bar2)




d4 <- data_prep(147)
f <- function(theta){
  res <- issf( 
    case ~ 
      state:cos_theta_pi:less_toxic_f + 
      state:cos_2theta:less_toxic_f + 
      state:sl:less_toxic_f +
      state:logsl:less_toxic_f +
      is_state1:sl:memdel(less_toxic, 
                    step_id, 
                    case, 
                    FUN, 
                    k = exp(theta)) + 
      strata(step_id),
    data = d4, 
    update = FALSE, 
    keep_data = TRUE
  ) %>% AIC()
  cat(sprintf("k = %s, AIC = %s \n", signif(exp(theta)), signif(res)))
  return(res)
}

f(-3)

optim2(
  -2, 
  method = "Brent", 
  fn = f, 
  lower = -5, 
  upper = 3
) -> out

out

memory_decay <- function(k, decay_prop = 0.05){
  # decay_prop = amount of memory left
  (1 / decay_prop - 1) / k
}

memory_decay(k = exp(-3.2)) / 10



d3 %>% 
  select(beta, rep_id, k1) %>% 
  as.data.frame() %>% 
  arrange(-abs(log(k1)))






