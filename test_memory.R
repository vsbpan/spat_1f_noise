source("spat1f/init_analysis.R")

data_prep <- function(ID, ref_data = get("ref_data", envir = globalenv())){
  ref_img <- fetch_trt_spec(ID, 
                            .ref_data = ref_data, 
                            quiet = FALSE)
  #cat(sprintf("Preparing data for repID = %s, beta = %s\n", ID))
  

  d4 <- fetch_events(ID) %>% # Get mask-R-CNN instances
    clean_events(ref_data = ref_data) %$% # Clean those instance
    move_seq(head_x, head_y, r_thresh = 0, inherit.theta = FALSE) # compute the step length and turn angl
  d4$state <- as.factor(viterbi(fit_HMM(as.moveData(d4))))
  d4 <- d4 %>%
    filter(!is.na(r)) %>%  # throw out steps where r is NA
    add_random_steps(n = 500L, # Simulate random available steps
                     sl_distr = fit_gamma(.$r), # Fit gamma step length dist
                     ta_distr = fit_genvonmises(.$theta_rel) # Generalized von Mises turn angle dist
    ) %>%
    flag_invalid_steps(remove = TRUE) %>% # Throw out any random step that is outside of the arena
    mutate(
      # Find the diet value at the simulated destination location
      less_toxic_start = 1 - read_value(x1, y1, c(1000, 1000),
                                        ref_img = ref_img), 
      less_toxic_end = 1 - read_value(x2, y2, c(1000, 1000),
                                        ref_img = ref_img), 
      less_toxic_start_f = as.factor(less_toxic_start),
      less_toxic_end_f = as.factor(less_toxic_end),
      is_state1 = ifelse(state == 1, 1, 0)
    ) %>%
    append_estimators(na_as_zero = TRUE) %>%   # Append correct step length and turn angle estimators
    rbind.fill(
      data.frame("step_id" = seq(-1000, 0), "less_toxic" = 1, "case" = TRUE)
    )
}

FUN <- function(t, k = 0.1){
  1 / (1 + t * k)
  # exp(-k * t)
}



out <- lapply(d3$rep_id, function(x){
  d4 <- data_prep(x)
  m <- issf( 
    case ~ 
      state:logsl:less_toxic_start_f +
      state:sl:less_toxic_start_f +
      is_state1:sl:less_toxic_start_f:foo +
      strata(step_id),
    data = d4 %>% 
      mutate(
        foo = memory(less_toxic_start, 
                     step_id, 
                     case,
                     weight_FUN = FUN, 
                     k = 0.1)) %>% 
      filter(
        step_id > 0
      ), 
    update = FALSE, 
    keep_data = FALSE
  )
  return(m)
})


out2 <- lapply(out, function(x){
  x %>% 
    coef(se = TRUE) %>% 
    map(.f = function(x){
      x[c(9,10)]
    }) %>% 
    do.call("rbind", .) %>% 
    as.data.frame() %>% 
    rename_all(function(x){
      gsub(":",".",x)
    }) %>% 
    as.matrix() %>% 
    flatten_mat_name() 
})

d5 <- d3 %>% 
  left_join(
    out2 %>% 
      bind_vec() %>% 
      cbind(rep_id = d3$rep_id)
  )

d5 %>% 
  mutate(
    lower = less_toxic_start_f0.sl.is_state1.foo__estimate - 
      2 * less_toxic_start_f0.sl.is_state1.foo__se, 
    upper = less_toxic_start_f0.sl.is_state1.foo__estimate + 
      2 * less_toxic_start_f0.sl.is_state1.foo__se
  ) %>% 
  mutate(
    sig = sign(lower) == sign(upper)
  ) %>% 
  ggplot(aes(x = var_trt, y = less_toxic_start_f0.sl.is_state1.foo__estimate, color = beta, group = beta)) + 
  geom_pointrange(stat = "summary", stroke = 2, linewidth = 3, color = "black", 
                  position = position_dodge(width = 0.5)) + 
  geom_pointrange(
    aes(ymin = less_toxic_start_f0.sl.is_state1.foo__estimate - 2 * less_toxic_start_f0.sl.is_state1.foo__se,
        ymax = less_toxic_start_f0.sl.is_state1.foo__estimate + 2 * less_toxic_start_f0.sl.is_state1.foo__se,
        alpha = sig),
    position = position_jitterdodge(jitter.width = 0.5)
  ) 


d5 %>% 
  mutate(
    lower = less_toxic_start_f1.sl.is_state1.foo__estimate - 
      2 * less_toxic_start_f1.sl.is_state1.foo__se, 
    upper = less_toxic_start_f1.sl.is_state1.foo__estimate + 
      2 * less_toxic_start_f1.sl.is_state1.foo__se
  ) %>% 
  mutate(
    sig = sign(lower) == sign(upper)
  ) %>% 
  ggplot(aes(x = var_trt, y = less_toxic_start_f1.sl.is_state1.foo__estimate, color = beta, group = beta)) + 
  geom_pointrange(stat = "summary", stroke = 2, linewidth = 3, color = "black", 
                  position = position_dodge(width = 0.5)) + 
  geom_pointrange(
    aes(ymin = less_toxic_start_f1.sl.is_state1.foo__estimate - 2 * less_toxic_start_f1.sl.is_state1.foo__se,
        ymax = less_toxic_start_f1.sl.is_state1.foo__estimate + 2 * less_toxic_start_f1.sl.is_state1.foo__se,
        alpha = sig),
    position = position_jitterdodge(jitter.width = 0.5)
  ) 


d5 %>% 
  mutate(
    lower = less_toxic_start_f1.sl.is_state1.foo__estimate - 
      2 * less_toxic_start_f1.sl.is_state1.foo__se, 
    upper = less_toxic_start_f1.sl.is_state1.foo__estimate + 
      2 * less_toxic_start_f1.sl.is_state1.foo__se
  ) %>% 
  mutate(
    sig = sign(lower) == sign(upper)
  ) %>% 
  filter(
    abs(less_toxic_start_f1.sl.is_state1.foo__estimate) < 1
  ) %>% 
  ggplot(aes(x = cat_pre_wt_log_scale, y = less_toxic_start_f1.sl.is_state1.foo__estimate, color = beta, group = beta)) + 
  geom_smooth(method = "lm") + 
  geom_pointrange(
    aes(ymin = less_toxic_start_f1.sl.is_state1.foo__estimate - 2 * less_toxic_start_f1.sl.is_state1.foo__se,
        ymax = less_toxic_start_f1.sl.is_state1.foo__estimate + 2 * less_toxic_start_f1.sl.is_state1.foo__se,
        alpha = sig),
    position = position_jitterdodge(jitter.width = 0.5)
  ) 





d4 <- data_prep(sample(d3$rep_id, 1))
m <- issf( 
  case ~ 
    state:logsl:less_toxic_start_f +
    state:sl:less_toxic_start_f +
    is_state1:sl:less_toxic_start_f:t + 
    is_state1:sl:less_toxic_start_f:foo +
    strata(step_id),
  data = d4 %>% 
    mutate(
      foo = memory(less_toxic_start, 
                   step_id, 
                   case,
                   weight_FUN = FUN, 
                   k = 0.01),
      t = step_id / 10 / 24
    ) %>% 
    filter(
      step_id > 0
    ), 
  update = FALSE, 
  keep_data = TRUE
);summary(m) 

memory_decay(k = 0.5) / 10

m <- issf( 
  case ~ 
    state:logsl:less_toxic_f +
    state:sl:less_toxic_f +
    is_state1:sl:less_toxic_f:foo +
    strata(step_id),
  data = d4 %>% 
    mutate(
      foo = memory(less_toxic, 
                   step_id, 
                   case,
                   weight_FUN = FUN, 
                   k = 0.5)
    ) %>% 
    filter(
      step_id > 0
    ), 
  update = FALSE, 
  keep_data = TRUE
);summary(m) 
d4 <- data_prep(55)




m2 <- issf( 
  case ~ 
    #state:cos_theta_pi:less_toxic_f + 
    #state:cos_2theta:less_toxic_f + 
    state:sl:less_toxic_f +
    state:logsl:less_toxic_f +
    is_state1:sl:less_toxic_f:foo +
    strata(step_id),
  data = d4 %>% 
    mutate(
      foo = memory(less_toxic, 
                   step_id, 
                   case,
                   weight_FUN = FUN, 
                   k = 0.1)
    ) %>% 
    filter(!is.na(foo)), 
  update = FALSE, 
  keep_data = TRUE
);summary(m2)


LRT(m, m2)



debug(bar2)




d4 <- data_prep(55)

f <- function(theta){
  res <- issf( 
    case ~ 
      state:sl:less_toxic_start_f +
      state:logsl:less_toxic_start_f +
      is_state1:sl:less_toxic_start:foo +
      strata(step_id),
    data = d4 %>% 
      mutate(
        foo = memory(less_toxic_start, 
                     step_id, 
                     case,
                     weight_FUN = FUN, 
                     k = exp(theta))
      ) %>% 
      filter(
        step_id > 0
      ), 
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
  lower = -8, 
  upper = 3
) -> out

out

memory_decay <- function(k, decay_prop = 0.05){
  # decay_prop = amount of memory left
  (1 / decay_prop - 1) / k
  
  #log(decay_prop) * -1 / k
}

memory_decay(k = exp(-3.2)) / 10
memory_decay(0.00684221) / 10


d3 %>% 
  select(beta, rep_id, k1) %>% 
  as.data.frame() %>% 
  arrange(-abs(log(k1)))






