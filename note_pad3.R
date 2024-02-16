m <- glmmTMB(
  RGR ~ 
    cat_pre_wt_log_scale + 
    I(cat_pre_wt_log_scale^2) +
    #cat_pre_wt_log_scale *
    #var_trt + beta + 
    scale(mean_toxic_conc) + 
    scale(log(area_herb+1)) + 
    scale(sl_mean_obs)  +
    scale(var_toxic_12) + 
    (1|session_id),
  data = d %>% 
    filter(
      mean_trt_numeric == 1 & var_trt != "constant"
    )
); summary(m)



z <- d %>% 
  filter(
    #var_trt != "constant"
    mean_trt_numeric == 1
  ) %>% 
  mutate(
    beta = ifelse(is.na(beta), "const", as.character(beta)),
    beta_red = ifelse(as.character(beta) == 5, 1, 0),
    beta_blue = ifelse(as.character(beta) == -5, 1, 0),
    beta_white = ifelse(as.character(beta) == 0, 1, 0),
    var_low = ifelse(var_trt == "var_low", 1, 0)
  )

m <- glmmTMB(
  RGR ~ 
    has_var:(beta_red * cat_pre_wt_log_scale + beta_white * cat_pre_wt_log_scale +
               var_high * cat_pre_wt_log_scale) +
    has_var/(cat_pre_wt_log_scale) + #cat_pre_wt_log_scale +
    #cat_pre_wt_log_scale + I(cat_pre_wt_log_scale^2) +  
    (1|session_id),
  data = d %>% 
    filter(
      #var_trt != "constant"
      mean_trt_numeric == 1
    ) %>% 
    mutate(
      beta = ifelse(is.na(beta), "const", as.character(beta)),
      beta_red = ifelse(as.character(beta) == 5, 1, 0),
      beta_blue = ifelse(as.character(beta) == -5, 1, 0),
      beta_white = ifelse(as.character(beta) == 0, 1, 0),
      var_high = ifelse(var_trt == "high_var", 1, 0)
    )
); summary(m)
plot_model(m, type = "eff", terms = c("cat_pre_wt_log_scale[all]", "beta"), show.data = TRUE)
check_model(m)
DHARMa::simulateResiduals(m) %>% plot()



plot_model(m, type = "eff", terms = c("has_var", "beta_red","cat_pre_wt_log_scale[-1,1]"))
plot_model(m, type = "eff", terms = c("cat_pre_wt_log_scale[all]", "beta_red", "has_var"))


m <- glmmTMB(
  RGR ~ 
    (var_trt + beta) * cat_pre_wt_log_scale +
    cat_pre_wt_log_scale + I(cat_pre_wt_log_scale^2) +  
    (1|session_id),
  data = d %>% 
    filter(
      var_trt != "constant" &
      mean_trt_numeric == 1
    ) %>% 
    mutate(
      beta = ifelse(is.na(beta), "const", as.character(beta)),
      beta_red = ifelse(as.character(beta) == 5, 1, 0),
      beta_blue = ifelse(as.character(beta) == -5, 1, 0),
      beta_white = ifelse(as.character(beta) == 0, 1, 0),
      var_high = ifelse(var_trt == "high_var", 1, 0)
    )
); summary(m)





d %>% 
  filter(!is.na(beta)) %>% 
  ggplot(aes(x = cat_pre_wt_log_scale, y = RGR, color = beta)) + 
  geom_point() + 
  geom_smooth()



m <- glmmTMB(
  mean_toxic_conc ~ 
    (beta + var_trt) * cat_pre_wt_log_scale  + 
    (1|session_id),
  data = d %>% 
    filter(
      var_trt != "constant"
    ) 
); summary(m)


m <- glmmTMB(
  mean_toxic ~ 
    var_trt + (beta) * cat_pre_wt_log_scale  + 
    (1|session_id),
  family = beta_family(),
  data = d %>% 
    filter(
      var_trt != "constant"
    ) %>% 
    mutate(
      mean_toxic = adjust_prop(mean_toxic, nudge.size = 0.01)
    )
); summary(m)




m <- glmmTMB(
  var_toxic_12 ~ 
    (beta + var_trt) + cat_pre_wt_log_scale  + 
    (1|session_id),
  data = d %>% 
    filter(
      var_trt != "constant"
    ) 
); summary(m)


car::Anova(m, type = "III")


m <- glmmTMB(
  sl_mean_obs ~ 
    var_trt + (beta) * cat_pre_wt_log_scale + #I(cat_pre_wt_log_scale^2) +
    (1|session_id),
  family = Gamma(link = "log"),
  data = d %>% 
    filter(
      var_trt != "constant" #& !rep_id %in% problem_ids
    ) %>% 
    mutate(
      #sl_kurt_obs = 6 / shape + 3
    ),
); summary(m) #DHARMa::simulateResiduals() %>% plot()

m <- glmmTMB(
  prob_move_obs ~ 
    (var_trt + beta) + cat_pre_wt_log_scale + #I(cat_pre_wt_log_scale^2) + 
    (1|session_id),
  data = d %>% 
    filter(
      var_trt != "constant"
    ) %>% 
    filter(prob_move_obs > 0),
  beta_family(),
); summary(m)

plot_model(m, type = "eff", terms = c("cat_pre_wt_log_scale", "var_trt"))

glmmTMB(
  log(area_herb+1) ~ 
    (var_trt + beta) * cat_pre_wt_log_scale + 
    (1|session_id),
  data = d %>% 
    filter(
      var_trt != "constant"
    ),
  gaussian(),
) %>% summary()


m <- glmmTMB(
  on_toxic ~ 
    (var_trt + beta * cat_pre_wt_log_scale) + 
    (1|session_id),
  data = d,
  family = beta_family(),
); summary(m)
car::Anova(m, type = "III")

m <- glmmTMB(
  toxic ~ 
    var_trt + (beta) * cat_pre_wt_log_scale,
  family = gaussian(),
  data = d %>% 
    filter(
      var_trt != "constant"
    ) %>% 
    mutate(
      session_id = as.factor(session_id)
    ),
); summary(m)


m <- glmmTMB(
  kappa1 ~ 
    var_trt + beta + cat_pre_wt_log_scale,
  family = gaussian(),
  data = d %>% 
    filter(
      var_trt != "constant"
    ) %>% 
    mutate(
      session_id = as.factor(session_id)
    ),
); summary(m)


ddist(make_genvonmises(0.3, 0.25)) %>% plot()


plot_model(m, type = "eff", terms = c("cat_pre_wt_log_scale", "beta"), show.data = TRUE)

m <- glmmTMB(
  ava_mean_toxin ~ 
    (var_trt + beta) * cat_pre_wt_log_scale + 
    (1|session_id),
  data = d %>% 
    filter(
      var_trt != "constant"
    ),
); summary(m)


m <- glmmTMB(
  ud_estimate ~ 
    beta + (var_trt) * cat_pre_wt_log_scale + 
    (1|session_id),
  family = Gamma(link = "log"),
  data = d %>% 
    filter(
      var_trt != "constant"
    ),
); summary(m)


d %>% 
  mutate(sl_mean = (shape * scale)) %>% 
  ggplot(aes(x = sl_mean / 1000 * 12, y = sl_mean_obs/1000*12)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") + 
  theme_bw(base_size = 15) +
  labs(x = "Selection free mean step length (cm)", y = "Observed mean step length (cm)")



d %>% 
  mutate(sl_kurt = 6/shape + 3) %>% 
  ggplot(aes(x = sl_kurt, y = sl_kurt_obs)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") + 
  theme_bw(base_size = 15) +
  labs(x = "Selection free kurt step length (cm)", y = "Observed kurt step length (cm)")


d %>% 
  mutate(sl_mean = (shape * scale)) %>% 
  ggplot(aes(x = sl_mean / 1000 * 12, y = sl_mean_obs/1000*12)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") 


d %>% 
  mutate(sl_mean = (shape * scale)) %>% 
  filter(kappa1 > 0 & kappa2 > 0) %>% 
  ggplot(aes(x = var_trt, color = beta, y = kappa2, fill = beta)) + 
  geom_pointrange(stat = "summary", 
                  position = position_dodge(width = 0.5),
                  color = "black") + 
  geom_point(position = position_jitterdodge(jitter.height = 0))



trans_mat <- matrix(c(0.03,0.03,
                      0.03*0.8,0.03*0.8), 
                    ncol = 2, nrow = 2, 
                    dimnames = list(
                      c("start_0","start_1"),
                      c("end_0","end_1")
                    ))


ref <- as.matrix(fetch_trt_spec(1))

t_max <- 10000



ids <- ref_data %>% 
  filter(!is.na(syn_id)) %>% 
  select(rep_id) %>% 
  unlist(use.names = FALSE)




w <- pb_par_lapply(
  rep(c(5,0,-5), 500),
  function(i, trans_mat){
    n <- 100
    o <- sim_mean_toxic(
      as.matrix(as.cimg(syn_spec(n, beta = i, plot = FALSE))), 
      trans_mat = trans_mat, 
      t_max = 150, 
      n = n) 
    cbind("beta" = i, o)
  },
  trans_mat = trans_mat,
  cores = 1, 
  inorder = FALSE
) %>% 
  do.call("rbind",.)


w %>% 
  ggplot(aes(x = as.factor(beta), y = (mean))) + 
  geom_point(position = "jitter", aes(color = as.factor(beta))) + 
  geom_pointrange(stat = "summary") + 
  facet_wrap(~start)


glmmTMB(mean ~ as.factor(beta), 
        data = w) %>% summary()



z<-pb_par_lapply(
  ids,
  function(i, ref_data, trans_mat){
    sim_mean_toxic(
      as.matrix(fetch_trt_spec(i, quiet = TRUE, .ref_data = ref_data)), 
      trans_mat = trans_mat, 
      t_max = 10000) 
  }, ref_data = ref_data, 
  trans_mat = trans_mat,
  cores = 8
) %>% 
  do.call("rbind", .)



z2 <- as.data.frame(z) %>% 
  cbind(
    "rep_id" = ids
  ) %>% 
  left_join(
    ref_data %>% 
      select(rep_id, beta),
    by = "rep_id"
  )


z2 %>% 
  ggplot(aes(x = beta, y = V1)) + 
  geom_point(position = "jitter") + 
  geom_pointrange(stat = "summary")
lm(V1~ as.factor(beta), z2) %>% summary()





sim_mean_toxic <- function(ref, trans_mat, t_max = 10000, n = 12){
  x <- y <- rep(NA_real_, t_max)
  x[1] <- sample.int(n, 1)
  y[1] <- sample.int(n, 1)
  
  
  for(t in 2:t_max){
    x[t] <- x[t-1] + (-1)^rbinom(n = 1, size = 1, prob = 0.5)
    y[t] <- y[t-1] + (-1)^rbinom(n = 1, size = 1, prob = 0.5)
    
    if(x[t] > n){
      x[t] <- x[t] - n
    }
    if(x[t] <= 0){
      x[t] <- n + x[t]
    }
    
    if(y[t] > n){
      y[t] <- y[t] - n
    }
    if(y[t] <= 0){
      y[t] <- n + y[t]
    }
    
    val_og <- ref[x[t-1], y[t-1]]
    val_new <- ref[x[t], y[t]]
    p_move <- trans_mat[val_og+1, val_new+1]
    
    if(rbinom(1, 1, p_move) == 0){
      x[t] <- x[t-1]
      y[t] <- y[t-1]
    }
  }
  
  out <- vapply(seq_along(x), 
         function(i){
           ref[x[i], y[i]]
         }, FUN.VALUE = numeric(1)) %>% 
    mean()
  return(data.frame("mean" = out, "start" = ref[x[1], y[1]]))
}








ids <- d %>% 
  filter(var_trt != "constant") %>% 
  filter(!is.na(camera_cutoff)) %>% 
  select(rep_id) %>% 
  filter(!rep_id %in% problem_ids) %>% 
  unlist(use.names = FALSE)












out <- pb_par_lapply(
  ids,
  function(i, ref_data){
    fetch_events(i) %>% 
      clean_events(ref_data = ref_data) %>% 
      insert_gaps() %>% 
      mutate(
        head_x = ifelse(score >=0.5, head_x, NA),
        head_y = ifelse(score >=0.5, head_y, NA)
      ) %$%
      move_seq(head_x, head_y, r_thresh = 0, inherit.theta = FALSE) %>% 
      as.data.frame() %>% 
      cbind("rep_id" = i) %>% 
      mutate(
        end_toxic = read_value(x2, y2, c(1000, 1000),
                           ref_img = fetch_trt_spec(i, .ref_data = ref_data, quiet = TRUE)),
        start_toxic = read_value(x2, y2, c(1000, 1000),
                               ref_img = fetch_trt_spec(i, .ref_data = ref_data, quiet = TRUE))
      )
  },
  ref_data = ref_data,
  cores = 8,
  inorder = FALSE
) %>% 
  do.call("rbind", .)





m <- glmmTMB(
  y ~ 
    (beta + var_trt) * cat_pre_wt_log_scale +
    (1|session_id), 
  family = gaussian(), 
  data = d %>% 
    mutate(
      toxic = toxic,
      toxic_time = `toxic:time`,
      y = toxic_time
    )
); summary(m)

plot_model(m, type = "eff", terms = c("cat_pre_wt_log_scale", "beta"))



m <- glmmTMB(
  r ~ 
  var_trt * cat_pre_wt_log_scale *
  time + 
    beta * cat_pre_wt_log_scale +
  (1|session_id) + (1 | rep_id) + ar1(step_id + 0 | rep_id), 
  family = Gamma(link = "log"), 
  data = 
    out %>% 
    left_join(d, by = "rep_id") %>% 
    filter(
      !rep_id %in% problem_ids
    ) %>% 
    mutate(
      time = step_id / 1200,
      step_id = as.factor(step_id)
    ) %>% 
    filter(
      !is.na(r)
    ) %>% 
    filter(
      !is.na(beta)
    )
); summary(m)







m <- glmmTMB(
  toxic ~ 
    (beta + var_trt) * cat_pre_wt_log_scale *
    scale(step_id) + 
    (1|session_id/rep_id), 
  family = binomial(), 
  data = 
    out %>% 
    left_join(d %>% 
                select(-toxic), by = "rep_id") %>% 
    filter(
      !rep_id %in% problem_ids
    )
); summary(m)

plot_model(m, type = "eff", terms = c("cat_pre_wt_log_scale", "beta"))

check_overdispersion(m)


plot_model(m, 
           type = "eff", 
           terms = c("cat_pre_wt_log_scale", "beta", "time[0, 0.5, 1]")) + 
  scale_y_continuous(trans = "log")








m <- glmmTMB(
  r ~ 
    (beta + var_trt) * cat_pre_wt_log_scale + 
    (1|session_id), 
  family = Gamma(link = "log"), 
  data = 
    out %>% 
    filter(step_id < 24 * 10) %>% 
    group_by(rep_id) %>% 
    summarise(r = mean(r, na.rm = TRUE)) %>% 
    left_join(d, by = "rep_id") %>% 
    filter(
      !rep_id %in% problem_ids
    ) 
); summary(m)





d %>% 
  filter(!is.na(camera_cutoff)) %>% 
  .$rep_id %>% 
  lapply(function(x){
    fetch_events(x) %>% 
      clean_events() %>% 
      filter(score > 0.9) %>% 
      summarise(p = sum(!is.na(head_x)) / unique(frames)) %>% 
      unlist()
  }) %>% 
  do.call("c",.) -> x


mean(x > 0.7)
summarise_vec(x)
hist(x)

d %>% filter(var_trt!= "constant") %>% filter(is.na(var_toxic_12) & !is.na(sl_mean_obs)) %>% View()




d2


subm_rgr <- glmmTMB(
  RGR ~ 
    cat_pre_wt_log_scale + 
    cat_pre_wt_log_scale_sq +
    mean_toxic_conc_scale + 
    area_herb_log_scale + 
    sl_mean_obs_log_scale +
    var_toxic_12_scale + 
    (1|session_id),
  data = d2
); summary(subm_rgr)

subm_var_toxic <- glmmTMB(
  var_toxic_12_scale ~ 
    beta + var_trt + cat_pre_wt_log_scale + on_toxic_logit_scale + ava_mean_toxin_scale +
    sl_mean_obs_log_scale + 
    (1|session_id),
  data = d2
); summary(subm_var_toxic)



subm_sl <- glmmTMB(
  sl_mean_obs_log_scale ~ 
    var_trt + beta * cat_pre_wt_log_scale + on_toxic_logit_scale + 
    (1|session_id),
  data = d2,
); summary(subm_sl)



subm_toxin_ingested <- glmmTMB(
  mean_toxic_conc_scale ~ 
    beta + var_trt * cat_pre_wt_log_scale + on_toxic_logit_scale + ava_mean_toxin_scale +
    (1|session_id),
  data = d2 
); summary(subm_toxin_ingested)


subm_on_toxic <- glmmTMB(
  on_toxic_logit_scale ~ 
    toxic_select_scale + ava_mean_toxin_scale + 
    (1|session_id),
  data = d2
); summary(subm_on_toxic)


subm_herb <- glmmTMB(
  area_herb_log_scale ~ 
    (var_trt + beta) + cat_pre_wt_log_scale + on_toxic_logit_scale + ava_mean_toxin_scale +
    (1|session_id),
  data = d2,
); summary(subm_herb)


car::Anova(m, type = "III")

subm_select <- glmmTMB(
  toxic_select_scale ~ 
    beta + var_trt*cat_pre_wt_log_scale + ava_mean_toxin_scale + 
    (1|session_id),
  data = d2,
); summary(subm_select)

contrast_by_pre_wt(subm_ava, "beta", type = ("trt"))


plot_model(subm_ava, type = "eff", terms = c("cat_pre_wt_log_scale", "beta"), show.data = TRUE)

subm_ava <- glmmTMB(
  ava_mean_toxin_scale ~ 
    var_trt + beta * cat_pre_wt_log_scale + 
    (1|session_id),
  data = d2,
); summary(subm_ava)




sem_fit <- psem(
    subm_rgr,
    subm_var_toxic,
    subm_on_toxic,
    subm_toxin_ingested,
    subm_select,
    subm_ava,
    subm_herb,
    subm_sl, 
    data = d2
)
sem_summary <- summary(sem_fit)

sem_summary

jutila.multigroup <- multigroup(utila, group = "grazed")

sem_summary <- suppressWarnings(suppressMessages(summary(sem_fit)))

piecewiseSEM::stdCoefs(list(subm_on_toxic), data = d2, standardize.type = "Menard.OE")

piecewiseSEM:::summary.psem




