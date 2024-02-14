m <- glmmTMB(
  RGR ~ 
    cat_pre_wt_log_scale + 
    #I(cat_pre_wt_log_scale^2) +
    cat_pre_wt_log_scale : (var_trt + beta) + 
    var_trt + beta + 
    scale(mean_toxic_conc) + 
    scale(log(area_herb+1)) + 
    scale(sl_mean_obs)  + 
    scale(var_toxic_12) + 
    (1|session_id),
  data = d %>% 
    filter(
      mean_trt_numeric == 1 & var_trt != "constant"
    ) %>% 
    filter(
      session_id != 1
    )
); summary(m)


m <- glmmTMB(
  RGR ~ 
    var_trt + beta + 
    (var_trt + beta) : cat_pre_wt_log_scale + 
    cat_pre_wt_log_scale + I(cat_pre_wt_log_scale^2) +  
    (1|session_id),
  data = d %>% 
    filter(
      var_trt != "constant"
    ) %>% 
    filter(
      !is.na(sl_mean_obs)
    )
); summary(m)
plot_model(m, type = "eff", terms = c("cat_pre_wt_log_scale[all]", "beta"), show.data = TRUE)
check_model(m)

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

car::Anova(m, type = "III")


m <- glmmTMB(
  sl_mean_obs ~ 
    var_trt + beta + cat_pre_wt_log_scale   + #I(cat_pre_wt_log_scale^2) +
    (1|session_id),
  family = Gamma(link = "log"),
  data = d %>% 
    filter(
      var_trt != "constant" #& !rep_id %in% problem_ids
    ) %>% 
    mutate(
      problem = ifelse(rep_id %in% problem_ids, 1, 0)
    ),
); summary(m) #DHARMa::simulateResiduals() %>% plot()

m <- glmmTMB(
  prob_move_obs ~ 
    (var_trt + beta) * cat_pre_wt_log_scale + #I(cat_pre_wt_log_scale^2) + 
    (1|session_id),
  data = d %>% 
    filter(
      var_trt != "constant"
    ),
  beta_family(),
); summary(m)

plot_model(m, type = "eff", terms = c("cat_pre_wt_log_scale", "beta"))

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
    var_trt + beta * cat_pre_wt_log_scale +
    (1|session_id),
  data = d %>% 
    filter(
      var_trt != "constant"
    ),
  family = beta_family(),
); summary(m)
car::Anova(m, type = "III")

m <- glmmTMB(
  toxic ~ 
    (var_trt + beta) + cat_pre_wt_log_scale + 
    (1|session_id),
  family = gaussian(),
  data = d %>% 
    filter(
      var_trt != "constant"
    ),
); summary(m)


plot_model(m, type = "eff", terms = c("cat_pre_wt_log_scale", "beta"), show.data = TRUE)

m <- glmmTMB(
  ava_mean_toxin ~ 
    (var_trt + beta) + cat_pre_wt_log_scale + 
    (1|session_id),
  data = d %>% 
    filter(
      var_trt != "constant"
    ),
); summary(m)


m <- glmmTMB(
  ud_estimate ~ 
    (var_trt + beta) + cat_pre_wt_log_scale + 
    (1|session_id),
  family = Gamma(link = "log"),
  data = d %>% 
    filter(
      var_trt != "constant"
    ),
); summary(m)


d %>% 
  mutate(sl_mean = (shape * scale)) %>% 
  filter(kappa1 > 0 & kappa2 > 0) %>% 
  ggplot(aes(x = var_trt, color = beta, y = log(sl_mean_obs), fill = beta)) + 
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









