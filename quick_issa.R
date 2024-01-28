source("helper_functions/init_analysis.R")
library(survival)
library(amt)

#mclogit::mclogit

ids <- ref_data %>% 
  filter(
    !is.na(camera_cutoff)
  ) %>% 
  filter(
    !is.na(syn_id)
  ) %>% 
  select(rep_id) %>% 
  unlist() %>% 
  unname()

ids <- ids[-92]


ids[92]

fit_list <- pb_par_lapply(
  ids, 
  function(i, ref_data){
    ID <- i
    d <- fetch_events(ID) %>% 
      insert_gaps() %$% 
      move_seq(head_x, head_y, r_thresh = 10, inherit.theta = FALSE) %>% 
      mutate(
        theta_for_fitting = ifelse(r_threshed == 0, NA, theta_rel)
      ) %>% 
      filter(r_threshed > 0) %>% 
      add_random_steps(n = 50,
                       sl_distr = amt::fit_distr(r_threshed, "gamma"),
                       ta_distr = amt::fit_distr(theta_for_fitting, "vonmises")
      ) %>% 
      flag_invalid_steps(remove = TRUE) %>% 
      mutate(
        moved = ifelse(r_threshed == 0, 0 , 1), 
        toxic = read_value(x2, y2, c(1000, 1000), 
                           ref_img = fetch_trt_spec(ID, .ref_data = ref_data)), 
        sl = r_threshed, 
        logsl = ifelse(moved == 1, log(r_threshed), 0)
      )
    
    out <- issf(
      case ~ 
        toxic + 
        (cos_theta + sl + logsl) + 
        strata(step_id),
      data = d 
    )
    
    return(out)
  }, 
  ref_data = ref_data,
  cores = 8, 
  inorder = TRUE
)


fit_list[[1]]

ID <- 100
d <- fetch_events(ID) %>% 
  insert_gaps() %$% 
  move_seq(head_x, head_y, r_thresh = 10, inherit.theta = TRUE) %>% 
  mutate(
    theta_for_fitting = ifelse(r_threshed == 0, NA, theta_rel)
  ) %>% 
  add_random_steps(n = 50,
    sl_distr = fit_zigamma2(r_threshed),
    ta_distr = fit_distr(theta_for_fitting, "vonmises")
  ) %>% 
  flag_invalid_steps(remove = TRUE) %>% 
  mutate(
    moved = ifelse(r_threshed == 0, 0 , 1), 
    toxic = read_value(x2, y2, c(1000, 1000), ref_img = fetch_trt_spec(ID)), 
    sl = r_threshed, 
    logsl = ifelse(moved == 1, log(r_threshed), 0)
  )

survival::clogit(
  case ~ 
    toxic + moved + 
    moved:(sin_theta + sl + logsl) + 
    strata(step_id),
  data = d 
)























z <- fit_list %>%
  lapply(function(x){
    c(coef(x$model), 
      #"p" = x$sl$params$p, 
      "shape" = x$sl_updated$params$shape, 
      "scale" = x$sl_updated$params$scale, 
      "kappa" = x$ta_updated$params$kappa)
  }) %>% 
  do.call("rbind", .) %>% 
  as.data.frame() %>% 
  cbind("rep_id" = ids)


z <- ref_data %>% 
  right_join(z, by = "rep_id")

z <- detection_report(z$rep_id) %>% 
  select(repID, n_keypoints, frames) %>% 
  mutate(
    prop_kp = n_keypoints / frames,
    rep_id = as.character(repID)
  ) %>% 
  select(-repID) %>% 
  right_join(z, by = "rep_id")

z %>% View()

z %>% 
  filter(prop_kp > 0.7) %>% 
  ggplot(aes(x = var_trt, y = toxic, group = beta, color = beta)) + 
  geom_pointrange(stat = "summary", color = "black", position = position_dodge(0.7)) + 
  geom_point(position = position_jitterdodge())


names(z)
z

z %>% 
  filter(prop_kp > 0.7) %>% 
  filter(scale > 0 & abs(scale) < 1000) %>% 
  ggplot(aes(x = log(cat_pre_wt), y = kappa, color = beta)) + 
  geom_point() +
  geom_smooth(method = "lm")


z %>% 
  filter(prop_kp > 0.7) %>% 
  filter(scale > 0 & abs(scale) < 1000) %>% 
  ggplot(aes(x = var_trt, y = kappa, group = beta, color = beta)) + 
  geom_pointrange(stat = "summary", color = "black", position = position_dodge(0.7)) + 
  geom_point(position = position_jitterdodge())






o$model
o$sl_updated
o$ta
o$ta_updated
reload()


iterate_random_steps <- function(start, issf_fit, n = 100){
  out.list <- vector(mode = "list", length = n+1)
  out.list[[1]] <- start
  ra<- random_numbers(issf_fit$ta_updated, 10^5)
  rr<- random_numbers(issf_fit$sl, 10^5)
  new_x2 <- start$x2
  new_y2 <- start$y2
  for(i in seq_len(n)){
    new <- TRUE
    while(new | !is_between(new_x2, c(0, 1000)) | 
                 !is_between(new_y2, c(0, 1000))
          ){
      out.list[[i+1]] <- out.list[[i]] %>% 
        add_random_steps(
          data = .,
          n = 1, 
          keep_obs = FALSE, 
          id = 1, 
          x_start = .$x2,
          y_start = .$y2,
          direction_start = .$theta_abs,
          sl_distr = NULL, 
          ta_distr = NULL, 
          ta_rand = ra,
          sl_rand = rr)
      new <- FALSE
      new_x2 <- out.list[[i+1]]$x2
      new_y2 <- out.list[[i+1]]$y2
    }
  }
  do.call("rbind.fill", out.list)
}

creat_start <- function(x, y, theta){
  data.frame("step_id" = 1, "x2" = x, "y2" = y, "theta_abs" = theta)
}

debug(iterate_random_steps)
z <- creat_start(500, 500, runif(1, 0, 2 * pi)) %>% 
  iterate_random_steps(start = ., issf_fit = o, n = 1000)






fetch_events(99) %$%
  move_seq(centroid_x, centroid_y, r_thresh = 10) %>% 
  plot_track(x1, y1)

z %>% 
  #.[500:1000,] %>% 
  ggplot(aes(x = x2, y = y2)) + 
  geom_path()  
  geom_count()



z %>% 
  #.[7000:8000,] %>% 
  ggplot(aes(x = x2, y = y2)) + 
  geom_density_2d_filled() + 
  geom_path()




debug(fit_cengamma)

reload()


dcengamma(rcengamma(10000, 1, 30, 10), 1, 30, 10, log = TRUE) %>% sum()



pgamma(
  5,shape =  1, scale = 400
)

thresh <- 10
test_d <- data.frame(
  "step_id" = 1:3000,
  "r" = rgamma(3000, shape = 1, scale = 250),
  "theta_rel" = random_numbers(make_genvonmises(0.5, 1), 3000)
) %>%
  mutate(
    r_threshed = ifelse(r > thresh, r, 0)
  ) %>% 
  add_random_steps(data = .,
                   n = 300L,
                   id = .$step_id,
                   x_start =  rep(NA, 3000),
                   y_start = rep(NA, 3000),
                   direction_start = rep(NA, 3000),
                   sl_distr = make_gamma(3, 50),
                   ta_distr = make_genvonmises(0.5, 1)) %>%
  #flag_invalid_steps(remove = TRUE) %>% 
  append_gamma_estimators() %>% 
  mutate(
    moved = ifelse(r > thresh, 1, 0),
  ) %>% 
  append_genvonmises_estimators(na_as_zero = TRUE) 



o <- issf(case ~ 
            moved:(cos_theta_pi + cos_2theta) + (sl + logsl) + 
            strata(step_id),
          data = test_d, 
          shape_estimator = "logsl",
          scale_estimator = "sl",
          kappa1_estimator = "moved:cos_theta_pi",
          kappa2_estimator = "moved:cos_2theta"
)

summary(test_d$r)
o
o$sl
o$sl_updated
o$ta_updated

debug(update_distr)
update_distr(o)


o$ta_updated
o$sl_updated

o$sl_updated

debug(issf)



test_d <-  fetch_events(90) %$% 
  move_seq(head_x, head_y, r_thresh = 10) %>% 
  add_random_steps(data = ., 
                   n = 100L,
                   id = .$step_id,
                   sl_distr = fit_gamma(.$r), 
                   ta_distr = fit_genvonmises(.$theta_rel)) %>% 
  flag_invalid_steps(remove = TRUE) %>% 
  mutate(
    invalid_step = ifelse(valid_step, 0, 1)
  ) %>% 
  append_estimators(na_as_zero = TRUE)



o <- issf(case ~ 
            moved:(cos_theta_pi + cos_2theta) + sl + logsl + 
            strata(step_id),
          data = test_d, 
          shape_estimator = "logsl",
          scale_estimator = "sl",
          kappa1_estimator = "moved:cos_theta_pi",
          kappa2_estimator = "moved:cos_2theta", 
          remove_invalid = TRUE
)



o$ta_updated
o$sl_updated


