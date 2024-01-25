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


z %>% 
  filter(scale < 0)

update_gamma
update_distr()

fit_list[[3]]$ta_updated

fit_list[[2]]$data$logsl

fit_list[[100]]$data %>% 
  issf(
    case ~ 
      toxic + 
      (cos_theta + sl + logsl) + 
      strata(step_id),
    data = .
  ) %>% .$ta_updated


source("helper_functions/init.R")


fit_list[[2]]$data %>% filter(case) %$% plot(x1, y1)


fetch_data_dict(2) %>% mask_area() %>% plot()

f<- forward_plot(fetch_data_dict(50), 100)

f()

fit_list[[43]]$data %>% 
  filter(case) %>% 
  .$theta_rel %>% 
  fit_bivonmises()

source("helper_functions/init.R")

rgenvonmises(1000, 2, 1) %>% hist(nclass=100)
dgenvonmises(x, 1, 0.5) %>% plot(x, .)
debug()


fit_genvonmises(rgenvonmises(10000, 0.3, 0.2), 
                init = c(1,0.5), parscale = c(1000, 1000))


circular::dgenvonmises(circular::circular(x), 
                       mu1 = circular::circular(0),
                       mu2 = circular::circular(pi-0.01), 
                       kappa1 = 0.1, 
                       kappa2 = 0) %>% 
  plot(x, .)






rbivonmises(100000, 100, 10) %>% 
  fit_bivonmises() 


make_vonmises_distr(1) %>% random_numbers(n = 1000) %>% hist()

ids[43]
fit_list[[43]]$data %>% filter(case) %>% .$theta_rel %>% hist()

fit_list[[40]]$data %>% filter(case) %>% 
  ggplot(aes(x = logsl, y = abs(theta_rel))) + 
  geom_point() + 
  geom_smooth(method = "lm")


source("helper_functions/init.R")

update_gamma

test_d <- data.frame(
      "step_id" = 1:3000,
      "r" = rgamma(3000, shape = 1, scale = 100),
      "theta_rel" = random_numbers(make_genvonmises(0.5, 1), 3000)
    ) %>%
      add_random_steps(n = 100L,
                       id = step_id,
                       x_start =  rep(NA, 3000),
                       y_start = rep(NA, 3000),
                       direction_start = rep(NA, 3000),
                       sl_distr = make_gamma_distr(1, 100),
                       ta_distr = make_genvonmises(1, 0)) %>%
      mutate(
        sl = r,
        logsl = log(r),
      ) %>% 
  append_genvonmises_estimators()

o <- issf(case ~ 
       (cos_theta_pi + cos_2theta + sl + logsl) + 
       strata(step_id),
     data = test_d)

o$ta_updated



o$ta_updated
o$sl_updated

make_bivonmises(0.2, 3) %>% random_numbers(1000) %>% hist()


fetch_events(5) %$% 
  move_seq(centroid_x, centroid_y, 20) %$%
  hist(theta_rel, nclass = 30)


f <- function(theta, k1, k2){
  exp(
    k1 * cos(theta)
  )  + 
  exp(k2 * cos(theta - pi))
}


x <- seq(-pi, pi, by = 0.01)
plot(
  x, 
  f(
    x, 
    3 * a,
    3.1 * a
  )
)

a <- 5



















