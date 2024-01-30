source("helper_functions/init_analysis.R")
library(survival)
library(amt)

# ids <- ref_data %>% 
#   filter(
#     !is.na(camera_cutoff)
#   ) %>% 
#   filter(
#     !is.na(syn_id)
#   ) %>% 
#   select(rep_id) %>% 
#   unlist() %>% 
#   unname()
# 
# 
# 
# ids <- ids[-c(17, 92, 115)]
# 
# 
# fit_list <- pb_par_lapply(
#   ids, 
#   function(i, ref_data){
#     # message(sprintf("%s \n", i))
#     ID <- i
#     d <- fetch_events(ID) %>% 
#       clean_events(ref_data = ref_data) %>% 
#       insert_gaps() %$% 
#       move_seq(head_x, head_y, r_thresh = 30, inherit.theta = FALSE) %>% 
#       add_random_steps(n = 100L,
#                        sl_distr = fit_gamma(.$r),
#                        ta_distr = fit_genvonmises(.$theta_rel)
#       ) %>% 
#       flag_invalid_steps(remove = TRUE) %>% 
#       mutate(
#         toxic = read_value(x2, y2, c(1000, 1000), 
#                            ref_img = fetch_trt_spec(ID, .ref_data = ref_data))
#       ) %>% 
#       append_estimators(na_as_zero = TRUE)
#     
#     out <- issf(
#       case ~ 
#         toxic + 
#         moved : (cos_theta_pi + cos_2theta) +
#         sl + logsl + 
#         strata(step_id),
#       data = d, 
#       scale_estimator = "sl", 
#       shape_estimator = "logsl", 
#       kappa1_estimator = "moved:cos_theta_pi",
#       kappa2_estimator = "moved:cos_2theta"
#     )
#     
#     return(out)
#   }, 
#   ref_data = ref_data,
#   cores = 8, 
#   inorder = TRUE
# )
# names(fit_list) <- ids
# 
# 
# ids2 <- ref_data %>% 
#   filter(
#     !is.na(camera_cutoff)
#   ) %>% 
#   filter(
#     var_trt == "constant"
#   ) %>% 
#   select(rep_id) %>% 
#   unlist() %>% 
#   unname()
# 
# ids2 <- ids2[-c(28)]
# 
# fit_list2 <- pb_par_lapply(
#   ids2, 
#   function(i, ref_data){
#     # message(sprintf("%s \n", i))
#     ID <- i
#     d <- fetch_events(ID) %>% 
#       clean_events(ref_data = ref_data) %>% 
#       insert_gaps() %$% 
#       move_seq(head_x, head_y, r_thresh = 30, inherit.theta = FALSE) %>% 
#       add_random_steps(n = 100L,
#                        sl_distr = fit_gamma(.$r),
#                        ta_distr = fit_genvonmises(.$theta_rel)
#       ) %>% 
#       flag_invalid_steps(remove = TRUE) %>% 
#       append_estimators(na_as_zero = TRUE)
#     
#     out <- issf(
#       case ~ 
#         moved : (cos_theta_pi + cos_2theta) +
#         sl + logsl + 
#         strata(step_id),
#       data = d, 
#       scale_estimator = "sl", 
#       shape_estimator = "logsl", 
#       kappa1_estimator = "moved:cos_theta_pi",
#       kappa2_estimator = "moved:cos_2theta"
#     )
#     
#     return(out)
#   }, 
#   ref_data = ref_data,
#   cores = 8, 
#   inorder = TRUE
# )
# names(fit_list2) <- ids2


#saveRDS(object = c(fit_list,fit_list2), "cleaned_data/issf_fit_list.rds")











z <- fit_list %>%
  lapply(function(x){
    c(coef(x$model), 
      "shape" = x$sl_updated$params$shape, 
      "scale" = x$sl_updated$params$scale, 
      "kappa1" = x$ta_updated$params$kappa1,
      "kappa2" = x$ta_updated$params$kappa2)
  }) %>% 
  do.call("rbind", .) %>% 
  as.data.frame() %>% 
  cbind("rep_id" = ids) %>% 
  filter(scale > 0)


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
  mutate(
    y = toxic
  ) %>% 
 # filter(abs(y) < 5) %>% 
  ggplot(aes(x = var_trt, y = y, group = beta, color = beta)) + 
  geom_pointrange(stat = "summary", color = "black", position = position_dodge(0.7)) + 
  geom_point(position = position_jitterdodge())

z %>% 
  filter(prop_kp > 0.7) %>% 
  mutate(
    y = toxic
  ) %>% 
  filter(abs(y) < 5) %>% 
  ggplot(aes(x = log(cat_pre_wt), y = y, color = var_trt)) + 
  geom_point() +
  geom_smooth(method = "lm")


z %>% 
  filter(prop_kp > 0.7) %>% 
  ggplot(aes(x = var_trt, y = as.numeric(estimate)*12*12, group = beta, color = beta)) + 
  geom_pointrange(stat = "summary", color = "black", position = position_dodge(0.7)) + 
  geom_point(position = position_jitterdodge(jitter.height = 0)) 
  scale_y_continuous(trans = "log10")

names(z)
z

z %>% 
  filter(prop_kp > 0.7) %>% 
  #filter(scale > 0 & abs(scale) < 1000) %>% 
  ggplot(aes(x = log(cat_pre_wt), y = as.numeric(estimate), color = beta)) + 
  geom_point() +
  geom_smooth(method = "lm")


z %>% 
  filter(prop_kp > 0.7) %>% 
  filter(scale > 0 & abs(scale) < 1000) %>% 
  ggplot(aes(x = var_trt, y = shape, group = beta, color = beta)) + 
  geom_pointrange(stat = "summary", color = "black", position = position_dodge(0.7)) + 
  geom_point(position = position_jitterdodge())







s <- z %>% 
  filter(prop_kp > 0.7) %>% 
  group_by(var_trt, beta) %>% 
  summarise(
    kappa1 = mean(kappa1),
    kappa2 = mean(kappa2),
    scale = mean(scale),
    shape = mean(shape)
  ) %>% 
  as.data.frame()


out <- lapply(1:nrow(s), function(i){
  make_genvonmises(s[i, "kappa1"],s[i, "kappa2"]) %>% 
    ddist(len = 21) %>% 
    cbind(rep_data.frame(s[i,], 21))
}) %>% 
  do.call("rbind", .)

out <- lapply(1:nrow(s), function(i){
  make_gamma(s[i, "shape"],s[i, "scale"]) %>% 
    ddist(len = 1000, x_max = 1000) %>% 
    cbind(rep_data.frame(s[i,], 100))
}) %>% 
  do.call("rbind", .)


ddist(make_gamma(1, 100))


out %>% 
  ggplot(aes(x = x, y = density, fill = beta)) + 
  geom_col(position = "dodge") + 
  coord_polar(theta = "x", start = pi) + 
  theme_bw(base_size = 15) + 
  labs(x = "Turn angle (radians)", y = "Density") + 
  facet_wrap(~var_trt)


out %>% 
  filter(x >0) %>% 
  ggplot(aes(x = x / 1000 * 12, y = density, color = beta)) + 
  geom_line( size = 2, aes(linetype = var_trt)) + 
  theme_bw(base_size = 15) + 
  labs(x = "Step length (cm)", y = "Density")
























ID <- 5
d <- fetch_events(ID) %>%
  filter(
    time <= filter(ref_data, rep_id == ID)$camera_cutoff
  ) %>% 
  insert_gaps() %$% 
  move_seq(head_x, head_y, r_thresh = 10, inherit.theta = FALSE) %>% 
  add_random_steps(n = 100L,
                   sl_distr = make_gamma(0.47, 944),#fit_gamma(.$r),
                   ta_distr = fit_genvonmises(.$theta_rel)
  ) %>% 
  flag_invalid_steps(remove = TRUE) %>% 
  mutate(
    toxic = read_value(x2, y2, c(1000, 1000),
                       ref_img = fetch_trt_spec(ID, .ref_data = ref_data))
  ) %>%
  append_estimators()

o <- issf(
  case ~ 
    toxic +  
    moved : (cos_theta_pi + cos_2theta) + 
    (sl + logsl) + 
    #x2 + y2 +
    #logsl : (cos_theta_pi + cos_2theta) + 
    #logsl : (cos_theta_pi:moved) +
    strata(step_id),
  data = d, 
  scale_estimator = "sl", 
  shape_estimator = "logsl", 
  kappa1_estimator = "moved:cos_theta_pi",
  kappa2_estimator = "moved:cos_2theta"
); summary(o)

o$ta_updated
o$sl_updated

z <- creat_start(500, 500, runif(1, 0, 2 * pi)) %>% 
  iterate_random_steps(start = ., 
                       issf_fit = o, n = 260, 
                       use_observed = TRUE, 
                       paired = TRUE)

z %>% 
  plot_track(x2, y2)

z %>% filter(scale < 0) %>% .$rep_id

fetch_events(147) %$%
  move_seq(head_x, head_y, r_thresh = 10) %>% 
  plot_track(x1, y1)


make_genvonmises(2, 2) %>% ddist(len = 1000) %>% plot()

o$sl_updated %>% ddist(len=1000, x_max = 500) %>% plot()









a2 <- fetch_repID() %>% 
  .[!.%in% ids & !grepl("trial|week", .)] %>% 
  pb_par_lapply(
  function(i){
    fetch_events(i) %$% ud_area(head_x, head_y, 0.95) /1000^2
  }, cores = 8
) %>% do.call("rbind", .)


a <- pb_par_lapply(
  ids, function(i){
    fetch_events(i) %$% ud_area(head_x, head_y, 0.95) /1000^2
  }, cores = 8
) %>% do.call("rbind", .)

z <- a %>% 
  as.data.frame() %>% 
  rename_all(.funs = function(x) {
    paste0(x,"50")
  }) %>% 
  cbind("rep_id" = ids) %>% 
  as.data.frame() %>% 
  right_join(z)

z

fetch_events(89) %>% 
  plot_track(head_x, head_y)

fetch_events(92) %$% 
  ud_area(head_x, head_y) / 1000^2 



names(z)
z <- z %>% filter(prop_kp > 0.5) %>% 
  mutate(
    estimate = as.numeric(estimate)
  )





library(glmmTMB)


z %>% 
  filter(prop_kp > 0.7) %>% 
  ggplot(aes(x = var_trt, y = as.numeric(estimate50)*12*12, group = beta, color = beta)) + 
  geom_pointrange(stat = "summary", color = "black", position = position_dodge(0.7)) + 
  geom_point(position = position_jitterdodge(jitter.height = 0)) 

m <- glmmTMB(log(estimate50) ~ 
          log(cat_pre_wt) * (var_trt) + beta + (1|session_id), 
        data = z); summary(m)

sjPlot::plot_model(m, type = "eff", terms = c("cat_pre_wt[all]", "var_trt"))
sjPlot::plot_model(m, type = "eff", terms = c("beta"))




ddist(make_genvonmises(0.1, 1)) %>% plot()






a2 <- a2 %>% 
  as.data.frame() %>% 
  cbind("rep_id" = fetch_repID() %>% 
          .[!.%in% ids & !grepl("trial|week", .)]) %>% 
  left_join(ref_data)

parse_trt <- function(x){
  as.numeric(gsub("mg/g", "", x))
}

a2 %>% 
  ggplot(aes(x = parse_trt(mean_trt), y = estimate)) + 
  geom_point() + 
  geom_smooth(method = "lm")


glmmTMB(
  log(estimate) ~ parse_trt(mean_trt) + log(cat_pre_wt) + (1|session_id), 
  data = a2
) %>% 
  summary()


