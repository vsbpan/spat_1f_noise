source("helper_functions/init.R")
source("helper_functions/init_analysis.R")

IDs <- fetch_repID()

# event_list <- IDs %>% 
#   pb_par_lapply(
#     function(x){
#       out <- fetch_events(x)
#       l <- fetch_data_dict(x)
#       trt_spec <- fetch_trt_spec(x, .ref_data = get("ref_data", parent.frame()))
#       out$head_high_trt <- read_value(
#         x = out$head_x, 
#         y = out$head_y, 
#         dim_xy = get_dim(l), 
#         ref_img = trt_spec)
#       out$head_high_trt <- read_value(
#         x = out$centroid_x, 
#         y = out$centroid_y, 
#         dim_xy = get_dim(l), 
#         ref_img = trt_spec)
#       out <- spat1f::insert_gaps(out)
#       return(out)
#     }, cores = 8, 
#     export_fun_only = FALSE
#   ) %>% 
#   append_name(paste0("rep",IDs))
# saveRDS(event_list, "cleaned_data/events_list.rds")
event_list <- readRDS("cleaned_data/events_list.rds") %>% 
  lapply(function(x) left_join(x, ref_data, by = "repID"))

move_list <- event_list %>% 
  lapply(function(d){
    cbind(
      spat1f::move_seq(x = d$head_x, y = d$head_y), # t to t + 1 r and theta
      d[-nrow(d),], # Use the meta data at t
      append_name(d[-1, "head_high_trt"], "head_high_trt_t1") # head_high_trt at t+1
      ) 
  })


move_list[[1]]

move_list %>% 
  do.call("rbind", .) %>% 
  filter(!is.na(theta_rel)) %>% 
  mutate(
    theta_rel = nearest_bin(theta_rel, seq(-2 * pi, 2 * pi, by = 1))
  ) %>% 
  group_by(beta, theta_rel) %>% 
  summarise(
    count = n()
  ) %>% 
  group_by(beta) %>% 
  mutate(
    tot = sum(count),
    p = count / tot
  ) %>% 
  ggplot(aes(x = theta_rel, y = p, fill = factor(beta))) +
  geom_col(position = "dodge")

d_move %>% 
  mutate(lr = log(r)) %>% 
  filter(!is.na(lr)) %>% 
  mutate(
    lr = nearest_bin(lr, seq(0, 8, by = 1))
  ) %>% 
  group_by(beta, lr) %>% 
  summarise(
    count = n()
  ) %>% 
  group_by(beta) %>% 
  mutate(
    tot = sum(count),
    p = count / tot
  ) %>% 
  ggplot(aes(x = lr, y = p, fill = factor(beta))) +
  geom_col(position = "dodge")


d_move2 <- d_move %>% 
  mutate(
    tz = ifelse(step_id < 24 * 6, "24h", "After"),
    r = r 
  )  %>% 
  group_by(
    rep_id, tz
  ) %>% 
  summarise(
    mean_lr = mean(log(r), na.rm = TRUE),
    var_lr = var(log(r), na.rm = TRUE),
    q10_lr = quantile(na.omit(log(r)), prob = 0.10),
    q90_lr = quantile(na.omit(log(r)), prob = 0.90),
    mean_theta_abs = mean(abs(theta_rel), na.rm = TRUE), 
    var_theta_abs  = var(abs(theta_rel), na.rm = TRUE)
  ) %>% 
  left_join(
    ref_dat %>% mutate(rep_id = as.character(rep_id)), by = "rep_id"
  )


d_move2 %>% 
  ggplot(aes(x = factor(var_trt), y = mean_lr)) + 
  geom_point(position = position_jitter(width = 0.2)) + 
  geom_pointrange(stat = "summary")


d_move2 %>% 
  ggplot(aes(x = factor(beta), y = q90_lr, color = var_trt)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6)) + 
  geom_pointrange(stat = "summary", position = position_dodge(0.6)) + 
  facet_wrap(~tz)

d_move2 %>% 
  ggplot(aes(x = factor(beta), y = mean_theta_abs, color = var_trt)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6)) + 
  geom_pointrange(stat = "summary", position = position_dodge(0.6))+ 
  facet_wrap(~tz)


d_move2 %>% 
  ggplot(aes(x = log(cat_pre_wt), y = mean_lr, color = factor(var_trt))) + 
  geom_point() + 
  geom_smooth(method = "lm")+ 
  facet_wrap(~tz)

a$r %>% log() %>% hist(nclass = 100)
b$r %>% log() %>% hist(nclass = 100)
a$theta_rel %>% hist(nclass= 100)
b$theta_rel %>% hist(nclass= 100)






d_event2 <- d_event %>% 
  mutate(
    tz = ifelse(rank < 24 * 6, "24h", "After")
  ) %>% 
  group_by(rep_id, tz) %>% 
  summarise(
    tot_sd = sqrt(var(head_x, na.rm = TRUE) + var(head_y, na.rm = TRUE))
  )%>% 
  left_join(
    ref_dat %>% mutate(rep_id = as.character(rep_id)), by = "rep_id"
  )


d_event2 %>% 
  ggplot(aes(x = factor(beta), y = tot_sd, color = var_trt)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6)) + 
  geom_pointrange(stat = "summary", position = position_dodge(0.6)) + 
  facet_wrap(~tz)


d_event2 %>% 
  ggplot(aes(x = log(cat_pre_wt), y = tot_sd, color = factor(beta))) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(
    aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~"))
  )+ 
  facet_wrap(~tz)
