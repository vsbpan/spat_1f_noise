source("helper_functions/init_analysis.R")

anchor_picker_app("misc_tests/prototype2.jpg")






IDs <- fetch_repID()

event_list <- IDs %>% 
  lapply(
    function(x){
      out <- fetch_events(x)
      l <- fetch_data_dict(x)
      trt_spec <- fetch_trt_spec(x)
      out$head_high_trt <- read_value(
        x = out$head_x, 
        y = out$head_y, 
        dim_xy = get_dim(l), 
        ref_img = trt_spec)
      out$head_high_trt <- read_value(
        x = out$centroid_x, 
        y = out$centroid_y, 
        dim_xy = get_dim(l), 
        ref_img = trt_spec)
      return(out)
    }
  ) %>% 
  append_name(IDs)






fetch_trt_spec(44)




fetch_events(36, append_detection_summary = F)




o2 <- lapply(z, function(x){
  d <- read_csv(sprintf("C:/R_Projects/spat_1f_noise/cleaned_data/events/rep%s.csv",x)) %>% suppressMessages()
  d 
})

names(o2) <- paste0("rep", z)
d_event <- lapply(seq_along(o2), function(x){
  cbind(o2[[x]], "rep_id" = gsub("rep","",names(o2)[x]))
}) %>% 
  do.call("rbind",.) %>%
  left_join(ref_dat %>% mutate(rep_id = as.character(rep_id)), by = "rep_id")



d_move<- lapply(seq_along(o), function(x){
  cbind(o[[x]], "rep_id" = gsub("rep","",names(o)[x]))
}) %>% 
  do.call("rbind",.) %>% 
  left_join(ref_dat %>% mutate(rep_id = as.character(rep_id)), by = "rep_id")



d_move %>% 
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
