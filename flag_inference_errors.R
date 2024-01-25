source("helper_functions/init.R")

# parse_inference

d <- fetch_events(93)

z <- fetch_data_dict(93)


rank2time(z, 926)

for (i in 1:50){
  plot(z, i, mask_col = "white", add = FALSE, main = i)
  cat(sprintf("Rank: %s \r", i))
  Sys.sleep(0.1)
}


f <- forward_plot(z, 800)

f()

# 803 - 805

get_mask(z, 22:25) %>% imappend("z") %>% mask_info()

plot(z, 800, mask_col = "white")

a <- get_mask(z, 22:25) %>% imappend("z") %>% mask_info()

f()


move_seq(a$centroid_x, a$centroid_y)$r > 2

moved_centroid <- move_seq(a$centroid_x, a$centroid_y)$r > 10
moved_head <- move_seq(b$head_x, b$head_y)$r > 10


moved_head * !moved_centroid 


b <- get_keypoints(z[22:25])

f <- forward_plot(z, 800)
f()

iou <- get_mask(z, 800:820) %>% 
  lapply(seq_len(length(.)-1), function(i, l){
    mask_IOU(l[[i]], l[[i+1]], use_C = FALSE)
  }, l = .) %>% 
  do.call("c", .)

move_IOU(z, c(1:10))



get_polygon(z[c(0,1)])



moved_head * !(iou < 0.5)



.detect_false_head_movement_stage1(z,d) %>% which()


names(d)


moved_centroid <- move_seq(a$centroid_x, a$centroid_y)$r > 10
moved_head <- move_seq(b$head_x, b$head_y)$r > 10


fetch_events(7) %>% View()

get_file_meta(z)


detect_false_head_movement(data_dict = fetch_data_dict(50), cores = 8)

get_file_meta(z)


fetch_data_dict(7) %>% plot(frame = 1, add = TRUE, kp_name = "head")
fetch_image(7, 1, transform = TRUE) %>% plot()


read_value

source("helper_functions/init_analysis.R")

fetch_trt_spec(114) %>% flip_xy() %>% plot()






plot_image_guide(fetch_trt_spec(114), transform = T)


d <- fetch_events(50)
z <- fetch_data_dict(50)


d$head_high <- read_value(d$head_x, d$head_y, dim_xy = get_dim(z), 
                          ref_img = fetch_trt_spec(50), 
                          transform = TRUE)

mean(d$head_high, na.rm = TRUE)


out <- pb_par_lapply(fetch_repID()[1:149], function(i, ref_data){
  d <- fetch_events(i)
  read_value(d$head_x, d$head_y, dim_xy = c(1000, 1000), 
             ref_img = fetch_trt_spec(i, quiet = TRUE, .ref_data = ref_data), 
             transform = TRUE)
}, ref_data = ref_data, cores = 8)


ids <- fetch_repID()[1:149]

z <- data.frame(
  "ids" = ids, 
  "h1" = out %>% 
    lapply(
      function(x) mean(x, na.rm = TRUE)
    ) %>% 
    do.call("c", .),
  "h2" = out2 %>% 
    lapply(
      function(x) mean(x, na.rm = TRUE)
    ) %>% 
    do.call("c", .)
)


z <- z %>% 
  mutate(
    ids = as.numeric(ids)
  ) %>% 
  arrange(ids)



d <- ref_data %>% 
  left_join(
    z %>% 
      select(ids, h1) %>% 
      mutate(rep_id = as.character(ids)), by = "rep_id"
  )










d %>% 
  mutate(size = ifelse(cat_pre_wt>0.02, "big", "small")) %>% 
  ggplot(aes(x = var_trt, y = h1, color = beta, group = beta)) + 
  geom_pointrange(stat = "summary", 
                  position = position_dodge(width = 0.7), 
                  color = "black") + 
  geom_point(position = position_jitterdodge())+ 
  facet_wrap(~size)





d %>% 
  mutate(size = ifelse(cat_pre_wt>0.01, "big", "small")) %>% 
  ggplot(aes(x = cat_pre_wt, y = h1, color = beta)) + 
  geom_point()+ 
  scale_x_continuous(trans = "log10") +
  geom_smooth(method = "lm")





source("helper_functions/init.R")

fetch_events(99) %$%
  move_seq(head_x, head_y, r_thresh = 1) %>% 
  mutate(
    moved = ifelse(r > 20, 1, 0)
  ) %>% 
  summarise(
    r = mean(r, na.rm = TRUE),
    moved = mean(moved, na.rm = TRUE), 
    n = n(),
    n_moved = n * moved
  )













