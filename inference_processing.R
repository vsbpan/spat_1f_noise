source("helper_functions/init.R")

inf_fn <- list.files("raw_data/inferences/", full.names = TRUE)

d <- read_csv(inf_fn[46])

l <- parse_inference(d)
d2 <- get_data(l)


d2 %>% 
  ggplot(aes(x = head_x, y = head_y)) + 
  geom_point()



trt_meta <- list.files("raw_data/trt_spectra_meta/", full.names = TRUE) %>% 
  lapply(function(x){
    invisible(read_csv(x))
  }) %>% 
  do.call("rbind",.)


trt_meta_iml <- trt_meta %>% 
  trt_meta_as_list()

trt_meta_iml$`syn_id__beta-5_id15` %>% 
  resize(size_x = get_dim(l)[1], size_y = get_dim(l)[2]) %>% 
  plot()

points(d2$head_x, d2$head_y, col = "red")

plot_image_guide(trt_meta_iml$`syn_id__beta-5_id15`)



d2$head_high_trt <- read_value(
           x = d2$head_x, 
           d2$head_y, 
           dim_xy = get_dim(l), 
           ref_img = fetch_trt_spec(get_repID(l), 
                                    ref_data = ref_dat, 
                                    trt_meta_iml = trt_meta_iml))

names(d2)

d2$file_path


d2 %>% 
  filter(size_px > 500) %>% 
  ggplot(aes(x = time, y = log(size_px))) + 
  geom_point()


d2 %>% 
  filter(size_px > 500) %>% 
  ggplot(aes(x = time, y = head_high_trt)) + 
  geom_point() + 
  geom_smooth()


d2 %>% 
  filter(size_px > 500) %>% 
  ggplot(aes(x = time, y = head_high_trt)) + 
  geom_point() + 
  geom_smooth()

d2 %>% 
  filter(!is.na(head_x))




d2 <- d2 %>% insert_gaps()


b <- move_seq(d2$head_x, d2$head_y) 




a$r %>% log() %>% hist(nclass = 100)
b$r %>% log() %>% hist(nclass = 100)
a$theta_rel %>% hist(nclass= 100)
b$theta_rel %>% hist(nclass= 100)









