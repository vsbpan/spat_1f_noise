library(tidyverse)
library(herbivar)
library(foreach)
library(doSNOW)
source("helper_functions/file_mgmt.R")
source("helper_functions/image_utils.R")
source("helper_functions/utils.R")
source("helper_functions/anchor_picker.R")
herbivar::pre_cmp_fun()




report_missing_photos("time_lapse_feed")



d <- read_table()
library(tidyverse)
library(glmmTMB)



d<-d %>% 
  mutate(
    start = as.POSIXct(paste(date_start, time_start, sep = " "), format = "%m-%d-%Y %H:%M"),
    end = as.POSIXct(paste(date_end_camera, time_end_camera, sep = " "), format = "%m-%d-%Y %H:%M")
  ) %>% 
  mutate(
    time_elapsed = as.numeric(as.difftime(end - start, units = "h"))
  ) %>% 
  mutate(
    RGR = log(cat_post_wt/cat_pre_wt) / time_elapsed
  )


d %>% 
  filter(
    session_id == 2 
  ) %>% 
  ggplot(aes(x = beta, y = RGR)) + 
  geom_point(position = position_jitter(height = 0, width = 0.2)) + 
  geom_pointrange(stat = "summary")

d %>% 
  #filter(var_trt = "constant") %>% 
  filter(mean_trt == "1 mg/g") %>% 
  ggplot(aes(x = var_trt, y = RGR)) + 
  geom_point(position = position_jitter(height = 0, width = 0.2)) + 
  geom_pointrange(stat = "summary")


d %>% 
  ggplot(aes(x = beta, y = cat_dead_cam_end)) + 
  geom_point(position = position_jitter(height = 0)) + 
  geom_pointrange(stat = "summary")
  

glmmTMB(
  RGR ~ 
    log(cat_pre_wt) + as.factor(beta) + var_trt + (1|session_id), 
  family = gaussian(), 
  data = d 
) %>% summary()


glmmTMB(
  cat_dead_cam_end ~ 
    log(cat_pre_wt) + var_trt + as.factor(beta) + (1|session_id), 
  family = binomial(), 
  data = d
) %>% summary()






