library(tidyverse)
library(herbivar)
library(foreach)
library(doSNOW)
source("spat1f/file_mgmt.R")
source("spat1f/image_utils.R")
source("spat1f/utils.R")
herbivar::pre_cmp_fun()

src_root <- "processed_feed"
trialID <- "rep21" # Set id
src_dir <- paste(src_root, trialID, sep = "/")
files <- list.files(src_dir, pattern = ".jpg") %>% files_reorder()
files_full_name <- paste(src_dir, files, sep = "/") %>% files_reorder()

vid <- load_video(files_full_name, cores = 6, thin.val = 5)


vid_cat_mask %>% 
  play(loop = TRUE)

undebug(detect_cat)

vid %>% 
  frame(900:910) %>% 
  detect_cat(clean = 4, lambda = 0.4, sat = 0, adjust = 1, thr = "otsu", cores = 1) %>% 
  imsplit("z") %>% 
  plot.imlist()

vid %>% 
  frame(900:910) %>% 
  imsplit("z") %>% 
  plot.imlist()

vid_cat_mask <- vid %>% detect_cat(lambda = 0.4, sat = 0, 
                                                      adjust = 1, thr = "otsu", 
                                                      cores = 6, clean = 4)

vid %>% 
  play(loop = TRUE)

vid_cat_mask %>% 
  frames(700:741) %>% 
  plot.imlist()


vid_cat_mask %>% 
  imsplit("z") %>% 
  lapply(function(x){
    sum(x)
  }) %>% do.call("c",.) %>% 
  plot()




