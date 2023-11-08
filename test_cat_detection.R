library(tidyverse)
library(herbivar)
library(foreach)
library(doSNOW)
source("helper_functions/file_mgmt.R")
source("helper_functions/image_utils.R")
herbivar::pre_cmp_fun()

src_root <- "processed_feed"
trialID <- "reptrial_5" # Set id
src_dir <- paste(src_root, trialID, sep = "/")
files <- list.files(src_dir, pattern = ".jpg") %>% files_reorder()
files_full_name <- paste(src_dir, files, sep = "/") %>% files_reorder()


vid <- lapply(seq_along(files_full_name),function(i){
  cat(sprintf("\rloading image %d", i))
  thin(fast_load_image(files_full_name[i], transform = FALSE), 5)
}) %>% imappend(axis = "z")


# vid <- load.video("reptrial_5.mp4", frames = 100, fps = 10, maxSize = 3)
# vid <- thin(vid, 5)

z <- 1:26





vid %>% 
  play(loop = TRUE)






vid_cat_mask <- vid %>% detect_cat(cores = 8)
vid %>% 
  play(loop = TRUE)

vid_cat_mask %>% 
  frames(700:741) %>% 
  plot.imlist()

vid %>% 
  frame(714) %>% 
  plot()

vid %>% 
  frame(711) %>% 
  detect_cat(w = c(0, 1)) %>% 
  plot()

vid_cat_mask %>% 
  play(loop = TRUE)

vid %>% 
  frame(711) %>% 
  color_index("all", plot = TRUE)

vid %>% 
  frame(711) %>% 
  color_index(index = c("BI","NG"), plot = FALSE) %>% 
  iml_prod(c(1, 1)) %>% 
  renorm(min = 0, max = 1) %>% 
  imagerExtra::SPE(0.1, s = 0.001, range = c(0,1)) %>% 
  threshold2(thr = 0.7, thr.exact = TRUE)



vid_cat_mask %>% 
  imsplit(axis = "z") %>% 
  lapply(function(x){
    sum(x, na.rm = TRUE)
  }) %>% 
  do.call("c", .) %>% 
  plot()

detect_jostle


vid %>% 
  frame(1:10) %>% 
  color_index(index = c("BI","NG"), plot = FALSE) %>% 
  iml_prod(c(5, 1)) %>% 
  get_gradient(axes = "z", scheme = -1) %>% 
  .[[1]] %>% 
  na_replace(0) %>% 
  imsplit("z") %>% 
  .[2:10] %>% 
  pb_par_lapply(function(x){
    imagerExtra::SPE(x, 0.1, s = 0.001, range = c(0,1))
  }) %>% 
  plot.imlist()

