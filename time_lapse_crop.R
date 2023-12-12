library(tidyverse)
library(herbivar)
library(foreach)
library(doSNOW)
rm(list = ls())
source("helper_functions/file_mgmt.R")
source("helper_functions/image_utils.R")
source("helper_functions/utils.R")
source("helper_functions/anchor_picker.R")
herbivar::pre_cmp_fun()

src_root <- "time_lapse_feed"
dest_root <- "processed_feed"

trialID <- list.dirs(src_root, full.names = FALSE, recursive = FALSE) %>% 
  .[!. %in% gsub(".rds","",list.files("raw_data/picked_anchors/"))] %>% 
  sample(., 1)
#trialID <- "rep16" # Set id
src_dir <- paste(src_root, trialID, sep = "/")
dest_dir <- paste(dest_root, trialID, sep = "/")
files <- list.files(src_dir, pattern = ".jpg") %>% files_reorder()
files_full_name <- paste(src_dir, files, sep = "/") %>% files_reorder()

# Initialize directory
init_dir(root_path = dest_root, dir_name = trialID)

# Read in test image and detect corners
img <- fast_load_image(files_full_name[1], transform = FALSE) 

pts <- anchor_picker_app(files_full_name[1], thin = 1, anchor_size = 5)
attr(pts, "indices") <- seq_along(files_full_name)
pts_list <- list(pts)
# 
# pts2 <- anchor_picker_app(files_full_name[max(attr(pts, "indices"))+1],
#                           thin = 1, anchor_size = 5)
# attr(pts2, "indices") <- (max(attr(pts, "indices"))+1):length(files_full_name)
# pts_list <- c(pts_list,list(pts2))
# # 
# 
# pts3 <- anchor_picker_app(files_full_name[max(attr(pts2, "indices"))+1],
#                           thin = 1, anchor_size = 5)
# attr(pts3, "indices") <- (max(attr(pts2, "indices"))+1):267 #length(files_full_name)
# pts_list <- c(pts_list,list(pts3))
# 
# pts4 <- anchor_picker_app(files_full_name[max(attr(pts3, "indices"))+1],
#                           thin = 1, anchor_size = 5)
# attr(pts4, "indices") <- (max(attr(pts3, "indices"))+1):length(files_full_name)
# pts_list <- c(pts_list,list(pts4))
# 
# pts5 <- anchor_picker_app(files_full_name[max(attr(pts4, "indices"))+1], 
#                           thin = 1, anchor_size = 5)
# attr(pts5, "indices") <- (max(attr(pts4, "indices"))+1):length(files_full_name)
# pts_list <- c(pts_list,list(pts5))

which_time(files_full_name, 68121)

plot(fast_load_image(files_full_name[500], transform = FALSE))
pts %>% 
  pt_list2df() %$% 
  points(x,y, col = c("green", "blue","blue","blue"), pch = 19, cex = 2)

pts2 %>% 
  pt_list2df() %$% 
  points(x,y, col = c("green", "blue","blue","blue"), pch = 19, cex = 2)






# Check reprojection quality
# img2 <- reproject_grid(fast_load_image(files_full_name[1], transform = FALSE), 
#                        init_pts = pts, 
#                        dest_size = 1000, 
#                        qc_plot = FALSE)
# plot(img2)

# If all is well, reproject for all images
# NOTE: The cropped image is mirrored by the x axis (looks correct with plot.cimg though bc the y-axis is flipped)
crop_raw_img(); save_anchors(trialID)






# 
# 
# 
# dest_root <- "processed_feed"
# list.dirs(dest_root, full.names = FALSE, recursive = FALSE)
# trialID <- "rep29"
# dest_dir <- paste(dest_root, trialID, sep = "/")
# 
# make_video(src_dir = dest_dir, file = trialID) # Make video
# 















































