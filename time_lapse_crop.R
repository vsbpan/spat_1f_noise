library(tidyverse)
library(herbivar)
library(foreach)
library(doSNOW)
source("helper_functions/file_mgmt.R")
source("helper_functions/image_utils.R")
source("helper_functions/utils.R")
source("helper_functions/anchor_picker.R")
herbivar::pre_cmp_fun()

src_root <- "time_lapse_feed"
dest_root <- "processed_feed"
trialID <- "rep21" # Set id
src_dir <- paste(src_root, trialID, sep = "/")
dest_dir <- paste(dest_root, trialID, sep = "/")
files <- list.files(src_dir, pattern = ".jpg")
files_full_name <- paste(src_dir, files, sep = "/")

# Initialize directory
init_dir(root_path = dest_root, dir_name = trialID)

# Read in test image and etect corners
img <- fast_load_image(files_full_name[1], transform = FALSE) 

pts <- anchor_picker_app(files_full_name[1], thin = 2, anchor_size = 5)
#undebug(detect_corners)

#pts <- detect_corners(img, qc_plot = TRUE, scale = 3, adjust = 1, list_of_list = TRUE)

plot(img)
pts %>% 
  pt_list2df() %$% 
  points(x,y, col = c("green", "blue","blue","blue"), pch = 19, cex = 2)



# Check reprojection quality
img2 <- reproject_grid(fast_load_image(files_full_name[1], transform = FALSE), 
                       init_pts = pts, 
                       dest_size = 1000, 
                       qc_plot = FALSE)
plot(img2)

# If all is well, reproject for all images
# NOTE: The cropped image is mirrored by the x axis (looks correct with plot.cimg though bc the y-axis is flipped)
crop_raw_img()

dest_root <- "processed_feed"
list.dirs(dest_root, full.names = FALSE, recursive = FALSE)
trialID <- "rep3"
dest_dir <- paste(dest_root, trialID, sep = "/")

make_video(src_dir = dest_dir, file = trialID) # Make video
















































