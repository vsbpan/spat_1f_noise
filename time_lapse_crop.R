library(tidyverse)
library(herbivar)
library(foreach)
library(doSNOW)
source("helper_functions/file_mgmt.R")
source("helper_functions/image_utils.R")
source("helper_functions/utils.R")
herbivar::pre_cmp_fun()

src_root <- "time_lapse_feed"
dest_root <- "processed_feed"
trialID <- "reptrial_2" # Set id
src_dir <- paste(src_root, trialID, sep = "/")
dest_dir <- paste(dest_root, trialID, sep = "/")
files <- list.files(src_dir, pattern = ".jpg")
files_full_name <- paste(src_dir, files, sep = "/")

# Initialize directory
init_dir(root_path = dest_root, dir_name = trialID)

# Read in test image and etect corners
img <- fast_load_image(files_full_name[1]) 
undebug(detect_corners)
pts <- detect_corners(img, qc_plot = TRUE, scale = 5, adjust = 1, list_of_list = TRUE)

plot(img)
list(
  list(184, 183), 
  list(285, 1140),
  list(1140, 15),
  list(1182, 1075)
) %>% pt_list2df() %$% 
  points(x,y, col = c("green", "blue","blue","blue"), pch = 19, cex = 2)

pts <- list( # trial 2
  list(184, 183), 
  list(285, 1140),
  list(1140, 125),
  list(1182, 1075)
)[c(2,1,3,4)]





# Check reprojection quality
img2 <- reproject_grid(fast_load_image(files_full_name[1]), 
                       init_pts = pts, 
                       dest_size = 1000, 
                       qc_plot = FALSE)
plot(img2)




# If all is well, reproject for all images
crop_raw_img()

dest_root <- "processed_feed"
trialID <- "reptrial_2"
dest_dir <- paste(dest_root, trialID, sep = "/")

make_video(src_dir = dest_dir, file = trialID) # Make video
















































