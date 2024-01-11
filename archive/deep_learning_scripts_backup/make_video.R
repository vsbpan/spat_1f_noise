library(tidyverse)
library(herbivar)
library(foreach)
library(doSNOW)
source("C:/R_projects/spat_1f_noise/helper_functions/file_mgmt.R")
source("C:/R_projects/spat_1f_noise/helper_functions/image_utils.R")
source("C:/R_projects/spat_1f_noise/helper_functions/utils.R")
herbivar::pre_cmp_fun()

dest_dir <- "detectron2/custom_data2/model5_vid_pred_15"
attach_order_file_name(dest_dir)
make_video(src_dir = dest_dir, file = "model5_vid_pred_15") # Make video



