source("helper_functions/init.R")


src_root <- "time_lapse_feed"
dest_root <- "processed_feed"

haves <- list.dirs(dest_root, full.names = FALSE, recursive = FALSE)
might_haves <- gsub(".rds","",list.files("raw_data/picked_anchors/"))

trialIDs <- might_haves[!might_haves %in% haves]
trialIDs <- trialIDs[!grepl("trial", trialIDs)]


for (trialID in trialIDs){
  src_dir <- paste(src_root, trialID, sep = "/")
  dest_dir <- paste(dest_root, trialID, sep = "/")
  files_full_name <- list.files(src_dir, pattern = ".jpg", full.names = TRUE) %>% 
    files_reorder()

  # Initialize directory
  init_dir(root_path = dest_root, dir_name = trialID)
  
  pts_list <- fetch_anchors(trialID)
  crop_raw_img(.pts_list = pts_list)
}

