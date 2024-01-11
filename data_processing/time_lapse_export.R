source("helper_functions/init.R")


source_dir <- "D:/time_lapse_export"
dest_dir_root <- "C:/R_Projects/spat_1f_noise/time_lapse_feed/"

# Move and rename files from flash drive. Need to enter repID and cameraID via prompt
export_time_lapse_img(source_dir, dest_dir_root)


# file.rename(
#   list.files("time_lapse_feed/rep138", full.names = T), 
#   shift_time_label("time_lapse_feed/rep138", 116779)
# )
report_missing_photos("time_lapse_feed")