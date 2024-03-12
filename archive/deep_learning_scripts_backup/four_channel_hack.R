library(tidyverse)
library(herbivar)
library(foreach)
library(doSNOW)
source("C:/R_projects/spat_1f_noise/spat1f/file_mgmt.R")
source("C:/R_projects/spat_1f_noise/spat1f/image_utils.R")
source("C:/R_projects/spat_1f_noise/spat1f/utils.R")
herbivar::pre_cmp_fun()

# 
# files <- list.files("C:/R_projects/deep_learning_playground/detectron2/custom_data2/datasets/caterpillar", full.names = TRUE)
# 
# candidate_files <- list.files(
#   "C:/R_projects/spat_1f_noise/processed_feed", 
#   full.names = TRUE, pattern = ".jpg", recursive = TRUE
# )
# 
# candidate_data <- file_meta(candidate_files)
# 
# file_data <- file_meta(files)

# 
# gsub("_diff","",files) %>%
#   data.frame("x" = .) %>%
#   group_by(x) %>%
#   tally() %>%
#   filter(n == 1) %>%
#   .$x -> x
# 
# file.remove(x)

files <- list.files("C:/R_projects/deep_learning_playground/test_data/rep15/", full.names = TRUE)



candidate_files <- files

candidate_data <- file_meta(candidate_files)

file_data <- file_meta(files)


tot <- nrow(file_data)
new_file_dir <- "C:/R_projects/deep_learning_playground/test_data/rep25"
for (i in seq_len(tot)){
  if(file_data$rank[i] < 11){
    next
  }
  img <- candidate_data %>% 
    filter(
      file_exists &
        repID == file_data$repID[i]
    ) %>% 
    filter(
      rank %in% c(as.numeric(file_data$rank[i]) + c(-10,0))
    ) %>% 
    .$file_path %>% 
    lapply(function(x){
      fast_load_image(x, transform = FALSE)
    }) %>% 
    as.imlist() %>% 
    imappend("z") %>% 
    get_gradient(.,axes = "z", scheme = 1) %>% 
    .[[1]] %>% 
    frame(1) %>% 
    imager::G()
  
  # img2 <- imappend(
  #   imlist = imlist(img, fast_load_image(file_data$file_path[i], transform = FALSE)), 
  #   axis = "c"
  # )
  
  new_name <- gsub(".jpg", "_diff.jpg", file_data$file_base[i])
  jpeg::writeJPEG(image = as.bmp(img), 
                  target = paste(new_file_dir, new_name, sep = "/"), quality = 1)
  cat(sprintf("Processing %s out of %s\r", i, tot))
}







