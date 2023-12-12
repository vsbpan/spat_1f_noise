library(tidyverse)




files <- list.dirs(
  "C:/R_projects/spat_1f_noise/processed_feed/", full.names = TRUE, recursive = FALSE
) %>% 
  .[!grepl("week|trial",.)] %>% 
  lapply(
  function(x){
    z <- list.files(x, full.name = TRUE)
    sample(z, 5)
  }
) %>% 
  do.call("c",.)

file_new <- gsub(".*/","C:/R_projects/deep_learning_playground/coco-annotator/datasets/caterpillar/",files)


file.copy(
  files, file_new
)


files <- list.files("C:/R_projects/deep_learning_playground/coco-annotator/datasets/caterpillar/", full.names = TRUE)


files <- files[gsub(".*rep|__cam.*","",files) %in% c(1:30)]

files_new <- gsub("caterpillar","small_caterpillar",files)

file.rename(files, files_new)



