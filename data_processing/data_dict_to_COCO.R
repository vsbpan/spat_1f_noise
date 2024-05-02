source("spat1f/init.R")


s <- list.files(
  "C:/R_projects/deep_learning_playground/coco-annotator/datasets/sample2", 
  recursive = FALSE, full.names = FALSE
)
s


coco_list <- lapply(
  fetch_repID(), 
  function(x){
    subset_COCO(
      as.Json(fetch_data_dict(x)), file_name = basename(s)
    )
  }
)

coco_list <- do.call("merge_COCO",coco_list)

coco_list <- coco_list %>% set_new_path("/datasets/sample2/")

export_COCO(coco_list, "annotations/sample2_modelv6_inference.json")

# 
# 
# a <- import_COCO("annotations/sample1_manual_annotation.json")
# a$images <- a$images %>% 
#   dplyr::select(all_of(names(coco_list$images))) %>% 
#   as.data.frame()
# a$categories <- a$categories %>% 
#   dplyr::select(all_of(names(coco_list$categories))) %>% 
#   as.data.frame()
# a$annotations <- a$annotations %>% 
#   dplyr::select(all_of(names(coco_list$annotations))) %>% 
#   as.data.frame()
# 
# 
# new_files <- s[!s %in% a$images$file_name]
# 
# 
# coco_final <- merge_COCO(subset_COCO(coco_list, new_files), a)
# 
# coco_final$categories <- coco_final$categories[1,]
# coco_final <- coco_final %>% set_new_path("/datasets/sample2/")
# 
# coco_final
# 
# export_COCO(coco_final, "annotations/modelv6_inference_with_manual.json")





