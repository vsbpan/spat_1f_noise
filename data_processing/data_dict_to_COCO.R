source("spat1f/init.R")



##a <- as.Json(fetch_data_dict(50))
#a <- a %>% set_new_path(path_root = "/datasets/rep50")

coco_list <- lapply(
  fetch_repID(), 
  function(x){
    subset_COCO(
      as.Json(fetch_data_dict(x)), file_name = basename(s)
    )
  }
)

coco_list <- do.call("merge_COCO",coco_list)

coco_list <- coco_list %>% set_new_path("/datasets/sample1/")

export_COCO(coco_list, "annotations/modelv5_inference.json")


s <- lapply(
  list.dirs("processed_feed", 
            recursive = FALSE, 
            full.names = TRUE),
  function(x){
    z <- list.files(x, full.names = TRUE, recursive = FALSE)
    n <- length(z)
    i <- round(seq(1, n, l = 20))
    z[i]
  }
) %>% 
  do.call("c", .)


s_new <- paste0("C:/R_projects/deep_learning_playground/coco-annotator/datasets/sample1/",
               basename(s))

#file.copy(s, s_new)


s <- list.files(
  "C:/R_projects/deep_learning_playground/coco-annotator/datasets/sample1", 
  recursive = FALSE, full.names = FALSE
)
s









