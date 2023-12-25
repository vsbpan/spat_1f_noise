source("helper_functions/init.R")

obj <- jsonlite::fromJSON("annotations/rep50-3.json")

jsonlite::write_json(obj, "annotations/rep50-4.json", pretty = TRUE)


obj <- as.Json(fetch_data_dict(50)) %>% 
  sample_COCO(100)


obj$images$path <- paste0(
  "/datasets/rep50/", 
  obj$images$file_name
)

export_COCO(obj, "annotations/rep50_inference2.json")

(a$images$path)


export_COCO(import_COCO(x = "annotations/rep50_COCO_annotator.json"), "annotations/rep50_COCO_annotator.json")




obj$images$events[[1]] %>% length()



db <- import_COCO("annotations/COCO_database.json")

db <- db$images %>% filter(dataset_id == "2")

dbdf <- data.frame(
  "id" = db$id,
  "file_name" = db$file_name
)



a$images

a <- as.Json(fetch_data_dict(50), db = dbdf)
a <- a %>% set_new_path(path_root = "datasets/rep50")











