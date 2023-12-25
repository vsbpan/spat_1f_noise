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
