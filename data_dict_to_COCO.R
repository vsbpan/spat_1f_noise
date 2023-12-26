source("helper_functions/init.R")

obj <- jsonlite::fromJSON("annotations/rep50-3.json")

jsonlite::write_json(obj, "annotations/rep50-4.json", pretty = TRUE)





a$images

a <- as.Json(fetch_data_dict(50))
a <- a %>% set_new_path(path_root = "/datasets/rep50")



export_COCO(a, "annotations/rep50_inference_test.json")
