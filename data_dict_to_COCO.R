source("helper_functions/init.R")

obj <- jsonlite::fromJSON("annotations/rep50-3.json")

jsonlite::write_json(obj, "annotations/rep50-4.json", pretty = TRUE)


export_COCO(as.Json(fetch_data_dict(50)) %>% 
              sample_COCO(1), "annotations/rep50-test.json")

(a$images$path)



obj$images$events[[1]] %>% length()
