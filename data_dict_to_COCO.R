source("helper_functions/init.R")

obj <- jsonlite::fromJSON("annotations/coco_test.json")

jsonlite::write_json(a, "abc_test2.json", pretty = TRUE)


(a$images$path)


