source("helper_functions/init.R")



d <- read_csv("raw_data/sample1_modelv6_validation_inference.csv") %>% 
  mutate(
    polygon = ifelse(type == "prediction" & score < 0.7, NA, polygon),
    keypoints = ifelse(type == "prediction" & score < 0.7, NA, keypoints)
  )
pred <- spat1f::parse_inference(d %>% 
                              filter(type == "prediction"), 
                            mode = "evaluate")

gt <- spat1f::parse_inference(d %>% 
                                  filter(type == "ground_truth"), 
                                mode = "evaluate")


keypoint_evaluator(pred, gt, k = c(100), keypoints = c("head"))
mask_evaluator(pred, gt)

detection_report(fetch_repID())






fetch_data_dict(5)



gt <- as.data_dict(import_COCO("annotations/sample1_manual_annotation.json"))
pred <- fetch_repID() %>% 
  lapply(
    function(k){
      fetch_data_dict(k)
    }
  ) %>% 
  do.call("c", .)

pred <- pred[which(get_file_meta(pred)$file_base %in% get_file_meta(gt)$file_base)]

gt <- gt[which(get_file_meta(gt)$file_base %in% get_file_meta(pred)$file_base)]

pred <- pred[match(get_file_meta(gt)$file_base, get_file_meta(pred)$file_base)]

keypoint_evaluator(pred, gt, k = c(100), keypoints = c("head"))
mask_evaluator(pred, gt)






