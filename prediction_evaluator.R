source("spat1f/init.R")



d <- read_csv("raw_data/sample1_modelv7_validation_inference.csv") %>% 
  mutate(
    polygon = ifelse(type == "prediction" & score < 0.9, NA, polygon),
    keypoints = ifelse(type == "prediction" & score < 0.9, NA, keypoints)
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




