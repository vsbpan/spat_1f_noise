source("helper_functions/init.R")



d <- read_csv("raw_data/sample1_modelv6_validation_inference.csv")
pred <- spat1f::parse_inference(d %>% 
                              filter(type == "prediction"), 
                            mode = "evaluate")

gt <- spat1f::parse_inference(d %>% 
                                  filter(type == "ground_truth"), 
                                mode = "evaluate")


keypoint_evaluator(pred, gt, k = c(100), keypoints = c("head"))
mask_evaluator(pred, gt)

detection_report(fetch_repID())
