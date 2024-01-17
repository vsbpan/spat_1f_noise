source("helper_functions/init.R")



d <- read_csv("sample1_top_model_inference.csv")
pred <- spat1f::parse_inference(d %>% 
                              filter(type == "prediction"), 
                            mode = "evaluate")

gt <- spat1f::parse_inference(d %>% 
                                  filter(type == "ground_truth"), 
                                mode = "evaluate")

z <- keypoint_evaluator(pred, gt, k = c(100), keypoints = c("head"))

i <- 100
plot_mask(gt, i)
plot_keypoint(gt, i, col = "red")
plot_keypoint(pred, i, col = "green")

mask_evaluator(pred, gt, size_range = c(1000, 100000))

plot_mask(gt, i)
plot_polygon(gt, i, col = "red")
plot_polygon(pred, i)







keypoint_evaluator(pred, gt, k = c(100), keypoints = c("head"))
mask_evaluator(pred, gt)
