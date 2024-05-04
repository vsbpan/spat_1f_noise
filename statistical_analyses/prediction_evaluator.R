source("spat1f/init.R")

##### testing data #####
# Load mask-RCNN inference
d <- read_csv("raw_data/master_modelv7_1_7999iter_test_inference.csv") %>% 
  mutate(
    polygon = ifelse(type == "prediction" & score < 0.9, NA, polygon),
    keypoints = ifelse(type == "prediction" & score < 0.9, NA, keypoints)
  )

# Parse inference and ground truths as data_dicts
pred <- spat1f::parse_inference(d %>% 
                              filter(type == "prediction"), 
                            mode = "evaluate")

gt <- spat1f::parse_inference(d %>% 
                                  filter(type == "ground_truth"), 
                                mode = "evaluate")



# Keypoint evaluation
keypoint_evaluator(pred, gt)

# Mask evalutation
mask_evaluator(pred, gt)


#### Training data #### 
# Load mask-RCNN inference
d <- read_csv("raw_data/master_modelv7_1_7999iter_train_inference.csv") %>% 
  mutate(
    polygon = ifelse(type == "prediction" & score < 0.9, NA, polygon),
    keypoints = ifelse(type == "prediction" & score < 0.9, NA, keypoints)
  )

# Parse inference and ground truths as data_dicts
pred <- spat1f::parse_inference(d %>% 
                                  filter(type == "prediction"), 
                                mode = "evaluate")

gt <- spat1f::parse_inference(d %>% 
                                filter(type == "ground_truth"), 
                              mode = "evaluate")



# Keypoint evaluation
keypoint_evaluator(pred, gt)

# Mask evalutation
mask_evaluator(pred, gt)



detection_report(fetch_repID())




