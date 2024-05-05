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


# Overall detection report from the Mask-R-CNN with score thresh of 0.3 and no error flagged
detection_report(fetch_repID())


# Percent of frames flagged as false head movement
out <- lapply(
  fetch_repID()[1:149], function(x){
    z <- fetch_events(x) %>% 
      clean_events(insert_gaps = FALSE, score_thresh = 0.9, keep_sus = TRUE)
    mean(z$is_sus)
  }
) %>% 
  do.call("c", .)

summarise_vec(out)


# Percent of frames flagged as implausible mask size
out <- lapply(
  fetch_repID()[1:149], function(x){
    z <- fetch_events(x) %>% 
      clean_events(insert_gaps = FALSE, score_thresh = 0.9, keep_sus = TRUE)
    mean(z$status %in% c("too small", "too big"))
  }
) %>% 
  do.call("c", .)

summarise_vec(out)


# Percent of frames flagged as false cluster
out <- lapply(
  fetch_repID()[1:149], function(x){
    z <- fetch_events(x) %>% 
      clean_events(insert_gaps = FALSE, score_thresh = 0.9, keep_sus = TRUE)
    mean(z$false_cluster)
  }
) %>% 
  do.call("c", .)

summarise_vec(out)




# Percent positive detection
out <- lapply(
  fetch_repID()[1:149], function(x){
    z <- fetch_events(x) %>% 
      clean_events(insert_gaps = TRUE, score_thresh = 0.9, keep_sus = FALSE)
    mean(!z$is_gap)
  }
) %>% 
  do.call("c", .)

summarise_vec(out)