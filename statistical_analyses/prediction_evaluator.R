source("spat1f/init_analysis.R")

##### testing data #####
# Load mask-RCNN inference
d <- read_csv("raw_data/master_modelv7_1_7999iter_test_inference.csv") %>% 
  mutate(
    polygon = ifelse(type == "prediction" & score < 0.7, NA, polygon),
    keypoints = ifelse(type == "prediction" & score < 0.7, NA, keypoints)
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
    polygon = ifelse(type == "prediction" & score < 0.7, NA, polygon),
    keypoints = ifelse(type == "prediction" & score < 0.7, NA, keypoints)
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
mask_evaluator(pred, gt, cores = 8)


# Overall detection report from the Mask-R-CNN with score thresh of 0.3 and no error flagged
detection_report(fetch_repID())



ids <- unname(unlist(id_list))[!unname(unlist(id_list)) %in% problem_ids]


# Percent of frames flagged as false head movement
# Percent of frames flagged as implausible mask size
# Percent of frames flagged as false cluster

out <- lapply(
  ids, function(x){
    z <- fetch_events(x, append_detection_summary = FALSE) %>% 
      clean_events(insert_gaps = FALSE, score_thresh = 0.7, keep_sus = TRUE)
    data.frame(
      "false_head" = mean(z$is_sus),
      "mask_error" = mean(z$status %in% c("too small", "too big")),
      "false_cluster" = mean(z$false_cluster)
    )
    
  }
) %>% 
  do.call("rbind", .)

apply(out, 2, summarise_vec, simplify = FALSE)



# Percent positive detection
# Total positive detection

out <- lapply(
  ids, function(x){
    z <- fetch_events(x, append_detection_summary = FALSE) %>% 
      clean_events(insert_gaps = TRUE, score_thresh = 0.7, keep_sus = FALSE)
    data.frame(
      "prop" = mean(!z$is_gap), 
      "sum" = sum(!z$is_gap)
    )
  }
) %>% 
  do.call("rbind", .)

apply(out, 2, summarise_vec, simplify = FALSE)
