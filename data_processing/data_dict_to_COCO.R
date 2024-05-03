source("spat1f/init.R")

#### Export mask-R-CNN predictions for a set of images in COCO Json format ####
# Read file names of sampled images
s <- list.files(
  "C:/R_projects/deep_learning_playground/coco-annotator/datasets/sample2", 
  recursive = FALSE, full.names = FALSE
)
s

# Find the model inference from stored data_dicts. 
coco_list <- lapply(
  fetch_repID(), 
  function(x){
    subset_COCO(
      as.Json(fetch_data_dict(x)), file_name = basename(s)
    )
  }
)

coco_list <- do.call("merge_COCO",coco_list)

coco_list <- coco_list %>% set_new_path("/datasets/sample2/")

# Export annotations
# export_COCO(coco_list, "annotations/sample2_modelv6_inference.json")




#### Clean and merge annotations exported from COCO Annotator ####

# Take the set of all images in a sample, wipe the annotations, and append any that is found in manual_COCO. Repeat for sample 1 and sample 2. 
# I annotated the images in three rounds. a set of ~400 images were used to train modelv5, which I used to do semi-automated/manual annotations for sample 1. Sample 1 was used to train modelv6, which I used for semi-automated/manual annotations for sample 2. The final dataset is what is used to train the final model, i.e. modelv7
a1 <- import_COCO("annotations/sample1_modelv5_inference.json")
b1 <- import_COCO("annotations/sample1_manual_annotation_updated.json")

z <- update_manual_COCO(a1, b1)
z$categories <- z$categories[1,]


a <- import_COCO("annotations/sample2_modelv6_inference.json")
b <- import_COCO("annotations/sample2_manual_annotation.json")

z2 <- update_manual_COCO(a, b)
z2$categories <- z2$categories[1,]

# Merge the two samples together
final_COCO <- merge_COCO(z, z2)

# Write the root path for where to look for the image files
final_COCO <- final_COCO %>% set_new_path("/datasets/master_sample/")

# Export master annotation
# export_COCO(final_COCO, path = "annotations/master_sample_manual_annotation.json")


#### Data splitting ####

final_COCO <- import_COCO("annotations/master_sample_manual_annotation.json")

final_COCO_split <- split_COCO(final_COCO, test = 0.1, val = 0.05)


# export_COCO(final_COCO_split$train, "annotations/coco_train.json")
# export_COCO(final_COCO_split$test, "annotations/coco_test.json")
# export_COCO(final_COCO_split$val, "annotations/coco_val.json")





