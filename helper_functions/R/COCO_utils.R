
# Set up generic
# Print method for COCO_Json
print.COCO_Json <- function(x){
  cat(sprintf("COCO annotation with %s images\n", nrow(x$images)))
  cat(sprintf("things: %s", x$categories["name"]))
}



# Convert data dict object to Json object
as.Json.data_dict <- function(x){
  
  format_set <- function(x){
    if(unique_len(x) == 1){
      n <- length(x)
      x <- unique(x)
    } else {
      n <- 1
    }
    
    trial_sets <- c("trial_1", "trial_2", "trial_3", "trial_4", "trial_5", "trial_6", 
                    "week2test", "week2test2", "week2test3", "week2test4", "week2test5")
    
    out <- vapply(x, function(x){
      ifelse(
        numbers_only(x), 
        sprintf("%02d%04d", 0, as.numeric(x)),
        ifelse(
          x %in% trial_sets, 
          sprintf("%02d%04d", 1, which(x %in% trial_sets)),
          sprintf("%02d%04d", 2, floor(runif(1, 0, 10000)))
        )
      )
    }, 
    character(1), 
    USE.NAMES = FALSE)
    
    return(
      rep(out, n)
    )
  }
  
  
  fmeta <- get_file_meta(x)
  
  img_id <- sprintf("1%s%05d", 
                    format_set(repID_clean(fmeta$repID)), 
                    as.numeric(fmeta$rank)
  )
  thing_id <- switch(unlist(unname(attr(x, "summary")["thing_class"])), 
                     "cat" = "1")
  
  
  kp <- lapply(map(x, "keypoints"), function(z){
    if(is.null(z)){
      return(NULL)
    } else {
      z <- rbind(t(z[,c(1,2)]), c(2,2,2))
      round(as.vector(z))
    }
  }) %>% 
    unname()
  
  seg <- map(x, "polygon") %>% 
    lapply(function(z){
      if(is.null(z)){
        return(NULL)
      } else {
        return(as.vector(t(z)))
      }
    })
  
  annotations <- data.frame(
    "category_id" = 1,
    "id" = as.numeric(sprintf("%s%s", img_id, thing_id)),
    "image_id" = as.numeric(img_id),
    "iscrowd" = FALSE,
    "num_keypoints" = ifelse(lapply(kp, is.null) %>% 
                               do.call("c",.), 
                             0,
                             3)
  )
  
  annotations$bbox <- unname(get_bbox(x)) 
  annotations$keypoints <- kp
  annotations$segmentation <- unname(seg)
  
  
  
  images <- data.frame(
    "file_name" = fmeta$file_base,
    "path" = fmeta$file_path, 
    "height" = get_dim(x)[2], 
    "width" = get_dim(x)[1],
    "id" = as.numeric(img_id)
  )
  
  categories <- data.frame(
    "color" = "#8186d5",
    "id" = thing_id,
    "name" = c("cat"),
    "supercategory" = ""
  )
  
  categories$keypoint_colors <- list(
    c("#ff0000", "#000", "#2bff00")
  )
  
  categories$keypoints <- list(
    c("head", "middle","tail")
  )
  
  categories$skeleton <- list(
    matrix(c(1,2,2,3), ncol = 2)
  )
  
  
  
  out <- list(
    "annotations" = annotations,
    "categories" = categories, 
    "images" = images,
    "info" = list(),
    "licenses" = list()
  )
  
  class(out) <- c("COCO_Json", "list")
  
  return(out)
}

as.Json <- function(x){
  UseMethod("as.Json")
}

# S3 method registration
registerS3method(genname = "as.Json", 
                 class = "data_dict", 
                 method = as.Json.data_dict)

# S3 method for print COCO_Json
registerS3method(genname = "print", 
                 class = "COCO_Json", 
                 method = print.COCO_Json)



import_COCO <- function(x){
  jsonlite::fromJSON(x)
}

export_COCO <- function(x, path){
  jsonlite::write_json(x, path, pretty = TRUE)
}

sample_COCO <- function(x, size){
  n <- nrow(x$images)
  s <- sample(seq_len(n), size, replace = FALSE)
  
  x$images <- x$images[s,]
  x$annotations <- x$annotations[s,]
  
  return(x)
}

split_COCO <- function(x, test, val = 0){
  n <- nrow(x$images)
  s_val <- sample(seq_len(n), round(n * val), replace = FALSE)
  s_test <- sample(seq_len(n)[-s_val], round(n * test), replace = FALSE)
  
  val <- x
  val$images <- x$images[s_val,]
  val$annotations <- x$annotations[s_val,]
  
  test <- x
  test$images <- x$images[s_test,]
  test$annotations <- x$annotations[s_test,]
  
  train <- x
  train$images <- x$images[-c(s_test, s_val),]
  train$annotations <- x$annotations[-c(s_test, s_val),]
  
  return(list(
    "train" = train, 
    "test" = test, 
    "val" = val
  ))
}

merge_COCO <- function(...){
  dots <- as.list(match.call())[-1]
  
  dots <- lapply(dots, function(x){
    eval(x)
  })
  
  categories <- lapply(dots, function(x){
    x$categories
  }) %>% 
    do.call("rbind",.) %>% 
    unique()
  
  annotations <- lapply(dots, function(x){
    x$annotations
  }) %>% 
    do.call("rbind",.) %>% 
    unique()
  
  images <- lapply(dots, function(x){
    x$images
  }) %>% 
    do.call("rbind",.) %>% 
    unique()
  
  
  out <- list(
    "annotations" = annotations,
    "categories" = categories, 
    "images" = images,
    "info" = list(),
    "licenses" = list()
  )
  
  class(out) <- c("COCO_Json", "list")
  
  return(out)
}




