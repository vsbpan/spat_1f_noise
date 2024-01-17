mask_evaluator <- function(prediction, ground_truth, IOU_thresh = c(0.5, 0.75, 0.9), size_range = NULL){
  
  if(!is.null(size_range)){
    mask_size_vec <- mask_area(prediction)
    subset_i <- is_between(mask_size_vec, size_range)
    prediction <- prediction[subset_i]
    ground_truth <- ground_truth[subset_i]
  }
  
  
  z <- lapply(seq_along(prediction), function(i){
    mask_IOU(get_mask(prediction, i)[[1]],
             get_mask(ground_truth, i)[[1]], 
             na.rm = TRUE)
  }) %>% 
    do.call("c",.)
  pred_exist <- lapply(get_polygon(prediction), function(x){
    !any(is_null_na(x))
  }) %>% 
    do.call("c",.)
  gt_exist <- lapply(get_polygon(ground_truth), function(x){
    !any(is_null_na(x))
  }) %>% 
    do.call("c",.)
  
  out <- .evaluator_calc(z, IOU_thresh, pred_exist, gt_exist, thresh_name = "IOU")
  
  print(out$result)
  invisible(out)
}


.evaluator_calc <- function(score, thresh, pred_exist, gt_exist, thresh_name = NULL){
  if(is.null(thresh_name)){
    thresh_name <- "score"
  }
  
  mat <- outer(thresh, score, 
               Vectorize(function(x, y){
                 y >= x
               }))
  
  
  TP <- as.vector(mat %*% gt_exist)
  TN <- as.vector((!mat) %*% (!gt_exist))
  FP <- as.vector((!mat) %*% pred_exist)
  FN <- as.vector((!mat) %*% gt_exist)
  
  n <- length(score)
  
  eval_res <- data.frame(
    "score" = thresh, 
    "TP" = TP,
    "FP" = FP,
    "TN" = TN,
    "FN" = FN,
    "n" = n
  )
  names(eval_res)[1] <- thresh_name
  
  eval_res$precision <- eval_res$TP / (eval_res$TP + eval_res$FP)
  eval_res$recall <- eval_res$TP / (eval_res$TP + eval_res$FN)
  eval_res$accuracy <- (eval_res$TP + eval_res$TN) / (eval_res$TP + eval_res$TN + eval_res$FP + eval_res$FN)
  
  out <- list("result" = eval_res, "score" = score)
  names(out)[2] <- thresh_name
  return(out)
}



keypoint_evaluator <- function(prediction, ground_truth, 
                               k = c(100, 1000, 100), # From visually assessing a sample whether this is reliable for distinguishing false positive and true positive.
                               keypoints = c("head", "middle", "tail"), 
                               OKS_thresh = c(0.5, 0.75, 0.9), 
                               size_range = NULL){
  
  keypoints <- match(match.arg(keypoints, several.ok = TRUE), c("head", "middle", "tail"))
  stopifnot(length(k) == length(keypoints))
  
  o <- order(keypoints)
  keypoints <- keypoints[o]
  k <- k[o]
  
  mask_size_vec <- mask_area(prediction)
  if(!is.null(size_range)){
    subset_i <- is_between(mask_size_vec, size_range)
    prediction <- prediction[subset_i]
    ground_truth <- ground_truth[subset_i]
  }
  
  
  
  z <- lapply(seq_along(prediction), function(i){
    kp1 <- map(prediction[i], "keypoints")[[1]]
    kp2 <- map(ground_truth[i], "keypoints")[[1]]
    s2 <- mask_size_vec[i] # Mask size
    
    
    
    pred_exist <- !is.null(kp1)
    gt_exist <- !is.null(kp2)
    
    if(!pred_exist | !gt_exist){
      oks <- 0
    } else {
      oks <- keypoint_OKS(
        kp1[keypoints,c(1,2), drop = FALSE], 
        kp2[keypoints,c(1,2), drop = FALSE], 
        s2 = s2, 
        k = k, 
        ground_truth_flag = kp1[keypoints,"score", drop = TRUE]
      )
    }
    
    c("oks" = oks, "pred_exist" = pred_exist, "gt_exist" = gt_exist)
    
  }) %>% 
    do.call("rbind", .) %>% 
    as.data.frame()
  
  out <- .evaluator_calc(z$ok, OKS_thresh, z$pred_exist, z$gt_exist, thresh_name = "OKS")
  
  print(out$result)
  invisible(out)
}


keypoint_OKS <- function(kp1, kp2, s2, k, ground_truth_flag){
  has_lab <- as.numeric(ground_truth_flag > 0)
  
  dist <- (kp1 - kp2)^2
  if(is.null(nrow(dist))){
    dist <- sum(dist)
  } else {
    dist <- rowSums(dist)
  }
  
  n <- sum(has_lab)
  
  if(n == 0){
    oks <- 0
  } else {
    oks <- sum(exp(-(dist)^2 / (s2 * 2 * k^2)) * has_lab) / n 
  }
  return(oks)
}


mask_union_area <- function(img, img2, na.rm = FALSE){
  stopifnot(all.equal(dim(img), dim(img2)))
  stopifnot(dim(img)[3] == 1)
  sum(as.pixset(img) | as.pixset(img2), na.rm = na.rm)
}


mask_intersection_area <- function(img, img2, na.rm = FALSE){
  stopifnot(all.equal(dim(img), dim(img2)))
  stopifnot(dim(img)[3] == 1)
  sum(as.pixset(img) & as.pixset(img2), na.rm = na.rm)
}


mask_IOU <- function(img, img2, na.rm = FALSE){
  mask_intersection_area(img, img2, na.rm) / mask_union_area(img, img2, na.rm)
}



bbox_area.default <- function(x, ...){
  if(is.null(x)){
    return(0)
  }
  z <- parse_bbox_vec(x)
  prod(abs(z[1,] - z[2,]))
}

bbox_area.data_dict <- function(x, ...){
  get_bbox(x) %>% 
    lapply(bbox_area.default) %>% 
    do.call("c",.) %>% 
    unname()
}


bbox_area <- function(x, ...){
  UseMethod("bbox_area")
}

registerS3method(genname = "bbox_area", 
                 class = "default", 
                 method = bbox_area.default)

registerS3method(genname = "bbox_area", 
                 class = "data_dict", 
                 method = bbox_area.data_dict)


mask_area.default <- function(x, ...){
  if(is.null(x)){
    return(0)
  }
  stopifnot(dim(x)[3] == 1)
  sum(as.pixset(x))
}

mask_area.data_dict <- function(x, ...){
    lapply(seq_along(x), function(i){
      mask_area.default(get_mask(pred, i)[[1]])
    }) %>% 
    do.call("c",.) %>% 
    unname()
}


mask_area <- function(x, ...){
  UseMethod("mask_area")
}

registerS3method(genname = "mask_area", 
                 class = "default", 
                 method = mask_area.default)

registerS3method(genname = "mask_area", 
                 class = "data_dict", 
                 method = mask_area.data_dict)



bbox_IOU <- function(bbox1, bbox2){
  
  x_left = max(bbox1[1], bbox2[1])
  y_top = max(bbox1[2], bbox2[2])
  x_right = min(bbox1[3], bbox2[3])
  y_bottom = min(bbox1[4], bbox2[4])
  
  if (x_right < x_left | y_bottom < y_top){
    return(0)
  }
  
  area_I <- (x_right - x_left) * (y_bottom - y_top)
  area_union <- bbox_area(bbox1) + bbox_area(bbox2) - area_I
  
  return(area_I / area_union)
}






