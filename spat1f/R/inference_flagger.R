# Compute the mask IOU, comparing with n frames ago
move_IOU <- function(data_dict, frames, n = 1, cores = 1, text = "Processing IOU of mask"){
  f <- function(i, data_dict){
    polyl <- get_polygon(data_dict[c(i,i-1)])
    if(length(polyl) !=2){
      return(NA_real_)
    } else {
      return(polygon_IOU(polyl[[1]], polyl[[2]])) 
    }
  }
  
  if(cores > 1){
    iou <- pb_par_lapply(frames, FUN = f, 
                         data_dict = data_dict, 
                         cores = cores, 
                         inorder = TRUE, 
                         export_fun_only = TRUE, 
                         loop_text = text)
  } else {
    iou <- lapply(frames, FUN = f, data_dict = data_dict) 
  }
  return(do.call("c", iou))
}

# Check if head moved based on if distance is greater than head_thresh
head_moved <- function(events, head_thresh = 10){
  head_r <- move_seq(events$head_x, events$head_y)$r
  moved_head <- head_r > head_thresh
  
  c(NA, moved_head) # No way of checking first frame
}

# Stage 1 of false head movement detection using head_moved
.false_move_stage1 <- function(events, 
                               target_moved,
                               cen_thresh = 10){
  
  cen_r <- move_seq(events$centroid_x, events$centroid_y)$r
  moved_centroid <- c(NA, cen_r > cen_thresh)
  
  as.logical(target_moved * !moved_centroid) 
}

# Stage 2 of false head movement detection using the more computationally expensive move_IOU()
.false_move_stage2 <- function(data_dict, target_moved, 
                               frames, iou_thresh = 0.5, n = 1, 
                               cores = 1, text = "Verifying suspicious frames:"){
  iou <- move_IOU(data_dict, frames, n = n, cores = cores, text = text)
  mask_moved <- c(iou < iou_thresh)
  
  frames[as.logical(target_moved[frames] * !(mask_moved))]
}


# Wrapper for detecting false head movement stage 1 & 2. 
detect_false_head_movement <- function(data_dict, events = NULL, 
                                       head_thresh = 10,
                                       cen_thresh = 10, 
                                       iou_thresh = 0.5, 
                                       n = 1, 
                                       cores = 1){
  if(is.null(events)){
    events <- get_data(data_dict)
  }
  
  if(head_thresh == "auto"){
    head_thresh <- median(sqrt(events$size_px), na.rm = TRUE) * 0.5
  }
  
  if(cen_thresh == "auto"){
    cen_thresh <- median(sqrt(events$size_px), na.rm = TRUE) * 0.5
  }
  
  
  moved_head <- head_moved(events = events, head_thresh = head_thresh)
  
  candidate_frames <- which(.false_move_stage1(
    events = events, 
    target_moved = moved_head, 
    cen_thresh = cen_thresh))
  
  sus_frames <- .false_move_stage2(data_dict, moved_head, candidate_frames,
                                   iou_thresh = iou_thresh, n = n, cores)
  
  return(sus_frames)
}


# Flag mask size Z score using a window of 31 
flag_mask_size <- function(events, w = 31){
  events %>% 
    mutate(
    size_z = roll_vapply(size_px, w = w, function(x){
      out <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
      if(length(out) == 0){
        return(NA)
      } else {
        return(out[length(out)/2+1])
      }
    })
  ) %>% 
    mutate(
      status = ifelse(
        is.na(size_z),
        "unclear",
        ifelse(size_z > 2, 
               "too big",
               ifelse(
                 size_z < -2, 
                 ifelse(out_of_frame, "out_of_frame", "too small"),
                 "no issue"
               ))
      )
    )
}

# Flag false cluster (long distance movement to the same location)
flag_false_cluster <- function(events, r_thresh = 100, bin_size = 50, r = 50, return_id = FALSE){
  events <- events %>% insert_gaps()
  move_df <- events %$% 
    move_seq(centroid_x, centroid_y)
  
  candidate_seq <- move_df %>% 
    filter(
      r > r_thresh
    ) %>% 
    mutate(
      cx1 = bin(x1, bin_size),
      cy1 = bin(y1, bin_size)
    ) %>% 
    group_by(
      cx1, cy1
    ) %>% 
    mutate(
      n = n()
    ) %>% 
    ungroup()
  
  z <- candidate_seq %>% 
    filter(n == max(n)) %>% 
    .[1,]
  
  
  i <- move_df %>% 
    mutate(
      sus = sqrt((x1 - z$cx1)^2 + (y1 - z$cy1)^2) < r
    ) %>% 
    filter(
      sus
    ) %>% 
    dplyr::select(step_id) %>% 
    unlist(use.names = FALSE)
  
  if(return_id){
    return(i)
  } else {
    return(
      events %>% 
        mutate(
          false_cluster = rank %in% i
        )
    )
  }
}




