# Computer the mask IOU, comparing with n frames ago
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

head_moved <- function(events, head_thresh = 10){
  head_r <- move_seq(events$head_x, events$head_y)$r
  moved_head <- head_r > head_thresh
  
  c(NA, moved_head) # No way of checking first frame
}

.false_move_stage1 <- function(events, 
                               target_moved,
                               cen_thresh = 10){
  
  cen_r <- move_seq(events$centroid_x, events$centroid_y)$r
  moved_centroid <- c(NA, cen_r > cen_thresh)
  
  as.logical(target_moved * !moved_centroid) 
}

.false_move_stage2 <- function(data_dict, target_moved, 
                               frames, iou_thresh = 0.5, n = 1, 
                               cores = 1, text = "Verifying suspicious frames:"){
  iou <- move_IOU(data_dict, frames, n = n, cores = cores, text = text)
  mask_moved <- c(iou < iou_thresh)
  
  frames[as.logical(target_moved[frames] * !(mask_moved))]
}

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
