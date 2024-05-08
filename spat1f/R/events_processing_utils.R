# Calculate turn angle with a sequence of absolute angles
turn_angle_calc <- function(theta){
  # theta in absolute angle in order of time steps
  where_NA <- is.na(theta)
  theta <- theta[!is.na(theta)]
  
  theta_rel <- c(0, (theta - c(theta[-1],NA))[-length(theta)])
  theta_rel <- theta_rel %% (2 * pi)
  theta_rel <-  ifelse(theta_rel < pi, 
                           theta_rel, # Right turns
                           (theta_rel %% pi) - pi) # Convert left turns into negative angle
  
  
  # Swap left right to match output of amt::steps()
  # step length passes check with identical()
  # turn angle passes check with error tolerance of .Machine$double.eps*4 radians
  theta_rel <- -theta_rel
  
  out <- rep(NA, length(where_NA))
  
  
  out[!where_NA] <- theta_rel # Insert calculated theta_rel back in place, assuming those with NA stay at the same place
  i <- (c(0,which(where_NA)) + 1) # Set the first turn angle of an uninterrupted sequence of movement vectors as NA 
  i <- i[i<length(where_NA)] # Remove if i is longer than the length of the target vector
  out[i] <- NA
  out
}

# Convert radians into degrees
rad2degree <- function(theta){
  theta / pi * 180
}

# Calculate step length and turn angle. step lengths smaller than r_thresh are considered none movement, causing the turn angle of that step to be NA
move_seq <- function(x,y, r_thresh = 1, inherit.theta = FALSE){ 
  stopifnot(length(x) == length(y))
  
  if(length(x) < 2){
    dim_names <- c(
      "step_id",
      "r",
      "r_threshed",
      "theta_abs",
      "theta_rel",
      "x1","x2", "y1", "y2"
    )
    return(
      data.frame(matrix(nrow = 0, 
                        ncol = length(dim_names), 
                        dimnames = list(NULL, 
                                        dim_names)))
    )
  }
  
  out <- lapply(seq_len(length(x) - 1), function(i){
    delta_x <- x[i] - x[i+1]
    delta_y <- y[i] - y[i+1]
    r <- sqrt((delta_x)^2 + (delta_y)^2)
    theta_abs <- atan2(delta_y, delta_x) # Angle with respect to the xy coordinate.
    c("r" = r, "r_threshed" = ifelse(r > r_thresh, r, 0), "theta_abs" = theta_abs)
  }) %>% bind_vec()
  
  out$theta_abs[out$r < r_thresh] <- NA # Treat those with very small steps as stand still
  
  out$theta_rel <-  turn_angle_calc(out$theta_abs) # Get turn angle with vector at t-1 as 0 degree reference
  n <- nrow(out)
  out <- data.frame("step_id" = seq_len(n), 
             out,
             "x1" = x[-n], 
             "x2" = x[-1], 
             "y1" = y[-n],
             "y2" = y[-1])
  
  if(inherit.theta){
    out <- out %>% 
      mutate(
        theta_abs = inherit_theta(theta_abs, r_threshed),
        theta_rel = turn_angle_calc(theta_abs)
      )
  }
  
  attr(out, "r_thresh") <- r_thresh
  return(out)
}


# Convert treatment metadata from data.frame to image list
trt_meta_as_list <- function(df){
  ufl_trt_iml <- lapply(seq_len(nrow(df)), function(i){
    x <- df[i,]
    unlist(x[,grepl("spec_", names(x))]) 
  }) %>% 
    lapply(
      image_unflatten
    )
  
  names(ufl_trt_iml) <- paste0("syn_id__", df$syn_id)
  ufl_trt_iml
}


# Insert gaps of NAs at positions where a photo is expected
insert_gaps <- function(df, frame_id = frame_id, expected_gap = 360){
  if("is_gap" %in% names(df)){
    df <- df %>% filter(!is_gap)
  }
  
  .expose_columns_interal()
  
  time <- c(file_time(frame_id), expected_gap * -1) # Append a negative time photo to ensure that there is always a beginning without NA at the same time
  time_grid <- 360 * (seq_len((diff(round(range(time) / expected_gap)) + 1)) - 1)
  o <- order(time)
  df <- df[o,]
  g <- ((diff(time[o]) - expected_gap) / expected_gap) %>% 
    round()
  g <- pmax(c(g, 0), 0)

  index <- lapply(seq_along(g), function(i, g){
    c(i, rep(1, g[i]))
  }, g = g) %>% 
    do.call("c",.)
  
  out <- df[index,][-1, ] # Drop the negative time photo
  rownames(out) <- NULL
  is_gap <- is.na(out$frame_id)
  out$is_gap <- is_gap
  out[is_gap,"frame_id"] <- paste0("gap", 
                                   seq_len(sum(g)), 
                                   "__s", 
                                   time_grid[which(is_gap)])
  return(out)
}

# Remove 'rep' from repID string for comparability reasons
repID_clean <- function(x){
  if(is.data_dict(x)){
   x <- get_repID(x) 
  }
  gsub("rep", "", x)
}

# Loop through each data_dict objects using the supplied repIDs vector and report detection rate
detection_report <- function(repIDs){
  repIDs <- repID_clean(repIDs)
  valid_ids <- fetch_repID()
  repIDs <- repIDs[repIDs %in% valid_ids]
  
  repIDs %>% 
    lapply(
      function(x){
        summary(fetch_data_dict(x))
      }
    ) %>% 
    do.call("rbind",.) %>% 
    mutate(
      mask = count_report(n_mask , frames),
      keypoints = count_report(n_keypoints, frames),
      bbox = count_report(n_bbox, frames)
    )
}


# Replace NA with theta from previous step if the step length is zero
inherit_theta <- function(theta, r){
  theta_inht <- theta[1]
  for(i in seq_along(theta)){
    if(is.na(theta[i])){
      if(!is.na(r[i]) && r[i] == 0){
        theta[i] <- theta_inht 
      }
    } else {
      theta_inht <- theta[i]
    }
  }
  return(theta)
}


# Clean events by throwing out frames beyond the camera cutoff, flagged errors, low scoring detections, and inserting gaps where needed. 
clean_events <- function(x, 
                         ref_data = get("ref_data", envir = globalenv()), 
                         score_thresh = 0.7,
                         insert_gaps = TRUE,
                         keep_sus = FALSE){
  repID <- unique(repID_clean(x$repID))
  stopifnot(length(repID) == 1)
  sus_frames <- fetch_sus_frames(repID)
  cut_off <- fetch_cutoff(repID, ref_data)
  
  x <- x %>%
    mutate(
      is_sus = ifelse(seq_len(nrow(x)) %in% sus_frames, TRUE, FALSE)
    ) %>% 
    filter(
      time <= cut_off
    )
  if(!keep_sus){
    x <- filter(x, !is_sus) %>% 
      filter(
        !status %in% c("too big", "too small") & !false_cluster
      )
  }
  
  if(!is.null(score_thresh)){
    x <- x %>% filter(score >= score_thresh)
  }
  
  if(insert_gaps){
    if(nrow(x) < 2){
      return(x)
    }
    x <- insert_gaps(x)
  }
  
  return(x)
}



