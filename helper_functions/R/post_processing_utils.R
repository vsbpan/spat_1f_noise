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
  out <- rep(NA, length(where_NA))
  
  
  out[!where_NA] <- theta_rel # Insert calculate theta_rel back in place, assuming those with NA stay at the same place
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
move_seq <- function(x,y, r_thresh = 1){ 
  stopifnot(length(x) == length(y))
  out <- lapply(seq_len(length(x) - 1), function(i){
    delta_x <- x[i] - x[i+1]
    delta_y <- y[i] - y[i+1]
    r <- sqrt((delta_x)^2 + (delta_y)^2)
    theta_abs <- atan2(delta_y, delta_x) # Angle with respect to the xy coordinate.
    c("r" = r, "theta_abs" = theta_abs)
  }) %>% bind_vec()
  
  out$theta_abs[out$r < r_thresh] <- NA # Treat those with very small steps as stand still
  
  out$theta_rel <-  turn_angle_calc(out$theta_abs) # Get turn angle with vector at t-1 as 0 degree reference
  out$step_id <- seq_len(nrow(out))
  out
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
# A bit buggy
insert_gaps <- function(df, expected_gap = 360){
  time <- file_time(df$frame_id)
  time_grid <- 360 * (seq_len((diff(round(range(time) / expected_gap)) +1)) - 1)
  o <- order(time)
  df <- df[o,]
  g <- ((diff(time[o]) - expected_gap) / expected_gap) %>% 
    round()
  g <- pmax(c(g, 0), 0)
  
  
  
  empty_row <- as.data.frame(matrix(rep(NA, ncol(df)), 
                                    nrow = 1, ncol = ncol(df), 
                                    dimnames = list(NULL, names(df))))
  
  out <- lapply(seq_along(g), function(i, g){
    rbind(df[i,], empty_row[rep(TRUE, g[i]),])
  }, g = g) %>% 
    do.call("rbind",.)
  
  rownames(out) <- NULL
  out[is.na(out$frame_id),"frame_id"] <- paste0("gap", 
                                                seq_len(sum(g)), 
                                                "__s", time_grid[is.na(out$frame_id)])
  out$is_gap <- grepl(out$frame_id,"gap")
  return(out)
}

# Remove 'rep' from repID string for comparability reasons
repID_clean <- function(x){
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



