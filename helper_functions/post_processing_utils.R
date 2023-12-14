
move_seq <- function(x,y){
  stopifnot(length(x) == length(y))
  out <- lapply(seq_len(length(x) - 1), function(i){
    delta_x <- x[i] - x[i+1]
    delta_y <- y[i] - y[i+1]
    r <- sqrt((delta_x)^2 + (delta_y)^2)
    theta_abs <- atan2(delta_y, delta_x)
    c("r" = r, "theta_abs" = theta_abs)
  }) %>% bind_vec()
  out$theta_rel <- c(NA, (out$theta_abs - c(out$theta_abs[-1],NA))[-(length(x)-1)])
  out$step_id <- seq_len(nrow(out))
  out
}


fetch_trt_meta <- function(path = "raw_data/trt_spectra_meta/master_trt_meta.csv"){
  suppressMessages(read_csv(path))
}

fetch_trt_spec <- function(repID, .ref_data = get("ref_data"), trt_meta_iml = NULL){
  
  repID <- repID_clean(repID) # Cleaning
  
  syn_id_matched <- as.vector(unlist(.ref_data[which(.ref_data$rep_id == repID), "syn_id"]))
  if(is.na(syn_id_matched)){
    warning(sprintf("Can't find syn_ID that maps to repID = %s", repID))
    return(NULL)
  }
  
  if(is.null(trt_meta_iml)){
    trt_meta_iml <- fetch_trt_meta() %>% 
      filter(syn_id == syn_id_matched) %>% 
      trt_meta_as_list()
  }
  
  cat(sprintf("Fetached synID '%s' for 'rep%s'\n", syn_id_matched, repID))
  index <- which(names(trt_meta_iml) == paste0("syn_id__",syn_id_matched))
  
  if(length(index) != 1){
    stop(
      sprintf(
        "Expected 1 matched spectrum, but %s found.", 
        length(index)
      )
    )
  }
  
  trt_meta_iml[[index]]
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

fetch_anchors <- function(repID, src_dir = "raw_data/picked_anchors/"){
  repID <- repID_clean(repID)
  path <- paste0(src_dir, "/rep",repID,".rds")
  
  if(file.exists(path)){
    return(readRDS(path))
  } else {
    warning(
      sprintf("repID: %s not found in %s", repID, src_dir)
    )
    return(NULL)
  }
}

fetch_events <- function(repID, append_detection_summary = TRUE, src_dir = "cleaned_data/events/"){
  repID <- repID_clean(repID)
  path <- paste0(src_dir, "/rep",repID,".csv")
  
  if(file.exists(path)){
    out <- suppressMessages(read_csv(path))
    if(append_detection_summary){
      s <- summary(
        fetch_data_dict(repID)) %>% 
        mutate(
          repID = paste0("rep", repID)
        ) %>% 
        dplyr::select(-camID)
      
      if(nrow(s) > 1){
        stop(sprintf("Expects 1 row, but got %s rows of data_dict summary"), nrow(s))
      }
      
      out <- out %>% 
        left_join(
          s, 
          by = "repID")
    }
    return(out)
  } else {
    warning(
      sprintf("repID: %s not found in %s", repID, src_dir)
    )
    return(NULL)
  }
}

fetch_data_dict <- function(repID, src_dir = "cleaned_data/data_dicts/"){
  repID <- repID_clean(repID)
  path <- paste0(src_dir, "/rep",repID,".rds")
  
  if(file.exists(path)){
    return(readRDS(path))
  } else {
    warning(
      sprintf("repID: %s not found in %s", repID, src_dir)
    )
    return(NULL)
  }
}


fetch_repID <- function(has = c("inference")){
  list.files("raw_data/inferences") %>% 
    gsub("rep|_inference.csv","",.)
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
  g <- c(g, 0)
  
  
  
  empty_row <- as.data.frame(matrix(rep(NA, ncol(df)), 
                                    nrow = 1, ncol = ncol(df), 
                                    dimnames = list(NULL, names(df))))
  
  out <- lapply(seq_along(g), function(i, g){
    rbind(d2[i,], empty_row[rep(TRUE, g[i]),])
  }, g = g) %>% 
    do.call("rbind",.)
  
  rownames(out) <- NULL
  out[is.na(out$frame_id),"frame_id"] <- paste0("gap", seq_len(sum(g)), "__s", time_grid[is.na(out$frame_id)])
  return(out)
}

repID_clean <- function(x){
  gsub("rep", "", x)
}

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



