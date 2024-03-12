# Fetch treatment spectra meta data stored on computer
fetch_trt_meta <- function(path = "raw_data/trt_spectra_meta/master_trt_meta.csv"){
  suppressMessages(read_csv(path))
}

# Fetch treatment spectra using repID. 
# For use in pb_par_lapply() or foreach(), set get("ref_data", envir = parent.frame())
fetch_trt_spec <- function(repID, .ref_data = get("ref_data", envir = globalenv()), 
                           trt_meta_iml = NULL, quiet = FALSE){
  
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
  if(!quiet){
    cat(sprintf("Fetached synID '%s' for 'rep%s'\n", syn_id_matched, repID))
  }
  
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


# Fetch store anchor list stored on computer using repID
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

# Fetch frame by frame data stored on computer using repID. Processed with get_data() on 'data_dict' object
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

# Fetch 'data_dict' objects (parsed form of RCNN inference data.frame) using repID
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

# Fetch cut off time for camera
fetch_cutoff <- function(repID, ref_data = get("ref_data", envir = globalenv())){
  repID <- repID_clean(repID)
 (ref_data[ref_data$rep_id == repID, ])$camera_cutoff
}


# Fetch suspicious frames where there is false head movement due to misidentification. 
fetch_sus_frames <- function(repID, sus_frame_list_path = "cleaned_data/sus_frames_list.rds"){
  repID <- repID_clean(repID)
  sus_list <- readRDS(sus_frame_list_path)
  map(sus_list, "sus_frames")[[which((do.call("c",map(sus_list, "rep_id")) == repID))]]
}

# Fetch all the repIDs with inference on the computer
fetch_repID <- function(has = c("inference", "processed", "raw", 
                                "events", "data_dict", "anchors")){
  has <- match.arg(has)
  
  if(has == "inference"){
    x <- list.files("raw_data/inferences") %>% 
      gsub("rep|_inference.csv","",.) 
  }
  if(has == "processed"){
    x <- list.dirs("processed_feed/", recursive = FALSE, full.names = FALSE) %>% 
      gsub("rep","",.)
  }
  if(has == "raw"){
    x <- list.dirs("time_lapse_feed/", recursive = FALSE, full.names = FALSE) %>% 
      gsub("rep","",.)
  }
  if(has == "events"){
    x <- list.files("cleaned_data/events") %>% gsub("rep|.csv","",.)
  }
  if(has == "data_dict"){
    x <- list.files("cleaned_data/data_dicts")%>% gsub("rep|.rds","",.)
  }
  if(has == "anchors"){
    x <- list.files("raw_data/picked_anchors")%>% gsub("rep|.rds","",.)
  }
  
  return(x)
}

# Fetch image with file rank and repID or `data_dict`
fetch_image <- function(x, rank = NULL, time = NULL, transform = TRUE){
  if(!is.data_dict(x)){
    x <- fetch_data_dict(x)
  }
  
  if(is.null(rank)){
    rank <- time2rank(x, time)
  }
  
  fm <- get_file_meta(x)
  
  fp <- fm[fm$rank == rank, "file_path"]
  stopifnot(length(fp) == 1)
  if(!file.exists(fp)){
    stop(sprintf("File does not exist on this machine: %s", basename(fp)))
  } else {
    fast_load_image(fp, transform = transform)
  }
}



# Method to get formatted keypoints from 'data_dict' objects
get_bbox <- function(x){
  stopifnot(is.data_dict(x))
  lapply(map(x, "bbox"), function(z){
    if(is.null(z)){
      return(NULL)
    } else {
      return(as.vector(t(z))) 
    }
  }) 
}


# Method to get formatted key points from 'data_dict' objects
get_keypoints <- function(x){
  stopifnot(is.data_dict(x))
  
  labs <- c("head", "middle", "tail")
  
  out <- lapply(1:3, function(i){
    map(x, "keypoints") %>% 
      lapply(function(x){
        o <- x[i,]
        if(is.null(o)){
          return(c("x" = NA, "y" = NA, "score" = NA))
        }
        o
      }) %>% 
      bind_vec(keep_row_names = TRUE) %>% 
      #cbind("label" = labs[i]) %>% 
      rename_all(.funs = function(x){
        paste0(labs[i],"_",x)
      }) 
  }) %>% 
    do.call("cbind",.)
  
  out <- cbind("frame_id" = rownames(out), out)
  rownames(out) <- NULL
  return(out)
}

# Method to get polygons 'data_dict' objects
get_polygon <- function(x){
  stopifnot(is.data_dict(x))
  map(x, "polygon")
}



# Method to get binary masks from from 'data_dict' objects
get_mask <- function(x, frames){
  stopifnot(is.data_dict(x))
  frames <- assert_frames(x, frames)
  
  x[frames] %>% 
    get_polygon() %>% 
    lapply(function(x){
      polygon2mask(x)
    }) %>% 
    as.imlist()
}

# Method to get summary of binary masks from from 'data_dict' objects
get_mask_summary <- function(x){
  stopifnot(is.data_dict(x))
  dim_xy <- get_dim(x)
  out <- x %>% 
    get_polygon() %>% 
    lapply(function(x){
      cbind(mask_info(x), "x_dim" = dim_xy[1], "y_dim" = dim_xy[2], 
            "out_of_frame" = out_of_frame(x))
    }) %>% 
    do.call("rbind",.)
  cbind("frame_id" = names(x),out)
}

# Method to get file meta data from from 'data_dict' objects
get_file_meta <- function(x){
  stopifnot(is.data_dict(x))
  out <- map(x, "file_meta") %>% 
    lapply(
      function(x){
        x
      }
    ) %>% 
    do.call("rbind",.)
  cbind("frame_id" = names(x),out)
}

# Method to get bbox confidence score from from 'data_dict' objects
get_score <- function(x){
  stopifnot(is.data_dict(x))
  map(x, "score") %>% 
    bind_vec(keep_row_names = FALSE, row_names_as_col = "frame_id") %>% 
    rename(score = V1)
}

# Method to get image dimension from from 'data_dict' objects
get_dim<- function(x){
  stopifnot(is.data_dict(x))
  as.numeric(strsplit(gsub("\\(|\\)","",attr(x,"summary")$dim), ",")[[1]])
}

# Method to get repID
get_repID <- function(x){
  stopifnot(is.data_dict(x))
  attr(x,"summary")$repID
}

# Method to get camID
get_camID <- function(x){
  stopifnot(is.data_dict(x))
  attr(x,"summary")$camID
}

# Wrapper function for various get_* methods on 'data_dict' objects. Joins the output as a single data.frame
get_data <- function(x, type = c("file_meta", "score", "keypoints", "mask_summary")){
  for (i in seq_along(type)){
    method <- switch(type[i],
                     "file_meta" = get_file_meta, 
                     "score" = get_score,
                     "keypoints" = get_keypoints,
                     "mask_summary" = get_mask_summary)
    if(i == 1){
      out <- method(x)
    } else {
      out <- full_join(out, method(x), by = "frame_id")
    } 
  }
  return(out)
}
