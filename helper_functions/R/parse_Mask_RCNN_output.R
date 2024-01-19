# Vectorized function that parse lists from python
parse_pylist <- function(x, coerce = as.numeric, simplify = TRUE){
  out <- gsub("\\[|\\]","",x) %>% 
    str_split(",") %>% 
    lapply(function(x){
      coerce(x)
    })
  
  if(simplify){
    out <- do.call("rbind",out)
    
    if(nrow(out) == 1){
      out <- as.vector(out)
    } 
  }
  
  return(out)
}

# Parse vector into n columns by looping through column index. Undo ravel() in python
split_n_steps <- function(x, n, name = paste0("x", seq_len(n))){
  if(all(is.na(x)) & length(x) == 1 || is.null(x)){
    #x <- rep(x, n)
    return(NULL)
  }
  
  index <- seq_along(x)
  
  out <- lapply(c(seq_len(n-1),0), function(ni){
    x[(index %% n == ni)]
  }) %>% 
    do.call("cbind",.)
  colnames(out) <- name
  
  return(out)
}

# Turn a vector of keypoint values into a matrix
parse_keypoints_vec <- function(x){
  out <- split_n_steps(x, 3, name = c("x","y","score"))
}

# Turn vector of bbox values into a matrix
parse_bbox_vec <- function(x){
  # First row is the top left corner
  # Second row is the bottom right corner
  split_n_steps(x,2, name = c("x","y"))
}

# Parse coco annotation format segmentation coordinates
parse_polygon_vec <- function(x){
  out <- split_n_steps(x, n = 2, name = c("x","y"))
  out <- unique(out)
  if(is.null(out)|| nrow(out) < 3){
    out <- NULL
  }
  return(out)
}

# Parse model info
parse_inference_info <- function(x){
  if(is.null(x) || is.na(x)){
    model_version <- NA
    model_inference_timestamp <- NA
  } else {
    model_version <- gsub("__.*", "", x)
    model_inference_timestamp <- gsub(".*__", "", x)
  }
  
  c("version" = model_version, "inf_time" = model_inference_timestamp)
}

# Take formatted detectron2 predictions from python and parse into 'data_dict' object 
parse_inference <- function(df, mode = c("rep","evaluate")){
  is_rep <- match.arg(mode) == "rep"
  
  
  if(is_rep){
    df <- df[order(file_time(df$file_name)),]
  }
  
  
  master_list <- list(
    "dim" = parse_pylist(df$image_size, simplify = FALSE),
    "file_meta" = df$file_name %>% 
      file_meta() %>% 
      lapply(seq_len(nrow(.)), function(i, df){
        df[i,]
      }, df = .),
    "bbox" = parse_pylist(df$bbox, simplify = FALSE) %>% 
      lapply(function(x){
        parse_bbox_vec(x)
      }),
    "score" = df$score %>% 
      lapply(function(x) {
        x
      }),
    "thing_class" = df$thing_class %>% 
      lapply(function(x) {
        x
      }),
    "keypoints" = parse_pylist(df$keypoints, simplify = FALSE) %>% 
      lapply(function(x){
        parse_keypoints_vec(x)
      }),
    "polygon" = parse_pylist(df$polygon, simplify = FALSE) %>% 
      lapply(function(x){
        parse_polygon_vec(x)
      }), 
    "inference_info" = df$inference_info %>% 
      lapply(function(x){
        parse_inference_info(x)
      })
  )
  
  out_list <- lapply(seq_len(nrow(df)), function(i){
    map(master_list, i)
  })
  names(out_list) <- paste0("frame",
                            file_rank(basename(df$file_name)),
                            "__s", 
                            file_time(basename(df$file_name))
  )
  
  
  
  
  
  out_list <- .data_dict_summary_calc(out_list, is_rep)
  
  return(out_list)
}


.data_dict_summary_calc <- function(out_list, is_rep){
  # Expects list for out_list not data_dict, but returns a data_dict
  stopifnot(!is.data_dict(out_list))
  empty_frames <- do.call("c",lapply(out_list, is.null))
  out_list <- out_list[!empty_frames]
  if(any(empty_frames)){
    message(sprintf("%s empty frames dropped.", sum(empty_frames)))
  }
  if(length(out_list) == 0){
    return(NULL)
  }
  
  fm <- do.call("rbind",map(out_list, "file_meta"))
  n <- length(out_list)
  n_bbox <- n - (map(out_list, "bbox") %>% lapply(is.null) %>% do.call("c",.) %>% sum())
  n_kp <- n - (map(out_list, "keypoints") %>% lapply(is.null) %>% do.call("c",.) %>% sum())
  n_mask <- n - (map(out_list, "polygon") %>% lapply(is.null) %>% do.call("c",.) %>% sum())
  
  if(is_rep){
    n_missing <- count_time_gaps(fm$time)
    n_dups <- count_time_dups(fm$time)
    run_time <- hms_format(diff(range(fm$time)))
    n_steps <- round(diff(range(fm$time))/ 360)
  } else {
    n_missing <- NA
    n_dups <- NA
    run_time <- NA
    n_steps <- NA
    
  }
  
  thing_class <- paste0(na.omit(do.call("c",unique(map(out_list, "thing_class")))), collapse = ",")
  inference_info <- do.call("rbind", map(out_list, "inference_info")) %>% unique()
  
  
  if(all(do.call("rbind", map(out_list, "dim")) %>% 
         apply(2,function(x){length(unique(x)) == 1}))){
    dim_formated <- paste0("(",paste0(out_list[[1]]$dim, collapse = ","),")")
  } else {
    dim_formated <- "Error: multiple dim detected."
  }
  
  
  repID <- paste0(null_to_NA(gsub("rep","",unique(fm$repID))), collapse = ", ")
  camID <- paste0(null_to_NA(gsub("cam","", unique(fm$camID))), collapse = ", ")
  
  
  attr(out_list, "is_rep") <- is_rep
  attr(out_list, "summary") <-  data.frame("repID" = repID, 
                                           "camID" = camID, 
                                           "frames" = n, 
                                           "missing" = n_missing,
                                           "duplicate" = n_dups,
                                           "time_steps" = n_steps,
                                           "duration" = run_time,
                                           "dim" = dim_formated,
                                           "thing_class" = thing_class,
                                           "n_bbox" = n_bbox,
                                           "n_keypoints" = n_kp,
                                           "n_mask" = n_mask,
                                           "version" = null_to_NA(inference_info[,1]),
                                           "inf_time" = null_to_NA(inference_info[,2]))
  
  class(out_list) <- c("data_dict", "list")
  return(out_list)
}



# S3 print method for data_dict objects
print.data_dict <- function(x, ...){
  m <- attr(x, "summary")
  
  n <- null_to_NA(m$frames)
  n_bbox <- null_to_NA(m$n_bbox)
  n_kp <- null_to_NA(m$n_keypoints)
  n_mask <- null_to_NA(m$n_mask)
  n_missing <- null_to_NA(m$missing)
  n_dups <- null_to_NA(m$duplicate)
  n_steps <- null_to_NA(m$time_steps)
  duration <- null_to_NA(m$duration)
  thing_class <- null_to_NA(m$thing_class)
  dim_formated <- null_to_NA(m$dim)
  repID <- null_to_NA(m$repID)
  camID <- null_to_NA(m$camID)
  version <- null_to_NA(m$version)
  inf_time <- null_to_NA(m$inf_time)
  
  cat(sprintf("rep_ID: %s\n\ncam_ID: %s\tdim: %s\n", 
              repID, 
              camID,
              dim_formated
  ))
  cat(sprintf("version: %s\tinference_time: %s\n", 
              version, 
              inf_time
  ))
  cat(sprintf("frames: %s\tmissing: %s\t duplicate: %s\t time_steps: %s \tduration: %s\n", 
              n,
              n_missing,
              n_dups,
              n_steps,
              duration
  ))
  cat(sprintf("thing_class: %s\t bbox: %s\t keypoints: %s\t polygon: %s\t\n", 
              thing_class,
              count_report(n_bbox, n), 
              count_report(n_kp, n),
              count_report(n_mask, n)
  ))
  
  return(
    invisible(
      m
    )
  )
}

# S3 print summary for data_dict objects. Fetches the 'summary' attribute
summary.data_dict <- function(x, ...){
  m <- attr(x, "summary")
  return(m)
}


# S3 '[]' for data_dict objects.
`[.data_dict` <- function(x, ...){
  args <- lapply(as.list(substitute(list(...)))[-1L], 
                 function(x){
                   eval(x, envir = parent.frame(3))
                 })
  class(x) <- c("list")
  out <- do.call("[", c(list(x), args))
  is_rep <- attr(x, "is_rep")
  if(is.null(is_rep)){
    is_rep <- TRUE
    warning("Outdated data_dict format. Setting 'is_rep' to TRUE.")
  }
  out <- .data_dict_summary_calc(out, is_rep)
  return(out)
}

# S3 methods registration
registerS3method(genname = "summary", 
                 class = "data_dict", 
                 method = summary.data_dict)


registerS3method(genname = "print", 
                 class = "data_dict", 
                 method = print.data_dict)

registerS3method(genname = "[", 
                 class = "data_dict", 
                 method = `[.data_dict`)



is.data_dict <- function(x, ...){
  inherits(x, "data_dict")
}
