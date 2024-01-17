# Find the centroid and size of a mask
mask_info <- function(vid){
  if(imager::spectrum(vid) > 1){
    vid <- vid[,,,1]
  }
  
  vid <- as.pixset(vid)
  
  if(imager::depth(vid) > 1){
    vid %>% 
      as.data.frame() %>% 
      group_by(z) %>% 
      summarise(
        centroid_x = mean(x),
        centroid_y = mean(y),
        size_px = n(),
        x_dim = nrow(vid),
        y_dim = ncol(vid)
      ) %>%
      as.data.frame()
  } else {
    vid %>% 
      as.data.frame() %>% 
      summarise(
        centroid_x = mean(x),
        centroid_y = mean(y),
        size_px = n(),
        x_dim = nrow(vid),
        y_dim = ncol(vid)
      ) %>%
      as.data.frame()
  }
}


# Read value at specified coordinates from reference image. Scales the x,y input with dim_xy of the target image so that it matches relatively to the ref_img
read_value <- function(x, y, dim_xy, ref_img){
  nx <- length(x)
  stopifnot(nx == length(y))
  
  if(is.null(ref_img)){
    return(rep(NA, nx))
  }
  
  stopifnot(
    imager::spectrum(ref_img) == 1 & imager::depth(ref_img) == 1
  )
  
  # As proportion of ref image
  x <- ceiling(x/dim_xy[1] * nrow(ref_img))
  y <- ceiling(y/dim_xy[2] * ncol(ref_img))
  
  return(
    vapply(seq_along(x), function(i){
      as.array(ref_img)[x[i], y[i], 1, 1]
    }, FUN.VALUE = numeric(1))
  )
}


# Apply mask_info() to all the files in src_dir
mask_summary <- function(src_dir, thin.val = 3, cores = 6){
  start_time <- Sys.time()
  
  f <- list.files(src_dir, pattern = ".jpg",full.names = TRUE, recursive = FALSE) %>% 
    files_reorder()
  
  out <- pb_par_lapply(
    f, 
    FUN = function(x, thin.val = thin.val){
      vid <- herbivar::thin(fast_load_image(x, transform = FALSE), thin.val)
      mask_data <- vid %>% mask_info()
      mask_data$time <- file_time(x)
      return(mask_data)
    },
    thin.val = thin.val,
    loop_text = "Reading mask",
    inorder = FALSE
  ) %>% 
    do.call("rbind",.)
  
  hms_runtime(as.numeric(Sys.time() - start_time, units = "secs"))
  
  return(out)
}

# Generate a dummy spectrum for testing
dummy_spec <- function(x){
  rbinom(144, 1, 0.5) %>% 
    image_unflatten()
}

# Generate a dummy binary mask for testing
dummy_mask <-function(z = 100, cores = 6){
  imager::imfill(x =1000, y = 1000, z = z) %>% 
    imsplit("z") %>% 
    pb_par_lapply(function(img){
      x <- runif(1, 1, nrow(img))
      y <- runif(1, 1, ncol(img))
      add_point(img, x, y, color = "white")
    }, cores = cores, inorder = FALSE) %>% 
    imappend("z")
}

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


# Turn polygon into a binary mask
# Kind of finiky and doesn't perform as well as graphics::polygon(), but it'll do. 
polygon2mask <- function(x,y = NULL, dim_xy = c(1000, 1000)){
  if(is.null(x)){
    mask <- imager::imfill(x = dim_xy[1], y = dim_xy[2], val = 0)
    return(mask)
  }
  
  if(is.null(y)){
    y <- x[,2]
    x <- x[,1]
  }
  
  x <- round(x)
  y <- round(y)
  
  hull <- spatstat.geom::convexhull.xy(x,y)
  
  if(is.null(hull)){
    mask <- imager::imfill(x = dim_xy[1], y = dim_xy[2], val = 0)
    return(mask)
  }
  
  mask <- hull %>% 
    spatstat.geom::as.mask(dimyx = c(diff(range(x)), diff(range(y))), 
                           xy = list(x = seq_len(dim_xy[1]), y = seq_len(dim_xy[2]))) %>% 
    spatstat.geom::as.array.im()
  
  dim(mask) <- c(dim(mask)[1:2], 1, dim(mask)[3])
  return(as.cimg(mask) %>% rotate_90())
}

# Turn a binary mask into a polygon
mask2polygon <- function(mask){
  l <- spatstat.geom::owin(mask = as.pixset(mask)[,,1,1]) %>% 
    spatstat.geom::as.polygonal() %>% 
    .$bdry %>% 
    .[[1]]
  return(do.call("cbind",l)[,c(2,1)])
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
  class(out_list) <- c("data_dict", "list")

  
  return(out_list)
}

.data_dict_summary_calc <- function(out_list, is_rep){
  
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
                                           "version" = inference_info[,1],
                                           "inf_time" = inference_info[,2])
  out_list
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
  class(out) <- c("data_dict", "list")
  out <- .data_dict_summary_calc(out, attr(x, "is_rep"))
  out
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
          return(NA)
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
  stopifnot(!missing(frames))
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
  out <- x %>% 
    get_polygon() %>% 
    lapply(function(x){
      cbind(mask_info(polygon2mask(x)), "out_of_frame" = out_of_frame(x))
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


is.data_dict <- function(x, ...){
  inherits(x, "data_dict")
}



