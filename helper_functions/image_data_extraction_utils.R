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

# Turn polygon into a binary mask
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
  return(as.cimg(mask))
}

# Turn a binary mask into a polygon
mask2polygon <- function(mask){
  l <- spatstat.geom::owin(mask = as.pixset(mask)[,,1,1]) %>% 
    spatstat.geom::as.polygonal() %>% 
    .$bdry %>% 
    .[[1]]
  return(do.call("cbind",l))
}





# Take formatted detectron2 predictions from python and parse into 'data_dict' object 
parse_inference <- function(df){
  df <- df[order(file_time(df$file_name)),]
  
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
  
  
  
  fm <- do.call("rbind",map(out_list, "file_meta"))
  n <- length(out_list)
  n_bbox <- n - (map(out_list, "bbox") %>% lapply(is.null) %>% do.call("c",.) %>% sum())
  n_kp <- n - (map(out_list, "keypoints") %>% lapply(is.null) %>% do.call("c",.) %>% sum())
  n_mask <- n - (map(out_list, "polygon") %>% lapply(is.null) %>% do.call("c",.) %>% sum())
  n_missing <- count_time_gaps(fm$time)
  run_time <- hms_format(diff(range(fm$time)))
  n_steps <- ceiling(diff(range(fm$time))/ 360)
  thing_class <- paste0(na.omit(do.call("c",unique(map(out_list, "thing_class")))), collapse = ",")
  
  if(all(do.call("rbind", map(out_list, "dim")) %>% 
         apply(2,function(x){length(unique(x)) == 1}))){
    dim_formated <- paste0("(",paste0(out_list[[1]]$dim, collapse = ","),")")
  } else {
    dim_formated <- "Error: multiple dim detected."
  }
  
  
  repID <- gsub("rep","",unique(fm$repID))
  camID <- gsub("cam","", unique(fm$camID))
  
  
  
  
  class(out_list) <- c("data_dict", "list")
  attr(out_list, "summary") <-  data.frame("repID" = repID, 
                                           "camID" = camID, 
                                           "frames" = n, 
                                           "missing" = n_missing,
                                           "time_steps" = n_steps,
                                           "duration" = run_time,
                                           "dim" = dim_formated,
                                           "thing_class" = thing_class,
                                           "n_bbox" = n_bbox,
                                           "n_keypoints" = n_kp,
                                           "n_mask" = n_mask)
  
  return(out_list)
}



# S3 print and summary methods for data_dict objects
print.data_dict <- function(x, ...){
  m <- attr(x, "summary")
  
  n <- m$frames
  n_bbox <- m$n_bbox
  n_kp <- m$n_keypoints
  n_mask <- m$n_mask
  n_missing <- m$missing
  n_steps <- m$time_steps
  duration <- m$duration
  thing_class <- m$thing_class
  dim_formated <- m$dim
  repID <- m$repID
  camID <- m$camID
  
  cat(sprintf("rep_ID: %s\n\ncam_ID: %s\tdim: %s\n", 
              repID, 
              camID,
              dim_formated
  ))
  cat(sprintf("frames: %s\tmissing: %s\t time_steps: %s \tduration: %s\n", 
              n,
              n_missing,
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


summary.data_dict <- function(x, ...){
  m <- attr(x, "summary")
  return(m)
}

registerS3method(genname = "summary", 
                 class = "data_dict", 
                 method = summary.data_dict, 
                 envir = asNamespace("herbivar"))


registerS3method(genname = "print", 
                 class = "data_dict", 
                 method = print.data_dict, 
                 envir = asNamespace("herbivar"))

# Method to get formatted key points from 'data_dict' objects
get_keypoints <- function(x){
  stopifnot(inherits(x, "data_dict"))
  
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

# Method to get binary masks from from 'data_dict' objects
get_mask <- function(x, frames){
  stopifnot(inherits(x, "data_dict"))
  x[frames] %>% 
    map("polygon") %>% 
    lapply(function(x){
      polygon2mask(x)
    }) %>% 
    as.imlist()
}

# Method to get summary of binary masks from from 'data_dict' objects
get_mask_summary <- function(x){
  stopifnot(inherits(x, "data_dict"))
  out <- x %>% 
    map("polygon") %>% 
    lapply(function(x){
      mask_info(polygon2mask(x))
    }) %>% 
    do.call("rbind",.)
  cbind("frame_id" = names(x),out)
}

# Method to get file meta data from from 'data_dict' objects
get_file_meta <- function(x){
  stopifnot(inherits(x, "data_dict"))
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
  stopifnot(inherits(x, "data_dict"))
  map(x, "score") %>% 
    bind_vec(keep_row_names = FALSE, row_names_as_col = "frame_id") %>% 
    rename(score = V1)
}

# Method to get image dimension from from 'data_dict' objects
get_dim<- function(x){
  stopifnot(inherits(x, "data_dict"))
  as.numeric(strsplit(gsub("\\(|\\)","",attr(x,"summary")$dim), ",")[[1]])
}

# Method to get repID
get_repID <- function(x){
  stopifnot(inherits(x, "data_dict"))
  attr(x,"summary")$repID
}

# Method to get camID
get_camID <- function(x){
  stopifnot(inherits(x, "data_dict"))
  attr(x,"summary")$camID
}

# Wrapper function for various get_* methods. Joins the output as a single data.frame
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

