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


# Read value at specified coordinates from reference image. Scales the x,y input with side_dim of the target image so that it matches relatively to the ref_img
read_value <- function(x, y, side_dim, ref_img){
  
  stopifnot(
    imager::spectrum(ref_img) == 1 & imager::depth(ref_img) == 1
  )
  
  # As proportion of ref image
  x <- ceiling(x/side_dim * nrow(ref_img))
  y <- ceiling(y/side_dim * ncol(ref_img))
  
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

parse_keypoints_vec <- function(x){
  split_n_steps(x, 3, name = c("x","y","score"))
}


parse_bbox_vec <- function(x){
  # First row is the top left corner
  # Second row is the bottom right corner
  split_n_steps(x,2, name = c("x","y"))
}

# Parse coco annotation format segmentation coordinates
parse_polygon_vec <- function(x){
  return(split_n_steps(x, n = 2, name = c("x","y")))
}

# Turn polygon into a binary mask
polygon2mask <- function(x,y = NULL, dim_xy = c(1000, 1000)){
  if(is.null(y)){
    y <- x[,2]
    x <- x[,1]
  }
  
  x <- round(x)
  y <- round(y)
  
  mask <-spatstat.geom::convexhull.xy(x,y) %>% 
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




# S3 print and summary methods for data_dict objects
print.data_dict <- function(x, ...){
  fm <- do.call("rbind",map(x, "file_meta"))
  n <- length(x)
  n_bbox <- n - (map(x, "bbox") %>% lapply(is.null) %>% do.call("c",.) %>% sum())
  n_kp <- n - (map(x, "keypoints") %>% lapply(is.null) %>% do.call("c",.) %>% sum())
  n_mask <- n - (map(x, "polygon") %>% lapply(is.null) %>% do.call("c",.) %>% sum())
  n_missing <- count_time_gaps(fm$time)
  run_time <- hms_format(diff(range(fm$time)))
  thing_class <- paste0(na.omit(do.call("c",unique(map(x, "thing_class")))), collapse = ",")
  dim_formated <- paste0("(",paste0(x[[1]]$dim, collapse = ","),")")
  repID <- gsub("rep","",unique(fm$repID))
  camID <- gsub("cam","", unique(fm$camID))

  cat(sprintf("rep_ID: %s\n\ncam_ID: %s\n", 
        repID, 
        camID
        ))
  cat(sprintf("frames: %s\tmissing: %s\t duration: %s\n", 
        n,
        n_missing,
        run_time
        ))
  cat(sprintf("dim: %s\n\nthing_class: %s\t bbox: %s\t keypoints: %s\t polygon: %s\t\n", 
              dim_formated,
              thing_class,
              count_report(n_bbox, n), 
              count_report(n_kp, n),
              count_report(n_mask, n)
              ))
  
  return(
    invisible(
      data.frame("repID" = repID, 
        "camID" = camID, 
        "frames" = n, 
        "missing" = n_missing,
        "duration" = run_time,
        "dim" = dim_formated,
        "thing_class" = thing_class,
        "n_bbox" = n_bbox,
        "n_kepoints" = n_kp,
        "n_mask" = n_mask
        )
    )
  )
}

registerS3method(genname = "summary", 
                 class = "data_dict", 
                 method = print.data_dict, 
                 envir = asNamespace("herbivar"))


registerS3method(genname = "print", 
                 class = "data_dict", 
                 method = print.data_dict, 
                 envir = asNamespace("herbivar"))


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
  names(out_list) <- gsub(".jpg", "",basename(df$file_name))
  
  class(out_list) <- c("data_dict", "list")
  return(out_list)
}

