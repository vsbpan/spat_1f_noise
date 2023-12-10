# Find the centroid and size of a mask
mask_info <- function(vid){
  if(imager::depth(vid) > 1){
    vid %>% 
      as.data.frame() %>% 
      dplyr::select(-cc) %>% 
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
      dplyr::select(-cc) %>% 
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

