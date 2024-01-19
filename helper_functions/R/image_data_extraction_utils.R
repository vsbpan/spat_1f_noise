# Find the centroid and size of a mask
mask_info <- function(x){
  if(is.data_dict(x)){
    n_frames <- length(x)
    x <- get_polygon(x)
    if(n_frames == 1){
      x <- x[[1]]
    }
    if(n_frames == 0){
      return(NULL)
    }
  }
  
  if(is.cimg(x) || is.pixset(x)){
    if(imager::spectrum(x) > 1){
      x <- x[,,,1]
    }
    
    x <- as.pixset(x)
    
    if(imager::depth(x) > 1){
      res <- x %>% 
        as.data.frame() %>% 
        group_by(z) %>% 
        summarise(
          centroid_x = mean(x),
          centroid_y = mean(y),
          size_px = n()
        ) %>%
        as.data.frame()
    } else {
      res <- x %>% 
        as.data.frame() %>% 
        summarise(
          centroid_x = mean(x),
          centroid_y = mean(y),
          size_px = n()
        ) %>%
        as.data.frame()
    }
  } else {
    if(is.matrix(x) || is.data.frame(x)){
      cen <- polygon_centroid(x)
      res <- data.frame(
        "centroid_x" = cen[1], 
        "centroid_y" = cen[2], 
        "size_px" = mask_area(x)
      )
    } else {
      
      if(is.list(x)){
        
        res <- lapply(x, function(xi){
          if(is.null(xi)){
            area <- NA
            cen <- c(NA, NA)
          } else {
            cen <- polygon_centroid(xi)
            area <- mask_area(xi)
          }
          
          data.frame(
            "centroid_x" = cen[1], 
            "centroid_y" = cen[2],
            "size_px" = area
          )
        }) %>% 
          do.call("rbind", .)
        res <- cbind(
          "z" = seq_len(nrow(res)), 
          res
        )
      } else {
        if(is.null(x)){
          res <- data.frame(
            "centroid_x" = NA, 
            "centroid_y" = NA,
            "size_px" = NA
          )
        } else {
          stop("Object type not supported.")
        }
      }
    }
  }
  
  return(res)
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


# # Apply mask_info() to all the files in src_dir
# mask_summary <- function(src_dir, thin.val = 3, cores = 6){
#   start_time <- Sys.time()
#   
#   f <- list.files(src_dir, pattern = ".jpg",full.names = TRUE, recursive = FALSE) %>% 
#     files_reorder()
#   
#   out <- pb_par_lapply(
#     f, 
#     FUN = function(x, thin.val = thin.val){
#       vid <- herbivar::thin(fast_load_image(x, transform = FALSE), thin.val)
#       mask_data <- vid %>% mask_info()
#       mask_data$time <- file_time(x)
#       return(mask_data)
#     },
#     thin.val = thin.val,
#     loop_text = "Reading mask",
#     inorder = FALSE
#   ) %>% 
#     do.call("rbind",.)
#   
#   hms_runtime(as.numeric(Sys.time() - start_time, units = "secs"))
#   
#   return(out)
# }

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



# Turn polygon into a binary mask
# Kind of finiky and doesn't perform as well as graphics::polygon(), but it'll do. 
polygon2mask <- function(x,y = NULL, dim_xy = c(1000, 1000), 
                         mini_mask = FALSE, raw_mat = FALSE){
  if(is.null(x)){
    mask <- imager::imfill(x = dim_xy[1], y = dim_xy[2], val = 0)
    if(raw_mat){
      return(NULL)
    }
    return(mask)
  }
  
  if(is.null(y)){
    y <- x[,2]
    x <- x[,1]
  }
  
  x <- round(x)
  y <- round(y)
  
  mask <- spatstat.geom::owin(poly = list(
    x = rev(x), y = rev(y), check = FALSE
  ))
  
  if(mini_mask){
    if(is.null(mask)){
      if(raw_mat){
        return(mask)
      } else {
        mask <- imager::imfill(x = diff(range(x)), y = diff(range(y)), val = 0)
      }
      return(mask)
    }
    mask <- mask %>% 
      spatstat.geom::as.mask(dimyx = c(diff(range(x)), diff(range(y))), 
                             xy = list(x = seq(min(x), max(x)), y = seq(min(y), max(y))))
    
    if(raw_mat){
      return(mask)
    } else {
      mask <- mask %>% spatstat.geom::as.array.im()
    }
    
  } else {
    if(is.null(mask)){
      mask <- imager::imfill(x = dim_xy[1], y = dim_xy[2], val = 0)
      return(mask)
    }
    mask <- mask %>% 
      spatstat.geom::as.mask(dimyx = c(diff(range(x)), diff(range(y))), 
                             xy = list(x = seq_len(dim_xy[1]), y = seq_len(dim_xy[2]))) %>% 
      spatstat.geom::as.array.im()
  }
  
  dim(mask) <- c(dim(mask)[1:2], 1, dim(mask)[3])
  return(as.cimg(mask) %>% rotate_90())
}

# Turn a binary mask into a polygon
mask2polygon <- function(mask){
  l <- spatstat.geom::owin(mask = as.pixset(mask)[,,1,1]) %>% 
    spatstat.geom::as.polygonal() %>% 
    .$bdry %>% 
    .[[1]]
  out <- do.call("cbind",l)[,c(2,1)]
  colnames(out) <- c("x", "y")
  return(out)
}

# Calculate mask centroid using polygon
polygon_centroid <- function(poly){
  o <- polygon2mask(poly, mini_mask = TRUE, raw_mat = TRUE)
  c(mean_wt(o$xcol, colSums(o$m)), mean_wt(o$yrow, rowSums(o$m)))
}



