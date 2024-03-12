
# Turns image into a single vector
image_flatten <- function(img){
  as.numeric(c(img))
}

# Restores flattened image into a cimg object
image_unflatten <- function(x){
  n <- length(x)
  ns <- sqrt(n)
  stopifnot(ns%%1 == 0)
  mat <- matrix(x, nrow = ns, ncol = ns)
  
  as.cimg(mat)
}

# Generate mapping function for quad to quad transformation (useful for perspective correction)
q2q_map <- function(plist, qlist, z = FALSE){
  p1x <- plist[[1]][[1]]
  p2x <- plist[[2]][[1]]
  p3x <- plist[[3]][[1]]
  p4x <- plist[[4]][[1]]
  p1y <- plist[[1]][[2]]
  p2y <- plist[[2]][[2]]
  p3y <- plist[[3]][[2]]
  p4y <- plist[[4]][[2]]
  
  
  q1x <- qlist[[1]][[1]]
  q2x <- qlist[[2]][[1]]
  q3x <- qlist[[3]][[1]]
  q4x <- qlist[[4]][[1]]
  q1y <- qlist[[1]][[2]]
  q2y <- qlist[[2]][[2]]
  q3y <- qlist[[3]][[2]]
  q4y <- qlist[[4]][[2]]
  
  
  mat <- c(
    p1x, p1y, 1, 0, 0, 0, -p1x*q1x, -p1y*q1x, 
    0, 0, 0, p1x, p1y, 1,  -p1x*q1y, -p1y*q1y,
    p2x, p2y, 1, 0, 0, 0, -p2x*q2x, -p2y*q2x, 
    0, 0, 0, p2x, p2y, 1,  -p2x*q2y, -p2y*q2y,
    p3x, p3y, 1, 0, 0, 0, -p3x*q3x, -p3y*q3x, 
    0, 0, 0, p3x, p3y, 1,  -p3x*q3y, -p3y*q3y,
    p4x, p4y, 1, 0, 0, 0, -p4x*q4x, -p4y*q4x, 
    0, 0, 0, p4x, p4y, 1,  -p4x*q4y, -p4y*q4y
  )
  
  A <- t(matrix(mat, nrow = 8, ncol = 8))
  B <- matrix(c(q1x, q1y, q2x, q2y, q3x, q3y, q4x, q4y), ncol = 1)
  
  X <- solve(A, B)
  
  if(z){
    map_fun <- function(x,y,z) {
      list(
        x = (X[1] * x + X[2] * y + X[3]) / (1 + X[7] * x + X[8] * y),
        y = (X[4] * x + X[5] * y + X[6]) / (1 + X[7] * x + X[8] * y),
        z = z
      )
    }
  } else {
    map_fun <- function(x,y) {
      list(
        x = (X[1] * x + X[2] * y + X[3]) / (1 + X[7] * x + X[8] * y),
        y = (X[4] * x + X[5] * y + X[6]) / (1 + X[7] * x + X[8] * y)
      )
    }
  }
  
  return(map_fun)
}

#' @title Apply quad to quad transformation on an image
#' @description
#' Apply quad to quad transformation on an image using \code{q2q_map()}. Useful for perspective correction. For more details see \code{imager::imwarp()}. 
#' @param img A cimg object
#' @param dest_pts,init_pts list of lists of four coordinates (x,y) setting the location of initial and destination of the transformation. If \code{NULL} (default), the four corners of the image will be selected. See examples. 
#' @param retain_size if \code{TRUE}, the returned image retains the same size. If \code{FALSE} (default), if the destination coordinates are outside of the original image, the image is enlarged. 
#' @param coordinates \code{"absolute"} or \code{"relative"} (default \code{"absolute"})
#' @param boundary boundary conditions: \code{"dirichlet", "neumann", "periodic"}. Default \code{"dirichlet"}.
#' @param interpolation \code{"nearest", "linear", "cubic"} (default \code{"nearest"})
#' @examples
#' img <- image_example() 
#' 
#' original_coords <- list(
#'     list(1,1),
#'     list(1, 443),
#'     list(810, 443),
#'     list(810, 1)
#'    )
#' trans_coords <- list(
#'       list(1,1), # Top left
#'       list(1, 200), # Bottom left
#'       list(500, 300), # Bottom right
#'       list(400, 1) # Top right
#'     )
#' 
#' # Apply perspective transformation
#' trans_img <- q2q_trans(
#'     img,   
#'     dest_pts = trans_coords,
#'     init_pts = original_coords
#'   )
#'  
#'  
#' # Undo perspective transformation
#'  backtrans_img <- q2q_trans(
#'     trans_img,   
#'     dest_pts = original_coords,
#'     init_pts = trans_coords
#'   )  
#'  
#'  
#'  plot(img)
#'  plot(backtrans_img) # Identical, but with lower resolution
#'  
#'  plot(trans_img) # Transformed image
#'   
#' 
q2q_trans <- function(img,
                      init_pts = NULL,
                      dest_pts = NULL,
                      retain_size = FALSE,
                      coordinates = c("absolute","relative"), 
                      interpolation = c("nearest", "linear", "cubic"),
                      boundary = c("dirichlet", "neumann", "periodic")){
  
  coordinates <- match.arg(coordinates)
  interpolation <- match.arg(interpolation)
  boundary <- match.arg(boundary)
  
  
  ny <- ncol(img)
  nx <- nrow(img)
  
  if(is.null(init_pts)){
    init_pts <- list(
      list(1,1),
      list(1, ny),
      list(nx,ny),
      list(nx, 1)
    )
  }
  
  if(is.null(dest_pts)){
    dest_pts <- list(
      list(1,1),
      list(1, ny),
      list(nx,ny),
      list(nx, 1)
    )
  }
  
  if(!retain_size){
    xmax <- max(do.call("c",purrr::map(dest_pts, 1)))
    ymax <- max(do.call("c",purrr::map(dest_pts, 2)))
    
    if(xmax > nx || ymax > ny){
      
      grow_x <- pmax(xmax - nx, 0)
      grow_y <- pmax(ymax - ny, 0)
      nspec <- dim(img)[4]
      
      x_edge <- as.cimg(array(rep(0, grow_x * ny * nspec), 
                              dim = c(grow_x, ny, 1, nspec)))
      y_edge <- as.cimg(array(rep(0, (nx + grow_x) * grow_y * nspec), 
                              dim = c((nx + grow_x), grow_y, 1, nspec)))
      
      
      img <- imager::imappend(
        list(
          imager::imappend(list(img, x_edge), "x"),
          y_edge), 
        "y")
      
    } 
  }
  
  out <- imager::imwarp(img, 
                        map = q2q_map(dest_pts, init_pts, z = dim(img)[3] > 1), # Swap q2q_map direction with backward algorithm
                        direction = "backward", 
                        coordinates = coordinates,
                        interpolation = interpolation, 
                        boundary = boundary)
  
  return(out)
}


# Corner detection engine
# detect_corners_engine <- function(img, list_of_list = TRUE, scale = 1){
#   pts <- image.CornerDetectionHarris:::detect_corners(img[,,1,1] * 255, 
#                                                       nx = nrow(img), 
#                                                       ny = ncol(img),
#                                                       Nselect = 4L, 
#                                                       strategy = 2L, 
#                                                       verbose = FALSE)
#   out <- data.frame("x" = as.numeric(pts$x) * scale, 
#                     "y" = as.numeric(pts$y) * scale) %>% 
#     arrange(x,y)
#   
#   if(nrow(out) != 4L){
#     stop("Failed to detect 4 points")
#   }
#   
#   if(list_of_list){
#     out <- lapply(seq.int(4), function(i){
#       as.list(out[i,])
#     })
#     return(out)
#   } else {
#     return(out)
#   }
# }

# Convert point list into a data.frame
pt_list2df<- function(x){
  data.frame("x" = do.call("c",map(x, 1)),
             "y" = do.call("c",map(x, 2)))
}

# Split image into different regions and returns the largest patch
split_max <- function(img){
  l <- split_connected(img)
  max_index <- lapply(l, function(x) {
    sum(x)
  }) %>% which.max()
  
  l[[max_index]]
}


# Clean up photo before corner detection
detect_corners <- function(img, qc_plot = FALSE, list_of_list = TRUE, 
                           scale = 5, lambda = 0.001, sat = 0.1, adjust = 1){
  # img2 <- color_index(img, index = "NG", plot = FALSE)$NG %>% 
  #   na_replace(0) %>% 
  #   threshold2(thr = 0.3, thr.exact = TRUE) %>%
  #   imager::bucketfill(x = 1, y = 1, z = 1, color = 1) %>% 
  #   suppressWarnings() %>% 
  #   imager::medianblur(n = 10) 
  img2 <- color_index(img, index = "NG", plot = FALSE)$NG %>% 
    na_replace(0) %>% 
    imagerExtra::SPE(lamda = lambda, s = sat, range = c(0,1)) %>% 
    threshold2(adjust = adjust) %>%
    thin(scale) %>% 
    medianblur(10) %>% 
    clean(10/scale^2) %>% 
    invert() %>% 
    fill(30/(scale * 2)) %>% 
    invert()
  
  
  # split_list <- img2 %>% 
  #   invert() %>% 
  #   threshold2(thr = 0, thr.exact = TRUE) %>% 
  #   split_connected()
  spl_max <- img2 %>% 
    threshold2() %>% 
    split_max()

  
  spl_max2 <- spl_max %>% 
    imager::bucketfill(x = 1, y = 1, z = 1, color = 1) %>% 
    suppressWarnings() %>% 
    invert() %>% 
    split_max()
  
  
  out <- detect_corners_engine(spl_max2, 
                               list_of_list = list_of_list, 
                               scale = scale)
  
  if(qc_plot){
    plot(img)
    if(!is.data.frame(out)){
      o <- pt_list2df(out)
    } else {
      o <- out
    }
    o %$% points(x,y, col = c("green", "blue","blue","blue"), pch = 19, cex = 3)
  }
  return(out)
}

# Perform q2q transformation on raw image, using corner detection 
# Can supply init_pts out side of function to reduce computation time
reproject_grid <- function(img, init_pts = NULL, dest_size = NULL, qc_plot = FALSE){
  if(is.null(init_pts)){
    init_pts <- detect_corners(img, qc_plot = qc_plot)
  }
  
  img_trans <- q2q_trans(
    img = img, 
    init_pts = init_pts,
    retain_size = FALSE
  ) 
  
  if(!is.null(dest_size)){
    img_trans <- img_trans %>% 
      imager::resize(size_x = dest_size,
                     size_y = dest_size, 
                     interpolation = 1)
    
    
  } 
  return(img_trans)
}

# Lighting correction for single source of light
# im_lighting_correction <- function(img){
#   d <- as.data.frame(img)
#   m <- sample_n(d,1e4) %>% lm(value ~ x*y,data=.) 
#   img-predict(m,d)
# }

# # Detect luster dust and apply shadow correction. Returns color inverted image
# detect_luster <- function(img, shadow_weight = 0.5){
#   luster_mask <- color_index(img, 
#                              index = c("HUE"), 
#                              plot = FALSE)[[1]] %>% 
#     imager::renorm(max = 1) %>% 
#     na_replace(0)
#   wt_mask <- color_index(img, index = "CI", plot = FALSE)[[1]] %>% 
#     renorm(min = 0.5, max = 1) %>% 
#     na_replace(1) # NAs classified as no shadows
#   wt_mask <- (shadow_weight / wt_mask) 
#   luster_mask_c <- imager::renorm(luster_mask + wt_mask, min = 0, max = 1)
#   invert(luster_mask_c)
# }

# Replace NA values with val
na_replace <- function(img, val){
  img[is.na(img)] <- val
  img
}

# Load image with JPEG 
fast_load_image <- function(path, transform = TRUE){
  bmp <- jpeg::readJPEG(path)
  if(is.na(dim(bmp)[3])){
    dim(bmp)[3] <- 1
  }
  if(transform){
    bmp <- bmp %>% aperm(c(2, 1, 3))
  }
  dim(bmp) <- c(dim(bmp)[1:2], 1, dim(bmp)[3])
  class(bmp) <- c("cimg", "imager_array", "numeric")
  bmp 
}

# Detect if a video has significant shift in pixel values
# detect_jostle <- function(vid, res = 10000, delta_thresh = 0.15, prop_px = 0.1){
#   if(!is.null(res) || !is.na(res)){
#     if(prod(dim(vid)[1:2]) > res){
#       side_len <- floor(sqrt(10000))
#       vid <- imager::resize(vid, size_x = side_len, size_y = side_len)  
#     }
#   }
#   
#   if(imager::spectrum(vid) > 1){
#     vid <- grayscale(vid)
#   }
#   
#   del <- get_gradient(vid, axes = "z", scheme = -1)[[1]] %>% 
#     imsplit(axis = "z")
#   
#   indices <- lapply(del, function(x){
#     mean(abs(c(x)) > delta_thresh) > prop_px
#   }) %>% 
#     do.call("c",.) %>% 
#     unname() %>% 
#     which()
#   return(indices)
# }


# Make video form a stack of images in a directory
make_video <- function(src_dir, file, fps = 10, extn = ".mp4"){
  f <- list.files(src_dir, full.names = TRUE)
  f2 <- paste(src_dir, gsub(".*_rank","image-",f), sep = "/")
  
  invisible(file.rename(
    f, 
    f2
  ))
  
  imager::make.video(fps = fps, 
                     pattern = "image-%d.jpg", 
                     dname = abs_path(src_dir), 
                     fname = abs_path(paste0(file,extn)), 
                     verbose = TRUE)
  
  invisible(file.rename(
    f2, 
    f
  ))
}

# Index the reproject points from specified index
choose_pts <- function(pts_list, index){
  i <- lapply(pts_list, function(x, index){
    index %in% attr(x, "indices")
  }, index = index) %>% 
    do.call("c",.) %>% 
    which()
  return(pts_list[[i]])
}

# Crop and reproject raw image as processed image
crop_raw_img <- function(
    indices = do.call("c", lapply(pts_list, function(x){attr(x, "indices")})),
    .pts_list = get("pts_list"),
    .files_full_name = get("files_full_name"), 
    .dest_dir = get("dest_dir"), 
    cores = parallel::detectCores(logical = FALSE) - 2
){
  
  start_time <- Sys.time()
  
  stopifnot(!any(duplicated(indices)))
  stopifnot(
    length(.files_full_name) == length(indices)
  )
  
  if(is.null(indices)){
    indices <- seq_along(.files_full_name)
  }
  
  new_file_names <- paste(.dest_dir, paste0("processed_", basename(.files_full_name)), sep = "/")
  file_exits <- gsub("_rank.*",".jpg",list.files(.dest_dir, full.names = TRUE))
  
  file_conflicts <- sum(new_file_names %in% file_exits)
  if(file_conflicts > 0){
    stop(sprintf("%s files already exist.", file_conflicts))
  }
  
  pb_par_lapply(indices, function(i, .pts_list, .files_full_name, .dest_dir){
    
    img2 <- reproject_grid(fast_load_image(.files_full_name[i], transform = FALSE), 
                           init_pts = choose_pts(.pts_list, i), 
                           dest_size = 1000, 
                           qc_plot = FALSE) %>% 
      mirror("x")
    dim(img2) <- dim(img2)[-3]
    jpeg::writeJPEG(img2,
                    paste(.dest_dir, paste0("processed_", 
                                            basename(.files_full_name)[i]), 
                          sep = "/"), 
                    quality = 1)
    
  }, 
  .pts_list = .pts_list,
  .files_full_name = .files_full_name, 
  .dest_dir = .dest_dir,
  cores = cores, 
  inorder = FALSE, 
  export_fun_only = TRUE) %>% 
    invisible()
  
  #message("\nInitializing parallel workers. . .")
  #cl <- makeCluster(cores, outfile = "")
  #registerDoSNOW(cl)
  
  # foreach(i = indices, 
  #         .export = ls(globalenv()),
  #         .combine = c, 
  #         .verbose = FALSE, 
  #         .final = invisible,
  #         .inorder = FALSE,
  #         .options.snow = list(
  #           progress = function(n) {
  #             cat(sprintf("\r Processing %d out of %d", n, length(indices)))
  #           }
  #         ),
  #         .packages = c("tidyverse", "herbivar")) %dopar% {
  #           
  #           img2 <- reproject_grid(fast_load_image(.files_full_name[i], transform = FALSE), 
  #                                  init_pts = choose_pts(pts_list, i), 
  #                                  dest_size = 1000, 
  #                                  qc_plot = FALSE) %>% 
  #             mirror("x")
  #           dim(img2) <- dim(img2)[-3]
  #           jpeg::writeJPEG(img2,
  #                           paste(.dest_dir, paste0("processed_", .files[i]), sep = "/"), 
  #                           quality = 1)
  #           #cat("\r Cropping", i, "of", length(.files))
  #         } 
  attach_order_file_name(.dest_dir)
  #stopCluster(cl)
  
  hms_runtime(as.numeric(Sys.time() - start_time, units = "secs"))
  message("\nDone!")
}





# # Detect caterpillar from a colored image using thresholding 
# detect_cat <- function(img, w = c(5,1), thr = "kmeans", adjust = 1.3, clean = 3,
#                        lambda = 0.1, sat = 0.1,  cores = 1){
#   start_time <- Sys.time()
#   if(TRUE){
#     out <- img %>% 
#       color_index(index = c("BI","NG"), plot = FALSE) %>% 
#       iml_prod(w) %>% 
#       renorm(min = 0, max = 1) %>% 
#       imsplit(axis = "z") %>% 
#       pb_par_lapply(function(x, lambda, sat, thr, adjust, clean){
#         imagerExtra::SPE(x, lambda, s = sat, range = c(0,1)) %>% 
#           threshold2(thr = thr, adjust = adjust) %>% 
#           clean(clean) %>% 
#           split_max()
#       }, 
#       clean = clean,
#       adjust = adjust,
#       thr = thr,
#       lambda = lambda, 
#       sat = sat,
#       cores = cores, 
#       inorder = TRUE) %>% 
#       as.imlist() %>% 
#       imappend("z")
#   }
#   hms_runtime(as.numeric(Sys.time() - start_time, units = "secs"))
#   return(out)
# }
# 
# 
# # Take two matrices and do element-wise multiplication with some weight
# iml_prod <- function(x, w = c(1,1)){
#   for(i in seq_along(w)){
#     if(i == 1){
#       if(w[i] > 0){
#         z <- x[[i]] * w[i]
#       } else {
#         z <- 1 / x[[i]] * w[i]
#       }
#       
#     } else {
#       if(w[i] > 0){
#         z <- z + x[[i]] * w[i]
#       } else {
#         z <- z + 1 / x[[i]] * w[i]
#         }
#       
#     }
#   }
#   z
# }

# Load a stack of images in a directory as video using parallelization. 
load_video <- function(file_paths, thin.val = 5, cores = 6){
  start_time <- Sys.time()
  
  out <- pb_par_lapply(file_paths, function(x, thin.val){
    herbivar::thin(fast_load_image(x, transform = FALSE), thin.val)
  }, cores = cores, 
  thin.val = thin.val,
  loop_text = "Loading image", 
  inorder = FALSE) %>% 
    imappend(axis = "z")
  
  hms_runtime(as.numeric(Sys.time() - start_time, units = "secs"))
  return(out)
}

# Reorder anchors for projection (top left, top right, bottom right, bottom left)
anchor_reorder <- function(x){
  is_ll <- (!inherits(x, "data.frame")) & is.list(x)
  if(is_ll){
    x <- pt_list2df(x)
  }
  o <- c(NA, NA, NA, NA)
  rs <- rowSums(x)
  min_o <- which.min(rs)
  max_o <- which.max(rs)
  o[1] <- min_o
  o[3] <- max_o
  
  o4_candidate <- which(x$x > x$y)
  o4_candidate <- o4_candidate[!o4_candidate %in% c(min_o, max_o)]
  o[4] <- o4_candidate
  o[2] <- c(1,2,3,4)[!c(1,2,3,4) %in% o]
  x <- x[o,]
  
  if(is_ll){
    x <- lapply(seq.int(4), function(i){
      as.list(x[i,])
    })
  }
  return(x)
}


# Reformat cimg as bitmap format for interface with jpeg package
as.bmp <- function(x){
  x <- as.cimg(x)
  dim(x) <- dim(x)[-3]
  x
}

# Add a point to an image by writing into the data array.
add_point <- function(img, x, y, r = ceiling(min(dim(img)[1:2])/80), color = "red"){
  if(imager::is.pixset(img) || imager::spectrum(img) < 3){
    img <- as.cimg_color(img)
  }
  
  y <- as.numeric(round(y))
  x <- as.numeric(round(x))
  
  color <- col2rgb(color)/255
  
  xmax <- as.numeric(pmin(x + r, nrow(img)))
  xmin <- as.numeric(pmax(x - r, 1))
  ymax <- as.numeric(pmin(y + r, ncol(img)))
  ymin <- as.numeric(pmax(y - r, 1))
  
  ptd <- expand.grid("xi" = seq(xmin, xmax, by = 1), "yi" = seq(ymin, ymax, by = 1))
  v <- sqrt((ptd$xi - x)^2 + (ptd$yi - y)^2) < r
  ptd <- ptd[v,]
  v <- (ptd$xi) + nrow(img) * (ptd$yi - 1) 
  
  img[,,,1][v] <- color[1,1]
  img[,,,2][v] <- color[2,1]
  img[,,,3][v] <- color[3,1]
  
  return(img)
}


# Check if object mask is out of frame (if polygon edge is within `tolerance` number of pixels of the boarders). Returns a boolean.  
out_of_frame <- function(poly, tolerance = 20, dim_xy = c(1000, 1000)){
  if(is.null(poly)){
    return(NA)
  }
  any(
    poly[,1] < (1 + tolerance) | 
      poly[,1] > (dim_xy[1] - tolerance)
  ) | 
    any(
      poly[,2] < (1 + tolerance) | 
        poly[,2] > (dim_xy[2] - tolerance)
    )
}

# Flip x and y axes
flip_xy <- function(img){
  f <- switch(class(img)[1], cimg = as.cimg, pixset = as.pixset, 
              array = as.array)
  return(f(aperm(img, c(2, 1, 3, 4))))
}

# Compute the HUE index of image with RGB
HUE <- function(img){
  color_index(img, "HUE", plot = FALSE)[[1]]
}

