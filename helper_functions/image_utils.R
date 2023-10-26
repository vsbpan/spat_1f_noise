
# Overlays grid cell labels on the treatment raster image
plot_image_guide <- function(img, col = "red", cex = 0.8, main = NULL){
  ny <- nrow(img)
  nx <- ncol(img)
  
  d <- expand.grid(
    "y" = seq_len(ny),
    "x" = seq_len(nx)
  ) %>% 
    mutate(
      label = paste0(LETTERS[y],x)
    ) %>% 
    mutate(
      y = y/ny*(ny-1) + 0.5,
      x = x/nx*(nx-1) + 0.5
    )
  
  plot(img, main = main)
  text(x = d$x, y = d$y, col = col, label = d$label, cex = cex)
}

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
detect_corners_engine <- function(img, list_of_list = TRUE){
  pts <- image.CornerDetectionHarris:::detect_corners(img[,,1,1] * 255, 
                                                      nx = nrow(img), 
                                                      ny = ncol(img),
                                                      Nselect = 4L, 
                                                      strategy = 2L, 
                                                      verbose = FALSE)
  out <- data.frame("x" = pts$x, "y" = pts$y) %>% 
    arrange(x,y)
  
  if(nrow(out) != 4L){
    stop("Failed to detect 4 points")
  }
  
  if(list_of_list){
    out <- lapply(seq.int(4), function(i){
      as.list(out[i,])
    })
    return(out)
  } else {
    return(out)
  }
}

# Clean up photo before corner detection
detect_corners <- function(img, qc_plot = FALSE, list_of_list = TRUE){
  img2 <- color_index(img, index = "NG", plot = FALSE)$NG %>% 
    na_replace(0) %>% 
    threshold2(thr = 0.3, thr.exact = TRUE) %>%
    imager::bucketfill(x = 1, y = 1, z = 1, color = 1) %>% 
    suppressWarnings() %>% 
    imager::medianblur(n = 10) 
  
  
  split_list <- img2 %>% 
    invert() %>% 
    threshold2(thr = 0, thr.exact = TRUE) %>% 
    split_connected()
    
  max_index <- lapply(split_list, function(x) {
    sum(x)
  }) %>% which.max()
    
  
  out <- detect_corners_engine(split_list[[max_index]], list_of_list = list_of_list)
  
  if(qc_plot){
    plot(img)
    if(!is.data.frame(out)){
      o <- do.call("rbind",out) %>% as.data.frame()
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
im_lighting_correction <- function(img){
  d <- as.data.frame(img)
  m <- sample_n(d,1e4) %>% lm(value ~ x*y,data=.) 
  img-predict(m,d)
}

# Detect lester dust and apply shadow correction. Returns color inverted image
detect_lester <- function(img, shadow_weight = 0.5){
  lester_mask <- color_index(img, 
                             index = c("HUE"), 
                             plot = FALSE)[[1]] %>% 
    imager::renorm(max = 1) %>% 
    na_replace(0)
  wt_mask <- color_index(img, index = "CI", plot = FALSE)[[1]] %>% 
    renorm(min = 0.5, max = 1) %>% 
    na_replace(1) # NAs classified as no shadows
  wt_mask <- (shadow_weight / wt_mask) 
  lester_mask_c <- imager::renorm(lester_mask + wt_mask, min = 0, max = 1)
  invert(lester_mask_c)
}

# Replace NA values with val
na_replace <- function(img, val){
  img[is.na(img)] <- val
  img
}

fast_load_image <- function(path){
  bmp <- jpeg::readJPEG(path) %>% aperm(c(2, 1, 3))
  dim(bmp) <- c(dim(bmp)[1:2], 1, dim(bmp)[3])
  class(bmp) <- c("cimg", "imager_array", "numeric")
  bmp 
}


detect_jostle <- function(vid, res = 10000, delta_thresh = 0.15, prop_px = 0.1){
  if(!is.null(res) || !is.na(res)){
    if(prod(dim(vid)[1:2]) > res)
      side_len <- floor(sqrt(10000))
    vid <- imager::resize(vid, size_x = side_len, size_y = side_len) 
  }
  
  if(imager::spectrum(vid) > 1){
    vid <- grayscale(vid)
  }
  
  del <- get_gradient(vid, axes = "z", scheme = -1)[[1]] %>% 
    imsplit(axis = "z")
  
  indices <- lapply(del, function(x){
    mean(abs(c(x)) > delta_thresh) > prop_px
  }) %>% 
    do.call("c",.) %>% 
    unname() %>% 
    which()
  return(indices)
}

### To Do ### 
## Add corner detection for artificial diet raster




