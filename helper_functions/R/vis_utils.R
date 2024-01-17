# Display binary mask using differnet methods
plot_mask <- function(x, frame, method = c("polygon", "mask"), col = "white", alpha = 1, ...){
  method <- match.arg(method)
  
  if(method == "mask"){
    x[frame] %>% get_mask(1) %>% .[[1]] %>% plot(...)
  } else {
    dims <- get_dim(x[frame])
    plot.cimg(
      imfill(dims[1], dims[2], z = 1, val = 0),
      rescale = FALSE
    )
    plot_polygon(x, frame, col = adjustcolor(col, alpha.f = 1), alpha = alpha, ...)
  }
  
}

# Display keypoints 
plot_keypoint <- function(data_dict, frame, x = "head_x", y = "head_y", col = "green", pch = 19, ...){
  d <- get_keypoints(data_dict[frame])
  points(x = d[,x], y = d[,y], col = col, pch = pch, ...)
}

# Display polygon from data_dict object. Looks better, but is slightly different to the result of polygon2mask()
plot_polygon <- function(x, frame, col = "blue", alpha = 0.5, ...){
  col <- adjustcolor(col = col, alpha.f = alpha)
  x[frame] %>% get_polygon() %>% .[[1]] %>% polygon(col = col, ...)
}
