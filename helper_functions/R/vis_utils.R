# Display binary mask using different methods. Results are very similar, but polygon looks tighter. Mask grows the boundaries a little bit. 
plot_mask <- function(data_dict, frame, method = c("polygon", "mask"), 
                      col = "white", alpha = 1, ...){
  method <- match.arg(method)
  frame <- assert_frame(data_dict, frame)
  stopifnot(is.data_dict(data_dict))
  
  if(method == "mask"){
    pred_mask <- data_dict[frame] %>% get_mask(1) %>% .[[1]]
    
    if(mask_area(pred_mask) == 0){
      message(sprintf("No %s detected.", "mask"))
    }
    pred_mask %>% plot(...)
  } else {
    dims <- get_dim(data_dict[frame])
    plot.cimg(
      imfill(dims[1], dims[2], z = 1, val = 0),
      rescale = FALSE, 
      ...
    )
    plot_polygon(data_dict, frame, col = adjustcolor(col, alpha.f = 1), alpha = alpha, ...)
  }
}

# Display keypoints 
plot_keypoint <- function(data_dict, frame, 
                          kp_name = c("head", "middle", "tail"), 
                          col = c("green", "yellow", "red"), pch = 19, ...){
  frame <- assert_frame(data_dict, frame)
  stopifnot(is.data_dict(data_dict))
  
  d <- get_keypoints(data_dict[frame])
  x_names <- paste0(kp_name, "_x")
  y_names <- paste0(kp_name, "_y")
  
  x <- d[,x_names]
  y <- d[,y_names]
  
  failed_pred <- is.na(x) | is.na(y)
  
  if(any(failed_pred)){
    message(sprintf(
      "No keypoint detection for %s",
      paste0(add_quote(kp_name[failed_pred]), collapse = ", ")
    ))
  }
  points(x = x, y = y, col = col, pch = pch, ...)
}

# Display polygon from data_dict object. Looks better, but is slightly different to the result of polygon2mask()
plot_polygon <- function(data_dict, frame, col = "blue", alpha = 0.5, ...){
  frame <- assert_frame(data_dict, frame)
  stopifnot(is.data_dict(data_dict))
  
  col <- adjustcolor(col = col, alpha.f = alpha)
  pred_polygon <- data_dict[frame] %>% get_polygon() %>% .[[1]]
  if(is.null(pred_polygon)){
    message(sprintf("No %s detected.", "mask"))
  } else {
    pred_polygon %>% polygon(col = col, ...)
  }
}

# S3 generic for data_dict object, calling bbox, polygon, and keypoint plot methods
plot.data_dict <- function(x, frame, add = FALSE, 
                           kp_col = c("green", "yellow", "red"), 
                           mask_col = "blue",
                           kp_name = c("head", "middle", "tail"),
                           bbox = FALSE,
                           boundary_col = "red",
                           ...){
  frame <- assert_frame(x, frame)
  kp_name <- match.arg(kp_name, several.ok = TRUE)
  if(add){
    plot_polygon(x, frame, col = mask_col, ...)
  } else {
    plot_mask(x, frame, col = mask_col, ...)
  }
  if(length(kp_col) != length(kp_name)){
    kp_col <- kp_col[match(kp_name, table = c("head", "middle", "tail"))] 
  }
  
  plot_keypoint(x, frame, col = kp_col, kp_name = kp_name, ...)
  
  if(bbox){
    plot_bbox(x, frame, col = boundary_col, ...)
  }
}


registerS3method(genname = "plot", 
                 class = "data_dict", 
                 method = plot.data_dict)


# Draw bounding box
draw_bbox <- function(coord, fill = "#00000000", col = "red", ...){
  rect(xleft = coord[1,1], xright = coord[2,1], ybottom = coord[1,2], ytop = coord[2,2], 
       col = fill, border = col, ...)
}

# Draw bounding box
plot_bbox <- function(data_dict, frame, fill = "#00000000", col = "red", ...){
  frame <- assert_frame(data_dict, frame)
  stopifnot(is.data_dict(data_dict))
  
  coord <- data_dict[frame] %>% get_bbox() %>% .[[1]] %>% parse_bbox_vec()
  
  if(is.null(coord)){
    message(sprintf("No %s detected.", "bounding box"))
  } else{
    draw_bbox(coord, fill = fill, col = col, ...) 
  }
}


# Overlays grid cell labels on the treatment raster image
plot_image_guide <- function(img, transform = FALSE, col = "red", cex = 0.8, main = NULL, mar = NULL, axes = TRUE, ...){
  
  if(!is.null(mar)){
    temp_mar <- par()$mar
    par(mar = mar)
  }
  
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
  
  if(transform){
    img <- flip_xy(img)
    x_temp <- d$x
    y_temp <- d$y
    d$y <- x_temp 
    d$x <- y_temp
  }
  plot.cimg(img, main = main, axes = axes, ...)
  text(x = d$x, y = d$y, col = col, label = d$label, cex = cex)
  
  if(!is.null(mar)){
    par(mar = temp_mar)
  }
}

# Store forward_plot() in an object, then run the object multiple times to move forward. 
forward_plot <- function(x,init_frame, mask_col = "white", ...){
  function(){
    plot(x, init_frame, 
         mask_col = mask_col, add = FALSE, 
         main = sprintf("Rank: %s", init_frame), 
         ...)
    init_frame <<- init_frame + 1
  }
}

# Visualize sequential animal track
plot_track <- function(data, x, y, type = c("track", "density")){
  .expose_columns_interal()
  type <- match.arg(type, several.ok = TRUE)
  
  g <- ggplot(data, aes(x = x, y = y))
  
  if("density" %in% type){
    g <- g + geom_density_2d_filled()
  }
  
  if("track" %in% type){
    g <- g + geom_path()
  }
  return(g)
}






