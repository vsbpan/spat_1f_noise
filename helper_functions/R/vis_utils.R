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


# Generate marginal effect data frame for plotting
marginal_effects <- function(model, terms, n = 300, ci = 0.95){
  
  predictor_frame <- insight::get_data(model)
  yname <- insight::find_response(model)
  predictor_frame <-predictor_frame[,!names(predictor_frame) %in% yname]
  
  rand_names <- insight::find_random(model)$random
  var_names <- names(predictor_frame)
  
  fam <- insight::get_family(model)
  if(!inherits(fam,"family")){
    linkinv <-function(x) x
  } else {
    linkinv <- fam$linkinv 
  }
  
  v <- lapply(terms, function(z){
    switch(as.character(grepl("\\[", z) && grepl("\\[", z)), 
           "TRUE" = as.numeric(unlist(strsplit(gsub(".*\\[|\\]","",z),","))),
           "FALSE" = NULL)
  })
  terms <- gsub("\\[.*","",terms)
  
  
  if(length(n) == 1){
    n <- rep(n, length(terms))
  }
  
  new_data <- expand.grid(lapply(seq_along(predictor_frame), function(i, d, n, v){
    x <- d[, i]
    
    if(names(d)[i] %in% terms){
      j <- which(terms %in% names(d)[i])
      
      if(is.null(v[[j]])){
        if(is.numeric(x)){
          x <- seq_interval(x, n[j])
        } else {
          x <- unique(x) 
        }
      } else {
        x <- v[[j]]
      }
      return(x) 
    } else {
      if(is.numeric(x)){
        x <- mean(x)
      } else {
        if(names(d)[i] %in% rand_names){
          x <- "foooooooooooooooo"
        } else {
          x <- unique(x)  
        }
      }
      return(x)
    }
  }, 
  d = as.data.frame(predictor_frame), 
  n = n, 
  v = v))
  
  names(new_data) <- var_names
  
  if(inherits(model, "clm")){
    pred <- suppressWarnings(predict(model, newdata = new_data, se = TRUE))
    
    new_data <- cbind(new_data, 
                      "yhat" = pred$fit, 
                      "se" = pred$se.fit)
    
    new_data <- new_data %>% 
      group_by(all_of(terms)) %>% 
      summarise_all(function(x){
        if(is.numeric(x)){
          mean(x)
        } else {
          NA
        }
      })
    
    new_data <- new_data %>% 
      gather(value = "yhat", key = "cat", paste("yhat",levels(model$y), sep = ".")) %>% 
      mutate(cat = factor(gsub(".*\\.","",cat), levels = levels(model$y), ordered = TRUE))
    names(new_data)[names(new_data) == "cat"] <- yname
    
  } else {
    pred <- suppressWarnings(predict(model, newdata = new_data, se = TRUE, type = "link"))
    
    new_data <- cbind(new_data, 
                      "yhat_link" = pred$fit, 
                      "se" = pred$se.fit)
    new_data$lower_link <- new_data$yhat_link + qnorm((1 - ci)/2) * new_data$se
    new_data$upper_link <- new_data$yhat_link + qnorm((1 - ci)/2, lower.tail = FALSE) * new_data$se
    
    new_data <- new_data %>% 
      mutate(se = se^2) %>% 
      group_by(across(terms)) %>% 
      summarise_all(function(x){
        if(is.numeric(x)){
          mean(x)
        } else {
          NA
        }
      }) %>% 
      mutate(se = sqrt(se)) 
    
    resp_inv <- insight::get_transformation(model)$inverse
    
    new_data$lower <- resp_inv(linkinv(new_data$lower_link))
    new_data$upper <- resp_inv(linkinv(new_data$upper_link))
    new_data$yhat <- resp_inv(linkinv(new_data$yhat_link))
  }
  
  return(new_data)
}

unscalelog <- function(logx){
  function(z) {
    exp(mean(logx, na.rm = TRUE) + sd(logx, na.rm = TRUE) * z)
  }
}



