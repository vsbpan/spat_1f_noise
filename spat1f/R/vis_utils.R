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
plot_track <- function(data, x, y, type = c("track", "density"), max_xy = c(1000, 1000)){
  .expose_columns_interal()
  type <- match.arg(type, several.ok = TRUE)
  
  g <- ggplot(data, aes(x = x, y = y)) +
    coord_cartesian(xlim = c(0, max_xy[1]), ylim = c(0, max_xy[2]))
  
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
  predictor_frame <-predictor_frame[,!names(predictor_frame) %in% yname, drop = FALSE]
  
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

# Handly function to back transform log scale variables
unscalelog <- function(logx){
  function(z) {
    exp(mean(logx, na.rm = TRUE) + sd(logx, na.rm = TRUE) * z)
  }
}


# Plot track and treatment spec together
plot_track_overlay.default <- function(events = NULL, repID = NULL, 
                               ref_data = get("ref_data", envir = globalenv()), 
                               score_thresh = 0.7,
                               plot_elements = c("ud","track"),
                               colored_track = c("none", "states","time"),
                               max_xy = c(1000, 1000),
                               trt_spec = "auto",
                               plot_title = sprintf("repID: %s", repID_clean(repID)), 
                               ...
                               ){
  colored_track <- match.arg(colored_track)
  plot_elements <- match.arg(plot_elements, several.ok = TRUE)
  
  if(is.null(repID)){
    repID <- na.omit(unique(events$repID))
  }
  
  stopifnot(length(repID) == 1)
  
  if(isTRUE(trt_spec == "auto")){
    trt_spec <- fetch_trt_spec(repID) 
  }
  if(is.null(trt_spec)){
    plot_d <- expand.grid(
      "dim1" = 1:12, 
      "dim2"= 1:12,
      "val" = 1
    )
  } else {
    plot_d <- trt_spec %>% 
      flip_xy() %>% 
      as.matrix() %>% 
      melt()
  }
  
  
  if(is.null(events)){
    events <- fetch_events(repID) %>%
      clean_events(insert_gaps = TRUE, 
                   score_thresh = score_thresh, 
                   ref_data = ref_data, 
                   keep_sus = FALSE)
  }
  
  if(colored_track == "states"){
    hmm_fit <- fit_HMM(as.moveData(move_seq(events$head_x, events$head_y)))
  }
  
  dim_x <- ncol(trt_spec) # swap col and row bc of flip_xy()
  dim_y <- nrow(trt_spec)
  
  events <- events %>% 
    mutate(
      head_x = head_x / max_xy[1] * dim_x + 0.5,
      head_y = head_y / max_xy[2] * dim_y + 0.5
    )
  
  g <- plot_d %>% 
    ggplot(aes(x = dim1, y = dim2)) + 
    geom_tile(
      fill = ifelse(plot_d$val == 1, "white", "grey"),
      alpha = 0.5)
  
  if("ud" %in% plot_elements){
    g <- g + 
      geom_density_2d_filled(
        data = events,
        aes(x = head_x, y = head_y), 
        inherit.aes = FALSE,
        alpha = 0.5, 
        show.legend = FALSE
      ) + 
      scale_fill_viridis_d(option = "mako")
  }
  
  if("track" %in% plot_elements){
    if(colored_track == "time"){
      g <- g +
        geom_path(
          data = events,
          aes(x = head_x, y = head_y, color = round(time / 360))
        ) + 
        scale_color_gradient2(midpoint = 600) + 
        labs(color = "time steps")
    } else if(colored_track == "states"){
      states <- as.character(viterbi(hmm_fit))
      g <- g +
        geom_path(
          data = events[-nrow(events),],
          aes(x = head_x, y = head_y, color = states, group = 1)
        ) + 
        labs(color = "states") + 
        scale_color_manual(values = .getPalette(unique_len(states)))
      
    } else {
      g <- g +
        geom_path(
          data = events,
          aes(x = head_x, y = head_y), color = "black"
        )
    }
    
      
  }
  g <- g +
    theme_void() + 
    geom_point(data = data.frame(), aes(x = 0.5, y = 0.5), alpha = 0) + # Force the edges to align
    geom_point(data = data.frame(), aes(x = dim_x + 0.5, y = dim_y + 0.5), alpha = 0) + # Force the edges to align
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) + 
    labs(title = plot_title) 
  
  return(g)
}


plot_track_overlay <- function(events, 
                               plot_elements = c("ud","track"),
                               colored_track = c("none", "states","time"),
                               max_xy = NULL,
                               trt_spec = NULL,
                               ...){
  UseMethod("plot_track_overlay")
}

registerS3method(genname = "plot_track_overlay", 
                 class = "default", 
                 method = plot_track_overlay.default)



plot_track_overlay.sim_steps <- function(events, 
                                         plot_elements = c("ud","track"),
                                         colored_track = c("none", "states","time"),
                                         max_xy = NULL,
                                         trt_spec = NULL,
                                         ...){
  
  if(is.null(max_xy)){
    max_xy <- attr(events, "max_xy")
  }
  
  if(is.null(trt_spec)){
    trt_spec <- attr(events, "ref_grid")
  }
  
  
  plot_track_overlay.default(
    events = events %>% 
      rename(head_x = x, 
             head_y = y),
    repID = "foo",
    plot_elements = plot_elements, 
    colored_track = colored_track,
    max_xy = max_xy,
    trt_spec = flip_xy(trt_spec), 
    plot_title = sprintf("Simulated tracks for %s steps", nrow(events) - 1),
    ...
  )
}


registerS3method(genname = "plot_track_overlay", 
                 class = "sim_steps", 
                 method = plot_track_overlay.sim_steps)



# Make log-log histogram
loghist <- function(x, nclass = 50, log.p = FALSE, log.x = TRUE, geom = c("line", "col"), draw_dist = NULL, ...){
  
  if(!log.p && length(geom) == 2){
    geom <- "col"
  } else {
    geom <- match.arg(geom)
  }
  
  if(log.x){
    p <- hist(log(x), plot = FALSE, nclass = nclass, ...)
    d <- data.frame(
      "x" = exp(p$mids),
      "p" = p$density
    )
  } else {
    p <- hist(x, plot = FALSE, nclass = nclass, ...)
    d <- data.frame(
      "x" = p$mids,
      "p" = p$density
    )
  }
  
  if(log.p){
    d <- d %>% 
      filter(
        p > 0
      )
  }
  
  g <- d %>% 
    ggplot(aes(x = x, y = p)) + 
    theme_bw(base_size = 15) + 
    labs(x = "x", y = "P(x)")
  
  if(log.x){
    g <- g + scale_x_continuous(trans = "log10", labels = fancy_scientificb)
  }
  
  if(log.p){
    g <- g + scale_y_continuous(trans = "log10")
  }
  
  if(geom == "line"){
    g <- g + geom_line()
  }
  
  if(geom == "col"){
    g <- g + geom_col()
  }
  
  if(!is.null(draw_dist)){
    if(inherits(draw_dist, "amt_distr")){
      draw_dist <- list(draw_dist) 
    }
    
    n <- length(draw_dist)
    den_data <- vector(mode = "list", length = n)
    
    for (i in seq_len(n)){
      dist_name <- draw_dist[[i]]$name
      
      den <- do.call(
        paste0("d",dist_name),
        c(
          list(d$x), 
          draw_dist[[i]]$params
        )
      )
      
      den[!is.finite(den)] <- NA_real_
      
      if(log.x){
        p <- den * d$x 
      } else {
        p <- den
      }
      
      den_data[[i]] <- data.frame(
        "x" = d$x, 
        "p" = p,
        "dist" = dist_name
      )
    }
    
    g <- g + 
      geom_line(
        data = do.call("rbind", den_data), 
        aes(color = dist),
        size = 1
      ) + 
      theme(legend.position = "top")
    
  }
  
  return(g)
}


fancy_linear <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = FALSE)
  # return this as an expression
  parse(text=l)
}

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # e+00 becomes 1
  l <- gsub("e\\+00", "", l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove prefactor 1
  l <- gsub("'1'e", "10^", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # remove plus
  l <- gsub("\\+", "", l)
  # return this as an expression
  parse(text=l)
}


fancy_scientificb <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove prefactor 1
  l <- gsub("'1'e", "10^", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # remove plus
  l <- gsub("\\+", "", l)
  # return this as an expression
  parse(text=l)
}

scientific_10_exp_labels <- scales::trans_format("log10", scales::math_format(10^.x) )
scientific_10_exp_breaks <- scales::trans_format("log10", function(x) 10^x )

scale_xy_log <- function(name = waiver(), 
                                axis = "xy",
                                trans = "log10", 
                                label = fancy_scientificb,  
                                ...){
  if(grepl("x", axis)){
    res <- scale_x_continuous(
      name, 
      trans = trans, 
      label = label,
      ...
    )
  }
  
  if(grepl("y", axis)){
    res <- scale_y_continuous(
      name, 
      trans = trans, 
      label = label,
      ...
    )
  }
  
  return(res)
}

scale_x_ta <- function(name = waiver(), ...){
  scale_x_continuous(
    name,
    breaks = c(-pi, -pi / 2, 0, pi / 2, pi), 
    labels = c(expression(-pi), expression(-pi/2), 0, expression(pi/2), expression(pi)), 
    ...)
}

