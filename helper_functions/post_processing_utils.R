
move_seq <- function(x,y){
  stopifnot(length(x) == length(y))
  out <- lapply(seq_len(length(x) - 1), function(i){
    delta_x <- x[i] - x[i+1]
    delta_y <- y[i] - y[i+1]
    r <- sqrt((delta_x)^2 + (delta_y)^2)
    theta_abs <- atan2(delta_y, delta_x)
    c("r" = r, "theta_abs" = theta_abs)
  }) %>% bind_vec()
  out$theta_rel <- c(NA, (out$theta_abs - c(out$theta_abs[-1],NA))[-(length(x)-1)])
  out
}



fetch_trt_spec <- function(rep_ID, ref_data, trt_meta_iml){
  sny_id <- ref_data[which(ref_data$rep_id == gsub("rep","",rep_ID)), "syn_id"]
  trt_meta_iml[[which(names(trt_meta_iml) == paste0("syn_id__",sny_id))]]
}



# Insert gaps of NAs at positions where a photo is expected
# A bit buggy
insert_gaps <- function(df, expected_gap = 360){
  time <- file_time(df$frame_id)
  time_grid <- 360 * (seq_len((diff(round(range(time) / expected_gap)) +1)) - 1)
  o <- order(time)
  df <- df[o,]
  g <- ((diff(time[o]) - expected_gap) / expected_gap) %>% 
    round()
  g <- c(g, 0)
  
  
  
  empty_row <- as.data.frame(matrix(rep(NA, ncol(df)), 
                                    nrow = 1, ncol = ncol(df), 
                                    dimnames = list(NULL, names(df))))
  
  out <- lapply(seq_along(g), function(i, g){
    rbind(d2[i,], empty_row[rep(TRUE, g[i]),])
  }, g = g) %>% 
    do.call("rbind",.)
  
  rownames(out) <- NULL
  out[is.na(out$frame_id),"frame_id"] <- paste0("gap", seq_len(sum(g)), "__s", time_grid[is.na(out$frame_id)])
  return(out)
}


