source("spat1f/init.R")



out2 <- fetch_repID() %>% 
  pb_par_lapply(
    function(x){
      #cat(sprintf("Testing rep %s\n", x))
      tryCatch(
        s <- detect_false_head_movement(fetch_data_dict(x),
                                        fetch_events(x), 
                                        head_thresh = "auto",
                                        cen_thresh = "auto",
                                        iou_thresh = 0.5,
                                        cores = 1),
        error = function(e){
          return("Error")
        }
      )
      
      out <- list("rep_id" = x, "sus_frames" = s)
      return(out)
    }, cores = 8,inorder = FALSE
  )


# saveRDS(out2, "cleaned_data/sus_frames_list.rds")














