source("spat1f/init_analysis.R")

plot_instance_overlay <- function(repID, rank, 
                                  mar = c(0.1,0.1,0.1,0.1), main = NULL, 
                                  trt_spec = resize(as.cimg_color(
                                    flip_xy(fetch_trt_spec(repID, quiet = TRUE))
                                    ), 
                                    size_x = 1000, 
                                    size_y = 1000
                                  ), 
                                  opacity = 0.08, ...){
  if(!is.null(mar)){
    temp_mar <- par()$mar
    par(mar = mar)
  }
  fetch_image(repID, rank = rank, type = "processed") %>% 
    imager::imdraw(trt_spec, 
                   opacity = opacity) %>% 
    plot(axes = FALSE, main = main)
  z <- fetch_data_dict(repID)[rank]
  time <- hm_format(get_file_meta(z)[,"time"])
  
  plot(z, add = TRUE, ...)
  text(1000 * 0.92, 1000 * 0.95, time, col = "white")
  
  if(!is.null(mar)){
    par(mar = temp_mar)
  }
}





repID <- 81
trt_spec <- resize(as.cimg_color(
  flip_xy(fetch_trt_spec(repID, quiet = TRUE))
), 
size_x = 1000, 
size_y = 1000
)



pb_par_lapply(
  seq_len(
    time2rank(fetch_data_dict(repID), fetch_cutoff(repID))
  ), 
  function(rankID, repID, trt_spec){
    fn <- sprintf("rep%s_rank%s", repID, rankID)
    jpeg(paste0("C:/R_projects/spat_1f_noise/invisible/test_video_overlay/",fn,".jpg"), width = 750, height = 750, units = "px", res = 250)
    plot_instance_overlay(repID, rankID, cex = 0.8, trt_spec = trt_spec)
    dev.off()
    return(invisible(NULL))
  }, 
  repID = repID, 
  trt_spec = trt_spec, 
  cores = 8, 
  inorder = FALSE
)

dest_dir <- "C:/R_projects/spat_1f_noise/invisible/test_video_overlay"
make_video(src_dir = dest_dir, 
           file = sprintf("rep%s_overlay", repID), fps = 10)



plot_track_overlay(repID = 45)
plot_track_overlay(repID = 46)
plot_track_overlay(repID = 50)
plot_track_overlay(repID = 36)

plot_track_overlay(repID = 144)
plot_track_overlay(repID = 141)
plot_track_overlay(repID = 138)
plot_track_overlay(repID = 132)