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


plot_track_overlay


repID <- 55
trt_spec <- resize(as.cimg_color(
  flip_xy(fetch_trt_spec(repID, quiet = TRUE))
), 
size_x = 1000, 
size_y = 1000
)

pb_par_lapply(
  1:100, 
  function(rankID, repID, trt_spec){
    fn <- sprintf("rep%s_rank%s", repID, rankID)
    jpeg(paste0("C:/R_projects/spat_1f_noise/invisible/test_video_overlay/",fn,".jpg"), width = 10, height = 10, units = "cm", res = 300)
    plot_instance_overlay(repID, rankID, cex = 0.8, trt_spec = trt_spec)
    dev.off()
    return(invisible(NULL))
  }, 
  repID = repID, 
  trt_spec = trt_spec, 
  cores = 1, 
  inorder = FALSE
)






