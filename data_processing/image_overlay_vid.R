trt_design_data <- read_csv("raw_data/trt_spectra_meta/Nov_20_week1_trt_spectra_meta.csv")

ufl_trt_iml <- lapply(seq_len(nrow(trt_design_data)), function(i){
  x <- trt_design_data[i,]
  unlist(x[,grepl("spec_", names(x))]) 
}) %>% 
  lapply(
    image_unflatten
  )

names(ufl_trt_iml) <- paste0("rep_id", trt_design_data$rep_id)




files <- list.files(dest_dir, full.names = TRUE)

img <- fast_load_image(files[1000], FALSE)

img_mask <- ufl_trt_iml$rep_id16 %>% 
  resize(size_x = dim(img)[1], size_y = dim(img)[2])
img_mask <- imappend(imlist(img_mask, img_mask, img_mask), axis = "c")



imager::imdraw(img, img_mask, opacity = 0.05) %>% 
  plot()

pb_par_lapply(
  files, 
  function(x){
    img <- fast_load_image(files[i], FALSE) # Fast load must be FALSE for the printed orientation to be up
    img_overlay <- imager::imdraw(img, img_mask, opacity = 0.05) 
    jpeg::writeJPEG(
      img_overlay %>% as.bmp(), 
      target = paste("C:/R_Projects/spat_1f_noise/invisible/test_video_overlay", 
                     file_base_name(files[i]), 
                     sep = "/"), 
      quality = 1
    )
  }, 
  cores = 8
) %>% invisible()


make_video(src_dir = "invisible/test_video_overlay", 
           file = paste0(trialID,"_overlay2")) 







