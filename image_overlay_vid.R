files <- list.files(dest_dir, full.names = TRUE)


img_mask <- ufl_trt_iml[[14]] %>% 
  resize(size_x = dim(img)[1], size_y = dim(img)[2])
img_mask <- imappend(imlist(img_mask, img_mask, img_mask), axis = "c")


pb_par_lapply(
  files, 
  function(x){
    img <- fast_load_image(files[i], FALSE)
    img_overlay <- imager::imdraw(img, img_mask, opacity = 0.05) 
    jpeg::writeJPEG(
      img_overlay %>% as.bmp(), 
      target = paste("C:/R_Projects/spat_1f_noise/misc_tests/test_video_overlay", 
                     file_base_name(files[i]), 
                     sep = "/"), 
      quality = 1
    )
  }, 
  cores = 8
)


make_video(src_dir = "misc_tests/test_video_overlay", 
           file = paste0(trialID,"_overlay2")) 









