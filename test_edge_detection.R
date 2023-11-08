library(tidyverse)
library(herbivar)
library(foreach)
library(doSNOW)
source("helper_functions/image_utils.R")


trt_design_data <- read_csv("raw_data/example_trt_spectra_meta.csv")

ufl_trt_iml <- lapply(seq_len(nrow(trt_design_data)), function(i){
  x <- trt_design_data[i,]
  unlist(x[,grepl("spec_", names(x))]) 
}) %>% 
  lapply(
    image_unflatten
  )

names(ufl_trt_iml) <- paste0("rep_id", trt_design_data$rep_id)

#################################################################################






img <- load.image("misc_tests/prototype3.jpg")

plot(img)

pts <- detect_corners(img, qc_plot = TRUE)



# plot(img)
# pts[[2]] <- list(x = 260, y = 1120)
# pts[[1]] <- list(x = 260, y = 180)

# do.call("rbind", pts) %>% as.data.frame() %$% points(x,y, col = "green", cex = 1.5)

img2 <- reproject_grid(img, init_pts = pts, dest_size = 1000, qc_plot = FALSE)
img2

plot(img2)


# img2 %>% 
#   immask(
#     ufl_trt_iml[[1]] %>% 
#       imager::resize(size_x = 1000, size_y = 1000, interpolation = 1) %>% 
#       as.pixset(), 
#     background = NA
#   ) %>% 
#   thin(3) %>% 
#   detect_luster() %>%
#   threshold2(thr = 0.65, thr.exact = TRUE) %>% 
#   herbivar::plot.pixset(col.na = "green")

plot(img2)
detect_luster(img2, shadow_weight = 0.7) %>%
  threshold2(thr = 0.55, thr.exact = TRUE) %>% 
  plot()

img3 <- detect_luster(img2) %>%
  threshold2(thr = 0.65, thr.exact = TRUE)


mask1 <- ufl_trt_iml[[1]] %>% 
  imager::resize(size_x = 1000, size_y = 1000, interpolation = 1) %>% 
  as.pixset()




imlist(
  img,
  img2,
  img3 %>% 
    as.cimg_color() %>% 
    immask(mask1, background = "blue"), 
  img3 %>% 
    as.cimg_color() %>% 
    immask(!mask1, background = "blue")
) %>% 
  plot.imlist(main.panel = c("Raw image", "Transformed diet image", 
                       "High toxin herbivory", "Low toxin herbivory"))









fast_load_image("misc_tests/prototype3.jpg")


plot(img2)

detect_luster(img2, shadow_weight = 0) %>% 
  imagerExtra::SPE(0.005) %>% 
  medianblur(n = 10) %>% 
  renorm(max = 1) %>% 
  threshold2(thr = 0.66, thr.exact = TRUE) %>% 
  plot()
  
  plot()
  threshold2(thr = 0.6, thr.exact = TRUE) %>% 
  plot()

  imagerExtra::ThresholdAdaptive(k = 0.01,range = c(0,1), windowsize = 101) %>% 
  plot()
?imagerExtra::ThresholdAdaptive
  
  
  
  
  
  
  
  
  
  
  

  
when_ago <- function(x, t_final, now = Sys.time()){
  now - (t_final - x)
}
  
when_ago(32434, 230276, as.POSIXct("2023-10-30 10:02:40"))






