library(tidyverse)
library(herbivar)



img <- load.image("misc_tests/test_cropped.jpg")




plot(img)


color_index(img)


indices <- c("VARI", "GLAI", "NB", "HUE","S")

img <- img %>% rotate_90()
color_index(img, index = indices)

best_index <- "HUE"

img_HUE <- color_index(img, "HUE", plot = FALSE)$HUE
img_thresh <- threshold2(img_HUE)





imlist(img, img_HUE, img_thresh %>% as.cimg) %>% 
  plot(main.panel = c("Raspberry pi image","HUE color index","Segmentation"), 
              interpolate = FALSE)







