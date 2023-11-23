#### Set up ####
## This script is the work flow for generating treatment spectra, saving them, and writing helper image guide.  

library(tidyverse)
library(herbivar)
# library(imager)
source("helper_functions/spectral_utils.R")
source("helper_functions/image_utils.R")
herbivar::pre_cmp_fun()


#### Generate spectra list by treatment groups ####
# Grid edge length
grid_nl <- 12

# Number of reps per trt
n <- 8

beta_0_iml <- lapply(1:n,function(x){
  z <- syn_spec(n = grid_nl, beta = 0, plot = FALSE)
}) %>% 
  as.imlist()
names(beta_0_iml) <- paste0("beta",0,"_id",1:n)


beta_n5_iml <- lapply(1:n,function(x){
  syn_spec(n = grid_nl, beta = -5, plot = FALSE)
}) %>% 
  as.imlist()
names(beta_n5_iml) <- paste0("beta",-5,"_id",1:n)


beta_5_iml <- lapply(1:n,function(x){
  syn_spec(n = grid_nl, beta = 5, plot = FALSE)
}) %>% 
  as.imlist()
names(beta_5_iml) <- paste0("beta",5,"_id",1:n)


#### Merge and shuffle spectra as rep order ####

beta_list <- c(beta_0_iml,beta_n5_iml,beta_5_iml)

beta_list <- beta_list[sample(seq_len(length(beta_list)), 
                 size = length(beta_list), 
                 replace = FALSE)]

##### Quality Checks ####

# Check if there is equal number of black and white tiles in each raster
((lapply(beta_list, function(x){
  sum(x)
}) %>% do.call("c",.)) == grid_nl^2/2) %>% 
  all()


#### Save generated spectra ####

trt_design_data <- lapply(beta_list, function(x){
  z <- image_flatten(x)
  
  nx <- ny <- sqrt(length(z))
  
  d <- expand.grid(
    "x" = seq_len(nx),
    "y" = seq_len(ny)
  ) %>% 
    mutate(
      label = paste0("spec_", LETTERS[y],x)
    )

  names(z) <- d$label
  
  return(z)
}) %>% 
  bind_vec(keep_row_names = FALSE, 
           row_names_as_col = "syn_id") %>% 
  mutate(
    beta = as.numeric(gsub("beta|_id.*","", syn_id)),
    rep_id = seq_len(nrow(.)),
    syn_date = as.character(Sys.time())
  )

trt_design_data

#write_csv(trt_design_data, "raw_data/Nov_20_week1_trt_spectra_meta.csv")




#### Read generated spectra and plot helper image guide ####


trt_design_data <- read_csv("raw_data/Nov_20_week1_trt_spectra_meta.csv")

ufl_trt_iml <- lapply(seq_len(nrow(trt_design_data)), function(i){
  x <- trt_design_data[i,]
  unlist(x[,grepl("spec_", names(x))]) 
  }) %>% 
  lapply(
    image_unflatten
  )

names(ufl_trt_iml) <- paste0("rep_id", trt_design_data$rep_id)



for(i in seq_along(ufl_trt_iml)){
  ufl_trt_iml[[i]] %>% 
    plot_image_guide(main = names(ufl_trt_iml)[i])
  Sys.sleep(1)
}


# for (i in seq_along(ufl_trt_iml)){
#   fn <- names(ufl_trt_iml)[i]
#   jpeg(paste0("templates_for_print/", fn, ".jpg"), width = 12, height = 14, units = "cm", res = 600)
#   ufl_trt_iml[[i]] %>% 
#     plot_image_guide(mar = c(0,0,1,0), axes = FALSE, main = fn)
#   dev.off()
# }


# Misc tests

# plot.imlist(beta_n5_iml, 
#             mar.panel=c(1,1,1,1), 
#             main.panel = NA)
# 
# 
# plot.imlist(beta_0_iml, 
#             mar.panel=c(1,1,1,1), 
#             main.panel = NA)
# 
# 
# plot.imlist(beta_5_iml, 
#             mar.panel=c(1,1,1,1), 
#             main.panel = NA)
# 
# 
# 
# img <- beta_0_iml[[1]]





