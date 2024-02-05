source("helper_functions/init_analysis.R")

# unname(unlist(id_list)) %>%
#   pb_par_lapply(function(x, ref_data){
#     m <- detect_herbivory(x, time = fetch_cutoff(x, ref_data), plot = FALSE)
#     jpeg::writeJPEG(image = as.bmp(m), quality = 1,
#                     target = sprintf("invisible/consumption_mask/%s.jpg", x))
#     invisible(NULL)
#   }, cores = 5, inorder = FALSE, ref_data = ref_data)


out_list <- list.files("invisible/consumption_mask/") %>%
  gsub(".jpg", "", .) %>%
  pb_par_lapply(
    function(i, ref_data){
      mask <- fast_load_image(sprintf("invisible/consumption_mask/%s.jpg", i),
                              transform = FALSE)
      a <- mask_area(mask)
      mt <- mask %>%
        as.pixset() %>%
        as.data.frame() %$%
        read_value(x, y,
                   dim_xy = c(1000, 1000),
                   ref_img = fetch_trt_spec(i, .ref_data = ref_data, quiet = TRUE)) %>%
        suppressWarnings() %>%
        mean(na.rm = TRUE) %>%
        NaN_to_NA()

      data.frame(
        "rep_id" = i,
        "area_herb" = a,
        "mean_toxic" = mt
        )
    }, cores = 8, inorder = FALSE, ref_data = ref_data
  )

z <- do.call("rbind",out_list)
write_csv(z, "cleaned_data/consumption_mask_derivative.csv")





detect_herbivory(33, plot = TRUE, offset = 0)
reload()
debug(detect_herbivory)

repID <- 1
fetch_image(fetch_data_dict(repID),time = fetch_cutoff(repID)) %>% 
  HUE() %>% 
  imagerExtra::SPE(lamda = 0.01,range = c(0,1)) %>% 
  threshold2(thr = "otsu") %>% 
  plot()








