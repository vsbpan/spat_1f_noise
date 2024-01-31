source("helper_functions/init_analysis.R")

unname(unlist(id_list)) %>% 
  pb_par_lapply(function(x, ref_data){
    m <- detect_herbivory(x, time = fetch_cutoff(x, ref_data), plot = FALSE)
    jpeg::writeJPEG(image = as.bmp(m), quality = 1, 
                    target = sprintf("invisible/consumption_mask/%s.jpg", x))
    invisible(NULL)
  }, cores = 5, inorder = FALSE, ref_data = ref_data)



