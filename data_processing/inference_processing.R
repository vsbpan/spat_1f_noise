source("spat1f/init.R")

IDs <- fetch_repID()

inf_fn <- list.files("raw_data/inferences/", full.names = TRUE)

inf_fn

list.files("cleaned_data/data_dicts/")
problem_ids_2 <- c("37", "65", "84", "87", "115", "140")
overwrite <- TRUE

# pb_par_lapply(
#   inf_fn,
#   function(x, overwrite){
# 
#     save_path <- paste0("cleaned_data/data_dicts/", gsub("_inference.*",".rds",basename(x)))
#     if(file.exists(save_path) && !overwrite){
#       return(NULL)
#     } else {
#       d <- read_csv(x)
#       l <- parse_inference(d)
#       saveRDS(l, save_path)
#     }
#   }, cores = 8, inorder = FALSE, export_fun_only = TRUE, overwrite = overwrite
# )

pb_par_lapply(fetch_repID("data_dict"),
    function(x, overwrite, problem_ids_2){

    save_path <- paste0("cleaned_data/events/rep",x, ".csv")
    if(file.exists(save_path) && !overwrite){
        return(NULL)
    } else {
      if(!x %in% problem_ids_2){
        return(NULL)
      }
      l <- fetch_data_dict(x)
      d <- get_data(l)
      d <- d %>%
        insert_gaps() %>%
        flag_mask_size() %>%
        filter(!is_gap)
      if(x %in% problem_ids_2){
        i <- d %>% 
          filter(score > 0.9) %>% 
          flag_false_cluster(bin_size = 60, r_thresh = 200, r = 60, return_id = TRUE)
        d <- d %>% 
          mutate(
            false_cluster = rank %in% i
          )
        
      } else {
        d <- d %>% 
          mutate(false_cluster = FALSE)
      }
      
      write_csv(d, save_path)
      invisible()
    }
  }, cores = 8, inorder = FALSE, export_fun_only = TRUE, 
  overwrite = overwrite, problem_ids_2 = problem_ids_2)





