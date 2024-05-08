source("spat1f/init_analysis.R")

inf_fn <- list.files("raw_data/inferences/", full.names = TRUE)

inf_fn
overwrite <- TRUE

pb_par_lapply(
  inf_fn,
  function(x, overwrite){

    save_path <- paste0("cleaned_data/data_dicts/", gsub("_inference.*",".rds",basename(x)))
    if(file.exists(save_path) && !overwrite){
      return(NULL)
    } else {
      d <- read_csv(x)
      l <- parse_inference(d)
      saveRDS(l, save_path)
    }
  }, cores = 8, inorder = FALSE, export_fun_only = TRUE, overwrite = overwrite
)


# has at least 1 false cluster
problem_ids_2 <- c("1","6", "9", "10", "20", "14", "22", "28", "61", "69", "18",
                   "106", "107", "11", "114", "119", "13", "130", "146", "84",
                   "15", "19", "2", "25", "4", "5", "77", "8", "80", "85",
                   "96", "98", "109", "136", "134", "87", "113", "115", "119", "140", 
                   "37", "60", "65")

pb_par_lapply(fetch_repID("data_dict"),
    function(x, overwrite, problem_ids_2){

    save_path <- paste0("cleaned_data/events/rep",x, ".csv")
    if(file.exists(save_path) && !overwrite){
        return(NULL)
    } else {
      # if(!x %in% problem_ids_2){
      #   return(NULL)
      # }
      l <- fetch_data_dict(x)
      d <- get_data(l)
      d <- d %>%
        insert_gaps() %>%
        flag_mask_size() %>%
        filter(!is_gap) %>%
        mutate(false_cluster = FALSE)
      
      if(x %in% problem_ids_2){
        for(k in 1){
          i <- d %>%
            filter(score >= 0.7) %>%
            flag_false_cluster(bin_size = 60, r_thresh = 200, 
                               r = 60, return_id = TRUE)
          
          d <- d %>%
            mutate(
              false_cluster = ifelse(rank %in% i, TRUE, false_cluster)
            )
        }

      } 
      write_csv(d, save_path)
      invisible()
    }
  }, cores = 6, inorder = FALSE, export_fun_only = TRUE,
  overwrite = overwrite, problem_ids_2 = problem_ids_2)








# Has at least 2 false clusters
problem_ids_3 <- c("6", "9", "10", "20", "14", "22", "28", "61", "69", "60",
                   "11","13", "15", "4", "85", "18", "65", "84",
                   "109", "136", "134", "87", "113", "115", "119", "140", "37")

pb_par_lapply(problem_ids_3,
              function(x, overwrite, problem_ids_3){

                save_path <- paste0("cleaned_data/events/rep",x, ".csv")
                if(file.exists(save_path) && !overwrite){
                  return(NULL)
                } else {
                  l <- fetch_data_dict(x)
                  d <- get_data(l)
                  d <- d %>%
                    insert_gaps() %>%
                    flag_mask_size() %>%
                    filter(!is_gap) %>%
                    mutate(false_cluster = FALSE)
                  
                  if(x %in% problem_ids_3){
                    for(k in 1:2){
                      i <- d %>%
                        filter(score >= 0.7) %>%
                        flag_false_cluster(bin_size = 60, r_thresh = 200, 
                                           r = 60, return_id = TRUE)
                      
                      d <- d %>%
                        mutate(
                          false_cluster = ifelse(rank %in% i, TRUE, false_cluster)
                        )
                    }
                  } 
                  write_csv(d, save_path)
                  invisible()
                }
              }, cores = 6, inorder = FALSE, export_fun_only = TRUE,
              overwrite = overwrite, problem_ids_3 = problem_ids_3)



# Has at least 3 false clusters
problem_ids_4 <- c("6", "9", "10", "14", "28", "11","13", "15", "4", "85", "136")
pb_par_lapply(problem_ids_4,
              function(x, overwrite, problem_ids_4){
                
                save_path <- paste0("cleaned_data/events/rep",x, ".csv")
                if(file.exists(save_path) && !overwrite){
                  return(NULL)
                } else {
                  l <- fetch_data_dict(x)
                  d <- get_data(l)
                  d <- d %>%
                    insert_gaps() %>%
                    flag_mask_size() %>%
                    filter(!is_gap) %>%
                    mutate(false_cluster = FALSE)
                  
                  if(x %in% problem_ids_4){
                    for(k in 1:3){
                      i <- d %>%
                        filter(score >= 0.7) %>%
                        flag_false_cluster(bin_size = 60, r_thresh = 200, 
                                           r = 60, return_id = TRUE)
                      
                      d <- d %>%
                        mutate(
                          false_cluster = ifelse(rank %in% i, TRUE, false_cluster)
                        )
                    }
                  } 
                  write_csv(d, save_path)
                  invisible()
                }
              }, cores = 6, inorder = FALSE, export_fun_only = TRUE,
              overwrite = overwrite, problem_ids_4 = problem_ids_4)




# Has at least 4 false clusters
problem_ids_5 <- c("85", "4", "15", "13", "11", "28", "6", "9", "10", "14")
pb_par_lapply(problem_ids_5,
              function(x, overwrite, problem_ids_5){
                
                save_path <- paste0("cleaned_data/events/rep",x, ".csv")
                if(file.exists(save_path) && !overwrite){
                  return(NULL)
                } else {
                  l <- fetch_data_dict(x)
                  d <- get_data(l)
                  d <- d %>%
                    insert_gaps() %>%
                    flag_mask_size() %>%
                    filter(!is_gap) %>%
                    mutate(false_cluster = FALSE)
                  
                  if(x %in% problem_ids_5){
                    for(k in 1:4){
                      i <- d %>%
                        filter(score >= 0.7) %>%
                        flag_false_cluster(bin_size = 60, r_thresh = 200, 
                                           r = 60, return_id = TRUE)
                      
                      d <- d %>%
                        mutate(
                          false_cluster = ifelse(rank %in% i, TRUE, false_cluster)
                        )
                    }
                  } 
                  write_csv(d, save_path)
                  invisible()
                }
              }, cores = 6, inorder = FALSE, export_fun_only = TRUE,
              overwrite = overwrite, problem_ids_5 = problem_ids_5)


# False cluster flagging failed. Throw these out
problem_ids_6 <- c("85", "4", "15", "13", "11", "28", "6", "9", "10", "14")





