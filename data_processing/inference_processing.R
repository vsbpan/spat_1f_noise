source("spat1f/init_analysis.R")

inf_fn <- list.files("raw_data/inferences/", full.names = TRUE)

inf_fn
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


# has at least 1 false cluster
problem_ids_2 <- c("6", "9", "10", "20", "18", "14", "22", 
                   "37", "28", "60", "61", "65", "69", "84", 
                   "87", "109", "113", "115", "134", "140", "136")

# pb_par_lapply(fetch_repID("data_dict"),
#     function(x, overwrite, problem_ids_2){
# 
#     save_path <- paste0("cleaned_data/events/rep",x, ".csv")
#     if(file.exists(save_path) && !overwrite){
#         return(NULL)
#     } else {
#       # if(!x %in% problem_ids_2){
#       #   return(NULL)
#       # }
#       l <- fetch_data_dict(x)
#       d <- get_data(l)
#       d <- d %>%
#         insert_gaps() %>%
#         flag_mask_size() %>%
#         filter(!is_gap)
#       if(x %in% problem_ids_2){
#         i <- d %>%
#           filter(score > 0.9) %>%
#           flag_false_cluster(bin_size = 60, r_thresh = 200, r = 60, return_id = TRUE)
#         d <- d %>%
#           mutate(
#             false_cluster = rank %in% i
#           )
# 
#       } else {
#         d <- d %>%
#           mutate(false_cluster = FALSE)
#       }
# 
#       write_csv(d, save_path)
#       invisible()
#     }
#   }, cores = 6, inorder = FALSE, export_fun_only = TRUE,
#   overwrite = overwrite, problem_ids_2 = problem_ids_2)








# Has at least 2 false clusters
problem_ids_3 <- c("9", "18", "28", "87", "113", "134")

# pb_par_lapply(problem_ids_3,
#               function(x, overwrite, problem_ids_3){
#                 
#                 save_path <- paste0("cleaned_data/events/rep",x, ".csv")
#                 if(file.exists(save_path) && !overwrite){
#                   return(NULL)
#                 } else {
#                   l <- fetch_data_dict(x)
#                   d <- get_data(l)
#                   d <- d %>%
#                     insert_gaps() %>%
#                     flag_mask_size() %>%
#                     filter(!is_gap)
#                   
#                   if(x %in% problem_ids_3){
#                     i <- d %>%
#                       filter(score > 0.9) %>%
#                       flag_false_cluster(bin_size = 60, r_thresh = 200, 
#                                          r = 60, return_id = TRUE)
#                     d <- d %>%
#                       mutate(
#                         false_cluster = rank %in% i
#                       )
#                     
#                     i <- d %>%
#                       filter(score > 0.9 & !false_cluster) %>%
#                       flag_false_cluster(bin_size = 60, r_thresh = 200, r = 60, return_id = TRUE)
#                     d <- d %>%
#                       mutate(
#                         false_cluster = ifelse(rank %in% i, TRUE, false_cluster)
#                       )
#                     
#                   }
#                   write_csv(d, save_path)
#                   invisible()
#                 }
#               }, cores = 1, inorder = FALSE, export_fun_only = TRUE,
#               overwrite = overwrite, problem_ids_3 = problem_ids_3)






# Has at least 3 false clusters
# Using manual removal
problem_ids_4 <- c("18", "28", "87")

# fetch_events(87) %>% 
#   mutate(
#     false_cluster = ifelse(
#       (is_between(head_x, c(920,935)) & is_between(head_y, c(550, 560))), 
#       TRUE,
#       false_cluster
#     )
#   ) %>% 
#   write_csv(paste0("cleaned_data/events/rep",87, ".csv"))

# fetch_events(28) %>%
#   mutate(
#     false_cluster = ifelse(
#       (is_between(head_x, c(940,980)) & is_between(head_y, c(800, 950))),
#       TRUE,
#       false_cluster
#     )
#   ) %>%
#   write_csv(paste0("cleaned_data/events/rep",28, ".csv"))

# fetch_events(18) %>%
#   mutate(
#     false_cluster = ifelse(
#       (is_between(head_x, c(940,965)) & is_between(head_y, c(766, 785))) | 
#         (is_between(head_x, c(945,985)) & is_between(head_y, c(860,910))),
#       TRUE,
#       false_cluster
#     )
#   ) %>%
#   write_csv(paste0("cleaned_data/events/rep",18, ".csv"))