source("helper_functions/init.R")

IDs <- fetch_repID()

inf_fn <- list.files("raw_data/inferences/", full.names = TRUE)

inf_fn

list.files("cleaned_data/data_dicts/")


# pb_par_lapply(
#   inf_fn,
#   function(x){
#     
#     save_path <- paste0("cleaned_data/data_dicts/", gsub("_inference.*",".rds",basename(x)))
#     if(file.exists(save_path)){
#       return(NULL)
#     } else {
#       d <- read_csv(x)
#       l <- parse_inference(d)
#       saveRDS(l, save_path)
#     }
#   }, cores = 8, inorder = FALSE, export_fun_only = TRUE
# )

# pb_par_lapply(list.files("cleaned_data/data_dicts/", full.names = TRUE),
#     function(x){
#       
#     save_path <- paste0("cleaned_data/events/", gsub(".rds",".csv",basename(x)))
#     if(file.exists(save_path)){
#         return(NULL)
#     } else {
#       l <- readRDS(x)
#       d <- get_data(l)
#       write_csv(d, save_path) 
#     }
#   }, cores = 8, inorder = FALSE, export_fun_only = TRUE)


