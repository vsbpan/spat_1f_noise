source("helper_functions/init_analysis.R")

IDs <- fetch_repID()

# inf_fn <- list.files("raw_data/inferences/", full.names = TRUE)
# 
# pb_par_lapply(
#   inf_fn, 
#   function(x){
#     d <- read_csv(x)
#     l <- parse_inference(d)
#     saveRDS(l, paste0("cleaned_data/data_dicts/", gsub("_inference.*",".rds",basename(x))))
#   }, cores = 6, inorder = FALSE, export_fun_only = TRUE
# )
# 
# pb_par_lapply(list.files("cleaned_data/data_dicts/", full.names = TRUE), 
#     function(x){
#     l <- readRDS(x)
#     d <- get_data(l)
#     write_csv(d, paste0("cleaned_data/events/", gsub(".rds",".csv",basename(x))))
#   }, cores = 6, inorder = FALSE, export_fun_only = TRUE)



