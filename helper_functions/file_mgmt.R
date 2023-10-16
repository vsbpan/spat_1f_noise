library(herbivar)
init_dir <- function(root_path, dir_name){
  stopifnot(length(dir_name) == 1)
  stopifnot(length(root_path) == 1)
  dest_dir <- paste(root_path, dir_name, sep = "/")
  dest_dir <- gsub("//", "/", dest_dir)
  
  
  if(dest_dir %in% list.dirs(root_path)){
    message("Directory already exists. Function terminated.")
    return(invisible(NULL))
  } else {
    message("Creating directory ", add_quote(dir_name)," ...")
    dir.create(dest_dir, showWarnings = TRUE)
    message("Done!")
    return(invisible(NULL))
  }
}

read_time <- function(x){
  date <- gsub("_.*", "", x)
  time <- gsub(".*_", "", x) %>% gsub("\\.", ":", .)
  time_formated <- as.POSIXct(paste(
    date, time, sep = " "
  ), format = c("%Y-%m-%d %H:%M:%S"))
  return(time_formated)
}


get_relative_time <- function(x, units = "secs"){
  stopifnot(all(is.finite.POSIXlt(x)))
  
  as.numeric(difftime(x, min(x), units = units))
}


export_time_lapse_img <- function(source_dir, dest_dir_root){
  repID <- paste0("rep",readline(prompt = "Enter replication ID:\t "))
  
  cameraID <- paste0("cam",readline(prompt = "Enter camera ID:\t "))
  
  init_dir(root_path = dest_dir_root, repID)
  dest_dir <- paste0(dest_dir_root, repID)
  
  move_files(source_dir, dest_dir, file_names = list.files(source_dir, pattern = ".jpg"))
  
  file_names_og <- list.files(dest_dir, full.names = TRUE)
  
  rel_time <- gsub(".*img|.jpg","",file_names_og) %>% 
    read_time() %>% 
    get_relative_time()
  
  file_names_new <- paste0(
    dest_dir, 
    "/",
    repID, 
    "__", 
    cameraID,
    "_s",
    rel_time,".jpg"
  )
  
  n <- sum(file.rename(from = file_names_og, to = file_names_new))
  message("Exported ", n, " files.")
}










