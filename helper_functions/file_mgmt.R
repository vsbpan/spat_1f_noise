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
  
  cat("Moving imges to local drive . . . \n")
  move_files_pb(source_dir, dest_dir, 
                file_names = list.files(source_dir, pattern = ".jpg"), 
                pb = TRUE)
  
  file_names_og <- list.files(dest_dir, full.names = TRUE)
  cat("\n Done!\n")
  
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
  
  cat("Renaming relative time . . . \n")
  n <- sum(file.rename(from = file_names_og, to = file_names_new))
  message("Exported ", n, " files.")
}




move_files_pb <- function (from_dir, to_dir, file_names, pb = FALSE){
  n <- length(file_names)
  for (i in seq_len(n)) {
    file.rename(paste(from_dir, file_names[i], sep = "/"), 
                paste(to_dir, file_names[i], sep = "/"))
    if(pb){
      cat("Moving file", i, "of", n  ,"\r")
    }
  }
}


abs_path <- function(x){
  paste(getwd(),x, sep = "/")
}


unrank_files <- function(src_dir){
  f <- list.files(src_dir, full.names = TRUE)
  f2 <- gsub("_rank[0-9].*\\.jpg",".jpg",f)
  invisible(file.rename(
    f, 
    f2
  ))
}


attach_order_file_name <- function(src_dir){
  f <- list.files(src_dir, full.names = TRUE)
  
  if(any(grepl("_rank", f))){
    warning(paste0(sum(grepl("_rank", f)), " files already has assigned rank. Unranking. . ."))
    unrank_files(src_dir)
    f <- list.files(src_dir, full.names = TRUE)
  }
  
  f <- f[gsub(".*_s|.jpg","",f) %>% 
           as.numeric() %>% 
           order()]
  
  for (i in seq_along(f)){
    fi <- f[i]
    file.rename(
      fi, 
      paste0(
        gsub(".jpg","",fi),
        "_rank",i,
        ".jpg")
    )
  }
}



files_reorder <- function(x){
  no_rank_present <- !any(grepl("_rank",x))
  stopifnot(all(grepl("_rank",x)) | no_rank_present)
  if(no_rank_present){
    o <- gsub(".*_s|.jpg","",x) %>% 
             as.numeric() %>% 
             order()
  } else{
    o <- order(as.numeric(gsub(".*_rank|.jpg","",x)))
  } 
  
  x[o]
}


rename_start_time <- function(x){
  time <- gsub(".*_s|.jpg|_rank[0-9].*","",x) %>% 
    as.numeric()
  
  
  file.rename(x,gsub_element_wise(time, time - min(time), x)) %>% invisible()
}

gsub_element_wise <- function(pattern, replacement, x, 
                              ignore.case = FALSE, perl = FALSE, 
                              fixed = FALSE, useBytes = FALSE){
  n <- length(x)
  if(n > 1){
    if(length(replacement) == 1){
      replacement <- rep(replacement, n)
    }
    if(length(pattern) == 1){
      pattern <- rep(pattern, n)
    }
  }

  lapply(seq_along(x), function(i){
    gsub(pattern[i], replacement[i], x[i])
  }) %>% 
    do.call("c",.)
}


paste_path <- function(x){
  gsub(("\\\\"),"/",readLines("clipboard"))
}


file_base_name <- function(x){
  gsub(".*/","",x)
}


file_root <- function(x){
  gsub_element_wise(file_base_name(x), "", x)
}


save_anchors <- function(trialID){
  dest <- paste0("raw_data/picked_anchors/",trialID,".rds")
  if(file.exists(dest)){
    choice <- menu(c("Yes", "No"), title = sprintf("'%s' already exists. Overwrite?", trialID))
    if(choice == 2){
      message("Exiting.")
    }
  }
  saveRDS(pts_list, file = dest)
  message(sprintf("'%s' anchors list saved!", trialID))
}



report_missing_photos <- function(root_dir){
  stopifnot(dir.exists(root_dir))
  dirs <- list.dirs(root_dir, recursive = FALSE)
  for (i in seq_along(dirs)){
    fs <- list.files(dirs[i])
    
    if(length(fs) < 2){
      next
    }
    v <- fs %>% 
      gsub(".*_s|.jpg|_rank.*","",.) %>% 
      as.numeric() %>% 
      sort() %>% 
      diff() %>% 
      sapply(function(x){
        z <- (x /360)
        !((z < 1.2) & (z > 0.8))
      }) %>% 
      sum()
    message(sprintf("%s has %s missing photos", dirs[i], v))
  }
}

