
# Set up generic
# Print method for COCO_Json
print.COCO_Json <- function(x){
  cat(sprintf("COCO annotation with %s images and %s annotations\n", 
              nrow(x$images),
              nrow(x$annotations)
              ))
  cat(sprintf("things: %s", paste0(x$categories[,"name"], collapse = ", ")))
}



# Convert data dict object to Json object
as.Json.data_dict <- function(x){
  
  format_set <- function(x){
    if(unique_len(x) == 1){
      n <- length(x)
      x <- unique(x)
    } else {
      n <- 1
    }
    
    trial_sets <- c("trial_1", "trial_2", "trial_3", "trial_4", "trial_5", "trial_6", 
                    "week2test", "week2test2", "week2test3", "week2test4", "week2test5")
    
    out <- vapply(x, function(x){
      ifelse(
        numbers_only(x), 
        sprintf("%01d%03d", 0, as.numeric(x)),
        ifelse(
          x %in% trial_sets, 
          sprintf("%01d%03d", 1, which(x %in% trial_sets)),
          sprintf("%01d%03d", 2, floor(runif(1, 0, 1000)))
        )
      )
    }, 
    character(1), 
    USE.NAMES = FALSE)
    
    return(
      rep(out, n)
    )
  }
  
  
  
  fmeta <- get_file_meta(x)
  n <- nrow(fmeta)
  #empty_df <- data.frame(row.names = seq_len(n))
  #empty_list <- vector(mode = "list", length = n)
  empty_list_list <- lapply(seq_len(n), function(x) vector(mode = "list", length = 0))
  
  img_id <- sprintf("1%s%04d",
                    format_set(repID_clean(fmeta$repID)),
                    as.numeric(fmeta$rank)
  )
  
  # id_lookup <- function(fn, db){
  #   vapply(
  #     fn, 
  #     function(x){
  #       db[(db$file_name == x), "id"]
  #     }, 
  #     numeric(1)
  #   ) %>% 
  #     unname()
  # }
  
  #img_id <- id_lookup(fmeta$file_base, db)
  #img_id <- seq_len(n)
  
  
  
  thing_id <- switch(unlist(unname(attr(x, "summary")["thing_class"])), 
                     "cat" = "1")
  
  
  kp <- lapply(map(x, "keypoints"), function(z){
    if(is.null(z)){
      return(NULL)
    } else {
      z <- rbind(t(z[,c(1,2)]), c(2,2,2))
      round(as.vector(z))
    }
  }) %>% 
    unname()
  
  seg <- map(x, "polygon") %>% 
    lapply(function(z){
      if(is.null(z)){
        return(
          vector(mode= "list", length = 0)
        )
      } else {
        return(
          list(as.vector(t(z)))
        )
      }
    })
  
  annotations <- data.frame(
    "category_id" = as.integer(1),
    "id" = as.integer(as.numeric(sprintf("%s%s", img_id, thing_id))),
    "image_id" = as.integer(img_id),
    "iscrowd" = FALSE,
    "num_keypoints" = as.integer(ifelse(lapply(kp, is.null) %>% 
                                          do.call("c",.), 
                                        0,
                                        3)),
    #"isbbox" = FALSE,
    "color" = "#f963b6"
    #"area" = as.integer(100)
  )
  
  #annotations$metadata <- empty_df
  annotations$bbox <- unname(get_bbox(x)) 
  annotations$keypoints <- kp
  annotations$segmentation <- unname(seg)
  
  #annotation_ord <- c("id", "image_id", "category_id", "segmentation", "bbox", "iscrowd", "isbbox", "color", "keypoints", "metadata", "num_keypoints")
  
  #annotations <- annotations[, annotation_ord]
  
  
  images <- data.frame(
    "file_name" = fmeta$file_base,
    "path" = fmeta$file_path, 
    "height" = as.integer(get_dim(x)[2]), 
    "width" = as.integer(get_dim(x)[1]),
    "id" = as.integer(img_id),
    #"regenerate_thumbnail" = FALSE,
    #"milliseconds" = as.integer(1),
    "deleted" = FALSE,
    "num_annotations" = as.integer(1)
    #"annotated" = FALSE
    #"dataset_id" = as.integer(1)
  )
  
  images$category_ids <- empty_list_list
  images$events <- empty_list_list
  images$annotating <- empty_list_list
  #images$metadata <- empty_df
  
  #img_lab_ord <- c("id", "category_ids", "path", "width", "height", "file_name", "annotated", "annotating", "num_annotations", "metadata", "deleted", "milliseconds", "events","regenerate_thumbnail")
  
  #images <- images[, img_lab_ord]
  
  
  categories <- data.frame(
    "color" = "#8186d5",
    "id" = as.integer(thing_id),
    "name" = c("cat"),
    "supercategory" = ""
  )
  
  categories$keypoint_colors <- list(
    c("#ff0000", "#000", "#2bff00")
  )
  
  categories$keypoints <- list(
    c("head", "middle","tail")
  )
  
  categories$skeleton <- list(
    matrix(c(1,2,2,3), ncol = 2)
  )
  
  
  
  out <- list(
    "images" = images,
    "categories" = categories, 
    "annotations" = annotations,
    "info" = list(),
    "licenses" = list()
  )
  
  class(out) <- c("COCO_Json", "list")
  
  return(out)
}

as.Json.list <- function(x, ...){
  class(x) <- c("COCO_Json", "list")
  x
}

as.Json <- function(x, ...){
  UseMethod("as.Json")
}

# S3 method registration
registerS3method(genname = "as.Json", 
                 class = "list", 
                 method = as.Json.list)

# S3 method registration
registerS3method(genname = "as.Json", 
                 class = "data_dict", 
                 method = as.Json.data_dict)

# S3 method for print COCO_Json
registerS3method(genname = "print", 
                 class = "COCO_Json", 
                 method = print.COCO_Json)


# Read in COCO format .json file and parse as COCO_jason object
import_COCO <- function(x){
  out <- jsonlite::fromJSON(x) %>% as.Json()
  
  
  i <- out$annotations[,"num_keypoints"]
  out$annotations[is.na(i),"num_keypoints"] <- 0
  
  out$images$num_annotations <- 1
  
  
  attr(out, "mtime") <- file.mtime(x) %>% 
    as.character() %>% 
    gsub("\\..*", "", .)
  attr(out, "src") <-  gsub("\\..*", "", basename(x))
  return(out)
}

# Export COCO_json as COCO format .json file
export_COCO <- function(x, path){
  stopifnot(is.COCO(x))
  jsonlite::write_json(x, path, pretty = TRUE)
}

# Sample COCO annotations
sample_COCO <- function(x, size){
  stopifnot(is.COCO(x))
  n <- nrow(x$images)
  s <- sample(seq_len(n), size, replace = FALSE)
  fns <- x$images$file_name[s,]
  
  subset_COCO(x, fns[s])
  return(x)
}

# Subset COCO annotation with file_name that is a subset of file_name. 
subset_COCO <- function(x, file_name){
  stopifnot(is.COCO(x))
  i <- x$images$file_name %in% file_name
  x$images <- x$images[i, ]
  ids <- x$images[,"id"]
  x$annotations <- x$annotations[x$annotations$image_id %in% ids, ]
  return(x)
}

# Split COCO annotations into test, train, and validation datasets and output list of COCO_json
split_COCO <- function(x, test, val = 0){
  stopifnot(is.COCO(x))
  n <- nrow(x$images)
  fns <- x$images$file_name
  s_val <- sample(seq_len(n), round(n * val), replace = FALSE)
  s_test <- sample(seq_len(n)[-s_val], round(n * test), replace = FALSE)
  s_train <- seq_len(n)[-c(s_val, s_test)]
  
  
  train <- subset_COCO(x, fns[s_train])
  test <- subset_COCO(x, fns[s_test])
  val <- subset_COCO(x, fns[s_val])
  
  return(list(
    "train" = train, 
    "test" = test, 
    "val" = val
  ))
}


# Merge COCO_json objects. Columns are filled with NA if missing. 
merge_COCO <- function(...){
  dots <- as.list(match.call())[-1]
  current_env <- parent.frame(n = 1)
  dots <- lapply(dots, function(x){
    x <- eval(x, envir = current_env)
    stopifnot(is.COCO(x))
    return(x)
  })
  
  nimg <- do.call("c", lapply(dots, function(x){
    nrow(x$images)
  }))
  dots <- dots[nimg > 0]
  
  categories <- lapply(seq_along(dots), function(i, x){
    o <- x[[i]]$categories
    rownames(o) <- paste0(i,"_",rownames(o))
    o
  }, x = dots) %>% 
    do.call("rbind",.) %>% 
    unique()
  
  
  annotations <- lapply(seq_along(dots), function(i, x){
    o <- x[[i]]$annotations
    if(nrow(o) == 0){
      return(NULL)
    }
    rownames(o) <- paste0(i,"_",rownames(o))
    o
  }, x = dots) %>% 
    do.call("rbind",.) %>% 
    unique()
  
  images <- lapply(seq_along(dots), function(i, x){
    o <- x[[i]]$images
    rownames(o) <- paste0(i,"_",rownames(o))
    o
  }, x = dots) %>% 
    do.call("rbind",.) %>% 
    unique()
  
  # n <- length(images$id)
  # images$id <- as.integer(seq_len(n))
  # annotations$id <- images$id
  # annotations$image_id <- images$id
  
  out <- list(
    "images" = images,
    "categories" = categories, 
    "annotations" = annotations,
    "info" = list(),
    "licenses" = list()
  )
  
  class(out) <- c("COCO_Json", "list")
  
  return(out)
}

# Set new root path for coco_Json object
set_new_path <- function(x, path_root){
  stopifnot(is.COCO(x))
  x$images$path <- paste0(path_root, "/", x$images$file_name)
  return(x)
}

# check if the object is COCO_Json
is.COCO <- function(x){
 inherits(x, "COCO_Json") 
}



# Reformat the target COCO_Json object such that only columns in the reference object are selected
force_format_COCO <- function(target_COCO, reference_COCO){
  stopifnot(is.COCO(target_COCO), is.COCO(reference_COCO))
  
  target_COCO$images <- target_COCO$images %>%
    dplyr::select(any_of(names(reference_COCO$images))) %>%
    as.data.frame()
  target_COCO$categories <- target_COCO$categories %>%
    dplyr::select(any_of(names(reference_COCO$categories))) %>%
    as.data.frame()
  target_COCO$annotations <- target_COCO$annotations %>%
    dplyr::select(any_of(names(reference_COCO$annotations))) %>%
    as.data.frame()
  return(target_COCO)
}

# Remove annotations from COCO_Json object
wipe_annotations_COCO <- function(x){
  stopifnot(is.COCO(x))
  
  #n <- length(x$annotations$bbox)
  #empty_list_list <- lapply(seq_len(n), function(x) vector(mode = "list", length = 0))
  
  x$images$num_annotations <- 0
  #x$annotations$num_keypoints <- 0
  #x$annotations$bbox <- empty_list_list
  #x$annotations$keypoints <- empty_list_list
  #x$annotations$segmentation <- empty_list_list
  nms <- names(x$annotations)
  x$annotations <- data.frame(matrix(ncol = length(nms), nrow = 0)) %>% 
    append_name(nms)
  return(x)
}


# Take a reference COCO_Json object, wipe the annotations, and append any that is found in manual_COCO. 
update_manual_COCO <- function(reference_COCO, manual_COCO){
  stopifnot(is.COCO(manual_COCO), is.COCO(reference_COCO))
  
  empty_COCO <- subset_COCO(reference_COCO, 
                            setdiff(reference_COCO$images$file_name,
                                    manual_COCO$images$file_name)
  )
  empty_COCO <- wipe_annotations_COCO(empty_COCO) 
  
  
  out <- merge_COCO(force_format_COCO(manual_COCO, reference_COCO), empty_COCO)
  return(out)
}
