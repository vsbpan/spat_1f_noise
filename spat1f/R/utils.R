# Format seconds into nice H:S:M format
hms_runtime <- function(x){
  h <- floor(x / 3600)
  m <- floor((x - h * 3600) / 60)
  s <- floor(x - h * 3600 - m * 60)
  cat(sprintf("\nRuntime %02d:%02d:%02d\n", h,m,s))
}

# Format seconds into nice H:M:S format 
hms_format <- function(x){
  h <- floor(x / 3600)
  m <- floor((x - h * 3600) / 60)
  s <- floor(x - h * 3600 - m * 60)
  return(sprintf("%02d:%02d:%02d", h,m,s))
}

# Format seconds into nice H:M format
hm_format <- function(x){
  h <- floor(x / 3600)
  m <- round((x - h * 3600) / 60)
  return(sprintf("%02d:%02d", h,m))
}

# Format counts in a nice way
count_report <- function(x,n){
  sprintf("%s (%s%s)", x, signif(x/n*100, digits = 2),"%")
}


# Find the frame time
frame_time <- function(x, fps){
  x <- x / fps # convert to seconds
  h <- floor(x / 3600)
  m <- floor((x - h * 3600) / 60)
  s <- floor(x - h * 3600 - m * 60)
  as.character(invisible(sprintf("%02d:%02d::%02d", h,m,s)))
}

# lapply() with progress bar and parallel support. 
pb_par_lapply <- function(x, FUN, cores = 1, ..., 
                          loop_text = "Processing",
                          silent = FALSE,
                          inorder = TRUE, 
                          export_fun_only = TRUE){

  has_clust <- inherits(cores, "cluster")
  
  if(!has_clust && (is.null(cores) || is.na(cores) || cores <= 1 || isFALSE(cores))){
    if(is.list(x)){
      indf <- function(x,i){
        x[[i]]
      } 
    } else {
      indf <- function(x,i) {
        x[i]
      }
    }
    
    n <- length(x)
    
    out <- lapply(seq_along(x), FUN = function(i){
      if(!silent){
        cat(sprintf("\r%s %d of %d",loop_text,i, n))
      }
      FUN(indf(x, i), ...)
    })
    cat("\n")
  } else {
    
    if(!has_clust){
      message(sprintf("\nInitializing %s parallel workers. . .", cores))
      
      check_CPU_request(cores, condition = "delay")
      
      cl <- makeCluster(cores, outfile = "")
      doSNOW::registerDoSNOW(cl)
    }
    
    # Remove spat1f package from list. foreach::`%dopar%` calls library(package) as some point, which would give an error
    pkg <- .packages()
    pkg <- pkg[pkg != "spat1f"]
    spat1f_path <- path.package("spat1f")
    load_spat1f <- load_all2
    
    
    environment(FUN) <- environment()
    indices <- seq_along(x)
    
    if(export_fun_only){
      # Export only functions in the global environment (faster)
      export <- do.call("c",lapply(ls(globalenv()), 
                                   function(x){
                                     switch(is.function(get(x)), x, NULL)
                                    }))
    } else {
      # Otherwise export other global variables as well
      export <- ls(globalenv())
    }
    
    out <- tryCatch(foreach(
      i = x, # Passing large list directly as elements to avoid memory overflow
      .export = export,
      .combine = c, 
      .verbose = FALSE,
      .inorder = inorder, 
      .final = function(x){
        if(!has_clust){
          message("\nClosing parallel workers. . .")
          stopCluster(cl)
          has_clust <- TRUE
        }
        invisible(x)
      },
      .options.snow = list(
        progress = function(n) {
          if(!silent){
            cat(sprintf("\r%s %d of %d",loop_text, n, length(indices)))
          }
        }
      ),
      .packages = pkg
    ) %dopar% {
      load_spat1f(path = spat1f_path, 
                  export_all = TRUE, 
                  quiet = TRUE)
      list(FUN(i, ...))
    }, error = function(e){
      if(!has_clust){
        message(e)
        message("\nClosing parallel workers. . .")
        stopCluster(cl)
      }
    })
  }
  
  return(out)
}

# read in data table from clipboard
read_table <- function(x){
  read.table("clipboard",sep = "\t", header = TRUE)
}

# Return a data.frame of CPU usage
cpu_report <- function(x){
  system2('powershell', 
          c('-Command', 
            'Get-WmiObject -Query \'select Name, PercentProcessorTime from Win32_PerfFormattedData_PerfOS_Processor\' | foreach-object { write-host "$($_.Name): $($_.PercentProcessorTime)" }; '), 
          stdout = TRUE) %>% 
    strsplit(split = " : ") %>% 
    do.call("rbind",.) %>% 
    as.data.frame() %>% 
    rename(
      "core" = V1, 
      "CPU_usage" = V2
    ) %>% 
    filter(core != "_Total") %>% 
    mutate_all(as.numeric) %>% 
    arrange(core)
}

# Check CPU resource usage and if the requested number of cores has usage above 50%, apply error condition. 
check_CPU_request <- function(cores, 
                              condition = c("stop", "warning", "delay")){
  
  av_cores_check <- function(food){
    cpu_report() %>% 
      filter(CPU_usage < 50) %>% 
      nrow()
  }
  
  av_cores <- av_cores_check()
  if(av_cores >= cores){
    return(invisible(NULL))
  } else {
    condition <- match.arg(condition)
    
    if(condition == "stop"){
      stop("Not enough available cores.")
    } else {
      if(condition == "warning"){
        warning("Not enough available cores.")
      } else {
        if(condition == "delay"){
          for (i in 1:100){ # timeout after 1500 seconds
            cat(
              sprintf("%s  Not eough cores at the moment. %s of %s requested cores fulfilled. Waiting . . . \r\n",
                      hms_format((i-1) * 15),
                      av_cores, 
                      cores)
              )
            Sys.sleep(15) # Recheck every 15 seconds
            av_cores <- av_cores_check()
            if(av_cores >= cores){
              break
            }
          }
          if(av_cores >= cores){
            return(invisible(NULL))
          } else {
            stop("Not enough available cores. Function timed out.")
          }
        }
      }
    }
    
  }
}


# Returns the same object with supplied names as attribute
append_name <- function(x, name){
  names(x) <- name
  x
}

# Convert NULL value to NA vector of length 'len'
null_to_NA <- function(x, len = 1){
  if(is.null(x)){
    return(rep(NA, len))
  } else {
    x
  }
}

# Returns TRUE if the string has only numbers
numbers_only <- function(x){
  !grepl("\\D", x)
}

# Returns TRUE if is odd
is.odd <- function(x){
  x %% 2 == 1
}

# Evaluate a function on a vector with a double sided window of size w. 'w' must be odd. 
roll_vapply <- function(x, w, FUN){
  
  stopifnot(is.odd(w))
  stopifnot(length(x) >= w)
  n <- length(x)
  i <- seq_len(n)
  out <- rep(NA, n)
  FUN <- match.fun(FUN)
  side_len <- (w - 1)/2
  for(k in seq(
    side_len + 1, 
    n - side_len + 1,
    by = 1
  )){
    indices <- seq(k - side_len, k + side_len, by = 1)
    out[k] <- FUN(x[indices])
  }
  
  return(out)
}

# Check if is null, else, return a vector of is.na(x)
is_null_na <- function(x){
  if(is.null(x)){
    return(TRUE)
  } else {
    return(is.na(x))
  }
}

# Is between range?
is_between <- function(x, range, inclusive = FALSE){
  if(inclusive){
    return(
      x >= range[1] & x <= range[2]
    )
  } else {
    return(
      x > range[1] & x < range[2]
    )
  }
}


# For checking frames are valid and if missing and length == 1, choose 1. Allow x be multiple integers
assert_frames <- function(x, frames){
  missing_frame <- missing(frames)
  
  if(missing_frame){
    if(length(x) == 1L){
      frames <- 1L
    } else {
      stop("Missing arguent 'frames' with no default")
    }
  } else {
    if(!is.numeric(frames) || any((frames %% 1) != 0)){
      stop(sprintf("Expecting 'frames' to be an integer."))
    }
  }
  
  return(frames)
}

# For checking frames are valid and if missing and length == 1, choose 1. x can only be an atomic integer
assert_frame <- function(x, frame){
  missing_frame <- missing(frame)
  
  if(missing_frame){
    if(length(x) == 1L){
      frame <- 1L
    } else {
      stop("Missing arguent 'frame' with no default")
    }
  } else {
    n <- length(frame)
    if(n != 1){
      stop(sprintf("Expecting 1 frame, but received %s", n))
    }
    if(!is.numeric(frame) || (frame %% 1) != 0){
      stop(sprintf("Expecting 'frame' to be an integer."))
    }
  }
  
  return(frame)
}

# Check if two values are identical. If difference is greater than thresh, return the difference.
check_identical <- function(x1,x2, thresh = .Machine$double.eps){
  if(identical(x1, x2)){
    return(TRUE)
  }
  z <- x1 - x2 
  cond <- z < thresh 
  if(!cond){
    return(z)
  } else {
    return(cond)
  }
}

# Function to be used in another function for exposing the columns of the first argument (data.frame or tibble) to the arguments. Useful for dplyr coding style. 
.expose_columns_interal <- function(){
  curr_env <- parent.frame(n = 1) # Parent function environment
  f <- formals(sys.function(sys.parent())) # Fetch formals of parent function
  mcall <- as.list(
    match.call(sys.function(sys.parent()), call = sys.call(sys.parent()))
  )[-c(1,2)] # Fetch parent function calls, dropping the first argument
  f[match(names(mcall), names(f))] <- mcall # replace default formals with user specified arguments
  
  l <- lapply(f[-1], function(x){ # Force evaluation of arguments in parent function env
    eval(x, envir = get(names(f)[1], envir = curr_env))
  })
  lapply(seq_along(l), function(i){ # Assign objects to parent env
    assign(names(l)[i], l[[i]], envir = curr_env)
  })
}


# Helper function for wipping custom functions in the global env
wipe_functions <- function(){
  as.character(lsf.str(pos = globalenv())) %>% 
    lapply(function(x){
      rm(list = x, envir = globalenv())
    })
  invisible()
}

# Export rbind.fill function from plyr
rbind.fill <- function(...){
  plyr::rbind.fill(...)
}

# Lazy wrapper for reloading spat1f if it is already loaded
reload <- function(dev = FALSE){
  if(dev){
    source(paste(path.package("spat1f"), "init_dev.R", sep = "/"))
  } else {
    source(paste(path.package("spat1f"), "init.R", sep = "/")) 
  }
}

# Recompile C++ code
recompile <- function(){
  # pkg <- "spat1f"
  # path <- pkgload::pkg_path(pkg)
  # dllname <- sprintf("%s.dll", pkg)
  # dllpath <- pkgload:::package_file("src", dllname, path = path)
  # library.dynam.unload(pkg, paste0(path,"/src"))
  # pkgload:::unload_dll(pkg)
  # file.remove(dllpath)
  # load_all2(path)
  pkgload::load_all(pkgload::pkg_path("spat1f"))
  
}

# Convert Nan to NA
NaN_to_NA <- function(x){
  ifelse(is.nan(x), NA, x)
}

# Handy function for parsing concentration
parse_conc <- function(x){
  as.numeric(gsub(" mg/g","", x))
}

# Set rowname as the first column
keep_rowname <- function(x){
  cbind("rn" = rownames(x), x)
}

# lapply() but with names(z) <- x
lapply_name <- function(x, FUN, ...){
  z <- lapply(x, FUN, ...)
  names(z) <- x
  z
}

# Wrapper for formatC()
sigfig <- function(x, digits = 2){
  formatC(x, digits = digits, format = "fg", flag = "#")
}

# Turn named matrix into a named vector
flatten_mat_name <- function(x){
  setNames(c(x),
           paste(rep(colnames(x), each = nrow(x)), rownames(x), sep="__")
  )
}

# Reverse name and value
reverse_names <- function(x){
  val <- x
  nms <- names(x)
  names(nms) <- val
  return(nms)
}


# Bind matrices together into an array along the z dimension
zbind <- function(...){
  mcall <- as.list(match.call(expand.dots = TRUE))[-1]
  
  array(c(...), 
        dim = c(
          dim(eval(mcall[[1]])),
          length(mcall)
        ), dimnames = c(
          dimnames(eval(mcall[[1]])),
          list(
            #do.call("c",lapply(mcall,deparse1))
          )
        ))
}

