# Format seconds into nice H:S:M format
hms_runtime <- function(x){
  h <- floor(x / 3600)
  m <- floor((x - h * 3600) / 60)
  s <- floor(x - h * 3600 - m * 60)
  cat(sprintf("\nRuntime %02d:%02d:%02d\n", h,m,s))
}

# Format seconds into nice H:S:M format 
hms_format <- function(x){
  h <- floor(x / 3600)
  m <- floor((x - h * 3600) / 60)
  s <- floor(x - h * 3600 - m * 60)
  return(sprintf("%02d:%02d:%02d", h,m,s))
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
                          inorder = TRUE, export_fun_only = TRUE){
  if(is.list(x)){
    indf <- function(x,i){
      x[[i]]
    } 
  } else {
    indf <- function(x,i) {
      x[i]
    }
  }
  
  has_clust <- inherits(cores, "cluster")
  
  if(!has_clust && (is.null(cores) || is.na(cores) || cores <= 1 || isFALSE(cores))){
    n <- length(x)
    out <- lapply(seq_along(x), FUN = function(i){
      cat(sprintf("\r%s %d of %d",loop_text,i, n))
      FUN(indf(x, i), ...)
    })
  } else {
    
    if(!has_clust){
      message(sprintf("\nInitializing %s parallel workers. . .", cores))
      
      check_CPU_request(cores, condition = "delay")
      
      cl <- makeCluster(cores, outfile = "")
      doSNOW::registerDoSNOW(cl)
    }
    
    
    indices <- seq_along(x)
    
    environment(FUN) <- environment()
    
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
    
    # Remove spat1f package from list. foreach::`%dopar%` calls library(package) as some point, which would give an error
    pkg <- .packages()
    pkg <- pkg[pkg != "spat1f"]
    
    out <- foreach(
      i = indices, 
      .export = export,
      .combine = c, 
      .verbose = FALSE,
      .inorder = inorder, 
      .final = invisible,
      .options.snow = list(
        progress = function(n) {
          cat(sprintf("\r%s %d of %d",loop_text, n, length(indices)))
        }
      ),
      .packages = 
    ) %dopar% {
      devtools::load_all(path = "helper_functions", 
                         export_all = TRUE, quiet = TRUE) # Load spat1f package
      
      list(FUN(indf(x, i), ...))
    }
    
    if(!has_clust){
      message("\nClosing parallel workers. . .")
      stopCluster(cl)
    }
    
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
              sprintf("%s  Not eough cores at the moment. %s of %s requested cores fulfilled. Waiting . . . \r",
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

