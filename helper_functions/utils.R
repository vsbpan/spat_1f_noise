# Format seconds into nice H:S:M format
hms_runtime <- function(x){
  h <- floor(x / 3600)
  m <- floor((x - h * 3600) / 60)
  s <- floor(x - h * 3600 - m * 60)
  cat(sprintf("\nRuntime %02d:%02d:%02d\n", h,m,s))
}


# Find the frame time
frame_time <- function(x, fps){
  x <- x / fps # convert to seconds
  h <- floor(x / 3600)
  m <- floor((x - h * 3600) / 60)
  s <- floor(x - h * 3600 - m * 60)
  as.character(invisible(sprintf("%02d:%02d::%02d", h,m,s)))
}

# lapply() with progressbar and parellel support. 
pb_par_lapply <- function(x, FUN, cores = 1, ..., 
                          loop_text = "Processing",
                          inorder = FALSE, export_fun_only = TRUE){
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
    
    out <- foreach(
      i = indices, 
      .export = ls(globalenv()),
      .combine = c, 
      .verbose = FALSE,
      .inorder = inorder, 
      .final = invisible,
      .options.snow = list(
        progress = function(n) {
          cat(sprintf("\r%s %d of %d",loop_text, n, length(indices)))
        }
      ),
      .packages = .packages()
    ) %dopar% {
      list(FUN(indf(x, i), ...))
    }
    
    if(!has_clust){
      message("\nClosing parallel workers. . .")
      stopCluster(cl)
    }
    
  }
  
  return(out)
}


read_table <- function(x){
  read.table("clipboard",sep = "\t", header = TRUE)
}

