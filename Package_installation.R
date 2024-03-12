# Installing packages and dependencies for this project

# Run this only when you first start working on this or when you need to update
# packages. This is just to help you get set up.

parse_dependencies <- function(description){
  l <- list(
    "imports" = strsplit(description$Imports, split = ",\n")[[1]],
   "depends" = strsplit(description$Depends, split = ",\n")[[1]],
   "suggests" = strsplit(description$Suggests, split = ",\n")[[1]]
  )
  
  l <- lapply(l, function(x){
    res <- gsub("\\(.*| ","",x)
    res[!res %in% "R"]
  })
  return(l)
}


check_installation_needed <- function(x, install = TRUE){
  owned_pkgs <- unname(installed.packages()[,1])
  x <- unlist(x)
  out <- x[!x %in% owned_pkgs]
  if(length(out) == 0){
    message("You are all set!")
    return(invisible())
  } else {
    if(install){
      message(sprintf("Trying to install missing packages: %s", 
                      paste0(out, collapse = ", ")))
      install.packages(out)
      return(invisible())
    } else {
      return(out)
    }
  }
}

# Read dependencies in the simulated package
pkg_discription <- packageDescription("spat1f", lib.loc = ".")
check_installation_needed(parse_dependencies(pkg_discription), install = TRUE)


# herbivar installation
# install.packages("devtools") # so you can install packages from GitHub
# install.packages("BiocManager") 
# BiocManager::install("EBImage")
# devtools::install_github("vsbpan/herbivar", build_vignettes = FALSE, dependencies = TRUE, force = TRUE)
