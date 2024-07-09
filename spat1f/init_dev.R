library(foreach)
library(doSNOW)
library(survival)

# Some custom functions in "./spat1f/R/*.R" are loaded as a simulated 'spat1f' package
.simulate_load_package <- function(path){
  #source(paste0(path,"/R/", "package_load.R"), local = TRUE)
  pkgload::load_all(path = path,
                    export_all = TRUE, 
                    quiet = TRUE, 
                    recompile = FALSE)
}

.simulate_load_package(path = "C:/R_projects/spat_1f_noise/spat1f")
