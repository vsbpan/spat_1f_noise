library(foreach)
library(doSNOW)
#library(amt)
library(survival)
#library(aster)

# Some custom functions in "./spat1f/R/*.R" are loaded as a simulated 'spat1f' package
devtools::load_all(path = "C:/R_projects/spat_1f_noise/spat1f", export_all = TRUE, quiet = TRUE)
