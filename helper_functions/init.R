library(foreach)
library(doSNOW)
#library(aster)

# Some custom functions in "./helper_functions/R/*.R" are loaded as a simulated 'spat1f' package
devtools::load_all(path = "helper_functions", export_all = TRUE, quiet = TRUE)
