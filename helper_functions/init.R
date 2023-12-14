library(foreach)
library(doSNOW)

# Some custom functions in "./helper_functions/R/*.R" are loaded as 'spat1f' package
devtools::load_all(path = "helper_functions", export_all = TRUE, quiet = TRUE)
