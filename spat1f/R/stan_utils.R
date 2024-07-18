
extract_posterior <- function(x){
  res <- rstan::summary(x, probs = c(0.025, 0.975))
  
  res <- res$summary %>% 
    as.data.frame() %>% 
    keep_rowname()
  names(res)[1] <- "variable"
  rownames(res) <- NULL
  res
}