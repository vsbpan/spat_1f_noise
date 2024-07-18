
extract_posterior <- function(x){
  res <- rstan::summary(x, probs = c(0.025, 0.975))
  
  res <- res$summary %>% 
    as.data.frame() %>% 
    keep_rowname()
  names(res)[1] <- "variable"
  rownames(res) <- NULL
  res
}

inherit_val <- function(x){
  curr_val <- x[1]
  for(i in seq_along(x)){
    xi <- x[i]
    if(!is.na(xi)){
      curr_val <- xi
    }
    x[i] <- curr_val
  }
  x
}