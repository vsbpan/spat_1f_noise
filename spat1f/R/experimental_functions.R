iterate_random_steps_states <- function(
                                  issf_fit = NULL,
                                  ta_sl_list = extract_ta_sl(issf_fit, updated = TRUE),
                                  start = make_start2(x = max_xy[1]/2, y = max_xy[2]/2),
                                  transition_mat = dummy_transition_mat(),
                                  n = 100,
                                  same_move = FALSE,
                                  max_xy = c(1000, 1000),
                                  ref_grid = dummy_spec(),
                                  rss_coef = 0, 
                                  keep_ref_grid = TRUE){
  stopifnot(nrow(start) == 1)
  
  
  ref_grid_flat <- c(ref_grid)
  dim_x <- nrow(ref_grid)
  dim_y <- ncol(ref_grid)
  
  
  if(is.null(ta_sl_list)){
    stop("'ta_sl_list' must not be NULL.")
  }
  
  stopifnot(
    inherits(ta_sl_list$ta[[1]], "ta_distr"),
    inherits(ta_sl_list$sl[[1]], "sl_distr")
  )
  
  ra <- rdist(ta_sl_list$ta, 1e5)
  rr <- rdist(ta_sl_list$sl, 1e5)
  
  ralen <- length(ra)
  rrlen <- length(rr)
  
  ntoxic <- ifelse(same_move, 1, 2)
  nstates <- ncol(transition_mat)
  
  if(ralen == 1 && rrlen == 1){
    if(!same_move){
      warning("'same_move' set to TRUE as only one sample of turn angles and step lengths are provided.") 
      same_move <- TRUE
    }
    ra <- rep(ra, nstates)
    rr <- rep(rr, nstates)
    
  } else {
    if(nstates * ntoxic != max(rrlen, ralen)){
      stop(sprintf("Length of step length and turn angle list (%s) not 1 or the same length as the number of unique state-toxic combinations (%s)", max(rrlen, ralen), nstates * ntoxic)) 
    }
  }
  
  if(ralen != rrlen){
    if(ralen > rrlen){
      fac <- ralen / rrlen
      if(fac %% 1 != 0){
        stop("Length of step length list is not a factor of the length of turn angle list.")
      }
      rr <- rep(rr, fac) # Reuse the list
    } else {
      fac <- rrlen / ralen
      if(fac %% 1 != 0){
        stop("Length of turn angle list is not a factor of the length of step length list.")
      }
      ra <- rep(ra, fac) # Reuse the list
    }
  }
  
  
  sim <- add_random_steps_iterate_statesC(
    n = n, 
    n_draws = 100L, 
    x_start = start$x,
    y_start = start$y, 
    direction_start = start$theta, 
    diet_start = start$on_toxic,
    state_start = start$state,
    sl_rand = rr, 
    ta_rand = ra, 
    rss_coef = rss_coef, 
    transition_mat = transition_mat,
    same_move = same_move,
    ref_grid_flat = ref_grid_flat, 
    max_x = max_xy[1],
    max_y = max_xy[2],
    dim_x = dim_x,
    dim_y = dim_y
  )
  
  res <- rbind.fill(start, sim)
  attr(res, "max_xy") <- max_xy
  attr(res, "ref_grid") <- ref_grid
  class(res) <- c("sim_steps", class(res))
  
  return(res)
}



dummy_transition_mat <- function(size = 2, sticky = TRUE){
  a <- rgamma(size, 10, scale = 0.03)
  res <- rdirichlet(size, a)
  
  if(sticky){
    concen_diag <- function(mat){
      swap <- function(matrixRow,x,y){
        indexY <- which(matrixRow == y)
        valX <- matrixRow[x]
        matrixRow[x] <- y
        matrixRow[indexY] <- valX
        return(matrixRow)
      }
      
      for(i in seq_len(nrow(mat))){
        rowI <- mat[i,]
        y <- max(rowI)
        mat[i,] <- swap(rowI, i, y)
      }
      return(mat)
    }
    res <- concen_diag(res)
  }
  return(res) 
}








