# Internal function that takes a fitted issf, use the updated sl and ta distributions to resample mean neighborhood quality
.get_ava_neighborhood_quality_engine <- function(model, ID, .ref_data, n = 100L){
  resample_ava_steps(model, n = 100L) %>% 
    mutate(
      toxic = read_value(x2, y2, c(1000, 1000),
                         ref_img = fetch_trt_spec(ID, .ref_data = .ref_data, quiet = TRUE))
    ) %>% 
    group_by(step_id) %>% 
    summarise(
      ava_neighborhood_qual = (1 - mean(toxic))
    ) %>% 
    summarise(
      ava_neighborhood_qual = mean(ava_neighborhood_qual)
    ) %>% 
    unlist() %>% 
    unname()
}


# Internal function that takes a repID and compute the average temporal variance over some time window
.get_temporal_var_engine <- function(repID, 
                                     .ref_data = get("ref_data", pos = globalenv()), 
                                     hours = 12L){
  x <- repID
  l_conc <- .ref_data %>% filter(rep_id == x) %>% .$low_diet_numeric
  h_conc <- .ref_data %>% filter(rep_id == x) %>% .$high_diet_numeric
  
  roll_vapply2 <- function(x, w, FUN){
    n <- length(x)
    
    if(n < w){
      return(rep(NA_real_, n))
    }
    
    i <- seq_len(n)
    out <- rep(NA, n)
    FUN <- match.fun(FUN)
    s <- seq.int(0, w - 1)
    for(k in seq_len(n - w + 1)){
      out[k] <- FUN(x[s + k])
    }
    return(out)
  }
  
  fetch_events(x, append_detection_summary = FALSE) %>% 
    clean_events(ref_data = .ref_data) %>% 
    mutate(
      toxic = read_value(head_x, head_y, 
                         ref_img = fetch_trt_spec(x, 
                                                  .ref_data = .ref_data, 
                                                  quiet = TRUE), 
                         c(1000,1000))
    ) %$%
    roll_vapply2(toxic, w = 10 * hours, FUN = function(xx){
      xx <- xx[!is.na(xx)]
      s <- xx == 0
      xx[s] <- l_conc
      xx[!s] <- h_conc
      
      var(xx, na.rm = TRUE)
    }) %>% 
    mean(na.rm = TRUE)
}




# Internal function that takes a repID and compute the average proportion time on toxic diet
.get_mean_on_toxic_engine <- function(repID, 
                                      .ref_data = get("ref_data", pos = globalenv())){
  x <- repID
  fetch_events(x, append_detection_summary = FALSE) %>% 
    clean_events(ref_data = .ref_data) %>% 
    mutate(
      toxic = read_value(head_x, head_y, 
                         ref_img = fetch_trt_spec(x, 
                                                  .ref_data = .ref_data, 
                                                  quiet = TRUE), 
                         c(1000,1000))
    ) %$%
    mean(toxic, na.rm = TRUE)
}


# Internal function that takes a fitted issf model and extract the observed step length and turn angle summaries
.get_obs_move_summary_engine <- function(repID, .ref_data, probes = c("mean", "Kurt")) {
  outer_named_vec <- function(x){
    rn <- x[,1, drop = TRUE]
    cn <- colnames(x)[-1]
    out <- unlist(x[,-1, drop = FALSE])
    names(out) <- tolower(mutate(expand.grid("rn" = rn, "cn" = cn), 
                                 nms = sprintf("%s_%s", rn, cn))$nms)
    return(out)
  }
  
  fetch_events(repID, append_detection_summary = FALSE) %>% 
    clean_events(ref_data = .ref_data) %$%
    move_seq(head_x , head_y) %>% 
    dplyr::select(r, theta_rel) %>% 
    gather(key = type, value = val) %>% 
    mutate(type = ifelse(type == "r", "sl", "ta")) %>% 
    group_by(type) %>% 
    reframe(
      data.frame(t(probe_dist(val, probes = probes, na.rm = TRUE)))
    ) %>% 
    outer_named_vec() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename_all(.funs = function(x){
      paste0(x, "_obs")
    })
}

# Internal function that takes a repID and compute the average proportion time in state 1
.get_prop_state1_engine <- function(repID, 
                                      .ref_data = get("ref_data", pos = globalenv())){
  x <- repID
  fetch_events(x, append_detection_summary = FALSE) %>% 
    clean_events(ref_data = .ref_data) %$% 
    move_seq(head_x, head_y) %>% 
    as.moveData() %>% 
    fit_HMM() %>% 
    stateProbs() %>% 
    .[,1,drop = TRUE] %>% 
    mean()
}



# Nice wrapper for `.get_ava_neighborhood_quality_engine()` with list as input
extract_ava_neighborhood_quality <- function(issf_fit_l, 
                                             .ref_data = get("ref_data", 
                                                             pos = globalenv()), 
                                             n = 100L){
  res <- lapply(
    seq_along(issf_fit_l),
    function(i, issf_fit_l, n, .ref_data){
      id <- names(issf_fit_l)[i]
      cat(sprintf("Processing repID = %s, %s of %s             \r", 
                  id, 
                  i,
                  length(issf_fit_l)))
      .get_ava_neighborhood_quality_engine(
        issf_fit_l[[i]], 
        id, 
        .ref_data = .ref_data,
        n = n
      )
    }, 
    issf_fit_l = issf_fit_l,
    n = n,
    .ref_data = .ref_data
  )
  
  out <- data.frame(
    "rep_id" = names(issf_fit_l),
    "ava_qual" = do.call("c", res)
  )
  cat("\n")
  return(out)
}



# Nice wrapper for `.get_temporal_var_engine()` with repIDs as input
extract_temporal_var <- function(repIDs, 
                                 .ref_data = get("ref_data", pos = globalenv()), 
                                 hours = 12L,
                                 cores = 1){
  v <- repIDs %>%
    pb_par_lapply(
      function(x, .ref_data, hours){
        cat(sprintf("\tProcessing repID = %s             \r", 
                    x))
        .get_temporal_var_engine(x, .ref_data = .ref_data, hours = hours)
      }, 
      .ref_data = .ref_data,
      hours = hours,
      cores = cores
    ) %>% 
    do.call("c", .)
  
  out <- data.frame(repIDs, v)
  names(out) <- c("rep_id", paste0("var_toxic_",hours))
  cat("\n")
  return(out)
}

# Nice wrapper for `.get_mean_on_toxic_engine()` with repIDs as input
extract_mean_on_toxic <- function(repIDs, 
                                  .ref_data = get("ref_data", pos = globalenv()),
                                  cores = 1){
  v <- repIDs %>%
    pb_par_lapply(
      function(x, .ref_data){
        cat(sprintf("\tProcessing repID = %s             \r", 
                    x))
        .get_mean_on_toxic_engine(x, .ref_data = .ref_data)
      }, 
      .ref_data = .ref_data,
      cores = cores
    ) %>% 
    do.call("c", .)
  
  out <- data.frame(repIDs, v)
  names(out) <- c("rep_id", "on_toxic")
  cat("\n")
  return(out)
}

# Extract model coefficients from fitted issf list
exract_model_coef <- function(issf_fit_l){
  res <- issf_fit_l %>%
    lapply(function(x){
      as.data.frame(t(
        c(
          coef(x$model),
          do.call("c",x$sl_updated[[1]]$params),
          do.call("c",x$ta_updated[[1]]$params)
        )
      ))
    }) %>%
    do.call("rbind.fill", .)
  
  cbind("rep_id" = names(issf_fit_l), res)
}

# Nice wrapper for `.get_obs_move_summary_engine()` with repIDs as input
extract_obs_move_summary <- function(repIDs, 
                                     .ref_data = get("ref_data", pos = globalenv()), 
                                     probes = c("mean", "Kurt"), 
                                     cores = 1){
  v <- seq_along(repIDs) %>%
    pb_par_lapply(
      function(i, repIDs, .ref_data, probes){
        id <- repIDs[i]
        cat(sprintf("\tProcessing repID = %s, %s of %s             \r", 
                    id, 
                    i,
                    length(repIDs)))
        .get_obs_move_summary_engine(id, 
                                     .ref_data = .ref_data, 
                                     probes = probes)
      }, 
      repIDs = repIDs,
      .ref_data = .ref_data,
      probes = probes,
      cores = cores
    ) %>% 
    do.call("rbind", .)
  
  out <- data.frame("rep_id" = repIDs, v)
  cat("\n")
  return(out)
}

# Nice wrapper for `.get_prop_state1_engine()` with repIDs as input
extract_prop_state1_summary <- function(repIDs, 
                                     .ref_data = get("ref_data", pos = globalenv()),
                                     cores = 1){
  v <- seq_along(repIDs) %>%
    pb_par_lapply(
      function(i, repIDs, .ref_data){
        id <- repIDs[i]
        cat(sprintf("\tProcessing repID = %s, %s of %s             \r", 
                    id, 
                    i,
                    length(repIDs)))
        .get_prop_state1_engine(id, .ref_data = .ref_data)
      }, 
      repIDs = repIDs,
      .ref_data = .ref_data,
      cores = cores
    ) %>% 
    do.call("rbind", .)
  
  out <- data.frame("rep_id" = repIDs, "prop_explore" = v)
  cat("\n")
  return(out)
}


