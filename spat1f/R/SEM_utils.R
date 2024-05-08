recursive_split <- function(l, sem_coef, target){
  lapply(
    l,
    function(x){
      final <- (x[[1]]$response == target)
      if(final){
        list(x)
      } else {
        split_df_list(sem_coef[sem_coef$predictor == x[[1]]$response,], append = x)
      } 
    }
  ) %>% 
    purrr::keep(
      function(x){
        length(x) > 0
      }
    ) %>% 
    unlist(recursive = FALSE)
}


split_df_list <- function(x, append = NULL){
  rownames(x) <- NULL
  x <- apply(x, 1, function(z){
    x <- as.data.frame(t(z))
    if(is.null(append) || length(append) == 0){
      return(list(x))
    } else {
      return(c(list(x),append))
    }
  }, simplify = FALSE)
  return(x)
}


drop_list <- function(x){
  n <- length(x)
  if(n == 1){
    x[[1]]
  } else {
    x
  }
}



drop_list_iter <- function(x){
  n <- length(x)
  while(n == 1){
    n <- length(x)
    x <- drop_list(x)
  }
  return(x)
}



model_list <- function(x){
  purrr::keep(
    x, 
    function(z){
      insight::is_model(z)
    }
  )
}


append_std_coef <- function(x, response, coef){
  i <- which(x$coefficients$Response == response & x$coefficients$Std.Estimate == "-")
  coef <- coef[match(x$coefficients$Predictor[i],names(coef))]
  x$coefficients$Std.Estimate[i] <- coef
  return(x)
}


append_std_coef2 <- function(x, response, coef){
  i <- which(x$Std.Estimate == "-" & x$Response == response)
  coef <- coef[match(x$Predictor[i],names(coef))]
  x$Std.Estimate[i] <- coef
  return(x)
}

fix_coef_name <- function(x){
  names(x)[names(x) == ""] <- "star"
  x
}

find_significant_pairs <- function(x, p = 0.05){
  x %>%
    fix_coef_name() %>% 
    mutate(
      Predictor = gsub(":cat_pre_wt_log_scale","",Predictor)
    ) %>% 
    filter(P.Value < p) %>% 
    mutate(
      res = paste0(Response, "~", Predictor)
    ) %>% 
    .$res %>% 
    unique()
}


SEM_pred_coef <- function(sem_fit, var, target, cat_size, og_set, exclude = NA, only = NA){
  
  sem_coef <- suppressWarnings(coefs(sem_fit)) %>% 
    as.data.frame()
  sem_coef <- fix_coef_name(sem_coef)
  mod_list_resp <- do.call("c",lapply(model_list(sem_fit), insight::find_response))
  
  
  sem_coef <- append_std_coef2(sem_coef, 
                               resp = "on_toxic",
                               round(obs_emp_std(
                                 model_list(sem_fit)[[which(mod_list_resp == "on_toxic")]]
                               ), digits = 4))
  sem_coef <- append_std_coef2(sem_coef, 
                               resp = "ava_qual",
                               round(obs_emp_std(
                                 model_list(sem_fit)[[which(mod_list_resp == "ava_qual")]]
                               ), digits = 4))

  sem_coef <- sem_coef %>%
    filter(
      !grepl("~~",Predictor) 
    ) %>% 
    mutate(
      Std.Estimate = as.numeric(Std.Estimate),
      Std.Error = as.numeric(Std.Error),
      std_est = Std.Estimate
    ) %>% 
    mutate(
      Estimate = ifelse(
        grepl(":cat_pre_wt_log_scale", Predictor),
        cat_size * Estimate,
        Estimate
      ),
      Predictor = gsub(":cat_pre_wt_log_scale","",Predictor),
      std_est = ifelse(Predictor %in% c("var_high","beta_red","beta_white","beta_numeric_scale"), 
                       Estimate,
                       std_est),
      response = Response,
      predictor = Predictor
    ) %>% 
    group_by(
      response, predictor
    ) %>% 
    summarise(
      std_est = sum(std_est)
    ) %>% 
    filter(
      paste0(response, "~", predictor) %in% og_set
    ) %>% 
    suppressMessages() %>% 
    ungroup() 
  
  if(!is.na(null_to_NA(exclude))){
    sem_coef <- sem_coef %>% 
      filter(
      !grepl(paste0(c("cat_pre_wt", as.character(exclude)),collapse = "|"),predictor)
    )
  }
  
  sem_coef <- as.data.frame(sem_coef)
  
  l <- sem_coef %>% 
    filter(predictor == var) %>% 
    as.data.frame() %>% 
    split_df_list()
  
  for(i in 1:10){
    l <- recursive_split(l, sem_coef = sem_coef, target = target)
  }
  
  if(!is.na(null_to_NA(only))){
    l <- keep(l, function(x){
      any(only %in% (path_vars(x)))
    })
  }
  
  l %>% 
    lapply(function(x){
      prod(as.numeric(do.call("rbind", x)$std_est))
    }) %>% 
    do.call("sum",.)
}


bootSEM <- function(x, nboot = 10, cores = 1, silent = FALSE){
  pb_par_lapply(seq_len(nboot), function(i, sem_fit){
    data <- sem_fit$data
    newdata <- slice_sample(data, replace = TRUE, n = nrow(data))
    lapply(
      model_list(sem_fit),
      function(x){
        update(x, data = newdata)
      }
    ) %>% 
      as.psem()
  }, sem_fit = x, cores = cores, inorder = FALSE, silent = silent)
}


path_vars <- function(l){
  unique(do.call("c",lapply(l, function(x){c(x$response, x$predictor)})))
}

# Observed-empirical approach to standardization from Grace et al. 2018
obs_emp_std <- function(model){
  pred <- predict(model, type = "response")
  R2 <- cor(insight::get_response(model), pred)^2
  sd_yhat <- sqrt(var(pred)/R2)
  beta <- fixef(model)$cond[-1]
  sd_x <- sqrt(apply(insight::get_modelmatrix(model)[, names(beta), drop = FALSE], 2, var))
  return(
    beta * sd_x /  sd_yhat
  )
}
