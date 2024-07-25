check_brms <- function(model, integer = NULL, plot = TRUE, asFactor = FALSE, ...) {
  if(is.null(integer)){
    integer <- switch(insight::get_family(model)$type, 
                      "real" = FALSE,
                      "int" = TRUE)
  }
  
  if(insight::get_family(model)$family %in% c("multinomial")){
    resp <- brms::posterior_predict(model, ndraws = 1000)
    y <- get_y(model)
    pred <- colMeans(brms::posterior_epred(model, ndraws = 1000, re.form = NA))
    
    for (i in seq_len(ncol(y))){
      
      DHARMa::createDHARMa(
        simulatedResponse = t(resp[,,i]),
        observedResponse = y[,i], 
        fittedPredictedResponse = pred[,i],
        integerResponse = FALSE) %>% 
        plot(asFactor = asFactor, title = dimnames(resp)[[3]][i], ...)
    }
    return(invisible(NULL))
  }
  
  
  if(insight::get_family(model)$family %in% c("categorical")){
    levels <- unique(get_y(model))
    y <- apply(levels,1,function(x){as.numeric(get_y(model) == x)})
    resp <- brms::posterior_predict(model, ndraws = 1000)
    resp <- array(do.call("c", lapply(levels,function(x){as.numeric(resp == x)})), 
                  dim = c(1000, nrow(y), length(levels)))
    pred <- colMeans(brms::posterior_epred(model, ndraws = 1000, re.form = NA))
    
    for (i in seq_len(ncol(y))){
      
      DHARMa::createDHARMa(
        simulatedResponse = t(resp[,,i]),
        observedResponse = y[,i], 
        fittedPredictedResponse = pred[,i],
        integerResponse = FALSE) %>% 
        plot(asFactor = asFactor, title = dimnames(resp)[[3]][i], ...)
    }
    return(invisible(NULL))
  }
  
  dharma.obj <- DHARMa::createDHARMa(
    simulatedResponse = t(brms::posterior_predict(model, ndraws = 1000)),
    observedResponse = get_y(model), 
    fittedPredictedResponse = colMeans(brms::posterior_epred(model, ndraws = 1000, re.form = NA)),
    integerResponse = integer)
  
  if (isTRUE(plot)) {
    plot(dharma.obj, asFactor = asFactor, ...)
  }
  
  invisible(dharma.obj)
  
}


# Get posterior draws
epred_draws <- function(model, newdata){
  if(missing(newdata)){
    newdata <- model$data
  }
  
  res <- brms::posterior_epred(model, 
                               newdata = newdata, 
                               allow_new_levels = TRUE, 
                               re_formula = NA) %>% 
    t()
  colnames(res) <- paste0("draw_", seq_len(ncol(res)))
  cbind(newdata,res) %>% 
    gather(key = draw, value = val, contains("draw_"))
}


pred_draws2 <- function(model, terms, along_n = 300, newdata2 = NULL, ...){
  predictors <- get_cleaned_newdata(model, terms = terms, n = along_n)
  
  if(!is.null(newdata2)){
    predictors <- cbind(predictors, newdata2)
  }
  
  epred <- posterior_epred(model, newdata = predictors, re_formula = NA, ...)
  
  foo <- function(x){
    data.frame(
      predictors,
      "y_hat" = t(x)
    )
  }
  if(length(dim(epred)) == 2){
    res <- foo(epred) %>% cbind("resp" = insight::find_response(model))
  } else {
    res <- lapply(
      seq_len(dim(epred)[3]), function(i){
        cbind(foo(epred[,,i]), "resp" = dimnames(epred)[[3]][i])
      }
    ) %>% 
      do.call("rbind", .)
  }
  
  res %>% 
    gather(key = draw, value = y_hat, contains("y_hat")) %>% 
    mutate(draw = gsub("y_hat\\.","draw_", draw))
}

get_cleaned_newdata <- function(model, terms = NULL, n = 300){
  
  predictor_frame <- insight::get_data(model)
  yname <- insight::find_response(model)
  predictor_frame <-predictor_frame[,!names(predictor_frame) %in% yname]
  
  rand_names <- insight::find_random(model)$random
  var_names <- names(predictor_frame)
  
  v <- lapply(terms, function(z){
    switch(as.character(grepl("\\[", z) && grepl("\\[", z)), 
           "TRUE" = as.numeric(unlist(strsplit(gsub(".*\\[|\\]","",z),","))),
           "FALSE" = NULL)
  })
  terms <- gsub("\\[.*","",terms)
  
  
  new_data <- expand.grid(lapply(seq_along(predictor_frame), function(i, d, n, v){
    x <- d[, i]
    
    if(names(d)[i] %in% terms){
      j <- which(terms %in% names(d)[i])
      
      if(is.null(v[[j]])){
        if(is.numeric(x)){
          x <- seq_interval(x, n[j])
        } else {
          x <- unique(x) 
        }
      } else {
        x <- v[[j]]
      }
      return(x) 
    } else {
      if(is.numeric(x)){
        x <- mean(x)
      } else {
        if(names(d)[i] %in% rand_names){
          x <- "foooooooooooooooo"
        } else {
          x <- unique(x)  
        }
      }
      return(x)
    }
  }, 
  d = as.data.frame(predictor_frame), 
  n = n, 
  v = v))
  
  names(new_data) <- names(predictor_frame)
  
  return(new_data)
}



