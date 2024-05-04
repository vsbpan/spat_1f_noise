source("spat1f/init.R")


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



d <- read_csv("raw_data/T ni conditioning study.csv") %>% 
  filter(Trial!= "58")
d <- d %>% 
  mutate_if(is.character, .funs = tolower) %>% 
  mutate(
    first_diet = ifelse(is.na(first_diet), "n", first_diet),
    second_diet = ifelse(is.na(first_diet), "n", first_diet),
    cat_wt_log = log(cat_wt),
    cat_wt_log_scale = as.numeric(scale(log(cat_wt)))
  )
d_long <- d %>% 
  gather(
    key = type, value = choice, c(first_diet,second_diet,end_of_day_diet,next_morning_choice)
  ) %>% 
  mutate(
    choice_h = ifelse(choice == "h", 1, 0),
    choice_l = ifelse(choice == "l", 1, 0),
    choice_n = ifelse(choice == "n", 1, 0),
    choice_c = ifelse(choice == "c", 1, 0)
  )





d %>% 
  ggplot(aes(x = log(cat_wt), color = condition_trt, y = choice_2h)) + 
  geom_point() + 
  facet_wrap(~center_block_trt) + 
  geom_smooth(method = "lm")




d %>% 
  ggplot(aes(x = bin(log(cat_wt), by = 2), color = condition_trt, y = first_diet)) + 
  geom_count(alpha = 0.5) + 
  facet_wrap(~center_block_trt)








d %>% 
  ggplot(aes(x = bin(log(cat_wt), by = 2), color = condition_trt, y = end_of_day_diet)) + 
  geom_count(alpha = 0.5) + 
  facet_wrap(~center_block_trt)


d %>% 
  ggplot(aes(x = bin(log(cat_wt), by = 2), color = condition_trt, y = next_morning_choice)) + 
  geom_count(alpha = 0.5) + 
  facet_wrap(~center_block_trt)


names(d)

library(brms)


d



d_long$y <- with(d_long, cbind(choice_h, choice_l, choice_n, choice_c))
names(d_long)

get_prior(
  choice ~
    type * (cat_wt_log_scale * 
              center_block_trt * condition_trt) + 
    (1|Trial),
  data = d_long,
  family = categorical(link = "logit")
)

m <- brm(
  choice ~
    cat_wt_log_scale * 
    center_block_trt * condition_trt,
  data = d_long %>% 
    filter(type == "first_diet"),
  family = categorical(link = "logit"),
  prior = c(
    set_prior("normal(0, 3)", class = "b"),
    set_prior("normal(0, 2.5)", class = "Intercept")
  ),
  cores = 6,
  iter = 4000,
  chains = 6,
  control = list(adapt_delta = 0.99)
)
#saveRDS(m, "invisible/fitted_models/first_diet_m.rds")

plot(m)
check_brms(m)


m2 <- brm(
  choice ~
    cat_wt_log_scale * 
    center_block_trt * condition_trt,
  data = d_long %>% 
    filter(type == "end_of_day_diet"),
  family = categorical(link = "logit"),
  prior = c(
    set_prior("normal(0, 3)", class = "b"),
    set_prior("normal(0, 2.5)", class = "Intercept")
  ),
  cores = 6,
  iter = 4000,
  chains = 6,
  control = list(adapt_delta = 0.99)
)
saveRDS(m2, "invisible/fitted_models/end_of_day_diet_m.rds")

m3 <- brm(
  choice ~
    cat_wt_log_scale * 
    center_block_trt * condition_trt,
  data = d_long %>% 
    filter(type == "next_morning_choice"),
  family = categorical(link = "logit"),
  prior = c(
    set_prior("normal(0, 3)", class = "b"),
    set_prior("normal(0, 2.5)", class = "Intercept")
  ),
  cores = 6,
  iter = 4000,
  chains = 6,
  control = list(adapt_delta = 0.99)
)
saveRDS(m3, "invisible/fitted_models/next_morning_choice_m.rds")


m4 <- brm(
  choice ~
    cat_wt_log_scale * 
    center_block_trt * condition_trt,
  data = d_long %>% 
    filter(type == "second_diet"),
  family = categorical(link = "logit"),
  prior = c(
    set_prior("normal(0, 3)", class = "b"),
    set_prior("normal(0, 2.5)", class = "Intercept")
  ),
  cores = 6,
  iter = 4000,
  chains = 6,
  control = list(adapt_delta = 0.99)
)
saveRDS(m4, "invisible/fitted_models/second_choice_m.rds")


marginal_effects(m, terms = c("cat_wt_log_scale", "center_block_trt","condition_trt"))


sjPlot::plot_model(m, type = "pred", 
                   terms = c("cat_wt_log_scale", "condition_trt"))

sjPlot::plot_model(m, type = "pred", 
                   terms = c("cat_wt_log_scale", "center_block_trt"))

sjPlot::plot_model(m, type = "pred", 
                   terms = c("condition_trt", "center_block_trt"))













get_mean_pred <- function(model, terms, along_n = 300, newdata2 = NULL){
  predictors <- get_cleaned_newdata(model, terms = terms, n = along_n)
  
  if(!is.null(newdata2)){
    predictors <- cbind(predictors, newdata2)
  }
  
  epred <- posterior_epred(model, newdata = predictors, re_formula = NA)
  
  foo <- function(x){
    data.frame(
      predictors,
      "y_hat" = matrixStats::colMeans2(x)
    )
  }
  if(length(dim(epred)) == 2){
    foo(epred) %>% cbind("resp" = insight::find_response(model))
  } else {
    lapply(
      seq_len(dim(epred)[3]), function(i){
        cbind(foo(epred[,,i]), "resp" = dimnames(epred)[[3]][i])
      }
    ) %>% 
      do.call("rbind", .)
  }
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












m1_pred <- get_mean_pred(m, 
              terms = c("cat_wt_log_scale", "center_block_trt","condition_trt"), 
              along_n = 100)

g <- m1_pred %>% 
  mutate(
    center_block_trt = paste0("center_block: ", toupper(center_block_trt)),
    condition_trt = paste0("condition: ", toupper(condition_trt))
  ) %>% 
  ggplot(aes(x = exp(cat_wt_log_scale), y = y_hat, fill = resp)) + 
  facet_wrap(~center_block_trt+condition_trt) + 
  geom_area() + 
  scale_fill_viridis_d() + 
  theme_minimal(base_size = 15) + 
  scale_x_continuous(trans = "log10") + 
  labs(x = "Cat weight (g)", y = "Proportion", fill = "Choice", subtitle = "First Choice") + 
  theme(legend.position = "right")



m2_pred <- get_mean_pred(m2, 
                         terms = c("cat_wt_log_scale", "center_block_trt","condition_trt"), 
                         along_n = 100)

g <- m2_pred %>% 
  mutate(
    center_block_trt = paste0("center_block: ", toupper(center_block_trt)),
    condition_trt = paste0("condition: ", toupper(condition_trt))
  ) %>% 
  mutate(
    resp = factor(resp, levels = c("h","l","n","c"))
  ) %>% 
  ggplot(aes(x = exp(cat_wt_log_scale), y = y_hat, fill = resp)) + 
  facet_wrap(~center_block_trt+condition_trt) + 
  geom_area() + 
  scale_fill_viridis_d() + 
  theme_minimal(base_size = 15) + 
  scale_x_continuous(trans = "log10") + 
  labs(x = "Cat weight (g)", y = "Proportion", fill = "Choice", subtitle = "End of Day Diet") + 
  theme(legend.position = "right")


ggsave("graphs/end_of_day_choice.png", g, dpi = 400, width = 7, height = 7, bg = "white")





m3_pred <- get_mean_pred(m3, 
                         terms = c("cat_wt_log_scale", "center_block_trt","condition_trt"), 
                         along_n = 100)

g <- m3_pred %>% 
  mutate(
    center_block_trt = paste0("center_block: ", toupper(center_block_trt)),
    condition_trt = paste0("condition: ", toupper(condition_trt))
  ) %>% 
  mutate(
    resp = factor(resp, levels = c("h","l","n","c"))
  ) %>% 
  ggplot(aes(x = exp(cat_wt_log_scale), y = y_hat, fill = resp)) + 
  facet_wrap(~center_block_trt+condition_trt) + 
  geom_area() + 
  scale_fill_viridis_d() + 
  theme_minimal(base_size = 15) + 
  scale_x_continuous(trans = "log10") + 
  labs(x = "Cat weight (g)", y = "Proportion", fill = "Choice", subtitle = "Next Day Diet") + 
  theme(legend.position = "right")



ggsave("graphs/next_day_choice.png", g, dpi = 400, width = 7, height = 7, bg = "white")




m4_pred <- get_mean_pred(m4, 
                         terms = c("cat_wt_log_scale", "center_block_trt","condition_trt"), 
                         along_n = 100)

g <- m4_pred %>% 
  mutate(
    center_block_trt = paste0("center_block: ", toupper(center_block_trt)),
    condition_trt = paste0("condition: ", toupper(condition_trt))
  ) %>% 
  mutate(
    resp = factor(resp, levels = c("h","l","n","c"))
  ) %>% 
  ggplot(aes(x = exp(cat_wt_log_scale), y = y_hat, fill = resp)) + 
  facet_wrap(~center_block_trt+condition_trt) + 
  geom_area() + 
  scale_fill_viridis_d() + 
  theme_minimal(base_size = 15) + 
  scale_x_continuous(trans = "log10") + 
  labs(x = "Cat weight (g)", y = "Proportion", fill = "Choice", subtitle = "Second Choice") + 
  theme(legend.position = "right")



ggsave("graphs/second_choice.png", g, dpi = 400, width = 7, height = 7, bg = "white")








