source("spat1f/init_analysis.R")

a <- detection_report(id_list$var) %>% 
  mutate(
    p = n_keypoints / frames
  ) %>% 
  filter(
    p > 0.8
  )

a <- a$repID %>% 
  lapply(
    function(x){
      p <- fetch_events(x) %>% 
        clean_events(keep_sus = FALSE) %>%
        filter(!out_of_frame & size_px > 100) %>% 
        insert_gaps() %>% 
        .$head_x %>% 
        is.na() %>% 
        mean()
      return(data.frame(
        "rep_id" = x,
        "p" = p
      ))
    }
  ) %>% 
  do.call("rbind", .)


b <- a %>% 
  filter((1-p) > 0.5)


ID <- sample(b$rep_id, 1)
fetch_events(ID) %>% 
  clean_events(keep_sus = TRUE) %>%
  filter(!out_of_frame & size_px > 100) %>% 
  insert_gaps() %>% 
  mutate(
    on_toxic = read_value(head_x, head_y, ref_img = fetch_trt_spec(ID)), 
    state = c(NA, viterbi(
      fit_HMM(as.moveData(move_seq(head_x, head_y)))
    ))
  ) %>% 
  ggplot(aes(x = factor(on_toxic), y = log(size_px), color = factor(on_toxic))) + 
  geom_point() + 
  geom_pointrange(stat = "summary", color = "black")





ID <- sample(b$rep_id, 1)
z <- fetch_events(ID) %>% 
  clean_events(keep_sus = TRUE) %>%
  filter(!out_of_frame & size_px > 100) %>% 
  insert_gaps() %>% 
  mutate(
    on_toxic = read_value(head_x, head_y, ref_img = fetch_trt_spec(ID)), 
    state = c(NA, viterbi(
      fit_HMM(as.moveData(move_seq(head_x, head_y)))
    ))
  )
z %>% 
  ggplot(aes(x = time, y = log(size_px), color = factor(on_toxic))) + 
  geom_point()















library(rstan)


standata <- within(list(), {
  y <- rnorm(100, 3)
  n <- length(y)
  sigma_irreg_sigma_prior = 1
})
model_file <- "test_stan.stan"
fit <- stan(file = model_file, data = standata,
            warmup = 1000, iter = 2000, chains = 2, verbose = TRUE)

fit






standata <- within(list(), {
  y <- as.numeric(na.omit(log(z$size_px)))
  n <- length(y)
  sigma_irreg_sigma_prior <- 1
  sigma_level_sigma_prior <- 1
  #n_mis <- sum(is.na(y))
  #y <- ifelse(is.na(y), Inf, y)
})
model_file <- "ssm.stan"
fit <- stan(file = model_file, data = standata,
            warmup = 1000, iter = 2000, chains = 2, verbose = TRUE)

rstan::get_posterior_mean(fit)
Rhat(sims = rstan::get_posterior_mean(fit))

extract(fit, "mu")$mu %>% 
  apply(2,Rhat)

extract_variable_matrix(example_draws(), "mu")






library(brms)
brm(x ~ ar(p = 1), data = ts(sample(c(rep(NA, 10), rnorm(100)), 100))) %>% 
  brms::stancode()


toyd <- data.frame(
  "t" = 1:100,
  "y" = sample(c(rep(NA, 10), rnorm(100)), 100)
)

get_prior(bf(y | mi() ~ mi(y_latent) + 0, family = gaussian()) + 
    bf(y_latent | mi() ~ ar(p = 1) + 0, family = gaussian()), 
    data = toyd %>% 
      mutate(
        y_latent = NA_real_
      )) 

brm_m <- brm(
  bf(y | mi() ~ mi(y_latent), family = gaussian()) + 
    bf(y_latent | mi() ~ ar(p = 1), family = gaussian()), 
  data = toyd %>% 
    mutate(
      y_latent = NA_real_
    ),
  prior = c(
    set_prior("normal(0,3)", class = "sigma", resp = "y"),
    set_prior("normal(0,3)", class = "sigma", resp = "ylatent"),
    set_prior("normal(0,3)", class = "b"),
    set_prior("normal(0,3)", class = "Intercept"),
    set_prior("unif(0,1)", class = "ar", resp = "ylatent")
  ), 
  cores = 2, 
  chains = 2
)

brm_m2 <- brm(
  bf(y | mi() ~ mi(y_latent) + 0, family = gaussian()) + 
    bf(y_latent | mi() ~ ar(p = 1), family = gaussian()), 
  data = data.frame(
    "y" = scale(log(z$size_px))
  ) %>% 
    mutate(
      y_latent = NA_real_
    ),
  prior = c(
    set_prior("normal(0,1)", class = "sigma", resp = "y"),
    set_prior("normal(0,1)", class = "sigma", resp = "ylatent"),
    set_prior("normal(1,1)", class = "b", coef = "miy_latent", resp = "y"),
    set_prior("normal(0,1)", class = "ar", resp = "ylatent")
  ), 
  cores = 2, 
  chains = 2
)

brms::standata(brm_m2) %>% str()
stancode(brm_m2)
posterior_interval(brm_m2) %>% 
  rownames() %>% 
  .[!grepl("[1-9][0-9]", .)]

brms::stancode(brm_m)






brms::standata(brm_m2) %>% str()

brms::standata(brm_m2)
foo <- function(r = 1, drift = 0, sigma_obs = 0.5, sigma_bm = 0.3, n_missing = 10){
  y <- c(1,rep(NA, 300))
  yobs <- rep(NA, length(y))
  for(i in seq_len(length(y) - 1)){
    y[i + 1] <- y[i] * r + rnorm(1, drift, sigma_bm)
    yobs[i + 1] <- rnorm(1, y[i+1], sigma_obs)
  }
  
  yobs[sample(seq_along(yobs), size = n_missing)] <- NA
  return(yobs[-1])
}
y <- foo(r = 1)
plot(y)
plot(foo(r = 1.005, drift = 0.001, sigma_obs = 0.2, sigma_bm = 0.1))

y <- as.numeric(scale(log(z$size_px), center = FALSE))
X <- cbind(z$on_toxic, as.numeric(z$state == 1))
X[1,] <- c(0, 1)
X[,1] <- inherit_val(X[,1])
X <- cbind(X, X[,1] * X[,2])

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

X


plot(y)

standata <- within(list(), {
  N_y <- length(y)
  Y_y <- ifelse(is.na(y), Inf, y)
  Jmi_y <- which(is.na(y))
  Nmi_y <- length(Jmi_y)
  lbmi_y <- -Inf
  ubmi_y <- Inf
  N_ylatent <- length(y)
  Y_ylatent <- rep(Inf, N_ylatent)
  Jmi_ylatent <-  which(is.infinite(Y_ylatent))
  Nmi_ylatent <- length(Jmi_ylatent)
  lbmi_ylatent <- -Inf
  ubmi_ylatent <- Inf
  N <- N_y
  Kar_ylatent <- 1
  J_lag_ylatent <- rep(1, N)
  X <- X
  N_predictors <- ncol(X)
})



model_file <- "stan/ssm_brm.stan"
fit <- rstan::stan(file = model_file, 
                   data = standata, 
                   warmup = 1000, 
                   iter = 2000, 
                   chains = 4, 
                   verbose = FALSE, 
                   cores = 4, 
                   control = list(max_treedepth = 15))


traceplot(fit, "Intercept_ylatent")


library(rstan)
library(tidybayes)

extract(fit) %>% str()

fit %>% 
  extract() %>% 
  lapply(function(x){
    apply(as.matrix(x), 2, rhat)
  })





extract_posterior <- function(x, FUN = posterior::summarise_draws){
  FUN <- match.fun(FUN)
  l <- extract(x)
  lapply(seq_along(l), function(i){
    param <- names(l)[i]
    FUN(as.matrix(l[[i]])) %>% 
      mutate(
        variable = sprintf("%s[%s]", param,gsub("\\.\\.\\.", "", variable))
      )
  }) %>% 
    do.call("rbind", .)
}




out <- extract_posterior(fit)

out %>% View()

out <- out %>% 
  mutate(
    t = as.numeric(gsub(".*\\[|]", "", variable)),
    variable_name = gsub("\\[.*", "", variable)
  ) 

out %>% 
  filter(
    variable_name == "Ymi_y"
  ) %>% 
  ggplot(aes(y = mean)) + 
  # geom_point(
  #   aes(x = which(is.na(y))[t])
  # ) + 
  geom_point(
    data = data.frame(
      t = seq_along(y),
      obs = y
    ),
    aes(y = obs, x = t), 
    color = "black"
  ) + 
  geom_line(
    data = out %>% 
      filter(
        variable_name == "Ymi_ylatent"
      ) %>% 
      cbind(
        "x" = X
      ),
    aes(x = t, color = factor(x.1), group = 1)
  ) + 
  coord_cartesian(ylim = c(0.8, 1.15))


out %>% 
  filter(!grepl("mi_", variable))

rstan::stan_plot(fit, c("beta_x_ylatent"))


