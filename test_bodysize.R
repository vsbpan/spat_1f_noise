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








foo <- function(r = 1, drift = 0, sigma_obs = 0.5, sigma_bm = 0.3, n_missing = 10, N = 1000){
  y <- c(1,rep(NA, N))
  yobs <- rep(NA, length(y))
  for(i in seq_len(length(y) - 1)){
    y[i + 1] <- y[i] * r + rnorm(1, drift, sigma_bm)
    yobs[i + 1] <- rnorm(1, y[i+1], sigma_obs)
  }
  
  yobs[sample(seq_along(yobs), size = n_missing)] <- NA
  return(yobs[-1])
}





y <- foo(r = 1.005, drift = 0.001, sigma_obs = 0.2, sigma_bm = 0.1)
X <- matrix(nrow = length(y), ncol = 0)
plot(y)
plot(foo(r = 1.005, drift = 0.001, sigma_obs = 0.2, sigma_bm = 0.1))


y <- as.numeric(scale(log(z$size_px), center = FALSE))
# X <- cbind(z$on_toxic, as.numeric(z$state == 1))
# X[1,] <- c(0, 1)
# X[,1] <- inherit_val(X[,1])
# X <- cbind(X, X[,1] * X[,2])


X <- cbind(z$on_toxic)
X[1,] <- c(0)
X[,1] <- inherit_val(X[,1])

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
  Kar_ylatent <- 1 # Number of lags
  J_lag_ylatent <- rep(1, N)
  X <- X
  N_predictors <- ncol(X)
  has_lag <- ifelse(Kar_ylatent > 0, 1, 0)
  sigma_y_prior <- 0.01
  sigma_ylatent_prior <- 0.01
})



model_file <- "stan/ssm_brm.stan"
fit2 <- rstan::stan(file = model_file, 
                   data = standata, 
                   warmup = 1500, 
                   iter = 2000, 
                   chains = 4, 
                   verbose = FALSE, 
                   cores = 4, 
                   control = list(max_treedepth = 15, adapt_delta = 0.95))


traceplot(fit2, "beta_x_ylatent")
traceplot(fit2, "sigma_y")


rstan::extract(fit2, "sigma_y")

as.array(fit2)[,,"sigma_y"] %>% 
  Rhat()

library(rstan)
library(tidybayes)

extract(fit) %>% str()





out <- extract_posterior(fit2)
out <- out %>% 
  mutate(
    t = as.numeric(gsub(".*\\[|]", "", variable)),
    variable_name = gsub("\\[.*", "", variable)
  ) 

out %>% 
  filter(!grepl("mi_", variable))





out %>% View()



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
    aes(x = t, group = 1)
  )
  coord_cartesian(ylim = c(0.8, 1.15))






rstan::stan_plot(fit, c("beta_x_ylatent"))

extract_posterior(fit1)





ID <- 55
z <- fetch_events(ID) %>% 
  clean_events() %>% 
  mutate(
    on_toxic = read_value(head_x, head_y, ref_img = fetch_trt_spec(ID)), 
    state = c(NA, viterbi(
      fit_HMM(as.moveData(move_seq(head_x, head_y)))
    ))
  )



z %>% 
  ggplot(aes(x = head_x, y = head_y)) + 
  geom_path()





