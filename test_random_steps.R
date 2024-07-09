source("spat1f/init_analysis.R")

ID <- 81

fetch_events(ID) %>% 
  clean_events() %$%
  move_seq(head_x, head_y) %>% 
  filter(!is.na(r)) %>% 
  add_random_steps(n = 100L,
                   sl_distr = fit_gamma(.$r),
                   ta_distr = fit_genvonmises(.$theta_rel)
  ) %>%
  flag_invalid_steps(remove = TRUE) %>%
  mutate(
    toxic = read_value(x2, y2, c(1000, 1000),
                       ref_img = fetch_trt_spec(ID, 
                                                .ref_data = ref_data, 
                                                quiet = TRUE)), 
    less_toxic = 1 - toxic
  ) %>%
  append_estimators(na_as_zero = TRUE) %>% 
  issf(
    case ~
      less_toxic +
      (cos_theta_pi + cos_2theta) +  
      (sl + logsl) + 
      strata(step_id),
    data = .) -> issf_fit




spec <- dummy_spec()
spec <- c(1:12 %% 2,
          (1:12+1) %% 2) %>% 
  rep(6) %>% 
  matrix(nrow = 12, ncol = 12) %>% 
  as.cimg()

fetch_trt_spec(5) %>% plot()

spec <- fetch_trt_spec(ID)

iterate_random_steps2(
  ta_sl_list = list(
    "sl" = list(
      make_gamma(1, 30),
      make_gamma(1, 10)
    ),
    "ta" = list(
      make_genvonmises(0.48, 0.267),
      make_genvonmises(0.48, 0.267)
    )
  ), 
  start = make_start2(0,500,500, 1), 
  n = 5000, 
  ref_grid = spec, 
  same_move = FALSE, 
  rss_coef = 0.5
) %>% 
  mutate(head_x = x, head_y = y) %>% 
  plot_track_overlay(
    repID = "Foo", 
    colored_track = "none", 
    trt_spec = flip_xy(spec), 
    plot_elements = "track"
  )

plot_track_overlay(
  repID = ID, 
  colored_track = FALSE, 
  trt_spec = spec, 
  plot_elements = "track"
)

reload()




ID <- sample(setdiff(fetch_repID(), problem_ids), size = 1)
x <- fetch_events(ID) %>% 
  clean_events() %$% 
  move_seq(head_x, head_y) %>% 
  .$r %>% 
  na.omit()

plot_track_overlay(repID = ID)

comp_dist(x)
x %>% 
  loghist(nclass = 50, log.p = FALSE, geom = "col", 
          draw_dist = list(fit_frechet(x), fit_invgamma(x), fit_lnorm(x), fit_gamma(x))
          )

extraDistr::r(10000, 5) %>% loghist(log.p = FALSE)


comp_dist(x, dist = c("frechet","invgamma","gamma","lnorm"))






survival_plot(x)







issf_fit <- list(
  "ta_updated" = list(
    make_genvonmises(kappa1 = 0.4867337, kappa2 = 0.2466265),
    make_genvonmises(kappa1 = 0.4867337, kappa2 = 0.2466265)
  ),
  "sl_updated" = list(
    make_gamma(shape = 0.5717441, scale = 50),
    make_gamma(shape = 0.5717441, scale = 30)
  )
)


# more toxic state 1
# less toxic state 1
# more toxic state 2
# less toxic state 2


spec <- rbinom(50^2, 1, 0.9) %>% image_unflatten()
spec <- syn_spec(beta = 3, plot = FALSE)
mat <- matrix(c(0.9, 0.1, 0.03, 0.97), nrow = 2, ncol = 2, byrow = TRUE)
mat <- matrix(1)

iterate_random_steps_states(ta_sl_list = list(
  "sl" = list(
    make_gamma(1, 150),
    make_gamma(1, 50),
    make_gamma(1, 6),
    make_gamma(1, 6)
  ),
  "ta" = list(
    make_unif(),
    make_unif(),
    make_genvonmises(0.48, 0.267),
    make_genvonmises(0.48, 0.267)
  )
), 
  n = 1000, 
  ref_grid = syn_spec(n = 100, beta = 2, plot = FALSE), 
  rss_coef = 0,
  transition_mat = dummy_transition_mat(), 
  max_xy = c(1000,1000)
) %>% 
  plot_track_overlay(
    colored_track = "none", 
    plot_elements = "track"
  )












iterate_random_steps_states(
  ta_sl_list = list(
    "sl" = list(
      make_gamma(1, 150),
      make_gamma(1, 5)
    ),
    "ta" = list(
      make_unif(),
      make_unif()
    )
  ), 
  n = 5000, 
  ref_grid = syn_spec(n = 100, beta = 1, plot = FALSE), 
  rss_coef = 0,
  transition_mat =  diag(1), 
  max_xy = c(5000,5000)
) %>% 
  plot_track_overlay(
    colored_track = "none", 
    plot_elements = "track"
  )



reload()




source("spat1f/init_analysis.R")

  






z <- iterate_random_steps_states(ta_sl_list = list(
  "sl" = list(
    make_gamma(1, 150),
    make_gamma(1, 50),
    make_gamma(1, 6),
    make_gamma(1, 6)
  ),
  "ta" = list(
    make_unif(),
    make_unif(),
    make_genvonmises(0.48, 0.267),
    make_genvonmises(0.48, 0.267)
  )
), 
n = 10000, 
ref_grid = syn_spec(n = 12, beta = 2, plot = FALSE), 
rss_coef = 0,
transition_mat = dummy_transition_mat(), 
max_xy = c(1000,1000)
)


z %>% 
  plot_track_overlay(plot_elements = "track")



cat_size <- c(1)
for (i in seq_along(z$on_toxic)){
  cat_size[i+1] <- grow_cat(cat_size[i], 
                            r = c(0.98, 0.98, 1.2,1.05), 
                            index = z$on_toxic[i] + 1 + (z$state[i] - 1) * 2)
}
plot(cat_size, log = "y")

grow_cat <- function(x0, r, index){
  x0 * r[index]
}















ID <- 44
fetch_events(ID) %>% 
  clean_events() %>% 
  mutate(on_toxic = read_value(head_x, head_y,ref_img = fetch_trt_spec(ID))) %>% 
  mutate(
    ratio = c(NA,diff(log(size_px))),
    avg_toxic = foo(on_toxic)
  ) %>% 
  group_by(
    avg_toxic
  ) %>% 
  do(as.data.frame(t(summarise_vec(.$ratio, na.rm = TRUE))))


foo <- function(x){
  (x + c(x[-1], NA))/2  
}


exp(0.00274 * 1000)



library(rEDM)

ID <- 44
d <- fetch_events(ID) %>% 
  clean_events()

d_meta <- d %>% 
  select(time, size_px)

d2 <- d %$% 
  move_seq(head_x, head_y) %>% 
  mutate(
    start_toxic = read_value(x1, y1, ref_img = fetch_trt_spec(ID)),
    on_toxic = read_value(x2, y2, ref_img = fetch_trt_spec(ID))
  )

d2$state <- moveHMM::viterbi(fit_HMM(
  d %$% move_seq(head_x, head_y) %>% 
    as.moveData()
))

d2 <- d2 %>% 
  cbind(d_meta[-1,])


d2$cat_size_log <- log(d2$size_px)
d2 <- d2 %>% 
  mutate(
    lr = log(r)
  )

EmbedDimension(dataFrame = d2, 
               lib = "1 500", 
               pred = "501 800", 
               columns = "lr",
               target = "lr")


PredictNonlinear(dataFrame = d2, 
                 lib = "1 500", 
                 pred = "501 800", 
                 columns = "lr", 
                 target = "lr", 
                 E = 7)



res <- SMap(dataFrame = d2,
     lib = "200 400", 
     pred = "401 600", 
     theta = 1,
     E = 2, 
     columns = "lr", 
     target = "lr", 
     embedded = FALSE)
rEDM::ComputeError(res$predictions$Observations, res$predictions$Predictions)





library(plotly)

d2 %>% 
  mutate(
    hour = round(time / 1)
  ) %>% 
  group_by(hour) %>% 
  summarise(
    lr = log(mean(r, na.rm = TRUE)),
    theta_rel = mean(theta_rel, na.rm = TRUE),
    on_toxic = mean(on_toxic, na.rm = TRUE),
    state = mean(state, na.rm = TRUE)
  ) %>% 
  .[-nrow(.),] %>% 
  mutate(
    lag_lr = c(lr[-c(1:10)], rep(NA,10)),
    lag_theta = c(theta_rel[-1], NA),
    lag2_lr = c(lr[-c(1,2)], NA, NA)
  ) %>% 
  .[-nrow(.),] %>% 
  .[-nrow(.),] ->w

w2 <- lapply(w$hour, function(x){
  cbind(w[w$hour <= x,], "frame" = x)
}) %>% 
  do.call("rbind", .)

w2 %>%
  plot_ly(
    x = ~lr, 
    y = ~lag_lr,
    frame = ~ frame
  ) %>% 
  add_markers() %>% 
  add_trace(name = 'trace 0', 
            mode = 'lines')

w %>% 
  ggplot(aes(x = lag_lr, y = lr)) + 
  geom_point() + 
  geom_path()


CCM(dataFrame = d2, 
    E = 6, 
    columns = "on_toxic", 
    target = "lr", 
    sample = 50, 
    tau = -6,
    libSizes = "50 600 300", 
    showPlot = TRUE)



w$lr %>% 
  acf(na.action = na.pass)
timeLag(df$sheep, technique = "ami", lag.max = 50)









vars <- c("r", "theta_rel", "on_toxic", "start_toxic", "cat_size_log", "state")
vars_pairs <- combn(vars, 2)
ccm_mat <- array(NA, dim = rep(length(vars), 2), dimnames = list(vars, vars))

for(i in 1:ncol(vars_pairs)){
  ccm_out <- CCM(dataFrame = d2, 
                 E = 6, 
                 columns = vars_pairs[1, i], 
                 target = vars_pairs[2, i], 
                 sample = 50, 
                 libSizes = paste(nrow(d) - 6, nrow(d)-6, 10, collapse = " "), 
                 showPlot = FALSE)
  outVars <- names(ccm_out)
  var_out <- unlist(strsplit(outVars[2], ":"))
  ccm_mat[var_out[2], var_out[1]] <- ccm_out[1,2]
  
  var_out <- unlist(strsplit(outVars[3], ":"))
  ccm_mat[var_out[2], var_out[1]] <- ccm_out[1,3]
}


diag(ccm_mat) <- 1

ccm_mat

cor_mat <- abs(cor(d2[,vars], use = "pair"))


ccm_mat2 <- ccm_mat
ccm_mat2[ccm_mat2 < cor_mat] <- 0


ccm_mat2[ccm_mat<0] <- 0

ccm_mat2







w <- fetch_events(100) %>% 
  clean_events()

w$state <- c(NA, viterbi(
  fit_HMM(as.moveData(move_seq(w$head_x, w$head_y)))
))

z <- lapply(seq_len(nrow(w)), function(d){
  x <- w$head_x
  y <- w$head_y
  
  x1 <- c(x[-seq_len(d)], rep(NA, d))
  y1 <- c(y[-seq_len(d)], rep(NA, d))
  sqd <- ((x - x1)^2 + (y - y1)^2)
  
  res <- tapply(sqd, w$state, function(x){
    mean(x, na.rm = TRUE)
  }) %>% 
    unname()
  
  c("msd1" = res[1], "msd2" = res[2], "msd" = mean(sqd, na.rm = TRUE))
}) %>% 
  do.call("rbind", .)


data.frame(
  "w" = seq_len(nrow(w)),
  z
) %>% 
  gather(
    key = state, value = msd, msd1:msd,
  ) %>% 
  mutate(
    msd = (msd),
    w = (w)
  ) %>% 
  mutate(
    msd_thr = ifelse((w) > 10, NA, msd)
  ) %>% 
  ggplot(aes(x = w, y = (msd), color = state)) +
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10() + 
  ggpmisc::stat_ma_line(
    aes(y = (msd_thr))
  ) + 
  ggpmisc::stat_ma_eq(
    eq.with.lhs = "italic(hat(y))~`=`~",
    ggpmisc::use_label(c("eq", "R2"))) 





foo <- function(ID, ref_data){
  w <- fetch_events(ID) %>% 
    clean_events(ref_data = ref_data)
  
  w$state <- c(NA, viterbi(
    fit_HMM(as.moveData(move_seq(w$head_x, w$head_y)))
  ))
  
  cut_off <- 30
  z <- lapply(1:cut_off, function(d){
    x <- w$head_x
    y <- w$head_y
    
    x1 <- c(x[-seq_len(d)], rep(NA, d))
    y1 <- c(y[-seq_len(d)], rep(NA, d))
    sqd <- ((x - x1)^2 + (y - y1)^2)
    
    res <- tapply(sqd, w$state, function(x){
      mean(x, na.rm = TRUE)
    }) %>% 
      unname()
    
    c("msd1" = log10(res[1]), "msd2" = log10(res[2]), "msd" = log10(mean(sqd, na.rm = TRUE)))
  }) %>% 
    do.call("rbind", .)
  
  res <- apply(z, 2, function(i){
    suppressMessages(lmodel2::lmodel2(log10(i) ~ log10(1:cut_off))$regression.results[2,3])
  })
  
  return(c(res, "rep_id" = ID))
}



a <- pb_par_lapply(i, 
                   FUN = function(ID, ref_data){
  tryCatch(foo(ID, ref_data), error = function(e){
    c(NA, NA, NA, ID)
  })
}, ref_data = ref_data, 
cores = 8, inorder = FALSE)


b <- ref_data %>% 
  left_join(bind_vec(a), by = "rep_id") %>% 
  gather(
    key = state, value = msd, msd1:msd
  ) 

b %>% 
  ggplot(aes(x = beta, y = as.numeric(msd), color = var_trt, group = var_trt)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5)) + 
  geom_pointrange(stat = "summary", position = position_dodge(0.5), color = "black") + 
  facet_wrap(~state) + 
  scale_y_log10()



b %>% 
  ggplot(aes(x = log(cat_pre_wt), y = as.numeric(msd), 
             color = beta)) +
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~state) + 
  scale_y_log10()
























