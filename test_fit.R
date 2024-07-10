source("spat1f/init.R")


test <- function(x){
  z <- mean(x) / se(x)
  p <- pnorm(abs(z), lower.tail = FALSE) * 2
  if(p < 0.001){
    p <- "< 0.001"
  } else {
    p <- paste0("= ", sigfig(p))
  }
  
  data.frame(Z = sigfig(z), P = p)
}



sim_fun <- function(){
  spec <- syn_spec(n = 12, beta = 0, plot = FALSE)
  d <- iterate_random_steps_states(
    ta_sl_list = list(
      "sl" = list(
        make_gamma(1, 100),
        make_gamma(1, 30)
      ),
      "ta" = list(
        make_unif(),
        make_unif()
      )
    ), 
    transition_mat = diag(1),
    start = make_start2(0,500,500, 1), 
    n = 4000, 
    ref_grid = spec, 
    same_move = FALSE, 
    rss_coef = 0 
  )
  
  # d <- d[seq_len(nrow(d)) %% 2 == 0,]
  #d <- d[sample.int(nrow(d), size = floor(nrow(d) * 0.33), replace = FALSE),]
  
  d %$% 
    move_seq(x, y) %>% 
    add_random_steps(n = 100, sl_distr = make_gamma(1, 30), 
                     ta_distr = make_vonmises(0)) %>% 
    flag_invalid_steps(remove_invalid = TRUE) %>% 
    mutate(
      less_toxic_start = 1 - (round(read_value(x1, y1, 
                                               ref_img = spec, transform = FALSE))),
      less_toxic_end = 1 - (round(read_value(x2, y2, 
                                             ref_img = spec, transform = FALSE)))
    ) %>% 
    append_estimators() %>% 
    issf(
      case ~ 
        less_toxic_end + 
        as.factor(less_toxic_end):sl
        as.factor(less_toxic_start):sl +
        as.factor(less_toxic_start):logsl +
        strata(step_id), 
      data = ., keep_data = FALSE, update = FALSE
    ) -> m
  
  
  z <- unname(coef(m))
  k <- -diff(log(unname(1/((1/m$sl$params$scale) - z[2:3]))))
  return(list(
    "model" = m,
    "coef" = c("rss" = z[1], "k" = k)
  ))
  # return("coef" = c("rss" = z[1], "k" = k))
}

sim_fun()

debug(sim_fun)


out6 <- pb_par_lapply(1:200, FUN = function(x, sim_fun){
  sim_fun()
}, sim_fun = sim_fun, 
cores = 8, 
inorder = FALSE)




out <- readRDS("foo.rds")

res <- do.call("rbind", out6) %>% 
  as.data.frame() %>% 
  gather(key = variable, value = val)

g6 <- res %>% 
  ggplot(aes(x = variable, y = exp(val))) + 
  tidybayes::stat_halfeye() + 
  geom_hline(aes(yintercept = (1)), color = "blue", linetype = "solid", size = 1) + 
  geom_hline(aes(yintercept = 3), color = "red", linetype = "dashed", size = 1) + 
  geom_text(
    data = res %>%
      group_by(variable) %>%
      do(test(.$val)),
    aes(y = exp(max(res$val)) * 1.12, label = sprintf("Z = %s, P %s", Z, P))
  ) +
  theme_bw(base_size = 15) + 
  scale_x_discrete(label = c("Arrestment", "Selection")) + 
  labs(x = "Variable", y = "Maximum likelihood estimate", subtitle = "0% drop out")


ggpubr::ggarrange(g2, g, g6, g5) %>% ggpubr::annotate_figure(top = "blue = true arrestment, red = true selection")







