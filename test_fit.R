source("spat1f/init_analysis.R")


sim_fun <- function(){
  spec <- syn_spec(n = 12,beta = 0, plot = FALSE)
  d <- iterate_random_steps2(
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
    start = make_start2(0,500,500, 1), 
    n = 2000, 
    ref_grid = spec, 
    same_move = FALSE, 
    rss_coef = 0
  )
  
  
  
  d %$% 
    move_seq(x, y) %>% 
    add_random_steps(n = 100, sl_distr = fit_gamma(.$r), ta_distr = fit_vonmises(.$theta_rel)) %>% 
    flag_invalid_steps(remove_invalid = TRUE) %>% 
    mutate(
      less_toxic_start = (round(read_value(x1, y1, ref_img = spec))),
      less_toxic_end = (round(read_value(x2, y2, ref_img = spec)))
    ) %>% 
    append_estimators() %>% 
    issf(
      case ~ 
        less_toxic_end + 
        as.factor(less_toxic_start):sl + 
        strata(step_id), 
      data = .
    ) -> m
  z <- unname(coef(m))
  k <- diff(log(unname(1/((1/m$sl$params$scale) - z[2:3]))))
  return(c("rss" = z[1], "k" = k))
}






out  <- pb_par_lapply(1:200, FUN = function(x, sim_fun){
  sim_fun()
}, sim_fun = sim_fun, 
cores = 8, 
inorder = FALSE)






