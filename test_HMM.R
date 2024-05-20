

ID <- 81
x <- fetch_events(ID) %>% 
  clean_events() %$% 
  move_seq(head_x, head_y) %>% 
  .$r %>% 
  na.omit()

plot_track_overlay(repID = ID)

x %>% 
  loghist(nclass = 50, log.p = FALSE, geom = "col", 
          draw_dist = list(fit_frechet(x), fit_invgamma(x), fit_lnorm(x), fit_gamma(x))
  )



d2 <- fetch_events(99) %>% 
  clean_events() %$% 
  move_seq(head_x, head_y) %>% 
  as.moveData(ID = "foo")

hmm_fit <- fit_HMM(d2)

x <- d2[viterbi(hmm_fit) == 2, "step"] %>% na.omit()

x %>% 
  loghist(nclass = 50, log.p = FALSE, geom = "col", 
          draw_dist = list(fit_lnorm(x), fit_gamma(x))
  )

comp_dist(x, dist = c("gamma", "lnorm"))


plot(hmm_fit, animals = 1)



