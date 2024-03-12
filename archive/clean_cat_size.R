source("spat1f/init_analysis.R")

d2 <- d %>% 
  mutate(
    repID = paste0("rep", repID)
  ) %>% 
  left_join(
    ref_data
  )

d2 <- d2 %>% 
  dplyr::select(-today) %>% 
  gather(
    key = period, value = slope, b1_500:b501_1000
  )

d2 <- d2 %>% 
  filter(
    slope < 0.0075
  )


d2 %>% 
  ggplot(
    aes(x = log(cat_pre_wt), y = slope, color = period)
  ) +
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(
    ~ cat_dead_cam_end + pupated_cam_end
  )


d2 %>% 
  ggplot(aes(x = var_trt, y = slope, group = period)) + 
  geom_pointrange(stat = "summary", position = position_dodge(width = 0.5)) + 
  geom_point(position = position_jitterdodge(
    jitter.width = 0.2, jitter.height = 0
  ), 
  aes(color = period)) + 
  facet_wrap(
    ~ cat_dead_cam_end + pupated_cam_end
  )


d <- data.frame(
  "repID" = 1:150, 
  "b1_500" = vapply(1:150, function(x){
    cat(sprintf("%s   \r", x))
    tryCatch(foo(x, nhours = 24, plot = FALSE, range = 1:500), 
             error = function(e){
               NA
             })
  }, numeric(1)),
  "b501_1000" = 
    vapply(1:150, function(x){
      cat(sprintf("%s   \r", x))
      tryCatch(foo(x, nhours = 24, plot = FALSE, range = 501:1000), 
               error = function(e){
                 NA
               })
    }, numeric(1))
)

debug(foo)


foo <- function(ID, nhours, plot = TRUE, range = 1:500){
  b <- fetch_data_dict(ID)
  #print(b)
  
  z <- ref_data %>% 
    filter(
      repID == paste0("rep", ID)
    ) %>% 
    dplyr::select(
      var_trt, beta, cat_pre_wt, cat_dead_cam_end, pupated_cam_end
    )
  
  test <- get_polygon(b) %>% 
    lapply(function(x){
      out_of_frame(x, tolerance = 10)
    }) %>% 
    do.call("c", .) %>% 
    unname()
  
  
  y <- fetch_events(ID)$size_px %>% 
    set_NA(index = test) %>% 
    roll_vapply(10 * nhours + 1, function(x){
      quantile(na.omit(x), probs = 0.9)
    })
  
  if(plot){
    plot(y, col = 
           ifelse(
             test, 
             "red",
             "black"
           ), main = sprintf(
             "var_trt: %s, beta: %s, cat_size: %s, dead: %s, pup: %s", 
             z$var_trt,
             z$beta, 
             ifelse(z$cat_pre_wt > 0.01, "big", "small"),
             z$cat_dead_cam_end,
             z$pupated_cam_end
           ))
  }
  
  y <- y[range]
  t <- seq_along(y)
  y2 <- y[y > 10]
  t <- t[y > 10]
  
  
  out <- tryCatch(
    coef(lm(log(y2) ~ t))[2], 
    error = function(e){
      NA
    }
  )
   
  
  if(sum(!(test | y == 0), na.rm = TRUE) < length(range) * 0.7){
    out <- NA
  }
  
  return(out)
}

set_NA <- function(x, index){
  x[index] <- NA
  x
}


















which_time(get_file_meta(b)$file_base, 14056)












