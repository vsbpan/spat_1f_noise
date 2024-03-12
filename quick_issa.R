source("spat1f/init_analysis.R")
library(survival)
library(amt)

# 
fit_list <- pb_par_lapply(
  unname(unlist(id_list)),
  function(i, ref_data){
    ID <- i
    d <- fetch_events(ID) %>%
      clean_events(ref_data = ref_data) %>%
      insert_gaps() %>% 
      mutate(
        head_x = ifelse(score >=0.9, head_x, NA),
        head_y = ifelse(score >=0.9, head_y, NA)
      ) %$%
      move_seq(head_x, head_y, r_thresh = 0, inherit.theta = FALSE) %>%
      filter(!is.na(r))

    if(length(d$r) < 10 | length(na.omit(d$theta_rel)) < 10){
      return(NULL)
    } else {
      d <- d %>%
        add_random_steps(n = 100L,
                         sl_distr = fit_gamma(.$r),
                         ta_distr = fit_genvonmises(.$theta_rel)
        ) %>%
        flag_invalid_steps(remove = TRUE) %>%
        mutate(
          start_toxic = ifelse(
            read_value(x1, y1, c(1000, 1000),
                       ref_img = fetch_trt_spec(ID, .ref_data = ref_data, quiet = TRUE)) == 1,
            "yes",
            "no"
          ),
          toxic = read_value(x2, y2, c(1000, 1000),
                             ref_img = fetch_trt_spec(ID, .ref_data = ref_data, quiet = TRUE))
        ) %>%
        append_estimators(na_as_zero = TRUE) %>% 
        mutate(
          time = step_id / 1000
        )
    }

    has_toxic <- !all(is.na(d$toxic))

    if(has_toxic){
      mod_form <- formula(
        case ~
          toxic +
          (cos_theta_pi + cos_2theta) +
          (sl + logsl) + 
          strata(step_id)
      )
    } else {
      mod_form <- formula(
        case ~
          (cos_theta_pi + cos_2theta) +
          (sl + logsl) +
          strata(step_id)
      )
    }

    out <- issf(
      mod_form,
      data = d,
      shape_estimator = c("logsl"),
      scale_estimator = c("sl"),
      kappa1_estimator = "cos_theta_pi",
      kappa2_estimator = "cos_2theta"
    )
    
    return(out)
  },
  ref_data = ref_data,
  cores = 8,
  inorder = TRUE
)
names(fit_list) <- unname(unlist(id_list))



#saveRDS(object = c(fit_list), "invisible/issf_fit_list.rds")

issf_fit_l <- readRDS("invisible/issf_fit_list.rds")

issf_fit_l <- issf_fit_l %>% 
  purrr:::keep(function(x){
  length(x) > 1
})


# extract_temporal_var <- function(issf_fit_l, hours = 12L){
#   v <- seq_along(issf_fit_l) %>% 
#     lapply(
#       function(i){
#         has_toxic <- !all(is.na(issf_fit_l[[i]]$data$toxic))
#         if(has_toxic){
#           l_conc <- ref_data %>% filter(rep_id == names(issf_fit_l)[i]) %>% .$low_diet_numeric
#           h_conc <- ref_data %>% filter(rep_id == names(issf_fit_l)[i]) %>% .$high_diet_numeric
#           
#           
#           out <- tryCatch(
#             issf_fit_l[[i]]$data %>% 
#               filter(case) %$% 
#               roll_vapply(toxic, w = 10 * hours + 1, FUN = function(xx){
#                 xx <- xx[!is.na(xx)]
#                 s <- xx == 0
#                 xx[s] <- l_conc
#                 xx[!s] <- h_conc
#                 
#                 var(xx, na.rm = TRUE)
#               }) %>% 
#               mean(na.rm = TRUE) ,
#             error = function(e){
#               NA
#             }
#           )
#         } else {
#           out <- 0
#         }
#         return(out)
#       }
#     ) %>% 
#     do.call("c", .)
#   
#   out <- data.frame(names(issf_fit_l), v)
#   names(out) <- c("rep_id", paste0("var_toxic_",hours))
#   
#   
#   return(out)
# }

# mix_means <- function(x){
#   x <- x[!is.na(x)]
#   cut <- cut_kmeans(x)
#   data.frame("lower_mean_cut" = mean(x[x<cut]), "cut" = cut, "upper_mean_cut" = mean(x[x>cut]))
# }







z <- issf_fit_l %>%
  lapply(function(x){
    c(
      coef(x$model),
      "shape" = x$sl_updated[[1]]$params$shape,
      "scale" = x$sl_updated[[1]]$params$scale,
      "kappa1" = x$ta_updated[[1]]$params$kappa1,
      "kappa2" = x$ta_updated[[1]]$params$kappa2)
  }) %>%
  lapply(function(x){
    as.data.frame(t(x))
  }) %>%
  do.call("rbind.fill", .) %>%
  cbind("rep_id" = names(issf_fit_l)) %>%
  cbind(
    lapply(
      issf_fit_l,
      function(x){
        has_toxic <- !all(is.na(x$data$toxic))

        if(has_toxic){
          on_toxic <- x$data %>%
            filter(case) %>%
            filter(!is.na(toxic)) %>%
            select(toxic) %>%
            colMeans()
          ava <- x$data %>% 
            filter(!case) %>% 
            filter(!is.na(toxic) & !is.na(start_toxic)) %>%
            group_by(step_id) %>% 
            summarise(
              "ava_mean_toxin" = mean(toxic),
              "ava_switch_toxin" = mean(
                ifelse(start_toxic == "yes", 
                       toxic == 0, 
                       toxic == 1)
              )
            ) %>% 
            dplyr::select(-step_id) %>% 
            summarise(
              "ava_mean_toxin" = mean(ava_mean_toxin, na.rm = TRUE),
              "ava_switch_toxin" = mean(ava_switch_toxin, na.rm = TRUE)
            ) %>% as.data.frame()
          
          out <- data.frame("on_toxic" = unname(on_toxic), ava)
          return(out)
        } else {
          return(
            data.frame("on_toxic" = NA, "ava_mean_toxin" = NA, "ava_switch_toxin" = NA)
          )
        }
      }
    ) %>%
      do.call("rbind", .),
    pb_par_lapply(issf_fit_l, function(x){
      x$data %>%
        filter(case) %$%
        ud_area(x2, y2) %>%
        t() %>%
        as.data.frame() %>%
        rename_all(
          function(x){
            paste0("ud_", x)
          }
        )
    }, cores = 1, inorder = TRUE, export_fun_only = TRUE) %>%
      do.call("rbind",.)
  )
z <- detection_report(z$rep_id) %>%
  select(repID, n_keypoints, frames) %>%
  mutate(
    prop_kp = n_keypoints / frames,
    rep_id = as.character(repID)
  ) %>%
  select(-repID) %>%
  right_join(z, by = "rep_id") 

w <- z %>%
  mutate(rep_id = as.character(rep_id)) %>%
  left_join(
    left_join(extract_temporal_var(issf_fit_l, 12), extract_temporal_var(issf_fit_l, 24))
  ) %>% 
  left_join(
    issf_fit_l %>%
      lapply(function(x){
        x$data %>%
          filter(case) %>%
          filter(!is.na(r)) %>%
          reframe(
            sl_mean_obs = mean(r),
            sl_kurt_obs = Kurt(r),
            prob_move_obs = mean(r > (sqrt(2) * 83)),
            n_valid_frames = sum(!is.na(r))
          )
      }) %>%
      do.call("rbind", .) %>%
      cbind("rep_id" = names(issf_fit_l)),
    by = "rep_id"
  )

w <- w %>% 
  mutate(
    n_keypoints  = n_valid_frames,
    prop_kp = n_valid_frames/frames
  )

# 
# 
# 
#write_csv(w, "cleaned_data/event_derivative.csv")



w %>%
  filter(scale > 0) %>% 
  #filter(prop_kp > 0.7) %>% 
  nrow()

w

w %>% 
  filter(scale > 0) %>% 
  filter(prop_kp > 0.5) %>% 
  ggplot(aes(x = prop_kp, y = (sl_mean_obs))) + 
  geom_point()




