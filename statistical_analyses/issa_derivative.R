source("spat1f/init_analysis.R")
library(survival)
library(amt)

# Fit a bunch of issf and store them in a list
fit_list <- pb_par_lapply(
  unname(unlist(id_list))[!unname(unlist(id_list)) %in% problem_ids],
  function(i, ref_data){
    ID <- i
    d <- fetch_events(ID) %>% # Get mask-R-CNN instances
      clean_events(ref_data = ref_data) %$% # Clean those instance
      move_seq(head_x, head_y, r_thresh = 0, inherit.theta = FALSE) %>% # compute the step length and turn angle
      filter(!is.na(r)) # throw out steps where r is NA
    
    # Throw out trials with too few observations to be useful
    if(length(d$r) < 30 | length(na.omit(d$theta_rel)) < 30){
      return(NULL)
    } else {
      d <- d %>%
        add_random_steps(n = 100L, # Simulate random available steps
                         sl_distr = fit_gamma(.$r), # Fit gamma step length dist
                         ta_distr = fit_genvonmises(.$theta_rel) # Generalized von Mises turn angle dist
        ) %>%
        flag_invalid_steps(remove = TRUE) %>% # Throw out any random step that is outside of the arena
        mutate(
          # Find the diet value at the simulated destination location
          toxic = read_value(x2, y2, c(1000, 1000),
                             ref_img = fetch_trt_spec(ID, 
                                                      .ref_data = ref_data, 
                                                      quiet = TRUE)), 
          less_toxic = 1 - toxic
        ) %>%
        append_estimators(na_as_zero = TRUE) # Append correct step length and turn angle estimators
    }
    
    has_toxic <- !all(is.na(d$toxic))
    
    if(has_toxic){
      mod_form <- formula(
        case ~
          less_toxic + # Habitat selection estimation
          (cos_theta_pi + cos_2theta) + # Turn angle update 
          (sl + logsl) + # step length update
          strata(step_id) # Stratify be step ID
      )
    } else {
      mod_form <- formula(
        case ~
          (cos_theta_pi + cos_2theta) +
          (sl + logsl) +
          strata(step_id)
      )
    }
    
    out <- issf( # Fit the issf
      mod_form,
      data = d,
      shape_estimator = c("logsl"),
      scale_estimator = c("sl"),
      kappa1_estimator = "cos_theta_pi",
      kappa2_estimator = "cos_2theta", 
      mulog_estimator = "logsl", 
      sdlog_estimator = "logslsq"
    )
    
    return(out)
  },
  ref_data = ref_data,
  cores = 8,
  inorder = TRUE
)
names(fit_list) <- unname(unlist(id_list))[!unname(unlist(id_list)) %in% problem_ids]

fit_list

saveRDS(object = c(fit_list), "invisible/issf_fit_list.rds")
rm("fit_list")

issf_fit_l <- readRDS("invisible/issf_fit_list.rds")

issf_fit_l <- issf_fit_l %>% 
  purrr:::keep(function(x){
    length(x) > 1 # filter out NULL
  })





i <- seq_along(issf_fit_l)

# Extract data from fitted models
res <- extract_ava_neighborhood_quality(issf_fit_l[i], 
                             .ref_data = ref_data, 
                             n = 100) %>% 
  left_join(
    extract_temporal_var(names(issf_fit_l)[i], 
                         hours = 12L, 
                         .ref_data = ref_data, 
                         cores = 6),
    by = "rep_id"
  ) %>% 
  left_join(
    exract_model_coef(issf_fit_l[i]),
    by = "rep_id"
  ) %>% 
  left_join(
    extract_mean_on_toxic(names(issf_fit_l)[i], .ref_data = ref_data, cores = 6),
    by = "rep_id"
  ) %>% 
  left_join(
    extract_obs_move_summary(names(issf_fit_l)[i], cores = 6),
    by = "rep_id"
  )

# write_csv(res, "cleaned_data/event_derivative.csv")

