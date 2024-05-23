source("spat1f/init_analysis.R")
library(survival)
library(amt)

#### Estimate arrestment with issa ####
# Fit a bunch of issf and store them in a list
issf_fit_l <- pb_par_lapply(
  unname(unlist(id_list$var))[!unname(unlist(id_list$var)) %in% problem_ids],
  function(i, ref_data){
    ID <- i
    d <- fetch_events(ID) %>% # Get mask-R-CNN instances
      clean_events(ref_data = ref_data) %$% # Clean those instance
      move_seq(head_x, head_y, r_thresh = 0, inherit.theta = FALSE) # compute the step length and turn angle
      
    
    # Throw out trials with too few observations to be useful
    if(length(na.omit(d$r)) < 30 | length(na.omit(d$theta_rel)) < 30){
      return(NULL)
    } else {
      # Append state column inferred from a fitted Hidden Markov Model
      d$state <- as.factor(viterbi(fit_HMM(as.moveData(d))))
      
      # Now enforce that each state must have >= 30 observations
      count <- d %>% 
        group_by(state) %>% 
        summarise(count = sum(!is.na(theta_rel))) %>% 
        .$count
      
      if(any(count < 30)){
        return(NULL)
      } 
      

      # Prepare data for issa
      d <- d %>%
        filter(!is.na(r)) %>%  # throw out steps where r is NA
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
          less_toxic = 1 - toxic,
          less_toxic_f = as.factor(less_toxic)
        ) %>%
        append_estimators(na_as_zero = TRUE) # Append correct step length and turn angle estimators
    }
    
    has_toxic <- !all(is.na(d$toxic))
    
    if(has_toxic){
      mod_form <- formula(
        case ~ 
          (less_toxic_f:state:cos_theta_pi + less_toxic_f:state:cos_2theta) + # Turn angle update 
          (state:less_toxic_f:sl + state:less_toxic_f:logsl) + # step length update
          strata(step_id) # Stratify be step ID
      )
    } else {
      mod_form <- formula(
        case ~
          (state:cos_theta_pi + state:cos_2theta) +
          (state:sl + state:logsl) +
          strata(step_id)
      )
    }
    
    out <- issf( # Fit the issf
      mod_form,
      data = d,
      sl_estimators = pick_default_estimators(
        "gamma", 
        list(
          c("less_toxic_f0", "less_toxic_f1"),
          c("state1", "state2")
        )
      ), 
      ta_estimators = pick_default_estimators(
        "genvonmises", 
        list(
          c("less_toxic_f0", "less_toxic_f1"),
          c("state1", "state2")
        )
      ), 
      keep_data = TRUE
    )
    
    return(out)
  },
  ref_data = ref_data,
  cores = 8,
  inorder = TRUE
)
names(issf_fit_l) <- unname(unlist(id_list$var))[!unname(unlist(id_list$var)) %in% problem_ids]


issf_fit_l <- issf_fit_l %>% 
  purrr:::keep(function(x){
    length(x) > 1 # filter out NULL
  })



i <- seq_along(issf_fit_l)

# Extract data from fitted models
res <- extract_temporal_var(names(issf_fit_l)[i], 
                         hours = 12L, 
                         .ref_data = ref_data, 
                         cores = 6) %>% 
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
  ) %>% 
  left_join(
    extract_prop_state1(names(issf_fit_l)[i], cores = 6),
    by = "rep_id"
  ) %>% 
  left_join(
    extract_ava_neighborhood_quality(issf_fit_l[i],
                                     .ref_data = ref_data,
                                     n = 100,
                                     cores = 6), 
    by = "rep_id"
  )

# write_csv(res, "cleaned_data/event_derivative_arrestment.csv")





#### Estimate selection strength in a different model ####
# Coefficients are not identifiable in a lot of reps if selection and arrestments are estimated in the same model. 

# Fit a bunch of issf and store them in a list
issf_fit_l <- pb_par_lapply(
  unname(unlist(id_list$var))[!unname(unlist(id_list$var)) %in% problem_ids],
  function(i, ref_data){
    ID <- i
    d <- fetch_events(ID) %>% # Get mask-R-CNN instances
      clean_events(ref_data = ref_data) %$% # Clean those instance
      move_seq(head_x, head_y, r_thresh = 0, inherit.theta = FALSE) # compute the step length and turn angle
    
    
    # Throw out trials with too few observations to be useful
    if(length(na.omit(d$r)) < 30 | length(na.omit(d$theta_rel)) < 30){
      return(NULL)
    } else {
      # Append state column inferred from a fitted Hidden Markov Model
      d$state <- as.factor(viterbi(fit_HMM(as.moveData(d))))
      
      # Now enforce that each state must have >= 30 observations
      count <- d %>% 
        group_by(state) %>% 
        summarise(count = sum(!is.na(theta_rel))) %>% 
        .$count
      
      if(any(count < 30)){
        return(NULL)
      } 
      
      
      # Prepare data for issa
      d <- d %>%
        filter(!is.na(r)) %>%  # throw out steps where r is NA
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
          less_toxic = 1 - toxic,
          less_toxic_f = as.factor(less_toxic)
        ) %>%
        append_estimators(na_as_zero = TRUE) # Append correct step length and turn angle estimators
    }
    
    has_toxic <- !all(is.na(d$toxic))
    
    if(has_toxic){
      mod_form <- formula(
        case ~ 
          (less_toxic_f:state:cos_theta_pi + less_toxic_f:state:cos_2theta) + # Turn angle update 
          (state:less_toxic_f:sl + state:less_toxic_f:logsl) + # step length update
          strata(step_id) # Stratify be step ID
      )
    } else {
      mod_form <- formula(
        case ~
          (state:cos_theta_pi + state:cos_2theta) +
          (state:sl + state:logsl) +
          strata(step_id)
      )
    }
    
    out <- issf( # Fit the issf
      mod_form,
      data = d,
      sl_estimators = pick_default_estimators(
        "gamma", 
        list(
          c("less_toxic_f0", "less_toxic_f1"),
          c("state1", "state2")
        )
      ), 
      ta_estimators = pick_default_estimators(
        "genvonmises", 
        list(
          c("less_toxic_f0", "less_toxic_f1"),
          c("state1", "state2")
        )
      ), 
      keep_data = TRUE
    )
    
    return(out)
  },
  ref_data = ref_data,
  cores = 8,
  inorder = TRUE
)
names(issf_fit_l) <- unname(unlist(id_list$var))[!unname(unlist(id_list$var)) %in% problem_ids]


issf_fit_l <- issf_fit_l %>% 
  purrr:::keep(function(x){
    length(x) > 1 # filter out NULL
  })



i <- seq_along(issf_fit_l)

# Extract data from fitted models
res <- extract_temporal_var(names(issf_fit_l)[i], 
                            hours = 12L, 
                            .ref_data = ref_data, 
                            cores = 6) %>% 
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
  ) %>% 
  left_join(
    extract_prop_state1(names(issf_fit_l)[i], cores = 6),
    by = "rep_id"
  ) %>% 
  left_join(
    extract_ava_neighborhood_quality(issf_fit_l[i],
                                     .ref_data = ref_data,
                                     n = 100,
                                     cores = 6), 
    by = "rep_id"
  )

# write_csv(res, "cleaned_data/event_derivative_select.csv")


#### Merge data as master file ####




o <- read_csv("cleaned_data/event_derivative_arrestment.csv", progress = FALSE) %>% 
  rename_all(.funs = function(x){
    o <- gsub("state","s",gsub("estimate", "est", gsub(":","\\.", x)))
    vapply(o, function(z){
      z <- str_split_1(z, pattern = "__")
      z2 <- z
      z2[1] <- z[length(z)]
      z2[length(z)] <- z[1]
      paste0(z2, collapse = ".")
    }, FUN.VALUE = character(1))
  }) %>% 
  left_join(
    read_csv("cleaned_data/event_derivative_select.csv", progress = FALSE) %>% 
      rename_all(.funs = function(x){
        o <- gsub("state","s",gsub("estimate", "est", gsub(":","\\.", x)))
        vapply(o, function(z){
          z <- str_split_1(z, pattern = "__")
          z2 <- z
          z2[1] <- z[length(z)]
          z2[length(z)] <- z[1]
          paste0(z2, collapse = ".")
        }, FUN.VALUE = character(1))
      }) %>% 
      mutate(
        sl_mean_pred1 = shape_1 * scale_1,
        sl_mean_pred2 = shape_2 * scale_2
      ) %>% 
      filter(
        sl_mean_pred2 < 2000) %>% 
      filter_at(vars(contains("shape_"), contains("scale_[0-9]")), 
                all_vars(. > 0)) %>% 
      filter(
        s1.less_toxic.se < 10 &
          s2.less_toxic.se < 10
      ) %>% 
      dplyr::select(rep_id, contains("less_toxic")),
    by = "rep_id"
  )

# write_csv(o, "cleaned_data/event_derivative_master.csv")