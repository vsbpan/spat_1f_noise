source("spat1f/init_analysis.R")
library(survival)
library(amt)

#### Estimate arrestment and selection with issa ####
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
          # Find the diet value at the origin location
          less_toxic_start = 1 - read_value(x1, y1, c(1000, 1000),
                             ref_img = fetch_trt_spec(ID, 
                                                      .ref_data = ref_data, 
                                                      quiet = TRUE)),
          less_toxic_end = 1 - read_value(x2, y2, c(1000, 1000),
                                     ref_img = fetch_trt_spec(ID, 
                                                              .ref_data = ref_data, 
                                                              quiet = TRUE)),
          less_toxic_start_f = as.factor(less_toxic_start)
        ) %>%
        append_estimators(na_as_zero = TRUE) # Append correct step length and turn angle estimators
    }
    
    has_toxic <- !all(is.na(d$less_toxic_start))
    
    if(has_toxic){
      mod_form <- formula(
        case ~ 
          state:less_toxic_end + 
          (less_toxic_start_f:state:cos_theta_pi + 
             less_toxic_start_f:state:cos_2theta) + # Turn angle update 
          (state:less_toxic_start_f:sl + 
             state:less_toxic_start_f:logsl) + # step length update
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
          c("state1", "state2"),
          c("less_toxic_start_f0", "less_toxic_start_f1")
        )
      ), 
      ta_estimators = pick_default_estimators(
        "genvonmises", 
        list(
          c("state1", "state2"),
          c("less_toxic_start_f0", "less_toxic_start_f1")
        )
      ), 
      keep_data = TRUE
    )
    
    return(out)
  },
  ref_data = ref_data,
  cores = 1,
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
    extract_model_coef(issf_fit_l[i]),
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
  )

# write_csv(res, "cleaned_data/event_derivative_joint_estimate.csv")


o <- read_csv("cleaned_data/event_derivative_joint_estimate.csv", progress = FALSE) %>% 
  rename_all(.funs = function(x){
    o <- gsub("state","s",gsub("estimate", "est", gsub(":","\\.", x)))
    vapply(o, function(z){
      z <- str_split_1(z, pattern = "__")
      z2 <- z
      z2[1] <- z[length(z)]
      z2[length(z)] <- z[1]
      paste0(z2, collapse = ".")
    }, FUN.VALUE = character(1))
  })

# write_csv(o, "cleaned_data/event_derivative_joint_estimate.csv")