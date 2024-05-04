library(ggpubr)
library(sjPlot)
library(performance)
library(tidybayes)
library(glmmTMB)
source("spat1f/init.R")



ref_data <- suppressMessages(read_csv("raw_data/1_f_noise_experiment data_May_03_2024.csv"))
today <- as.character("2024-01-31")#Sys.Date()
problem_ids <- c("2","6","7", "9", "10","11", "20","28", "85", "113", "146")


ref_data <- ref_data %>% 
  filter(!is.na(cat_dead_cam_end)) %>% 
  mutate(
    exp_start = as.POSIXct(paste(date_start, time_start, sep = " "), format = "%m-%d-%Y %H:%M"),
    exp_end = as.POSIXct(paste(date_end_camera, time_end_camera, sep = " "), format = "%m-%d-%Y %H:%M"),
    pupation_date = as.POSIXct(pupation_date, format = "%m-%d-%Y"),
    death_date = as.POSIXct(death_date, format = "%m-%d-%Y"),
    eclosure_date = as.POSIXct(eclosure_date, format = "%m-%d-%Y"),
    assemble_date = as.POSIXct(assemble_date, format = "%m-%d-%Y"),
    date_start = as.POSIXct(date_start, format = "%m-%d-%Y"),
    date_end_camera = as.POSIXct(date_end_camera, format = "%m-%d-%Y"),
    today = as.POSIXlt(today, format = "%Y-%m-%d")
  ) %>% 
  mutate(
    exp_time_elapsed = as.numeric(as.difftime(exp_end - exp_start, units = "h")),
    diet_age = as.numeric(as.difftime(date_start - assemble_date, units = "h")),
    surv_time = ifelse(
      is.na(death_date),
      ifelse(
        error == 0, 
        NA,
        as.numeric(as.difftime(today - date_start, units = "d"))
      ),
      as.numeric(as.difftime(death_date - date_start, units = "d"))
    )
  ) %>% 
  mutate(
    RGR = log(cat_post_wt/cat_pre_wt) / exp_time_elapsed,
    dead = ifelse(
      is.na(death_date),
      ifelse(
        error == 0,
        0,
        NA
      ),
      1
    ),
    pupated = ifelse(
      is.na(pupation_date),
      ifelse(
        is.na(death_date),
        NA,
        0
      ),
      1
    ),
    eclosed = ifelse(
      is.na(eclosure_date),
      ifelse(
        is.na(death_date),
        NA,
        0
      ),
      1
    )
  ) %>% 
  mutate(
    time_to_pupation = ifelse(
      is.na(pupated),
      NA,
      ifelse(
        pupated == 1, 
        as.numeric(as.difftime(pupation_date - date_start, units = "d")),
        ifelse(
          dead == 1,
          NA,
          as.numeric(as.difftime(today - date_start, units = "d"))
        )
      )
    ), 
    time_to_eclosure = ifelse(
      is.na(eclosed),
      NA,
      ifelse(
        pupated == 1, 
        as.numeric(as.difftime(eclosure_date - date_start, units = "d")),
        ifelse(
          dead == 1,
          NA,
          as.numeric(as.difftime(today - date_start, units = "d"))
        )
      )
    ), 
    adult_time = ifelse(
      is.na(eclosed), 
      NA, 
      ifelse(
        eclosed == 0, 
        NA, # or 0
        as.numeric(as.difftime(death_date - eclosure_date, units = "d"))
      )
    ),
    pupation_time = ifelse(
      is.na(pupated) | pupated == 0, 
      NA, 
      ifelse(
        eclosed == 0, 
        NA, # or 0
        as.numeric(as.difftime(eclosure_date - pupation_date, units = "d"))
      )
    )
  ) %>% 
  mutate(
    beta = as.factor(ifelse(is.na(beta), "const", beta)), 
    rep_id = as.character(rep_id),
    session_id = gsub("\\.[0-9]","", as.character(session_id)),
    repID = paste0("rep",rep_id)
  ) %>% 
  left_join(
    read_csv("cleaned_data/consumption_mask_derivative.csv") %>% 
      suppressMessages() %>% 
      mutate(rep_id = as.character(rep_id)), 
    by = "rep_id" 
  ) %>% 
  left_join(
    read_csv("cleaned_data/event_derivative.csv") %>% 
      suppressMessages() %>% 
      mutate(
        rep_id = as.character(rep_id)
        ) %>% 
      filter(
        !rep_id %in% problem_ids
      ),
    by = "rep_id" 
    ) %>% 
  mutate(
    mean_trt_numeric = parse_conc(mean_trt),
    low_diet_numeric = parse_conc(low_diet),
    high_diet_numeric = parse_conc(high_diet)
  ) %>% 
  mutate(
    mean_toxic_conc = ifelse(
      var_trt == "constant",
      mean_trt_numeric,
      (mean_toxic) * high_diet_numeric + (1 - mean_toxic) * low_diet_numeric
    ),
    on_toxic_conc = ifelse(
      var_trt == "constant",
      mean_trt_numeric,
      (on_toxic) * high_diet_numeric + (1 - on_toxic) * low_diet_numeric
    )
  ) %>% 
  mutate(
    observed_dead = ifelse(pupated == 1 & eclosed == 0, 0, 1), # !Cats that died at pupal stage
    surv_time = ifelse(observed_dead == 1, 
                       surv_time, 
                       # censor at pupation date bc not sure when they died. 
                       as.numeric(as.difftime(pupation_date - date_start, units = "d"))
    )
  ) 



ref_data


id_list <- list("var" = ref_data %>% 
       filter(!is.na(camera_cutoff)) %>% 
       filter(var_trt != "constant") %>% 
       dplyr::select(rep_id) %>% 
       unlist() %>% 
       unname(), 
     "const" = ref_data %>% 
       filter(!is.na(camera_cutoff)) %>% 
       filter(var_trt == "constant") %>% 
       dplyr::select(rep_id) %>% 
       unlist() %>% 
       unname())

cat(sprintf("Data last updated on %s", today))
rm("today")

