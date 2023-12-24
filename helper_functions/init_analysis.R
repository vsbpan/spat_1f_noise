library(ggpubr)
library(sjPlot)
library(performance)
library(tidybayes)
library(glmmTMB)
#source("helper_functions/init.R")



ref_data <- read_csv("raw_data/1_f_noise_experiment data_Dec_24_2023.csv")

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
    today = as.POSIXlt(Sys.Date(), format = "%Y-%m-%d")
  ) %>% 
  mutate(
    exp_time_elapsed = as.numeric(as.difftime(exp_end - exp_start, units = "h")),
    diet_age = as.numeric(as.difftime(date_start - assemble_date, units = "h")),
    surv_time = ifelse(
      is.na(death_date),
      ifelse(
        error == 0, 
        NA,
        as.numeric(as.difftime(today - date_start, units = "h"))
      ),
      as.numeric(as.difftime(death_date - date_start, units = "h"))
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
    pupation_time = ifelse(
      is.na(pupated),
      NA,
      ifelse(
        pupated == 1, 
        as.numeric(as.difftime(pupation_date - date_start, units = "d")),
        as.numeric(as.difftime(today - date_start, units = "d"))
      )
    ), 
    eclosure_time = ifelse(
      is.na(eclosed),
      NA,
      ifelse(
        pupated == 1, 
        as.numeric(as.difftime(eclosure_date - date_start, units = "d")),
        as.numeric(as.difftime(today - date_start, units = "d"))
      )
    )
  ) %>% 
  mutate(
    beta = as.factor(beta), 
    rep_id = as.character(rep_id),
    session_id = as.character(session_id),
    repID = paste0("rep",rep_id)
  )




