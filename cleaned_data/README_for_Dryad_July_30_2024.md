# Large Differences in Herbivore Performance Emerge from Simple Herbivore Behaviors and Fine-Scale Spatial Heterogeneity in Phytochemistry

## Table of Contents
1. [Overview](#Overview)
2. [Description](#Description)
3. [Data](#Data)



## Overview <a name="Overview"></a>

This repository hosts code and data used in the manuscript titled *Large Differences in Herbivore Performance Emerge from Simple Herbivore Behaviors and Fine-Scale Spatial Heterogeneity in Phytochemistry*. The version history of the code, intermediary data products, and demos can be found on GitHub (https://github.com/vsbpan/spat_1f_noise and https://github.com/vsbpan/herbivar). These repositories will be updated as the project progresses beyond the scope of this manuscript. Some of the important code and data files from these repositories at the conclusion of this present manuscript's analyses are archived here for reproducibility purposes. 


## Description <a name="Description"></a>

* `spat_1f_noise.zip`: Contains the raw data, intermediary data products (Mask-R-CNN predictions, movement tracks, etc), code for data processing, code for analysis, and a custom R package *spat1f*. See the `README.md` file in the zip folder for a description of the files in the folder. Dependent packages and their versions are reported in section 6 of the file `spat_1f_noise/statistical_analyses/spat1f_MS1_final_analysis.html`. You can also run `spat_1f_noise/Package_installation.R` to install package dependencies automatically. The code and data for the two supplemental experiments (conditioning experiment, Appendix 14 & constant toxin experiment, Appendix 13) are also included in this zipped folder. 
* `herbivar_0.2.1.tar.gz` A custom R package which the R package *spat1f* is dependent on. You can run `install.packages("herbivar_0.2.1.tar.gz", repos = NULL, type = "source")` to install this version of the package. 
* `ref_data_for_Dryad_July_11_2024.csv` is a spreadsheet of cleaned data used in the statistical analyses. 


## Data <a name="Data"></a>

* `ref_data_for_Dryad_July_11_2024.csv` Final experimental data (some of the variables were not analyzed in manuscript 1)
    * **rep_id**: caterpillar ID / replication ID
    * **beta**: the spectral exponent used to generate the treatment spectrum
    * **syn_id**: the treatment spectrum ID that maps to `raw_data/trt_spectra_meta/master_trt_meta.csv`
    * **cam_id**: camera ID
    * **arena_id**: arena ID
    * **mean_trt**: mean xanthotoxin dose in the diet 
    * **var_trt**: variation treatment
    * **low_diet**: xanthotoxin dose in the low dose diet
    * **high_diet**: xanthotoxin dose in the high dose diet
    * **assemble_date**: diet landscape assembly date in mm-dd-YYYY
    * **date_start**: date of experiment start (when the caterpillar is placed on treatment diet) in mm-dd-YYYY
    * **time_start**: time of experiment start in HH:MM
    * **cat_pre_wt**: caterpillar pre-experimental weight in grams
    * **date_end_camera**: date when the experiment ended in mm-dd-YYYY. 
    * **time_end_camera**: time when the experiment ended in HH:MM. 
    * **cat_post_wt**: caterpillar weight at the time when the experiment ended. 
    * **cat_dead_cam_end**: 1 (TRUE) or 0 (FALSE) the caterpillar died within five days. 
    * **pupation_date**: date when caterpillar begun to spin for eventual pupation in mm-dd-YYYY
    * **pupal_weight**: the weight of the melanized pupa in grams
    * **pupated_cam_end**: 1 (TRUE) or 0 (FALSE) the caterpillar begun to spin within five days.
    * **eclosure_date**: date when the moth emerged from the pupa in mm-dd-YYYY
    * **sex**: male or female moth
    * **death_date**: date when the caterpillar died in mm-dd-YYYY
    * **deformed_adult**: 1 (TRUE) or 0 (FALSE) whether the adult had deformed wings or body
    * **notes**: misc notes
    * **temperature**: temperature of the room at which the experiment was done in Celsius. 
    * **session_id**: experimental session ID
    * **error**: 1 (TRUE) or 0 (FALSE) a fatal error occurred for the replication that warrants its removal from analyses. Reason written for when TRUE. 
    * **premature_camera_end**: whether the camera stopped taking pictures before the end of the experiment. 
    * **camera_cutoff**: the time (seconds) limit after which the time lapse photos are no longer relevant (e.g. because the caterpillar died)
    * **surv_time**: how many days the caterpillar was seen alive
    * **RGR**: hourly relative growth rate
    * **pupated**: 1 (TRUE) or 0 (FALSE) the caterpillar pupated
    * **eclosed**: 1 (TRUE) or 0 (FALSE) the caterpillar eclosed
    * **time_to_pupation**: days to pupation
    * **time_to_eclosure**: days to eclosure
    * **adult_time**: days the adult was alive
    * **pupation_time**: days the caterpillar spent in the pupal stage
    * **repID**: another name for rep_id
    * **area_herb**: the number of pixels worth of diet the caterpillar consumed
    * **mean_toxic_conc**: the weighed toxin concentration consumed by the caterpillar in mg/g
    * **var_toxic_12**: the average variance in toxin concentration experienced by the caterpillar over 12 hours in mg^2/g^2.
    * **shape_***: the shape parameter of the step length kernel under different conditions (1: Exploration state & on more toxic diet; 2: Resting/feeding state & on more toxic diet; 3: Exploration state & on less toxic diet; 4: Resting/feeding state & on less toxic diet)
    * **scale_***: the scale parameter of the step length kernel under different conditions (1: Exploration state & on more toxic diet; 2: Resting/feeding state & on more toxic diet; 3: Exploration state & on less toxic diet; 4: Resting/feeding state & on less toxic diet)
    * **kappa1_***: the kappa1 parameter of the turn angle kernel under different conditions (1: Exploration state & on more toxic diet; 2: Resting/feeding state & on more toxic diet; 3: Exploration state & on less toxic diet; 4: Resting/feeding state & on less toxic diet)
    * **kappa2_***: the kappa2 parameter of the turn angle kernel under different conditions (1: Exploration state & on more toxic diet; 2: Resting/feeding state & on more toxic diet; 3: Exploration state & on less toxic diet; 4: Resting/feeding state & on less toxic diet)
    * **on_toxic**: the proportion of time the caterpillar spends on the more toxic diet. 
    * **sl_mean_obs**: the average Euclidean displacement of the caterpillar every six minutes in pixels
    * **ta_mean_obs**: the average turn angle of the caterpillar every six minutes in radians
    * **sl_kurt_obs**: the kurtosis of the Euclidean displacement of the caterpillar every six minutes in pixels
    * **ta_kurt_obs**: the kurtosis of the turn angle of the caterpillar every six minutes in radians
    * **prop_explore**: the proportion of time the caterpillar is observed exploring
    * **s_.less_toxic.est**: the maximum likelihood estimate of the relatove selection strength (RSS, aka. immigration) towards the less toxic diet compared to the more toxic diet (in log odds). (s1: exploration state; s2: feeding/resting state).
    * **s_.less_toxic_start.se**: the standard error of the maximum likelihood estimate of the relative selection strength (RSS, aka. immigration) towards the less toxic diet compared to the more toxic diet (in log odds). (s1: exploration state; s2: feeding/resting state).
    * **observed_dead**: 1 (TRUE) or 0 (FALSE) the caterpillar was physically observed to die.


    
    
    
    
