# 1/f Noise Project

## Table of Contents
1. [Overview](#Overview)
2. [Installation](#Installation)
3. [Description](#Description)
    * [Code library](#CodeLibrary)
    * [File structure](#FileStructure)
    * [Data](#Data)



## Overview <a name="Overview"></a>

This repo hosts code for image processing and data analysis for an INTERN funded project, examining the consequences of spatial autocorrelation and variance in fine-scale plant toxin distribution on caterpillar outcomes. 


![**Manuscript 1 methods figure**](demo/methods_figure.jpg)

Summary of the data processing pipeline that derives measures of herbivore movement and feeding behavior that may explain differential herbivore performance. Lighter squares represent diet tiles that contain a higher concentration of xanthotoxins. A layer of red luster dust was sprayed on top of the diet to delineate eaten and uneaten diet. 




## Installation <a name="Installation"></a>

Run `Package_installation.R` to install repository dependencies. For the dependency *herbivar*, see  [herbivar](https://github.com/vsbpan/herbivar) GitHub page on installation instructions. 


## Description <a name="Description"></a>

### Code library <a name="CodeLibrary"></a>

The custom code written for this project are bundled as a pseudo simulated package *spat1f*. It can be imported using `source(spat1f/init.R)`, or `soure(spat1f/init_analysis.R)` to get the cleaned experimental metadata in the global environment as well.  


* `spat1f/` The root directory of the *spat1f* package. `.R` files in this directory are not loaded in the package. 
    * `spat1f/R/` Contains the R source code 
    * `spat1f/src/` Contains the C++ source code
    * `DESCRIPTION` Package description file
    * `NAMESPACE` Namespace file

* Key objects
    * `cimg`: RGB image tensor array from *imager* that represents an image in R. Can be used with functions from *imager* and some functions in *spat1f*.
    * `pixset`: Binary image tensor array from *imager* that represents binary masks in R. Can be used with functions from *imager* and some functions in *spat1f*. 
    * `imlist`: A list of `cimg` or `pixset` objects. 
    * `data_dict`: Various Mask-R-CNN outputs and metadata for a collection of images. Has functions for importing, conversion to `COCO_Json` and `events`, data extraction, summary, subsetting, merging, visualization, evaluation, and validation. 
    * `COCO_Json`: A list of data.frames and lists that represents a COCO annotation file. Has functions for importing, exporting, merging, randomization, subsetting, and summary. 
    * `issf_fit`: Fitted integrated step selection function. Has functions for summary and step simulation.
    * `events`: Not an actual class, but used throughout. A data.frame of location and associated metadata of each caterpillar at each time point. Has functions for validation, data extraction, summary, and visualization. Usually paired with `clean_events()` to keep only valid caterpillar identifications. 
    * `amt_distr` No S3 generics, but used throughout. A list of parameter values and metadata for a parametric distribution. Has functions for creation, update, finding the density, random number generation, and fitting.
    * `moveHMM` A class of model object from the *moveHMM* package. Some methods are exported to *spat1f* so that user won't need to load the *moveHMM* library, which conflicts with *tidyverse*. 


### File structure <a name="FileStructure"></a>

* `annotations/` Houses `.json` files of image annotations in COCO format. This is mostly used for back up purposes. 
    * `coco_test.json` COCO annotation file for testing data
    * `coco_train.json` COCO annotation file for training data
    * `coco_val.json` COCO annotation file for validation data
    * `master_sample_manual_annotation.json` Full manually validated COCO annotation file
    * `sample1_manual_annotation_updated.json` Batch 1 of manually validated COCO annotation file
    * `sample2_manual_annotation.json` Batch 2 of manually validated COCO annotation file

* `archive/` Houses misc scripts that are no longer in use or is there only as a backup. 

* `cleaned_data/` Houses results of expensive computations. 
    * `cleaned_data/events/` Houses `.csv` files on the location and associated metadata at each time point for each caterpillar. Called by `fetch_events()`. Named after each trial ID. 
    * `cleaned_data/data_dicts/` Houses `data_dict` class objects that represents a cleaned and formatted version of the Mask-R-CNN output of each frame. `data_dict` objects are at the center of all data analyses and are supported by various S3 methods. The `.rds` files are called by `fetch_data_dict()`. Named after each trial ID. 
    * `consumption_mask_derivative.csv` Results from end of experiment herbivory masks. 
    * `event_derivative.csv` Results from movement tracks and ISSF fits. 
    * `SEM_sim*.csv` Simulation results from `SEM_pred_coef()` for each bootstrap refit of the data. Basically the results of total effects through certain paths or after removing certain nodes. 
    * `sus_frames_list.rds` A list of frames for each rep for which they are flagged as having false head movements. 

* `data_processing/` Houses various scripts used to complete certain computer vision / data processing goals. 
    * `crop_raw_image_auto.R` Read saved anchors and for those images that have yet to be cropped, crop them. 
    * `data_dict_to_COCO.R` Write `data_dict` objects as `.json` files in COCO format. Compatible with *COCO annotator* and *detectron2*.  
    * `detect_herbivory.R` Compute binary mask and associated summary statistics of herbivory from processed images. 
    * `flag_inference_errors.R` Flag potential errors in Mask-R-CNN predictions.  
    * `generate_spec.R` Generate treatment spectra and save them on the computer for later reference.  
    * `make_video.R` Overlay treatment spectra and/or Mask-R-CNN predictions on processed image stack and stitch together as a video using *FFmpeg*. 
    * `inference_processing.R` Turn raw inference from *detectron2* into `data_dict` objects and events. Flag implausible mask size and repeated long-distance jumps. 
    * `time_lapse_crop.R` Launches the anchor picker app, saves anchors, and crop selected sets of raw images. 
    * `time_lapse_export.R` Code for exporting and renaming image files from USB drive containing raw image files from raspberry pi. 

* `statistical_analyses/`
    * `SEM_simulation.R` Bootstrap refit of SEM and subsequent calculation of the total effects of selected paths. 
    * `manuscript_plots.R` Code for generating some plots in manuscript 1.
    * `prediction_evaluator.R` Code for evaluating the performance of Mask-R-CNN, computing sample size, and error flagging rate  
    * `conditional_experiment.R` Some preliminary analysis code for a conditioning experiment not in manuscript 1
    * `constant_trt_supplemental_analysis.R` Code for analyses and graph in appendix 10 of manuscript 1 (constant toxin concentration experiment).
    * `issa_derivative.R` Fit ISSF to each movement track then extract key proximal variables from the tracks
    * `spat1f_MS1_final_analysis.html` Knitted output of `spat1f_MS1_final_analysis.Rmd`
    * `spat1f_MS1_final_analysis.Rmd` Code for performing the main analyses included in manuscript 1

* `demo/` Houses images for populating this README file. 

* `deep_learning/` Houses code for performing deep learning
    * `detectron2/` Houses some modified skeleton code of detectron2 used to train the Mask-R-CNN
    * `modelv7_config.ipynb` Jupyter notebook of code and model configurations used to train, predict from, and format the output from the Mask-R-CNN

* `graphs/` Houses graphs generated for this project. Due to the folder size, it is not uploaded to GitHub. 

* `misc_tests/` Houses prototype images used to develop the image processing capabilities of this repository. Now defunct. 

* `notes/` Houses some messy notes on the rationale and details on some of the analysis / processing done. For internal use only. 

* `processed_feed/` Houses sub-directories named after the trial ID (e.g. `rep50/`) containing the output of cropped time lapse images. The naming convention is as such: `processed_repID_camID_sTimeInSeconds_rankN.jpg` e.g. `processed_rep3_cam3_s4325_rank13.jpg`. This directory is not uploaded to GitHub due to its size. 

* `raw_data/` Houses raw data and some intermediary data products. 
    * `raw_data/inferences/` Houses `.csv` outputs from Mask R-CNN predictions generated in python. Named after each trial ID.  
    * `raw_data/picked_anchors/` Houses lists of lists of anchors picked from the shiny app (`anchor_picker_app()`) used for perspective correction and cropping of the raw time lapse images. Named after each trial ID
    * `raw_data/spectra_mta/`Houses treatment spectra metadata in `.csv` format generated from `data_processing/generate_spec.R`
    * `master_modelv7_1_7999iter_test_inference.csv` Mask-R-CNN inference on testing dataset along with ground truths
    * `master_modelv7_1_7999iter_train_inference.csv` Mask-R-CNN inference on training dataset along with ground truths
    * `1_f_noise_experiment data_DATE.csv` Experimental metadata, including herbivore performance outcomes saved at different dates. The file with the latest date is the one used in all analyses. 
    * `T ni conditioning study.csv` *T. ni* conditioning experiment data. Not part of manuscript 1.

* `time_lapse_feed/` Houses sub-directories named after the trial ID (e.g. `rep50/`) containing raw time lapse images. The naming convention is as such: `repID_camID_sTimeInSeconds.jpg` e.g. `rep3_cam3_s4325.jpg`. This directory is not uploaded to GitHub due to its size. 

* `video/` Houses videos generated by *FFmpeg* that stitches together the time lapse image stack. This directory is not uploaded to GitHub due to its size. 

* `invisible/` Houses various intermediary data products too large to be uploaded to GitHub. 

* `Package_installation.R` Nice code for helping you install dependencies of *spat1f*

* `README.md` This README documentation


### Data <a name="Data"></a>


* `raw_data/1_f_noise_experiment data_May_03_2024.csv` Final experimental data (some of the variables were not analyzed in manuscript 1)
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
    * **camera_cutoff**: the time (seconds) limit after which the time lapse photos are no longer relevant (e.g. because the caterpillar died).


* `raw_data/trt_spectra_meta/master_trt_meta.csv` treatment spectra metadata
    * **syn_id**	the treatment spectrum ID
    * **syn_date** system time when the spectrum was generated
    * **spec_A1 -- spec_L12** whether the spectrum should have a high (1) or low (0) xanthotoxin diet at the corresponding coordinate. 

<img src="demo/rep_id49.jpg" alt="Example of printed template for diet landscape assembly" style="width:400px"/>


* `raw_data/T ni conditioning study.csv` *T. ni* conditioning experiment data. Not part of manuscript 1.
    * **Trial** Trial ID
    * **cat_wt** Caterpillar weight at the start of the choice assay in grams
    * **condition_trt** Was the cat reared on 0 mg/g (L) or 2 mg/g (H) concentration diet for 7 days? (H/L)
    * **center_block_trt** Whether there is a center block with 2 mg/g xanthotoxin at which the cat starts (Y/N)
    * **first_diet**	Is the first chosen diet a high or low diet (2 mg/g or 0 mg / g)? NA if no first choice is made (H/L/NA)
    * **choice_2h**	choice by 2 hours. 1 yes 0 no
    * **second_diet**	Is the second chosen diet a high or low diet (2 mg/g or 0 mg / g)? NA if no second choice is made (H/L/NA)
    * **end_of_day_diet**	What diet is the caterpillar on at the end of day (H = high, low = low, c = center, n = matrix)
    * **next_morning_choice** What diet is the caterpillar on next day (H = high, low = low, c = center, n = matrix)	
    * **notes** Notes
    


