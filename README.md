## cchfv-vaccine-trial-multicountry

A repository of code and data used for modelling transmission of CCHFV between livestock and spillover into humans in three endemic areas (Herat in Afghanistan, several provinces in Turkey and South Africa). The code includes routines for calibration, projection, and for analysis of vaccine trials in these settings, as presented in the published reference by Vesga and Metras 2022 et al ([here](https://www.medrxiv.org/content/10.1101/2022.06.09.22276201v1)). *Note: this is a working repository, this means code is likely to change over time*

**Quick start guide**

First, set local path in R to GitHub directory, e.g.: `setwd("~/Documents/GitHub/`cchfv-vaccine-trial-multicountry`")`

Calibration procedures can be run from `scripts/main_calibration_SA.r` , `scripts/main_calibration_AFG.r` and `scripts/main_calibration_TKY.r`

These scripts call several data loading scripts (e.g., `R/setup_data.r`, `R/setup_model.r`) and the core code for ordinary differential equations executed in Odin (e.g, `odin_model_vaccination_AFG.r`)

Other post-calibration analysis can be called from the R files in scripts folder, such as `scripts/main_trial_cohort.r` , or plotting results with `scripts/combine_plots_vaccines,r`