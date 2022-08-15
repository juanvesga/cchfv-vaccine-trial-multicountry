# 12/08/2022
# JF Vesga - London School of Hygiene and Tropical Medicine
# This script plots model fits to data as seen in manuscript appendix
# Figs S1-S3

graphics.off()
rm(list = ls())
# Load pacman package and necessary packages
library("pacman")
p_load(
  Hmisc, reshape2, zoo, plyr, ggplot2, gmodels, data.table, devtools, fitR,
  coda, mvtnorm, ISOweek, lubridate, dplyr, zoo, stringr, here, psych, splines,
  data.table, broom, psych, matrixStats, profvis, odin, FreqProf
)


# Part 0. select country -------------------------------------------------------
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

country <- "TKY" # country : AFG; TKY; SA
driver <- "sat_def" # Envir drivers: sat_def; soil_t; rel_hum; ndvi
runs <- "saved" # saved or simulate
samps <- 500 # Samples on posterior distribution


# Part 1. Import set up: Load data and model parameters and functions ---------
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (country == "AFG") {
  source(here("R", "odin_models_vaccination_AFG.R"))
  source(here("R", "setup_data_AFG.R"))
  source(here("R", "setup_model_AFG.R"))
  fig_file <- c("FigS1.tiff")
} else

if (country == "SA") {
  source(here("R", "odin_models_vaccination_SA.R"))
  source(here("R", "setup_data_SA.R"))
  source(here("R", "setup_model_SA.R"))
  fig_file <- c("FigS2.tiff")
} else

if (country == "TKY") {
  source(here("R", "odin_models_vaccination_TKY.R"))
  source(here("R", "setup_data_TKY.R"))
  source(here("R", "setup_model_TKY.R"))
  fig_file <- c("FigS3.tiff")
}

source(here("R", "plot_fits_function.R"))

# Load posterior samples
if (country == "AFG") {
  posteriors <- read.table(
    here("output", country, paste("posterior_", driver, ".txt", sep = "")),
    sep = "\t"
  )

  nsim <- samps

  # Sample
  rand <- posteriors[sample(nrow(posteriors), nsim), ]

  # Pass parameters
  theta <- rand %>%
    select(
      A,
      D_lact_liv,
      D_imm_liv,
      F_risk,
      O_factor,
      clin_frac
    )
} else {
  posteriors <- read.table(
    here("output", country, paste("posterior_", driver, ".txt", sep = "")),
    sep = "\t"
  )

  nsim <- samps

  # Sample the posterior
  rand <- posteriors[sample(nrow(posteriors), nsim), ]

  # Pass parameters
  theta <- rand %>%
    select(
      A,
      D_lact_liv,
      D_inf_L,
      D_imm_liv,
      F_risk,
      O_factor,
      RRreport,
      clin_frac
    )
}

########################################################################################
# Part 2. Run model instances from sampled paarmeters or load pre-saved results
########################################################################################

if (runs == "simulate") {
  sim.sir <- get_output_sir(params, theta)

  params$prev_inf_liv_all <- sim.sir$prev_inf_liv

  sim.seir <- get_output_seir(params, theta, 100)

  sim <- c(sim.sir, sim.seir)

  saveRDS(sim, file = here("output", country, paste("simulations_", driver, ".rds", sep = "")))
} else {
  sim <- readRDS(here("output", country, paste("simulations_", driver, ".rds", sep = "")))
}

########################################################################################
# Part 3. Plot
########################################################################################

plot_fits_function(sim, observations, country)


########################################################################################
# Part 4. Save plots
########################################################################################

ggsave(path = here("output"), filename = fig_file, plot = last_plot(), device = "tiff")
