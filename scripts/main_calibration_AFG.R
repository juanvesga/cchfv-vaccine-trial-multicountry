# 12/08/2022
# JF Vesga - London School of Hygiene and Tropical Medicine
# This script runs an MCMC algorithm (fitR) using Afghanistan specific data
# This script calls several function files and routines in folder R.
# Matrices of posterior parameter values and LLK are saved to be used later 
# to run simulations by sampling from the posterior

rm(list = ls())
library("pacman") # Load pacman package
p_load(
  Hmisc, reshape2, zoo, plyr, ggplot2, gmodels, data.table, devtools, fitR,
  coda, mvtnorm, ISOweek, lubridate, dplyr, zoo, stringr, here, psych, splines,
  data.table, broom, psych, matrixStats, profvis, odin
)
########################################################################################
# Part 0. select country and environmental driver of tick activity
########################################################################################
country <- "AFG" # leave as AFG for Afghanistan
driver <- "sat_def" # Leave as sat_def for saturation deficit


########################################################################################
# Part 1. Import set up: Load data and model parameters and functions
########################################################################################
source(here("R", "setup_data_AFG.R"))
source(here("R", "setup_model_AFG.R"))
source(here("R", "odin_models_vaccination_AFG.R"))



# Set MCMC important parameters

chain <- "chain1.csv"
n_iterations <- 20000
start_from_chain <- 0
n.particles <- 100
n.replicates <- 500

###########################################################################################################
# PART 2. INFERENCE: Model: prior, posterior and inference
###########################################################################################################

# Initialize parameter values

theta <- data.frame(
  A = 0.3249264, # driving temperature dependent force of infection
  D_lact_liv = 8.015156, # mnths of waning collostrum immunity
  D_imm_liv = 58.21934, # mnths of waning acquired immunity
  F_risk = 0.3230339, # risk for farmers
  O_factor = 0.4401415, # factor reducing transmission among "others"
  clin_frac = 0.20      # Fraction of incidence developing clinical symptoms
) 


# 2.1. MY PRIOR function
###########################################################################################################
my_prior <- function(theta) {
  ## uniform prior on A
  log.prior.A <- dunif(theta[["A"]], min = 0, max = 100, log = TRUE)

  ## uniform prior on D_lact_liv
  log.prior.D_lact_liv <- dunif(theta[["D_lact_liv"]], min = 2, max = 12, log = TRUE)

  ## uniform prior on D_imm_liv
  log.prior.D_imm_liv <- dunif(theta[["D_imm_liv"]], min = 36, max = 120, log = TRUE)

  ## uniform prior on risk to farmers infected
  log.prior.F_risk <- dunif(theta[["F_risk"]], min = 0, max = 10, log = TRUE) # 0.01 / 5
  ## uniform prior on risk to others: scalar for risk to other occupations
  log.prior.O_factor <- dunif(theta[["O_factor"]], min = 0, max = 1, log = TRUE)

  log.prior.clin_frac <- dunif(theta[["clin_frac"]], min = 0.06, max = 0.9, log = TRUE)


  return(log.prior.A +
    log.prior.D_lact_liv +
    log.prior.D_imm_liv +
    log.prior.F_risk +
    log.prior.O_factor +
    log.prior.clin_frac)
}
my_prior(theta)

# 2.2. MY POSTERIOR
###########################################################################################################
my_posterior <- function(theta) {
  
  theta_in <- bind_rows(theta)
  # 1. Run SIR deterministic model and get output
  #########################################################################################################

  sim.sir <- get_output_sir(params, theta_in)

  # 2. Estimate LogLikelihood and output posterior density
  #########################################################################################################

  Llk_weights <- c(
    length(observations$prev_liv_age_pos_IgG),
    length(observations$prev_liv_all_pos_IgG),
    length(observations$cases_human_mo),
    length(observations$cases_human_yr),
    length(observations$deaths_human_yr),
    length(observations$prev_farmer_denoma),
    length(observations$prev_other_denoma)
  )


  loglik_liv_prev_age <- sum(dbinom(
    x = observations$prev_liv_age_pos_IgG, size = observations$prev_liv_age_denom,
    prob = sim.sir$l_prev_age, log = TRUE
  ), na.rm = TRUE) 

  loglik_liv_prev_all <- dbinom(
    x = observations$prev_liv_all_pos_IgG, size = observations$prev_liv_all_denom,
    prob = sim.sir$l_prev_all, log = TRUE
  ) / Llk_weights[2]


  log.likelihood <- loglik_liv_prev_age + loglik_liv_prev_all

  log.prior <- my_prior(theta)

  log.posterior_livestock <- log.prior + log.likelihood

  # 3. Run SEIR stochastic model and get output
  #########################################################################################################

  # Prepare parameters

  params$theta <- theta_in

  # Infectious prevalence in livestock
  risk_livestock <- data.matrix(sim.sir$prev_inf_liv)
  
  params$prev_inf_liv_all <- risk_livestock 

  sim <- get_output_seir(params, theta_in, n.replicates)

  loglik_human_reported_month <- sum(dpois(
    x = round(observations$cases_human_mo),
    lambda = as.numeric(sim$h_incidence_reported_month[observations$index_mo_cases]),
    log = TRUE
  ),
  na.rm = TRUE
  ) / Llk_weights[3]

  loglik_human_reported_year <- sum(dpois(
    x = round(observations$cases_human_yr),
    lambda = as.numeric(sim$h_incidence_reported_year[observations$index_yr_cases]),
    log = TRUE
  ),
  na.rm = TRUE
  ) / Llk_weights[4]


  loglik_human_fatality_year <- sum(dpois(
    x = round(observations$deaths_human_yr),
    lambda = as.numeric(sim$h_fatality_year[observations$index_yr_deaths]),
    log = TRUE
  ),
  na.rm = TRUE
  ) / Llk_weights[5]


  loglik_prev_farmer <- dbinom(
    x = observations$prev_farmer_Posa, size = observations$prev_farmer_denoma,
    prob = sim$h_prev_farmer, log = TRUE
  ) / Llk_weights[6]

  loglik_prev_other <- dbinom(
    x = observations$prev_other_Posa, size = observations$prev_other_denoma,
    prob = sim$h_prev_other, log = TRUE
  ) / Llk_weights[7]


  return(log.density = log.posterior_livestock +
    loglik_human_reported_month +
    loglik_human_reported_year * 0.5 +
    loglik_human_fatality_year * 0.5 +
    loglik_prev_farmer +
    loglik_prev_other)
}

my_posterior(theta)

###########################################################################################################
#  PART 3. Inference MCMC-MH
###########################################################################################################

if (start_from_chain == 1) {
  trace1 <- read.csv(here("output", country, "trace_chain1_saturation_deficit.csv"))
  cov1 <- read.csv(here("output", country, "cov_mat_chain1_saturation_deficit.csv"))

  init.theta <- as.numeric(trace1[nrow(trace1), 2:(ncol(trace1) - 1)])

  covarianceMat <- data.matrix(cov1)[, -1]
  rownames(covarianceMat) <- colnames(covarianceMat)
  names(init.theta) <- rownames(covarianceMat)
} else {
  init.theta <- c(
    A = 0.34, 
    D_lact_liv = 7, 
    D_imm_liv = 70, 
    F_risk = 0.54, 
    O_factor = 0.08,
    clin_frac = 0.2
  )

  covarianceMat <- NULL
}


proposal.sd <- init.theta / 12



n.iterations <- n_iterations
print.info.every <- 100

limits <- list(
  lower = c(
    A = 0.1,
    D_lact_liv = 2, 
    D_imm_liv = 36, 
    F_risk = 0,
    O_factor = 0,
    clin_frac = 0.06
  ), 

  upper = c(
    A = 5,
    D_lact_liv = 12, 
    D_imm_liv = 120, 
    F_risk = 10,
    O_factor = 1,
    clin_frac = 0.9
  )
) 

adapt.size.start <- 500
adapt.size.cooling <- 0.999
adapt.shape.start <- 500
acceptance.rate.weight <- NULL
acceptance.window <- NULL
covmat <- NULL 

# 3.1. Chains
##############################################################################################################
start_time <- Sys.time()


res_mcmc <- mcmcMH(
  target = my_posterior,
  init.theta = init.theta,
  proposal.sd = proposal.sd,
  n.iterations = n.iterations,
  limits = limits,
  adapt.size.start = adapt.size.start,
  adapt.size.cooling = adapt.size.cooling,
  adapt.shape.start = adapt.shape.start,
  print.info.every = print.info.every,
  covmat = covmat,
  verbose = FALSE
)
end_time <- Sys.time()
timelag <- end_time - start_time
##############################################################################################################
# 3.2. Export results ########################################################################################
##############################################################################################################
vc <- res_mcmc$trace
trace <- apply(vc, 2, unlist)
trace <- data.frame(trace)
acc_rate <- res_mcmc$acceptance.rate
cov_mat <- res_mcmc$covmat.empirical

write.csv(trace, here("output", country, paste("trace_", chain, sep = "")))
write.csv(cov_mat, here("output", country, paste("cov_mat_", chain, sep = "")))
write.csv(acc_rate, here("output", country, paste("acc_rate_", chain, sep = "")))

#########################################################################################################
### END OF CODE
#########################################################################################################
