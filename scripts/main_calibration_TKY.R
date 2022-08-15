# 12/08/2022
# JF Vesga - London School of Hygiene and Tropical Medicine
# This script runs an MCMC algorithm (fitR) using Turkey specific data
# This script calls several function files and routines in folder R.
# Matrices of posterior parameter values and LLK are saved to be used later 
# to run simulations by sampling from the posterior

rm(list = ls())
library("pacman") # Load pacman package
p_load(
  Hmisc, reshape2, zoo, plyr, ggplot2, gmodels, data.table, devtools, fitR,
  coda, mvtnorm, ISOweek, lubridate, dplyr, zoo, stringr, here, psych, splines,
  data.table, broom, psych, matrixStats, profvis, odin, neldermead
)
##############################################################################################################
# Import set up: Load data and model parameters and functions ########################################################################################
##############################################################################################################
country <- "TKY"
driver <- "rel_hum" # soil_t, sat_def, rel_hum, NDVI
runtype <- "MCMC" # MCMC or simplex or lhs
start_set <- "MLE" # random or MLE or MLE_LHS


source(here("R", "setup_data_TKY.R"))
source(here("R", "setup_model_TKY.R"))
source(here("R", "odin_models_vaccination_TKY.R"))


##############################################################################################################
# MCMC parameters ########################################################################################
##############################################################################################################


chain <- "chain1"
n_iterations <- 3000
start_from_chain <- 0
n.particles <- 100
n.replicates <- 500

###########################################################################################################
# PART 2. INFERENCE: Model: prior, posterior and inference
###########################################################################################################

theta <- data.frame(
  A = 0.3249264, # driving temperature dependent force of infection
  D_inf_L = 7 / 30,
  D_lact_liv = 8.015156, # mnths of waning collostrum immunity
  D_imm_liv = 58.21934, # mnths of waning acquired immunity
  F_risk = 6e-2, # risk for farmers
  O_factor = 0.4401415, # Fcator for otehrs
  RRreport = 0.57,
  clin_frac = 0.12
) # Clinical fraction

names <- colnames(theta)
# 2.1. MY PRIOR

###########################################################################################################
my_prior <- function(theta) {
  ## uniform prior on A
  log.prior.A <- dunif(theta[["A"]], min = 1, max = 10, log = TRUE)

  log.prior.D_inf_L <- dunif(theta[["D_inf_L"]], min = 5 / 30, max = 9 / 30, log = TRUE)

  ## uniform prior on D_lact_liv
  log.prior.D_lact_liv <- dunif(theta[["D_lact_liv"]], min = 1, max = 18, log = TRUE)

  ## uniform prior on D_imm_liv
  log.prior.D_imm_liv <- dunif(theta[["D_imm_liv"]], min = 24, max = 60, log = TRUE)

  ## uniform prior on risk to farmers infected
  log.prior.F_risk <- dunif(theta[["F_risk"]], min = 10e-3 * 1e3, max = 10e-2 * 1e3, log = TRUE) # 0.01 / 5
  ## uniform prior on risk to others: scalar for risk to other occupations
  log.prior.O_factor <- dunif(theta[["O_factor"]], min = 0, max = 1, log = TRUE)

  log.prior.RRreport <- dunif(theta[["RRreport"]], min = 0.5, max = 1, log = TRUE)

  log.prior.clin_frac <- dunif(theta[["clin_frac"]], min = 0.1, max = 0.4, log = TRUE)



  return(log.prior.A +
    log.prior.D_inf_L +
    log.prior.D_lact_liv +
    log.prior.D_imm_liv +
    log.prior.F_risk +
    log.prior.O_factor +
    log.prior.RRreport +
    log.prior.clin_frac)
}
my_prior(theta)

# 2.2. MY POSTERIOR
###########################################################################################################
my_posterior <- function(theta, label = names) {
  signllk <- 1
  if (size(theta)[2] == 1) {
    signllk <- -1
    theta <- as.data.frame(t(theta))
    colnames(theta) <- label
  }



  theta_in <- bind_rows(theta)

  theta_in$F_risk <- theta_in$F_risk / 1e3
  # 1. Run SIR deterministic model and get output
  #########################################################################################################

  sim.sir <- get_output_sir(params, theta_in)

  # 2. Estimate LogLikelihood and output posterior density
  #########################################################################################################

  Llk_weights <- c(
    length(observations$prev_liv_age_pos_IgG),
    length(observations$prev_liv_all_pos_IgG),
    length(observations$fitcases_human_mo),
    length(observations$prev_farmer_denoma),
    length(observations$prev_other_denoma)
  )



  loglik_liv_prev_age <- sum(dbinom(
    x = observations$prev_liv_age_pos_IgG, size = observations$prev_liv_age_denom,
    prob = sim.sir$l_prev_age4g, log = TRUE
  ), na.rm = FALSE) / Llk_weights[1]

  if (is.na(loglik_liv_prev_age)) {
    loglik_liv_prev_age <- Inf
  }

  loglik_liv_prev_all <- sum(dbinom(
    x = observations$prev_liv_all_pos_IgG, size = observations$prev_liv_all_denom,
    prob = sim.sir$l_prev_all, log = TRUE
  ), na.rm = FALSE) / Llk_weights[2]

  if (is.na(loglik_liv_prev_all)) {
    loglik_liv_prev_all <- Inf
  }


  log.prior <- my_prior(theta)

  log.posterior_livestock <- loglik_liv_prev_all

  # 3. Run SEIR stochastic model and get output
  #########################################################################################################

  # Prepare parameters
  params$theta <- theta_in

  # Infectious prevalence in livestock
  risk_livestock <- data.matrix(sim.sir$prev_inf_liv)

  params$prev_inf_liv_all <- risk_livestock

  sim <- get_output_seir(params, theta_in, n.replicates)


  loglik_human_reported_month <- sum(dpois(
    x = round(observations$fitcases_human_mo),
    lambda = as.numeric(sim$h_incidence_reported_month[observations$index_mo_fitcases]),
    log = TRUE
  ), na.rm = FALSE) / (Llk_weights[3])
  if (is.na(loglik_human_reported_month)) {
    loglik_human_reported_month <- Inf
  }



  loglik_prev_farmer <- sum(dbinom(
    x = observations$prev_farmer_Posa, size = observations$prev_farmer_denoma,
    prob = sim$h_prev_farmer, log = TRUE
  ), na.rm = TRUE) / Llk_weights[4]

  if (is.na(loglik_prev_farmer)) {
    loglik_prev_farmer <- Inf
  }

  loglik_prev_other <- sum(dbinom(
    x = observations$prev_other_Posa, size = observations$prev_other_denoma,
    prob = sim$h_prev_other, log = TRUE
  ), na.rm = TRUE) / Llk_weights[5]

  if (is.na(loglik_prev_other)) {
    loglik_prev_other <- Inf
  }


  return(log.density = signllk * (log.prior +
    log.posterior_livestock +
    loglik_human_reported_month +
    loglik_prev_farmer +
    loglik_prev_other))
}


###########################################################################################################
#  PART 3. Inference MCMC-MH
###########################################################################################################
limits <- list(
  lower = c(
    A = 1,
    D_inf_L = 5 / 30,
    D_lact_liv = 1,
    D_imm_liv = 24,
    F_risk = 10e-3 * 1e3,
    O_factor = 0,
    RRreport = 0.5,
    clin_frac = 0.1
  ),
  upper = c(
    A = 10,
    D_inf_L = 9 / 30,
    D_lact_liv = 18,
    D_imm_liv = 60,
    F_risk = 10e-2 * 1e3,
    O_factor = 1,
    RRreport = 1,
    clin_frac = 0.4
  )
)

if (start_from_chain == 1) {
  trace1 <- read.csv(here("output", country, paste("trace_", chain, "_", driver, ".csv", sep = "")))
  cov1 <- read.csv(here("output", country, paste("cov_mat_", chain, "_", driver, ".csv", sep = "")))
  init.theta <- as.numeric(trace1[nrow(trace1), 2:(ncol(trace1) - 1)])

  covarianceMat <- data.matrix(cov1)[, -1]
  rownames(covarianceMat) <- colnames(covarianceMat)
  names(init.theta) <- rownames(covarianceMat)
} else if (start_set == "MLE") {
  x <- read.csv(here("output", country, paste("MLE_", driver, ".csv", sep = "")))
  x <- x[, 2]

  init.theta <- c(
    A = x[1],
    D_inf_L = x[2],
    D_lact_liv = x[3],
    D_imm_liv = x[4],
    F_risk = x[5],
    O_factor = x[6],
    RRreport = x[7],
    clin_frac = x[8]
  )


  covarianceMat <- NULL
} else if (start_set == "MLE_LHS") {
  x <- read.csv(here("output", country, paste("MLE_LHS_", driver, ".csv", sep = "")))
  x <- x[, 2]

  init.theta <- c(
    A = x[1],
    D_inf_L = x[2],
    D_lact_liv = x[3],
    D_imm_liv = x[4],
    F_risk = x[5],
    O_factor = x[6],
    RRreport = x[7],
    clin_frac = x[8]
  )

  covarianceMat <- NULL
} else if (start_set == "random") {
  init.theta <- limits$lower
  for (ii in 1:length(theta)) {
    init.theta[ii] <- runif(1, min = limits$lower[ii], max = limits$upper[ii])
  }
}


proposal.sd <- init.theta * 0.2

n.iterations <- n_iterations
print.info.every <- 50



if (runtype == "simplex") {
  fminsearch <- neldermead::fminsearch
  opt <- optimset(TolX = 1.e-2, Display = "iter", MaxFunEvals = 500)
  x1 <- fminsearch(
    fun = my_posterior,
    x0 = init.theta,
    options = opt
  )
  x <- x1$optbase$xopt

  write.csv(x, here("output", country, paste("MLE_", driver, ".csv", sep = "")))
} else if (runtype == "lhs") {
  require(lhs)
  set.seed(1976)

  grid <- randomLHS(n_iterations, length(limits$lower))
  thetas <- sapply(1:length(limits$lower), function(z) {
    qunif(
      grid[, z],
      limits$lower[z], limits$upper[z]
    )
  })
  colnames(thetas) <- names
  llk <- seq(1:n_iterations)

  for (jj in 1:size(grid, 1)) {
    llk[jj] <- my_posterior(thetas[jj, ])
  }


  final <- cbind(thetas, llk)

  write.csv(trace, here("output", country, paste("LHS", "_", driver, ".csv", sep = "")))


  id <- which(llk == max(llk))
  MLE <- final[id, ]
  MLE <- MLE[1:size(grid, 2)]

  write.csv(MLE, here("output", country, paste("MLE_LHS_", driver, ".csv", sep = "")))
} else {
  adapt.size.start <- 1000
  adapt.size.cooling <- 0.999
  adapt.shape.start <- 500
  acceptance.rate.weight <- NULL
  acceptance.window <- NULL
  # max.scaling.sd <- 5
  covmat <- covarianceMat

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
  # Export results ########################################################################################
  ##############################################################################################################
  vc <- res_mcmc$trace
  trace <- apply(vc, 2, unlist)
  trace <- data.frame(trace)
  acc_rate <- res_mcmc$acceptance.rate
  cov_mat <- res_mcmc$covmat.empirical

  write.csv(trace, here("output", country, paste("trace_", chain, "_", driver, ".csv", sep = "")))
  write.csv(cov_mat, here("output", country, paste("cov_mat_", chain, "_", driver, ".csv", sep = "")))
  write.csv(acc_rate, here("output", country, paste("acc_rate_", chain, "_", driver, ".csv", sep = "")))
}
##############################################################################################################
# End of code ########################################################################################
##############################################################################################################
