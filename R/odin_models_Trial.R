# 12/08/2022
# JF Vesga - London School of Hygiene and Tropical Medicine
# This script contains the simplified Odin compartmental models for CCHF in 
#livestock (SIR)and among humans (SEIR), used for running cohort trial models
# Also wrapper functions to format and package results
# to be used elsewhere

library(odin)

## Stochastic SEIR for spillover from livestock into humans - Cohort
seir_cohort_stoc <- odin::odin(
  {
    ## Core equations for transitions between compartments:
    update(F_S) <- F_S - newInfections_F + seed_F
    update(F_E) <- F_E + newInfections_F - newInfectious_F
    update(F_I) <- F_I + newInfectious_F - newRecovered_F - newfatal_F
    update(F_R) <- F_R + newRecovered_F

    # Infection in Others
    update(O_S) <- O_S - newInfections_O + seed_O
    update(O_E) <- O_E + newInfections_O - newInfectious_O
    update(O_I) <- O_I + newInfectious_O - newRecovered_O - newfatal_O
    update(O_R) <- O_R + newRecovered_O

    output(F_incidence) <- newInfectious_F
    output(O_incidence) <- newInfectious_O
    output(F_fatality) <- newfatal_F
    output(O_fatality) <- newfatal_O


    # Initial conditions
    initial(F_S) <- states0[1]
    initial(F_E) <- states0[2]
    initial(F_I) <- states0[3]
    initial(F_R) <- states0[4]
    initial(O_S) <- states0[5]
    initial(O_E) <- states0[6]
    initial(O_I) <- states0[7]
    initial(O_R) <- states0[8]

    ################### HUMANS
    ProbaInf_F <- 1 - exp(-beta_farmer * risk_livestock)
    ProbaInf_O <- 1 - exp(-beta_other * risk_livestock)
    p_EtoI <- 1 - exp(-time_to_infous)
    p_ItoR <- 1 - exp(-time_immune_human)
    p[1] <- 1 - CFR * clin_frac
    p[2] <- CFR * clin_frac


    ## Draws from binomial distributions for numbers changing between
    ## compartments:

    # farmers S--> E: Prob of infection
    newInfections_F <- round(rbinom(F_S, ProbaInf_F))

    # E-->I
    newInfectious_F <- round(F_E * p_EtoI)
    # I --> R
    noutI_F <- round(rbinom(F_I, p_ItoR))
    outcome_F[] <- rmultinom(noutI_F, p)
    newRecovered_F <- outcome_F[1]
    newfatal_F <- outcome_F[2]


    # OTHERS S--> E: Prob of infection
    newInfections_O <- round(rbinom(O_S, ProbaInf_O))

    # E-->I
    newInfectious_O <- round(O_E * p_EtoI)
    # I --> R
    noutI_O <- round(rbinom(O_I, (p_ItoR)))
    outcome_O[] <- rmultinom(noutI_O, p)
    newRecovered_O <- outcome_O[1]
    newfatal_O <- outcome_O[2]



    # User defined input

    # Transition rates and user defined
    risk_livestock <- interpolate(risk_livestock_t, risk_livestock_y, "linear")
    seed_F <- interpolate(risk_livestock_t, n_seed_f_y, "linear")
    seed_O <- interpolate(risk_livestock_t, n_seed_o_y, "linear")

    beta_farmer <- params[1]
    beta_other <- params[2]
    time_to_infous <- params[3]
    time_immune_human <- params[4]
    CFR <- params[5]
    clin_frac <- params[6]


    states0[] <- user()
    risk_livestock_t[] <- user()
    risk_livestock_y[] <- user()
    n_seed_f_y[] <- user()
    n_seed_o_y[] <- user()
    params[] <- user()

    # boilerplate from odin
    dim(states0) <- user()
    dim(risk_livestock_t) <- user()
    dim(risk_livestock_y) <- user()
    dim(n_seed_f_y) <- user()
    dim(n_seed_o_y) <- user()
    dim(params) <- user()
    dim(p) <- 2
    dim(outcome_F) <- 2
    dim(outcome_O) <- 2
  },
  verbose = FALSE
)


# Wrapping function to get parameters , orgaanize input, create model and run SEIR
run_seir_cohort <- function(init.seir, params, times, n.replicates = 1) {
  risk_livestock_t <- seq(1, params$nt)
  risk_livestock_y <- params$prev_inf_liv


  if (params$country == "AFG") {
    seir.params <- c(
      params$theta[["F_risk"]],
      params$theta[["F_risk"]] * params$theta[["O_factor"]],
      params$time_to_infous,
      params$time_immune_human,
      params$CFR,
      params$theta[["clin_frac"]],
      params$vacc_eff
    )
  } else {
    seir.params <- c(
      params$theta[["F_risk"]],
      params$theta[["F_risk"]] * params$theta[["O_factor"]],
      params$time_to_infous,
      params$time_immune_human,
      params$CFR,
      params$clin_frac,
      params$vacc_eff
    )
  }

  # Vaccination vector
  seed_f <- round(params$pop_coverage_f)
  seed_o <- round(params$pop_coverage_o)


  n_seed_f_y <- seq(1, params$nt) * 0
  n_seed_o_y <- seq(1, params$nt) * 0
  n_seed_f_y[params$start_vax] <- seed_f
  n_seed_o_y[params$start_vax] <- seed_o


  # Generate a stoch model instance with the input
  seir.instance <- seir_cohort_stoc$new(
    states0 = init.seir,
    risk_livestock_t = risk_livestock_t,
    risk_livestock_y = risk_livestock_y,
    n_seed_f_y = n_seed_f_y,
    n_seed_o_y = n_seed_o_y,
    params = seir.params
  )


  # Run the model
  t <- times

  if (n.replicates > 1) {
    out.seir <- (replicate(100, seir.instance$run(t)[, -1]))

    # Find mean of 100 replicates of stoch model
    out <- rowMeans(out.seir, na.rm = TRUE, dims = 2)
  } else {
    out.seir <- seir.instance$run(t)
    out <- out.seir
  }



  return(out)
}

# Call model function and organize/format model outputs SEIR
get_output_cohort <- function(params, theta, n.replicates = 1) {
  nruns <- nrow(theta)
  mo_length <- params$nt
  yr_length <- length(unique(as.Date(temp_month$month, format = "%Y")))

  # Allocate memory
  incidence_farmer <- matrix(NA, nrow = nruns, ncol = mo_length)
  incidence_other <- matrix(NA, nrow = nruns, ncol = mo_length)
  clincases_farmer <- matrix(NA, nrow = nruns, ncol = mo_length)
  clincases_other <- matrix(NA, nrow = nruns, ncol = mo_length)
  fatalities_farmer <- matrix(NA, nrow = nruns, ncol = mo_length)
  fatalities_other <- matrix(NA, nrow = nruns, ncol = mo_length)
  susceptible_farmer <- matrix(NA, nrow = nruns, ncol = mo_length)
  susceptible_other <- matrix(NA, nrow = nruns, ncol = mo_length)
  exposed_farmer <- matrix(NA, nrow = nruns, ncol = mo_length)
  exposed_other <- matrix(NA, nrow = nruns, ncol = mo_length)
  #######  Prepare for stochastic model
  init.seir <- c(
    F_S = 0,
    F_E = 0,
    F_I = 0,
    F_R = 0,
    O_S = 0,
    O_E = 0,
    O_I = 0,
    O_R = 0
  )

  for (jj in 1:nruns) {
    if (params$country == "AFG") {
      params$theta <- c(
        F_risk = theta$F_risk[jj] / params$risk_scale,
        O_factor = theta$O_factor[jj],
        clin_frac = theta$clin_frac[jj]
      )
    } else {
      # Pass calibrated parameters
      params$theta <- c(
        F_risk = theta$F_risk[jj] / params$risk_scale,
        O_factor = theta$O_factor[jj],
        RRreport = theta$RRreport[jj]
      )
    }

    # Risk of transmission in other occupations
    params$risk_O <- params$theta[["O_factor"]] * params$theta[["F_risk"]]


    params$prev_inf_liv <- params$prev_inf_liv_all[jj, ]
    t <- seq(1, params$nt)
    out <- as.data.frame(run_seir_cohort(init.seir, params, t, n.replicates))

    #------------------------------------ 
    # Process and format model output
    #------------------------------------ 


    incidence_farmer[jj, ] <- out$F_incidence
    incidence_other[jj, ] <- out$O_incidence
    clincases_farmer[jj, ] <- out$F_incidence * params$clin_frac
    clincases_other[jj, ] <- out$O_incidence * params$clin_frac
    fatalities_farmer[jj, ] <- out$F_fatality
    fatalities_other[jj, ] <- out$O_fatality
    susceptible_farmer[jj, ] <- out$F_S
    susceptible_other[jj, ] <- out$O_S
    exposed_farmer[jj, ] <- out$F_E * (1 - params$clin_frac)
    exposed_other[jj, ] <- out$O_E * (1 - params$clin_frac)
  }


  sim <- list(
    incidence_farmer = incidence_farmer,
    incidence_other = incidence_other,
    clincases_farmer = clincases_farmer,
    clincases_other = clincases_other,
    fatalities_farmer = fatalities_farmer,
    fatalities_other = fatalities_other,
    susceptible_farmer = susceptible_farmer,
    susceptible_other = susceptible_other,
    exposed_farmer = exposed_farmer,
    exposed_other = exposed_other
  )


  return(sim)
}


# Run a trial cohort
run_cohort <- function(params, theta) 
  {
  sim.sir <- get_output_sir(params, theta)
  
  risk_livestock <- data.matrix(sim.sir$prev_inf_liv)
  params$prev_inf_liv_all <- t(replicate(nrow(theta), colMeans(risk_livestock)))
  # Run SEIR
  sim.seir <- get_output_cohort(params, theta, 500)
  
  sims <- c(sim.sir, sim.seir)
  
  return(sims)
}



##########
## End Code
###########
