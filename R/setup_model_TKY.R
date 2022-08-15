# 12/08/2022
# JF Vesga - London School of Hygiene and Tropical Medicine
# This script setsup all the necessary parameteres for the model


###########################################################################################################
# Part 1:  MODEL PARAMETERS
###########################################################################################################

params <- list(
  nt = nrow(temp_month),
  ## ------------------- Livestock parameters

  # Number of livestock age groups
  ya = 5,

  # Population size N animals based on FAO report from region in 2008
  N_liv = 596386 * 0.5985738, # (USE ANIMAL/HUMAN RATIO FROM AFGHANISTAN)TurkStat 2014 Sheep + Cattle


  # livestock population age structure  and death rates
  # proportion of each age-group in the population
  pop_st = c(0.4918, 0.2499, 0.1270, 0.0646, 0.0667),
  mu_1 = 0.002643826 * 30,

  # daily death rate
  deathd = c(0.07618946, 0.07436794, 0.07466935, 0.07447458, 0.07470417),
  imm_t0 = c(0.29, 0.48, 0.8, 0.87, 0.87) * 0.1,

  # Duration of infectiousness in livestock (in days)
  D_inf_L = 7 / 30,
  D_lact_liv = 6,
  D_imm_liv = (5 * 12),

  ## ------------------- Human parameters

  N_hum = 596386, # Rural population in Erzican, Erzurum, Gumashane,Sivas & Tokat

  # Farmer labour 2011 file:///C:/Users/JuanVesga/Downloads/Turkish_agriculture_at_a_glance.pdf
  prop_o = 1 - (6143000 / 21101000),


  # Proportion of farmers immune at time 0. from data
  prop_F_R = 0.15,

  # Proportion of other occupations immune at time 0. from data
  prop_O_R = 0.08,

  # Reporting fraction. fraction of human infections reported by the surveillance
  RR = 0.75,
  RRreport = 0.75, # Increse of reporting over the yeras (2013 to 2019)

  # Latent period in humans in days (to calculate the rate at which human move from latent to infectious)
  D_lat_H = 4 / 30,

  # Infectious period in humans in days (to calculate the rate at which human move from infectious to immune)
  D_inf_H = 9 / 30,

  # Duration immunity in humans
  D_imm_H = 50 * 12,

  # case fatality rate
  CFR = 0.054, # Reported in Turkey by MoH

  # Clinical fraction
  clin_frac = 0.3, # 1-0.88, # Turkey data (Bodur et al.)

  # mortality and birth rates in human for constant pop size
  # life span in years  between 2008 and 2014 https://data.worldbank.org/indicator/SP.DYN.LE00.IN?locations=AF
  L = 61.5,
  b_d = (1 / (61.5 * 12)), # daily birth and death rate for constant pop

  sp_hz = 10, # spline hazard multiplier. This value will be multiplied by spline


  ## baseline vaccination paramts
  start_vax = 70,
  stop_vax = 73,
  vacc_eff = 0.9,
  time_imm = 0.5,
  pop_coverage_f = 0,
  pop_coverage_o = 0,
  p_vacc_livestock = 0,
  liv_vacc_yearly = FALSE
)

# Scale correction for risk
params$risk_scale <- 1

# Rate at which livestock and humans becomes infectious and immune
params$recover <- 1 / params$D_inf_L # rate at which infectious livestock become immune
params$time_to_infous <- 1 / params$D_lat_H # rate at which latent human become infectious
params$time_immune_human <- 1 / params$D_inf_H # rate at which infectious human become immune
params$time_susceptible_human <- 1 / params$D_imm_H # rate at which humans become susceptible again
params$time_passimm_loss_livestock <- 1 / params$D_lact_liv # rate of moving from passive immnity at birth to susceptible
params$time_susceptible_livestock <- 1 / params$D_imm_liv # rate of moving from passive immnity at birth to susceptible
params$time_protection <- 1 / params$time_imm
# Added parameters defined with others
params$indices_infL <- c(6, 7, 8, 9, 10) - 1
params$indices_livestock <- seq(1, 15) - 1


# total population of other occupations
params$N_O <- round(params$N_hum * params$prop_o, digit = 0)
# 17767
# total farming population
params$N_F <- params$N_hum - params$N_O
params$Birth_F <- params$N_F * params$b_d
params$Birth_O <- params$N_O * params$b_d
params$Birth <- sum(params$pop_st * params$deathd)
params$Ageing <- (1 - params$deathd) / 12
