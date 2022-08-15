# 12/08/2022
# JF Vesga - London School of Hygiene and Tropical Medicine
# This script takes a sample of the posterior distribution of the calibrated model
# to run simulated trial cohorts exploring different levels of sample size and
# Vaccine efficacy. After the procedures, plots are produced and saved
# Note: to produce manuscript plots use scripts plot_fig1_map and plot_figs2_4 

rm(list = ls())
library("pacman") # Load pacman package
p_load(
  Hmisc, stats, reshape2, zoo, plyr, ggplot2, gmodels, data.table, devtools, fitR,
  coda, mvtnorm, ISOweek, lubridate, dplyr, zoo, stringr, here, psych, splines,
  data.table, broom, psych, matrixStats, profvis, odin, foreach, doParallel,
  pwr, metR, scales
)
########################################################################################
# Part 0. select country, driver and new simulation or saved runs
########################################################################################

country <- "AFG"
driver <- "sat_def"
runs <- "saved"

labc <- "Afghanistan"
if (country == "TKY") {
  labc <- "Turkey"
} else if (country == "SA") {
  labc <- "South Africa"
}


########################################################################################
# Part 1. Import set up: Load data and model parameters and functions
########################################################################################
source(here("R", "odin_models_Trial.R"))
source(here("R", "vaccine_functions.R"))

if (country == "AFG") {
  source(here("R", "odin_models_vaccination_AFG.R"))
  source(here("R", "setup_data_AFG.R"))
  source(here("R", "setup_model_AFG.R"))

  vax_t <- "2008-05-01"
} else if (country == "SA") {
  source(here("R", "odin_models_vaccination_SA.R"))
  source(here("R", "setup_data_SA.R"))
  source(here("R", "setup_model_SA.R"))

  vax_t <- "2008-09-15"
} else if (country == "TKY") {
  source(here("R", "odin_models_vaccination_TKY.R"))
  source(here("R", "setup_data_TKY.R"))
  source(here("R", "setup_model_TKY.R"))

  vax_t <- "2008-03-15"
}


# Load posterior parameters from calibration
posteriors <- read.table(here("output", 
                              country, 
                              paste("posterior_", driver, ".txt", sep = "")
                              , sep = ""))

# Sample from posterior
nsim <- 200

rand <- posteriors[sample(nrow(posteriors), nsim), ]

# Pass parameters to structure
if (country == "AFG") {
  theta <- rand %>%
    select(A, 
           D_lact_liv, 
           D_imm_liv, 
           F_risk, 
           O_factor, 
           clin_frac)

  theta_median <- data.frame(
    A = median(theta$A),
    D_lact_liv = median(theta$D_lact_liv),
    D_imm_liv = median(theta$D_imm_liv),
    F_risk = median(theta$F_risk),
    O_factor = median(theta$O_factor),
    clin_frac = median(theta$clin_frac)
  )
} else if ((country == "TKY")) {
  theta <- rand %>%
    select(A,
           D_inf_L,
           D_lact_liv,
           D_imm_liv,
           F_risk,
           O_factor,
           RRreport,
           clin_frac)

  theta_median <- data.frame(
    A = median(theta$A),
    D_inf_L = median(theta$D_inf_L),
    D_lact_liv = median(theta$D_lact_liv),
    D_imm_liv = median(theta$D_imm_liv),
    F_risk = median(theta$F_risk),
    O_factor = median(theta$O_factor),
    RRreport = median(theta$RRreport),
    clin_frac = median(theta$clin_frac)
  )
} else if ((country == "SA")) {
  theta <- rand %>%
    select(A,
           D_inf_L,
           D_lact_liv,
           D_imm_liv, 
           F_risk,
           O_factor,
           RRreport,
           clin_frac)

  theta_median <- data.frame(
    A = median(theta$A),
    D_inf_L = median(theta$D_inf_L),
    D_lact_liv = median(theta$D_lact_liv),
    D_imm_liv = median(theta$D_imm_liv),
    F_risk = median(theta$F_risk),
    O_factor = median(theta$O_factor),
    RRreport = median(theta$RRreport),
    clin_frac = median(theta$clin_frac)
  )
}

# Months Time vector
mo <- as.Date(seq(as.Date(observations$start_date), 
                  by = "month", length.out = params$nt))

# Find best vaccination start: tvax was selected visually to match start of high peak
tvax <- which(mo == vax_t)


# Baseline conditions
params$start_vax <- tvax
params$stop_vax <- tvax + 3
params$vacc_eff <- 0.9
params$time_protection <- 1 / 0.5
params$pop_coverage_f <- 0
params$pop_coverage_o <- 0
params$p_vacc_livestock <- 0
params$liv_vacc_yearly <- FALSE
params$country <- country


#########################################
### Run cohort
#########################################
seeds <- seq(1000, 100000, 1000)
nruns <- length(seeds)

if (runs == "simulate") {

  # Create memory arrays
  inc <- array(NA, c(nsim, nruns, params$nt))
  clin <- array(NA, c(nsim, nruns, params$nt))
  fatal <- array(NA, c(nsim, nruns, params$nt))
  cuminc <- array(NA, c(nsim, nruns, params$nt))
  cumclin <- array(NA, c(nsim, nruns, params$nt))
  cumfatal <- array(NA, c(nsim, nruns, params$nt))
  ARU_clin <- array(NA, c(nsim, nruns, params$nt))
  ARU_inc <- array(NA, c(nsim, nruns, params$nt))
  ARU_fat <- array(NA, c(nsim, nruns, params$nt))

  sus <- array(NA, c(nsim, nruns, params$nt))
  expo <- array(NA, c(nsim, nruns, params$nt))

  for (hh in 1:nruns) {
    print(hh)
    params$pop_coverage_f <- seeds[hh]

    # Run simulation
    sims <- run_cohort(params, theta)

    # Collect results in arrays
    inc[, hh, ] <- sims$incidence_farmer
    clin[, hh, ] <- sims$clincases_farmer
    sus[, hh, ] <- sims$susceptible_farmer
    expo[, hh, ] <- sims$exposed_farmer
    fatal[, hh, ] <- sims$fatalities_farmer
    cuminc[, hh, ] <- t(apply(sims$incidence_farmer, 1, cumsum))
    cumclin[, hh, ] <- t(apply(sims$clincases_farmer, 1, cumsum))
    cumfatal[, hh, ] <- t(apply(sims$fatalities_farmer, 1, cumsum))
    ARU_clin[, hh, ] <- cumclin[, hh, ] / sims$susceptible_farmer[hh, tvax + 1]
    ARU_inc[, hh, ] <- cuminc[, hh, ] / sims$susceptible_farmer[hh, tvax + 1]
    ARU_fat[, hh, ] <- cumfatal[, hh, ] / sims$susceptible_farmer[hh, tvax + 1]
  }

  # Passs into list structure
  cohort <- list()
  cohort$inc <- inc
  cohort$clin <- clin
  cohort$fatal <- fatal
  cohort$cuminc <- cuminc
  cohort$cumclin <- cumclin
  cohort$cumfatal <- cumfatal
  cohort$ARU_clin <- ARU_clin
  cohort$ARU_inc <- ARU_inc
  cohort$ARU_fat <- ARU_fat

  cohort$sus <- sus
  cohort$expo <- expo

  saveRDS(cohort, 
          file = here("output", country,
                      paste("cohort_Sept", "_", driver, ".rds", sep = "")))
} else {
  # if running from saved results
  cohort <- readRDS(here("output", country,
                         paste("cohort_Sept", "_", driver, ".rds", sep = "")))
}

#####################
# Attack rates at 6 months for infection, clinical case and fatality
#####################

# Clinical case AR

id <- which(is.na(cohort$ARU_clin))
AR <- cohort$ARU_clin
AR[id] <- 0
qtls <- as.data.frame(rowQuantiles(t(AR[, 1, ]),
  probs = c(0.025, 0.5, 0.975)
))
qtls$x <- mo

# Infection case AR
id <- which(is.na(cohort$ARU_inc))
ARi <- cohort$ARU_inc
ARi[id] <- 0
qtlsi <- as.data.frame(rowQuantiles(t(ARi[, 1, ]),
  probs = c(0.025, 0.5, 0.975)
))
qtlsi$x <- mo


# Fatality case AR
id <- which(is.na(cohort$ARU_fat))
ARf <- cohort$ARU_fat
ARf[id] <- 0
qtlsf <- as.data.frame(rowQuantiles(t(ARf[, 1, ]),
  probs = c(0.025, 0.5, 0.975)
))
qtlsf$x <- mo

# Define colors for groups
colors <- c("Infection" = "lightblue",
            "Case" = "Purple",
            "Fatal" = "black")


# Attack rates plots
attackR <- ggplot() +
  geom_line(data = qtlsi, aes(x = x, y = `50%`, col = "Infection"), lwd = 1) +
  geom_ribbon(data = qtlsi, aes(x = x, ymin = `2.5%`, ymax = `97.5%`), fill = "lightblue", alpha = 0.5) +
  geom_line(data = qtls, aes(x = x, y = `50%`, col = "Case"), lwd = 1) +
  geom_ribbon(data = qtls, aes(x = x, ymin = `2.5%`, ymax = `97.5%`), fill = "Purple", alpha = 0.2) +
  geom_line(data = qtlsf, aes(x = x, y = `50%`, col = "black"), lwd = 1) +
  geom_ribbon(data = qtlsf, aes(x = x, ymin = `2.5%`, ymax = `97.5%`), fill = "grey", alpha = 0.6) +
  labs(title = labc, x = " ", y = "CCHFV Attack Rate", color = "") +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(
    axis.text = element_text(colour = "black", size = 10, face = "bold"),
    axis.title = element_text(colour = "black", size = 12, face = "bold"),
    title = element_text(colour = "black", size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "bottom",
    legend.text = element_text(colour = "black", size = 10, face = "bold")
  ) +
  scale_x_date(
    date_breaks = "2 month",
    date_labels = "%b-%Y",
    limits = as.Date(c(mo[tvax], mo[tvax + 6]))
  ) +
  ylim(0, max(qtlsi$`97.5%`[tvax + 6]) * 1.1)


save(attackR, file = here("output", country, "attackRFig.rdata"))



########################
# Months to reach endpoints analysis
########################

AR <- mean(cohort$ARU_clin[, nruns, tvax + 6])
AR_i <- mean(cohort$ARU_inc[, nruns, tvax + 6])
AR_f <- mean(cohort$ARU_fat[, nruns, tvax + 6])

# Define VE for H0 and H1
VE_0 <- 0.3 # H0
VE_1 <- 0.6 # H1

# Estimate necessary endpoints (n) by type of endpoint
# Clinical
p0 <- AR * (1 - VE_0) / (AR * (1 - VE_0) + AR)
p1 <- AR * (1 - VE_1) / (AR * (1 - VE_1) + AR)

# Infection
p0_i <- AR_i * (1 - VE_0) / (AR_i * (1 - VE_0) + AR_i)
p1_i <- AR_i * (1 - VE_1) / (AR_i * (1 - VE_1) + AR_i)

# Fatality
p0_f <- AR_f * (1 - VE_0) / (AR_f * (1 - VE_0) + AR_f)
p1_f <- AR_f * (1 - VE_1) / (AR_f * (1 - VE_1) + AR_f)


# Estimate effective size
eff.size <- ES.h(p0, p1) #  same as doing this -> 2*asin(sqrt(p0))-2*asin(sqrt(p))

# Run power test , one-sided
N <- pwr.p.test(h = eff.size, sig.level = 0.025, power = 0.90, alternative = "greater")

# endpoints
endpoints <- ceiling(N$n)


# Estimate time taken in months to reach necesary endpoints
samplesize <- matrix(NA, nrow = nruns)
ttoend <- matrix(NA, nrow = nruns, ncol = nsim)
ttoendi <- matrix(NA, nrow = nruns, ncol = nsim)
ttoendf <- matrix(NA, nrow = nruns, ncol = nsim)


for (cc in 1:nruns) {
  samplesize[cc] <- seeds[cc]
  for (ss in 1:nsim) {
    id <- which(cohort$cumclin[ss, cc, ] >= endpoints)
    idi <- which(cohort$cuminc[ss, cc, ] >= endpoints)
    idf <- which(cohort$cumfatal[ss, cc, ] >= endpoints)



    if (length(id) > 0) {
      ttoend[cc, ss] <- min(id) - tvax
    } else {
      ttoend[cc, ss] <- params$nt - tvax
    }

    if (length(idi) > 0) {
      ttoendi[cc, ss] <- min(idi) - tvax
    } else {
      ttoendi[cc, ss] <- params$nt - tvax
    }
    if (length(idf) > 0) {
      ttoendf[cc, ss] <- min(idf) - tvax
    } else {
      ttoendf[cc, ss] <- params$nt - tvax
    }
  }
}


####################
# Plot time to reach endpoints anaysis
####################
dat <- as.data.frame(cbind(ttoend, samplesize))
datm <- melt(dat, id = "V201")
names(datm) <- c("sampleSize", "sim", "months")

dati <- as.data.frame(cbind(ttoendi, samplesize))
datmi <- melt(dati, id = "V201")
names(datmi) <- c("sampleSize", "sim", "months")

datf <- as.data.frame(cbind(ttoendf, samplesize))
datmf <- melt(datf, id = "V201")
names(datmf) <- c("sampleSize", "sim", "months")


qtls <- as.data.frame(rowQuantiles((ttoend),
  probs = c(0.025, 0.5, 0.975)
))
qtls$x <- samplesize

qtlsi <- as.data.frame(rowQuantiles((ttoendi),
  probs = c(0.025, 0.5, 0.975)
))
qtlsi$x <- samplesize

qtlsf <- as.data.frame(rowQuantiles((ttoendf),
  probs = c(0.025, 0.5, 0.975)
))
qtlsf$x <- samplesize


colors <- c("Infection" = "orange", "Case" = "blue")

qtls_plot <- ggplot() +
  geom_ribbon(data = qtlsi, aes(x = x, ymin = `2.5%`, ymax = `97.5%`), fill = "orange", alpha = 0.5) +
  geom_line(data = qtlsi, aes(x = x, y = `50%`, col = "Infection"), lwd = 1)
if (length(unique(qtls$`50%`)) > 1) {
  qtls_plot <- qtls_plot +
    geom_line(data = qtls, aes(x = x, y = `50%`, col = "Case"), lwd = 1)
}

qtls_plot <- qtls_plot +
  geom_ribbon(data = qtls, aes(x = x, ymin = `2.5%`, ymax = `97.5%`), fill = "blue", alpha = 0.2)

qtls_plot <- qtls_plot +
  labs(title = "", x = "Sample Size", y = "Months to reach\n 150 end points", color = "") +
  scale_x_continuous(breaks = seq(0, max(seeds), 20000), labels = scales::comma) +
  scale_y_continuous(breaks = seq(0, 170, 10)) + # , limits = c(0,170))+
  scale_color_manual(values = colors) +
  theme_minimal() +
  theme(
    axis.text = element_text(colour = "black", size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "bottom",
    legend.text = element_text(colour = "black", size = 10, face = "bold")
  ) +
  geom_hline(yintercept = 6, linetype = "dotted", color = "black", lwd = 1.2)

qtls_plot

save(qtls_plot, file = here("output", country, "endPFig.rdata"))


########################
## Patient-months analysis
########################

# Create memory arrays
pmo <- (apply(cohort$sus, c(1, 2), FUN = cumsum))
cases <- (apply(cohort$cumclin, c(1, 2), FUN = cumsum))
icases <- (apply(cohort$cuminc, c(1, 2), FUN = cumsum))
fcases <- (apply(cohort$cumfatal, c(1, 2), FUN = cumsum))

# Vector with different number of endpoint thresholds (VE=50, VE=60, VE=90)
endpoints_vec <- seq(5, 400, 5)

nt <- length(endpoints_vec)

pmtoend <- matrix(NA, nrow = nt, ncol = nsim)
pmtoendi <- matrix(NA, nrow = nt, ncol = nsim)
pmtoendf <- matrix(NA, nrow = nt, ncol = nsim)

for (cc in 1:nt) {
  for (ss in 1:nsim) {
    id <- which(cohort$cumclin[ss, , tvax + 6] >= endpoints_vec[cc])
    idi <- which(cohort$cuminc[ss, , tvax + 6] >= endpoints_vec[cc])
    idf <- which(cohort$cumfatal[ss, , tvax + 6] >= endpoints_vec[cc])

    if (length(id) > 0) {
      pmtoend[cc, ss] <- pmo[tvax + 6, ss, min(id)]
    } else {
      pmtoend[cc, ss] <- max(pmo[tvax + 6, ss, ])
    }

    if (length(idi) > 0) {
      pmtoendi[cc, ss] <- pmo[tvax + 6, ss, min(idi)]
    } else {
      pmtoendi[cc, ss] <- max(pmo[tvax + 6, ss, ])
    }
    if (length(idf) > 0) {
      pmtoendf[cc, ss] <- pmo[tvax + 6, ss, min(idi)]
    } else {
      pmtoendf[cc, ss] <- max(pmo[tvax + 6, ss, ])
    }
  }
}

# Get pervcentiles for plots
qtls <- as.data.frame(rowQuantiles((pmtoend),
  probs = c(0.025, 0.5, 0.975)
))
qtls$x <- endpoints_vec
iqtls <- as.data.frame(rowQuantiles((pmtoendi),
  probs = c(0.025, 0.5, 0.975)
))
iqtls$x <- endpoints_vec

# Define colors
colors <- c("Infection" = "orange", "Case" = "blue")

qtls_plot <- ggplot() +
  geom_ribbon(data = qtls, aes(x = x, ymin = `2.5%`, ymax = `97.5%`), fill = "blue", alpha = 0.5) +
  geom_line(data = qtls, aes(x = x, y = `50%`, col = "Case"), lwd = 1) +
  geom_ribbon(data = iqtls, aes(x = x, ymin = `2.5%`, ymax = `97.5%`), fill = "orange", alpha = 0.5) +
  geom_line(data = iqtls, aes(x = x, y = `50%`, col = "Infection"), lwd = 1) +
  labs(title = "", x = "Cases", y = "Person-months follow-up", color = "") +
  scale_x_continuous(breaks = seq(0, 400, 20)) +
  scale_y_continuous(breaks = seq(0, 6e5, 5e4), labels = scales::comma) + # , limits = c(0,600000))+
  scale_color_manual(values = colors) +
  theme_minimal() +
  theme(
    axis.text = element_text(colour = "black", size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "bottom",
    legend.text = element_text(colour = "black", size = 10, face = "bold")
  ) +
  geom_vline(xintercept = 18, color = "grey35", lwd = 0.5) +
  geom_vline(xintercept = 150, color = "red", lwd = 0.5) +
  geom_vline(xintercept = 399, color = "darkgreen", lwd = 0.5)


qtls_plot <- qtls_plot +
  annotate(
    geom = "text",
    x = 18,
    y = 620000,
    label = "VE ~90% ",
    color = "grey35",
    size = 4,
    fontface = "bold",
    hjust = -0.2
  ) +
  annotate(
    geom = "text",
    x = 150,
    y = 620000,
    label = "VE ~60% ",
    color = "red",
    size = 4,
    fontface = "bold",
    hjust = -0.2
  ) +
  annotate(
    geom = "text",
    x = 380,
    y = 620000,
    label = "VE ~50% ",
    color = "darkgreen",
    size = 4,
    fontface = "bold",
    hjust = 0.5
  ) +
  coord_cartesian(ylim = c(0, 6.2e5), clip = "off")



# save plot
save(qtls_plot, file = here("output", country, "pmonth_Fig.rdata"))

##################################################################
# Cllin
id <- which(is.na(cohort$cumclin))
AR <- cohort$cumclin
AR[id] <- 0
qtls <- as.data.frame(rowQuantiles(t(AR[, nruns, ]),
  probs = c(0.025, 0.5, 0.975)
))
qtls$x <- seq(0, length(mo) - 1, 1) - tvax

# Inf
id <- which(is.na(cohort$cuminc))
ARi <- cohort$cuminc
ARi[id] <- 0
qtlsi <- as.data.frame(rowQuantiles(t(ARi[, nruns, ]),
  probs = c(0.025, 0.5, 0.975)
))
qtlsi$x <- seq(0, length(mo) - 1, 1) - tvax

# fatal
id <- which(is.na(cohort$cumfatal))
ARf <- cohort$cumfatal
ARf[id] <- 0
qtlsf <- as.data.frame(rowQuantiles(t(ARf[, nruns, ]),
  probs = c(0.025, 0.5, 0.975)
))
qtlsf$x <- seq(0, length(mo) - 1, 1) - tvax

colors <- c("Infection" = "lightblue", "Case" = "Purple", "Fatal" = "black")


cumcases <- ggplot() +
  geom_line(data = qtlsi, aes(x = x, y = `50%`, col = "Infection"), lwd = 1) +
  geom_ribbon(data = qtlsi, aes(x = x, ymin = `2.5%`, ymax = `97.5%`), fill = "lightblue", alpha = 0.5) +
  geom_line(data = qtls, aes(x = x, y = `50%`, col = "Case"), lwd = 1) +
  geom_ribbon(data = qtls, aes(x = x, ymin = `2.5%`, ymax = `97.5%`), fill = "Purple", alpha = 0.2) +
  geom_line(data = qtlsf, aes(x = x, y = `50%`, col = "black"), lwd = 1) +
  geom_ribbon(data = qtlsf, aes(x = x, ymin = `2.5%`, ymax = `97.5%`), fill = "grey", alpha = 0.6) +
  labs(title = labc, x = "Months", y = "Cumulative cases", color = "") +
  geom_hline(yintercept = 150, linetype = "dashed", color = "grey", lwd = 1) +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(
    axis.text = element_text(colour = "black", size = 10, face = "bold"),
    axis.title = element_text(colour = "black", size = 12, face = "bold"),
    title = element_text(colour = "black", size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    # axis.text.x=element_text(angle=60, hjust=1),
    legend.position = "bottom",
    legend.text = element_text(colour = "black", size = 10, face = "bold")
  ) +
  scale_x_continuous(
    breaks = seq(0, length(mo) - 1, 5) - tvax,
    limits = c(0, 36)
  ) +
  ylim(0, max(qtlsi$`97.5%`[tvax + 36]) * 1.1)


save(cumcases, file = here("output", country, "cumCasesFig.rdata"))

#################################################################
## Countour plots of coverage vs sample size vs VE
#################################################################

# Attack rate vector
AR <- mean(cohort$ARU_clin[, nruns, tvax + 6])

# vector for VE to examine
ve <- seq(0.4, 0.9, 0.05)
nscen <- length(ve)

# Empty arrays
samplesize <- matrix(NA, nrow = nruns, ncol = nscen)
ttoend <- array(NA, c(nruns, nsim, nscen))
ttoendi <- array(NA, c(nruns, nsim, nscen))
ttoendf <- array(NA, c(nruns, nsim, nscen))
endps <- matrix(NA, nrow = nscen)


for (vv in 1:nscen) {
  VE_0 <- 0.3
  VE_1 <- ve[vv]

  p0 <- AR * (1 - VE_0) / (AR * (1 - VE_0) + AR)
  p1 <- AR * (1 - VE_1) / (AR * (1 - VE_1) + AR)

  # Re estimate endpoints according to VE H1
  eff.size <- ES.h(p0, p1)
  N <- pwr.p.test(h = eff.size, sig.level = 0.025, power = 0.90, alternative = "greater")
  endpoints <- ceiling(N$n)
  endps[vv] <- endpoints

  for (cc in 1:nruns) {
    samplesize[cc, vv] <- seeds[cc]
    for (ss in 1:nsim) {
      id <- which(cohort$cumclin[ss, cc, ] >= endpoints)
      idi <- which(cohort$cuminc[ss, cc, ] >= endpoints)
      idf <- which(cohort$cumfatal[ss, cc, ] >= endpoints)

      if (length(id) > 0) {
        ttoend[cc, ss, vv] <- min(id) - tvax
      } else {
        ttoend[cc, ss, vv] <- params$nt - tvax
      }

      if (length(idi) > 0) {
        ttoendi[cc, ss, vv] <- min(idi) - tvax
      } else {
        ttoendi[cc, ss, vv] <- params$nt - tvax
      }
      if (length(idf) > 0) {
        ttoendf[cc, ss, vv] <- min(idf) - tvax
      } else {
        ttoendf[cc, ss, vv] <- params$nt - tvax
      }
    }
  }
}


## Countour plot with clinical cases as endpoint
r <- as.data.frame(apply(ttoend, c(1, 3), quantile, .5, names = FALSE))

r$sample <- seq(1000, 100000, 1000)
colnames(r) <- c(seq(0.4, 0.9, 0.05), "Sample")
r_melt <- reshape2::melt(r, id = "Sample")
class(r_melt$Sample)
class(r_melt$variable)
class(r_melt$value)

r_melt$variable <- as.numeric(levels(r_melt$variable))[r_melt$variable] * 100

breaks_co <- c(4, 6, 12, 36, 60)

heat_cases <- ggplot(r_melt, aes(x = variable, y = Sample)) +
  geom_raster(aes(fill = value)) +
  geom_contour(aes(z = value), breaks = breaks_co, colour = "white", size = 0.2, alpha = 1) +
  geom_text_contour(aes(z = value),
    breaks = breaks_co, colour = "black",
    label.placer = label_placer_random(), skip = 0, stroke = 0.1
  ) +
  labs(
    title = labc,
    x = "Target vaccine efficacy (%)",
    y = "Sample size",
    fill = "Follow-up (months)"
  ) +
  scale_fill_gradientn(
    colors = c("orange", "purple", "blue"),
    breaks = c(6, 36, 60, 120, max(r[, 1:11])),
    limits = range(r_melt$value)
  ) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 10, face = "bold"),
    legend.position = "right", panel.background = element_blank(),
    title = element_text(colour = "black", size = 12, face = "bold"),
    axis.text = element_text(colour = "black", size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 7), legend.key = element_blank()
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(40, 90, 10))

save(heat_cases, file = here("output", country, "heat_cases.rdata"))



## Countour plot with infection cases as endpoint

r <- as.data.frame(apply(ttoendi, c(1, 3), quantile, .5, names = FALSE))

r$sample <- seq(1000, 100000, 1000)
colnames(r) <- c(seq(0.4, 0.9, 0.05), "Sample")
r_melt <- reshape2::melt(r, id = "Sample")
class(r_melt$Sample)
class(r_melt$variable)
class(r_melt$value)

r_melt$variable <- as.numeric(levels(r_melt$variable))[r_melt$variable] * 100

breaks_co <- c(4, 6, 12, 36, 60)

heat_inf <- ggplot(r_melt, aes(x = variable, y = Sample)) +
  geom_raster(aes(fill = value)) +
  geom_contour(aes(z = value), breaks = breaks_co, colour = "white", size = 0.2, alpha = 1) +
  geom_text_contour(aes(z = value),
    breaks = breaks_co, colour = "black",
    label.placer = label_placer_random(), skip = 0, stroke = 0.1
  ) +
  labs(
    title = labc,
    x = "Target vaccine efficacy (%)",
    y = "Sample size",
    fill = "Follow-up (months)"
  ) +
  scale_fill_gradientn(
    colors = c("seagreen3", "blue", "midnightblue"),
    breaks = c(6, 36, 60, 120, max(r[, 1:11])),
    limits = range(r_melt$value)
  ) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 10, face = "bold"),
    legend.position = "right", panel.background = element_blank(),
    title = element_text(colour = "black", size = 12, face = "bold"),
    axis.text = element_text(colour = "black", size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 7), legend.key = element_blank()
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(40, 90, 10))

save(heat_inf, file = here("output", country, "heat_inf.rdata"))


###########
## end of code
###########
