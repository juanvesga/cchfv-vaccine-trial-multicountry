# 12/08/2022
# JF Vesga - London School of Hygiene and Tropical Medicine
# This script loads, organise, re-format and clean data for South Africa

## 1. Human data incidence
##############################################################################

# Start dates
temp_start_date <- as.Date("2000-01-15") # start data of data
end_date <- as.Date("2021-12-15")
temp_start_mo <- strftime(temp_start_date, format = "%Y-%m")
end_mo <- strftime(end_date, format = "%Y-%m")
fwork_start <- as.Date("2017-05-01")
fwork_end <- as.Date("2017-11-30")
fwork_mid <- as.Date("2017-07-15")

fwork_starty <- as.Date("2017-05-01")
fwork_endy <- as.Date("2018-05-31")
fwork_midy <- as.Date("2017-11-15")
format_date <- "%Y-%m-%d"

# Load Human incidence data
human_inc <- read.csv(here("data", country, "cchfSA_human_monthly2000_2021.csv"),
  header = TRUE
)


human_inc$time <- 1 + (interval(temp_start_date, human_inc$date) %/% months(1))

human_month <- human_inc %>%
  mutate(cases = ifelse(cases == 0, NA, as.numeric(paste(cases))))


year_vec <- seq(2000, 2021, by = 1)



## 2. Environmental drivers data
##############################################################################

# Load drivers data
drivers <- read.csv(here("data", country, "drivers.csv"), header = TRUE, sep = ",")
temp_month <- drivers

# Set in month calendar format
temp_month$day <- as.Date(drivers$month, format = "%d/%m/%Y")
temp_month$month <- strftime(temp_month$day, format = "%Y-%m")

# Pass to independent structures
soil_t <- drivers$soil
rh <- drivers$rh
air_temp <- drivers$air2m

# Calculate stauration deficit from RH and air T

sat_def <- (1 - rh) * 4.9463 * exp(0.0621 * (air_temp))
df <- data.frame(x = air_temp, y = sat_def)

# Fit lienar model to find relation between airT and Sat-def (for setting caps)
fit3 <- lm(y ~ poly(x, 3, raw = TRUE), data = df)

# Load NDVI data
NDVI <- read.csv(here("data", country, "ndvi.csv"), header = TRUE, sep = ",")
ndvi <- NDVI$ndvi * 100


# 3. Import human prevalence by occupation farming and other occupations
#########################################################################

# define the dates of the field work

field_work_start <- (interval(temp_start_date, fwork_start) %/% months(1))
field_work_end <- (interval(temp_start_date, fwork_end) %/% months(1))
field_work_mid <- (interval(temp_start_date, fwork_mid) %/% months(1))

field_work_mid_date <- (temp_start_date %m+% months(field_work_mid))

prev_other <- read.csv(here("data", country, "others.csv"), header = TRUE)
prev_other$low_ci <- round(binconf(prev_other$Posa, prev_other$denoma)[2], digit = 4)
prev_other$up_ci <- round(binconf(prev_other$Posa, prev_other$denoma)[3], digit = 4)

prev_farmer <- read.csv(here("data", country, "farmers.csv"), header = TRUE)
prev_farmer$low_ci <- round(binconf(prev_farmer$Posa, prev_farmer$denoma)[2], digit = 4)
prev_farmer$up_ci <- round(binconf(prev_farmer$Posa, prev_farmer$denoma)[3], digit = 4)

# 4. Livestock age-stratified prevalence
###############################################################################
prev_liv <- read.csv(here("data", country, "livestock_seroprev.csv"), header = TRUE)
prev_liv$low_ci <- round(binconf(prev_liv$pos_igg, prev_liv$denom)[, 2], digit = 4)
prev_liv$up_ci <- round(binconf(prev_liv$pos_igg, prev_liv$denom)[, 3], digit = 4)

prev_liv_all <- read.csv(here("data", country, "livestock_seroprev_wo_age.csv"), header = TRUE)
prev_liv_all$low_ci <- round(binconf(prev_liv_all$pos_igg, prev_liv_all$denom)[, 2], digit = 4)
prev_liv_all$up_ci <- round(binconf(prev_liv_all$pos_igg, prev_liv_all$denom)[, 3], digit = 4)


# 5. Get all data in a list structure
###############################################################################

# Index for monthly human cases
tmp <- which(!is.na(human_month$cases))

observations <- list(
  start_date = temp_start_date,
  start_month = temp_start_mo,
  index_mo_cases = human_month$time[tmp],
  cases_human_mo = human_month$cases[tmp],
  prev_liv_age_pos_IgG = prev_liv$pos_igg,
  prev_liv_age_denom = prev_liv$denom,
  prev_liv_age = prev_liv$prop,
  prev_liv_age_low_ci = prev_liv$low_ci,
  prev_liv_age_up_ci = prev_liv$up_ci,
  prev_liv_all_pos_IgG = prev_liv_all$pos_igg,
  prev_liv_all_denom = prev_liv_all$denom,
  prev_liv_all = prev_liv_all$prop,
  prev_liv_all_low_ci = prev_liv_all$low_ci,
  prev_liv_all_up_ci = prev_liv_all$up_ci,
  prev_farmer_Posa = prev_farmer$Posa,
  prev_farmer_denoma = prev_farmer$denoma,
  prev_farmer = prev_farmer$prev,
  prev_farmer_low_ci = prev_farmer$low_ci,
  prev_farmer_up_ci = prev_farmer$up_ci,
  prev_other_Posa = prev_other$Posa,
  prev_other_denoma = prev_other$denoma,
  prev_other = prev_other$prev,
  prev_other_low_ci = prev_other$low_ci,
  prev_other_up_ci = prev_other$up_ci
)

tmp <- data.frame(
  time = observations$index_mo_cases,
  obs = observations$cases_human_mo
)

data1_human <- tmp[order(tmp$time), ]

data2_human <- data.frame(
  time = c(field_work_start:field_work_end),
  num_f = observations$prev_farmer_Posa,
  deno_f = observations$prev_farmer_denoma,
  num_o = observations$prev_other_Posa,
  deno_o = observations$prev_other_denoma
)

# . 6 Functions for setting temperature cap in tick activity according to driver
#############################################################################

# Find FOI factor according to Soil temperature
temp_cap <- 30
temp_threshold <- 12
if (driver == "sat_def") {
  environment_factor <- sat_def

  temp_threshold <- as.numeric(predict(fit3, data.frame(x = 12)))

  temp_cap <- as.numeric(predict(fit3, data.frame(x = 30)))

  envrdriver_foi_func <- function(y, temp_foi_factor, min_limit = min(sat_def)) {
    if (y >= temp_threshold) {
      x <- (temp_foi_factor * (min(y, temp_cap) - min_limit))
    } else if (y < temp_threshold) {
      x <- temp_foi_factor * 1
    }

    return(x)
  }
} else if (driver == "soil_t") {
  environment_factor <- soil_t

  envrdriver_foi_func <- function(y, temp_foi_factor, min_limit = min(soil_t)) {
    if (y >= temp_threshold) {
      x <- (temp_foi_factor * (min(y, temp_cap) - min_limit))
    } else if (y < temp_threshold) {
      x <- temp_foi_factor * 1
    }


    return(x)
  }
} else if (driver == "ndvi") {
  environment_factor <- ndvi

  envrdriver_foi_func <- function(y, temp_foi_factor, min_limit = min(ndvi)) {
    x <- temp_foi_factor * (y - min_limit)

    return(x)
  }
} else if (driver == "rel_hum") {
  environment_factor <- (1 - rh) * 100

  envrdriver_foi_func <- function(y, temp_foi_factor, min_limit = min(ndvi)) {
    x <- temp_foi_factor * (y - min_limit)

    return(x)
  }
}

#########
# End of Code
#########
