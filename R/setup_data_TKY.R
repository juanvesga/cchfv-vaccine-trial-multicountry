# 12/08/2022
# JF Vesga - London School of Hygiene and Tropical Medicine
# This script loads, organise, re-format and clean data for Turkey

## 1. Human data incidence
##############################################################################

# Start dates
temp_start_date <- as.Date("2004-01-15") # start data of data
end_date <- as.Date("2021-12-01")
temp_start_mo <- strftime(temp_start_date, format = "%Y-%m")
end_mo <- strftime(end_date, format = "%Y-%m")
format_date <- "%Y-%m-%d"

# Human incidence data
human_inc <- read.csv(here("data", country, "cchfTKY_human_monthly2004_2017.csv"), header = TRUE)
human_inc$cases <- round(human_inc$cases * 1)
human_inc$time <- 1 + (interval(temp_start_date, human_inc$date) %/% months(1))
human_month <- human_inc %>%
  mutate(cases = ifelse(cases == 0, NA, as.numeric(paste(cases))))


human_month_fit <- human_month %>%
  filter((date > "2005-01-01" &
    date < "2006-01-01") | (date > "2016-01-01" &
    date < "2018-01-01"))

year_vec <- seq(2004, 2021, by = 1)

## 2. Environmental drivers data
##############################################################################

drivers <- read.csv(here("data", country, "drivers.csv"), header = TRUE, sep = ",")
temp_month <- drivers

temp_month$day <- as.Date(drivers$month, format = "%d/%m/%Y")
temp_month$month <- strftime(temp_month$day, format = "%Y-%m")

soil_t <- drivers$soil
rh <- drivers$rh
air_temp <- drivers$air2m


# Estimate stauration deficit
sat_def <- (1 - rh) * 4.9463 * exp(0.0621 * (air_temp))
df <- data.frame(x = air_temp, y = sat_def)

# Fit lienar model to find relation between airT and Sat-def (for setting caps)
fit3 <- lm(y ~ poly(x, 3, raw = TRUE), data = df)

# Load NDVI
NDVI <- read.csv(here("data", country, "ndvi.csv"), header = TRUE, sep = ",")
ndvi <- NDVI$ndvi * 100



# 3. Import human prevalence by occupation farming and other occupations
##############################################################################
prev_other <- read.csv(here("data", country, "others.csv"), header = TRUE)
prev_other$low_ci <- round(binconf(prev_other$Posa, prev_other$denoma)[, 2], digit = 4)
prev_other$up_ci <- round(binconf(prev_other$Posa, prev_other$denoma)[, 3], digit = 4)

prev_farmer <- read.csv(here("data", country, "farmers.csv"), header = TRUE)
prev_farmer$low_ci <- round(binconf(prev_farmer$Posa, prev_farmer$denoma)[, 2], digit = 4)
prev_farmer$up_ci <- round(binconf(prev_farmer$Posa, prev_farmer$denoma)[, 3], digit = 4)

fstartdateH <- as.Date(prev_farmer$fstart, format = "%d/%m/%Y")
fenddateH <- as.Date(prev_farmer$fend, format = "%d/%m/%Y")

# 4. Livestock age-stratified prevalence from 2009 (5 age-groups)
#############################################################################
prev_liv <- read.csv(here("data", country, "livestock_seroprev.csv"), header = TRUE)
prev_liv$low_ci <- round(binconf(prev_liv$pos_igg, prev_liv$denom)[, 2], digit = 4)
prev_liv$up_ci <- round(binconf(prev_liv$pos_igg, prev_liv$denom)[, 3], digit = 4)

prev_liv_all <- read.csv(here("data", country, "livestock_seroprev_wo_age.csv"), header = TRUE)
prev_liv_all$low_ci <- round(binconf(prev_liv_all$pos_igg, prev_liv_all$denom)[, 2], digit = 4)
prev_liv_all$up_ci <- round(binconf(prev_liv_all$pos_igg, prev_liv_all$denom)[, 3], digit = 4)

fstartdateL <- as.Date(prev_liv_all$fstart, format = "%d/%m/%Y")
fenddateL <- as.Date(prev_liv_all$fend, format = "%d/%m/%Y")



# define the dates of the field work

field_work_startH<-interval(temp_start_date,fstartdateH)%/% months(1) 
field_work_endH <- interval(temp_start_date,fenddateH) %/% months(1) 

field_work_startLall<-interval(temp_start_date,fstartdateL)%/% months(1) 
field_work_endLall <- interval(temp_start_date,fenddateL) %/% months(1) 

field_work_startL<-field_work_startLall[1]
field_work_endL<-field_work_endLall[1]

# 5. Get all data in a list structure
###############################################################################

# Index for monthly human cases
tmp <- which(!is.na(human_month$cases))
tmpf <- which(!is.na(human_month_fit$cases))


observations <- list(
  start_date = temp_start_date,
  start_month = temp_start_mo,
  index_mo_cases = human_month$time[tmp],
  cases_human_mo = human_month$cases[tmp],
  index_mo_fitcases = human_month_fit$time[tmpf],
  fitcases_human_mo = human_month_fit$cases[tmpf],
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