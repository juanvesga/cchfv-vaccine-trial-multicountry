# 12/08/2022
# JF Vesga - London School of Hygiene and Tropical Medicine
# This script loads, organise, re-format and clean data for Afghanistan

## 1. Human data incidence
##############################################################################

temp_start_date <- as.Date("2008-04-01") # start data of data
end_date <- as.Date("2008-10-31")
temp_start_mo <- strftime(temp_start_date, format = "%Y-%m")
end_mo <- strftime(end_date, format = "%Y-%m")
fwork_start <- as.Date("2009-08-30")
fwork_end <- as.Date("2009-08-15")
fwork_mid <- as.Date("2009-08-31")

fwork_starty <- as.Date("2008-09-01")
fwork_endy <- as.Date("2009-07-31")
fwork_midy <- as.Date("2009-02-28")
format_date <- "%d/%m/%Y"


# Get weekly data 2008
human_inc1 <- read.csv(here("data", country, "cchf_human_weekly2008.csv"),
  header = TRUE
)

human_inc1$date <- as.Date(human_inc1$date, format = "%d/%m/%Y")

human_inc1$month <- strftime(human_inc1$date, format = "%Y-%m")

# Aggregtae by month
human_inc1 <- human_inc1 %>%
  group_by(month) %>%
  summarise(cases = sum(obs))

human_inc1 <- human_inc1[human_inc1$month >= temp_start_mo, ]

human_inc1$time <- seq(1, NROW(human_inc1$cases), 1)

# Get monthly data
human_inc2 <- read.csv(here("data", country, "cchf_human_monthly2017_2018.csv"),
  header = TRUE
)
human_inc2$date <- as.Date(human_inc2$date, format = "%d/%m/%Y")

human_inc2$month <- strftime(human_inc2$date, format = "%Y-%m")

human_inc2$time <- 1 + (interval(temp_start_date, human_inc2$date) %/% months(1))

human_inc2 <- human_inc2 %>%
  rename(cases = obs) %>%
  select(month, cases, time)

# Final monthly cases data
human_month <- rbind(human_inc1, human_inc2)
human_month <- human_month %>%
  mutate(cases = ifelse(cases == 0, NA, as.numeric(paste(cases))))

# Get Yearly data
human_inc3 <- read.csv(here("data", country, "cchf_human_yearly2009_2018.csv"),
  header = TRUE
)
human_inc3$date <- as.Date(human_inc3$date, format = "%d/%m/%Y")
human_inc3$year <- strftime(human_inc3$date, format = "%Y")
human_year <- human_inc3 %>%
  rename(cases = obs)

human_year <- human_year %>%
  mutate(cases = ifelse(cases == 0, NA, as.numeric(paste(cases))))
human_year <- human_year[human_year$year >= 2008, ]

years_fit <- human_year$year


# Get Yearly deaths data
human_mrt <- read.csv(here("data", country, "cchf_human_deaths_yearly2008_2018.csv"),
  header = TRUE
)
human_mrt$date <- as.Date(human_mrt$date, format = "%d/%m/%Y")
human_mrt$year <- strftime(human_mrt$date, format = "%Y")
humanD_year <- human_mrt %>%
  rename(deaths = obs)

humanD_year <- humanD_year %>%
  mutate(deaths = ifelse(deaths == 0, NA, as.numeric(paste(deaths)))) %>%
  select(date, year, deaths)

humanD_year <- humanD_year[humanD_year$year >= 2008, ]

yearsD_fit <- humanD_year$year

year_vec <- seq(2008, 2019, by = 1)

## 2. Environmental drivers data
##############################################################################

temp_day <- read.csv(here("data", country, "soil_temp.csv"),
  header = TRUE, sep = ","
)

temp_day$date <- as.Date(temp_day$vector_days, format = format_date)

temp_day <- temp_day[temp_day$date >= temp_start_date, ]

temp_day$month <- strftime(temp_day$date, format = "%Y-%m")

temp_month <- temp_day %>%
  group_by(month) %>%
  summarise(soil_t = mean(Temperature))
temp_month$time <- c(1:nrow(temp_month))

soil_t <- temp_month$soil_t


## Relative humidity
rel_hum <- read.csv(here("data", country, "rel_hum.csv"),
  header = TRUE, sep = ","
)
rel_hum$date <- strftime(rel_hum$vector_month, format = "%Y-%m")
id <- rel_hum$date %in% temp_day$month
rel_hum <- rel_hum[id, ]

rh <- rel_hum$rh


## Air temp
airT <- read.csv(here("data", country, "air_temp.csv"),
  header = TRUE, sep = ","
)
airT$date <- strftime(airT$vector_month, format = "%Y-%m")
id <- airT$date %in% temp_day$month
airT <- airT[id, ]

air_temp <- airT$Temperature



# Estimate saturation deficit
sat_def <- (1 - rh) * 4.9463 * exp(0.0621 * (air_temp))
df <- data.frame(x = air_temp, y = sat_def)

# Linear func to relation satdef temp
fit3 <- lm(y ~ poly(x, 3, raw = TRUE), data = df)


# Load NDVI
NDVI <- read.csv(here("data", country, "ndvi.csv"), header = TRUE, sep = ",")
id <- NDVI$vector_month %in% temp_day$month
NDVI <- NDVI[id, ]
ndvi <- NDVI$ndvi * 100


# 3. Import human prevalence by occupation farming and other occupations (from CCHF 2009 dataset)
#################################################################################

# define the dates of the field work August 2009
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

# 4. Livestock age-stratified prevalence from 2009 (5 age-groups)
##############################################################################
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
  index_yr_cases = which(year_vec %in% years_fit),
  index_yr_deaths = which(year_vec %in% yearsD_fit),
  cases_human_yr = human_year$cases,
  deaths_human_yr = humanD_year$deaths,
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
  time = c(observations$index_mo_cases, observations$index_yr_cases * 12),
  obs = c(observations$cases_human_mo, observations$cases_human_yr)
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
#############
## End of code
#############