#GHG perturbation model
# Model follows:
#Frolking et al. 2006 J. Geophys. Res.
# Neubauer & Megonigal 2015 Ecosystems
#Code by: Shizhou Ma
# Version Date: October 2025

# This script applies a greenhouse gas perturbation model to wetland CH4 and CO2
# flux inputs and estimates:
#   1. Instantaneous radiative forcing
#   2. Cumulative radiative forcing
#   3. Net radiative forcing savings relative to baseline conditions

# Required packages to run the code:
# tidyverse V2.0.0

# Set your won file path here
input_path <- "FILLIN"
#the input data unit is C-CO2, C-CH4
GHGinputssp126 <- read.csv(file.path(input_path, "GHGpertubinputssp126.csv"))
GHGinputssp245 <- read.csv(file.path(input_path, "GHGpertubinputssp245.csv"))
GHGinputssp370 <- read.csv(file.path(input_path, "GHGpertubinputssp370.csv"))
GHGinputssp585 <- read.csv(file.path(input_path, "GHGpertubinputssp585.csv"))
#change the unit from kg/ha to g/m2
datasets <- c("GHGinputssp126",
              "GHGinputssp245",
              "GHGinputssp370",
              "GHGinputssp585")

# columns to rescale
cols_to_scale <- c("meanch4", "sdch4", "meanco2", "sdco2")

for (nm in datasets) {
  # pull the data.frame
  df <- get(nm)
  
  # multiply those columns by 0.1
  df[ , cols_to_scale] <- df[ , cols_to_scale] * 0.1
  
  # put it back into your workspace
  assign(nm, df)
}

GHGpertubationconversion_dataset <- function(data, type = c("ins", "cu"), runs = 1000, stepsPerYear = 5) {
  # Conversion factors
  co2m <- 44/12
  ch4m <- 16/12
  
  # parse a decade string into start and end years
  parse_decade <- function(decade_string) {
    parts <- strsplit(decade_string, ":")[[1]]
    # Convert each part into a year number
    parse_year <- function(x) {
      x <- trimws(x)
      if (grepl("BC", x, ignore.case = TRUE)) {
        return(-as.numeric(gsub("BC", "", x, ignore.case = TRUE)))
      } else if (grepl("AD", x, ignore.case = TRUE)) {
        return(as.numeric(gsub("AD", "", x, ignore.case = TRUE)))
      } else {
        return(as.numeric(x))
      }
    }
    start_year <- parse_year(parts[1])
    end_year   <- parse_year(parts[2])
    return(c(start_year, end_year))
  }
  
  # Process the dataset to extract the time intervals and number of years per interval.
  n_intervals <- nrow(data)
  year_counts <- numeric(n_intervals)
  start_end   <- matrix(NA, nrow = n_intervals, ncol = 2)
  
  for (i in 1:n_intervals) {
    yrs <- parse_decade(data$decade[i])
    start_end[i, ] <- yrs
    year_counts[i] <- yrs[2] - yrs[1] + 1  # inclusive of both start and end
  }
  total_years <- sum(year_counts)
  total_steps <- total_years * stepsPerYear
  
  # Create an overall time vector spanning the simulation.
  overall_start <- start_end[1, 1]
  overall_end   <- start_end[n_intervals, 2]
  time <- seq(overall_start, overall_end, length.out = total_steps)
  
  # Preallocate matrices to store outputs for each run.
  net_out       <- matrix(NA, nrow = total_steps, ncol = runs)
  ch4_out       <- matrix(NA, nrow = total_steps, ncol = runs)
  co2_out       <- matrix(NA, nrow = total_steps, ncol = runs)
  totalbase_out <- matrix(NA, nrow = total_steps, ncol = runs)
  savings_out   <- matrix(NA, nrow = total_steps, ncol = runs)
  temp_out      <- matrix(NA, nrow = total_steps, ncol = runs)
  tempbase_out  <- matrix(NA, nrow = total_steps, ncol = runs)
  
  # Constants and Parameters (same as in the original code)
  pool1lt <- 0.2173; pool2lt <- 0.224; pool3lt <- 0.2824; pool4lt <- 0.2763
  radeff.co2 <- 1.7509e-15  # radiative efficiency for CO2 [W m-2 (kg CO2)-1]
  radeff.ch4 <- 1.2758e-13  # radiative efficiency for CH4 [W m-2 (kg CH4)-1]
  radeff.temp <- 0.813      # temperature conversion factor [W m-2 K-1]
  ie.ch4 <- 1.65            # indirect effect multiplier for CH4
  ch4decay <- 0.98400034
  co2pool1decay <- 1
  co2pool2decay <- 0.999493029
  co2pool3decay <- 0.9945415
  co2pool4decay <- 0.95459472
  
  # For the baseline simulation, assume pre-conversion conditions from the first row.
  pre_co2   <- data$meanco2[1]
  pre_co2sd <- data$sdco2[1]
  pre_ch4   <- data$meanch4[1]
  pre_ch4sd <- data$sdch4[1]
  
  # Monte Carlo simulation runs
  for (run in 1:runs) {
    # Build annual time series for the perturbed simulation from the dataset.
    co2_annual <- numeric(0)
    ch4_annual <- numeric(0)
    for (j in 1:n_intervals) {
      n_years <- year_counts[j]
      co2_vals <- rnorm(n_years,
                        mean = data$meanco2[j] * co2m,
                        sd   = data$sdco2[j] * co2m)
      ch4_vals <- rnorm(n_years,
                        mean = data$meanch4[j] * ch4m,
                        sd   = data$sdch4[j] * ch4m)
      co2_annual <- c(co2_annual, co2_vals)
      ch4_annual <- c(ch4_annual, ch4_vals)
    }
    
    # Baseline simulation: use pre conversion conditions over the entire simulation period.
    co2_annual_base <- rnorm(total_years, mean = pre_co2 * co2m, sd = pre_co2sd * co2m)
    ch4_annual_base <- rnorm(total_years, mean = pre_ch4 * ch4m, sd = pre_ch4sd * ch4m)
    
    # Convert annual values into subannual time steps.
    co2_simulated_step      <- rep(co2_annual / stepsPerYear, each = stepsPerYear)
    ch4_simulated_step      <- rep(ch4_annual / stepsPerYear, each = stepsPerYear)
    co2_simulated_stepbase  <- rep(co2_annual_base / stepsPerYear, each = stepsPerYear)
    ch4_simulated_stepbase  <- rep(ch4_annual_base / stepsPerYear, each = stepsPerYear)
    
    # Baseline simulation computation (using pre-conversion conditions)
    ch4vecbase <- numeric(total_steps)
    ch4vecbase[1] <- ch4_simulated_stepbase[1]
    for (t in 2:total_steps) {
      ch4vecbase[t] <- ch4_simulated_stepbase[t] + ch4vecbase[t-1] * ch4decay
    }
    ch4invfinalbase <- ch4vecbase
    ch4oxfbase <- c(0, (1 - ch4decay) * 44/16 * ch4invfinalbase)
    ch4oxffinalbase <- ch4oxfbase[1:total_steps]
    
    co2vec1base <- numeric(total_steps)
    co2vec2base <- numeric(total_steps)
    co2vec3base <- numeric(total_steps)
    co2vec4base <- numeric(total_steps)
    co2_simulated_stepfbase <- co2_simulated_stepbase - ch4_simulated_stepbase * 44/16
    for (t in 2:total_steps) {
      co2vec1base[t] <- (co2_simulated_stepfbase[t] + ch4oxffinalbase[t]) * pool1lt +
        co2vec1base[t-1] * co2pool1decay
      co2vec2base[t] <- (co2_simulated_stepfbase[t] + ch4oxffinalbase[t]) * pool2lt +
        co2vec2base[t-1] * co2pool2decay
      co2vec3base[t] <- (co2_simulated_stepfbase[t] + ch4oxffinalbase[t]) * pool3lt +
        co2vec3base[t-1] * co2pool3decay
      co2vec4base[t] <- (co2_simulated_stepfbase[t] + ch4oxffinalbase[t]) * pool4lt +
        co2vec4base[t-1] * co2pool4decay
    }
    co2invfinalbase <- co2vec1base + co2vec2base + co2vec3base + co2vec4base
    
    if (type[1] == "ins") {
      ch4RFbase <- ie.ch4 * ch4invfinalbase * radeff.ch4 * 1e-3
      co2RFbase <- radeff.co2 * co2invfinalbase * 1e-3
      netRFbase <- ch4RFbase + co2RFbase
      tempbase  <- netRFbase * radeff.temp
    } else if (type[1] == "cu") {
      ch4RFbase   <- ie.ch4 * ch4invfinalbase * radeff.ch4 * 1e-3
      co2RFbase   <- radeff.co2 * co2invfinalbase * 1e-3
      ch4RFstepbase <- ch4RFbase / stepsPerYear
      co2RFstepbase <- co2RFbase / stepsPerYear
      ch4RFbase   <- cumsum(ch4RFstepbase)
      co2RFbase   <- cumsum(co2RFstepbase)
      netRFbase   <- ch4RFbase + co2RFbase
      tempbase    <- netRFbase * radeff.temp
    }
    
    # Perturbed simulation using dataset parameters
    ch4vec <- numeric(total_steps)
    ch4vec[1] <- ch4_simulated_step[1]
    for (t in 2:total_steps) {
      ch4vec[t] <- ch4_simulated_step[t] + ch4vec[t-1] * ch4decay
    }
    ch4invfinal <- ch4vec
    ch4oxf <- c(0, (1 - ch4decay) * 44/16 * ch4invfinal)
    ch4oxffinal <- ch4oxf[1:total_steps]
    
    co2vec1 <- numeric(total_steps)
    co2vec2 <- numeric(total_steps)
    co2vec3 <- numeric(total_steps)
    co2vec4 <- numeric(total_steps)
    co2_simulated_stepf <- co2_simulated_step - ch4_simulated_step * 44/16
    for (t in 2:total_steps) {
      co2vec1[t] <- (co2_simulated_stepf[t] + ch4oxffinal[t]) * pool1lt +
        co2vec1[t-1] * co2pool1decay
      co2vec2[t] <- (co2_simulated_stepf[t] + ch4oxffinal[t]) * pool2lt +
        co2vec2[t-1] * co2pool2decay
      co2vec3[t] <- (co2_simulated_stepf[t] + ch4oxffinal[t]) * pool3lt +
        co2vec3[t-1] * co2pool3decay
      co2vec4[t] <- (co2_simulated_stepf[t] + ch4oxffinal[t]) * pool4lt +
        co2vec4[t-1] * co2pool4decay
    }
    co2invfinal <- co2vec1 + co2vec2 + co2vec3 + co2vec4
    
    if (type[1] == "ins") {
      ch4RF <- ie.ch4 * ch4invfinal * radeff.ch4 * 1e-3
      co2RF <- radeff.co2 * co2invfinal * 1e-3
      netRF <- ch4RF + co2RF
      temp  <- netRF * radeff.temp
    } else if (type[1] == "cu") {
      ch4RF   <- ie.ch4 * ch4invfinal * radeff.ch4 * 1e-3
      co2RF   <- radeff.co2 * co2invfinal * 1e-3
      ch4RFstep <- ch4RF / stepsPerYear
      co2RFstep <- co2RF / stepsPerYear
      ch4RF   <- cumsum(ch4RFstep)
      co2RF   <- cumsum(co2RFstep)
      netRF   <- ch4RF + co2RF
      temp    <- netRF * radeff.temp
    }
    
    net_savings <- netRF - netRFbase
    
    # Save results from this run
    net_out[, run]       <- netRF
    ch4_out[, run]       <- ch4RF
    co2_out[, run]       <- co2RF
    totalbase_out[, run] <- netRFbase
    savings_out[, run]   <- net_savings
    tempbase_out[, run]  <- tempbase
    temp_out[, run]      <- temp
  }
  
  # Calculate summary statistics across all runs
  output_df <- data.frame(
    time = time,
    net_RF       = apply(net_out, 1, mean),
    net_sd         = apply(net_out, 1, sd),
    co2_cum        = apply(co2_out, 1, mean),
    co2_sd         = apply(co2_out, 1, sd),
    ch4_cum        = apply(ch4_out, 1, mean),
    ch4_sd         = apply(ch4_out, 1, sd),
    totalbase_cum  = apply(totalbase_out, 1, mean),
    totalbase_sd   = apply(totalbase_out, 1, sd),
    savings_cum    = apply(savings_out, 1, mean),
    savings_sd     = apply(savings_out, 1, sd),
    temp_mean      = apply(temp_out, 1, mean),
    temp_sd        = apply(temp_out, 1, sd),
    tempbase_mean  = apply(tempbase_out, 1, mean),
    tempbase_sd    = apply(tempbase_out, 1, sd)
  )
  
  return(output_df)
}
#Run the built function with rondom forest estimated wetland GHG fluxes
result126 <- GHGpertubationconversion_dataset(GHGinputssp126, type = "ins", runs = 100)
result245 <- GHGpertubationconversion_dataset(GHGinputssp245, type = "ins", runs = 100)
result370 <- GHGpertubationconversion_dataset(GHGinputssp370, type = "ins", runs = 100)
result585 <- GHGpertubationconversion_dataset(GHGinputssp585, type = "ins", runs = 100)

result126c <- GHGpertubationconversion_dataset(GHGinputssp126, type = "cu", runs = 100)
result245c <- GHGpertubationconversion_dataset(GHGinputssp245, type = "cu", runs = 100)
result370c <- GHGpertubationconversion_dataset(GHGinputssp370, type = "cu", runs = 100)
result585c <- GHGpertubationconversion_dataset(GHGinputssp585, type = "cu", runs = 100)

# For time <= 0, take only one dataset (here pertubssp126) 
common_paleo <- result126 %>%
  filter(time <= 1900) %>%
  mutate(scenario = "paleo")

common_hist <- result126 %>%
  filter(time>1900, time <= 2020) %>%
  mutate(scenario = "hist")

# For time > 0, add the appropriate scenario label for each dataset
ssp126_post <- result126 %>%
  filter(time > 2020) %>%
  mutate(scenario = "ssp126")

ssp245_post <- result245 %>%
  filter(time > 2020) %>%
  mutate(scenario = "ssp245")

ssp370_post <- result370 %>%
  filter(time > 2020) %>%
  mutate(scenario = "ssp370")

ssp585_post <- result585 %>%
  filter(time > 2020) %>%
  mutate(scenario = "ssp585")

# Combine the datasets together
GHGpertubmerged<- bind_rows(common_paleo, common_hist, ssp126_post, ssp245_post, ssp370_post, ssp585_post)

