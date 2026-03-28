# GWP* CO2-equivalent conversion
# Code by: Shizhou Ma
# Version Date: October 2025

# This script combines wetland CH4 and CO2 sequestration model outputs and
# converts CH4 fluxes into CO2-equivalent values using the GWP* approach.
# It generates:
#   1. Combined CH4 and CO2 sequestration datasets for historical and future scenarios
#   2. Area-weighted annual CH4 and CO2 sequestration summaries
#   3. GWP*-based CO2-equivalent budgets
#
# Required packages to run the code:
# tidyverse V2.0.0, tidymodels V1.3.0, ranger V0.17.0,
# lubridate V1.9.4, vip V0.4.1, ggpubr V0.6.0, DALEXtra V2.2.0,
# doParallel V1.0.17, ggtext V0.1.2, patchwork V1.3.0, ggpmisc V0.6.0


#Fill in your own file path
setwd_path <- "FILL-IN"
# CH4 input files
ch4paleo  <- read.csv(paste0(setwd_path, "CH4randomforestpredictpre1984.csv"))
ch4       <- read.csv(paste0(setwd_path, "CH4randomforestpredictcontemporary.csv"))
ch4ssp126 <- read.csv(paste0(setwd_path, "CH4randomforestpredictssp126.csv"))
ch4ssp245 <- read.csv(paste0(setwd_path, "CH4randomforestpredictssp245.csv"))
ch4ssp370 <- read.csv(paste0(setwd_path, "CH4randomforestpredictssp370.csv"))
ch4ssp585 <- read.csv(paste0(setwd_path, "CH4randomforestpredictssp585.csv"))

# CO2 input files
cseqpaleo  <- read.csv(paste0(setwd_path, "CO2randomforestpredictpre1984.csv"))
cseq       <- read.csv(paste0(setwd_path, "CO2randomforestpredictcontemporary.csv"))
cseqssp126 <- read.csv(paste0(setwd_path, "CO2randomforestpredictssp126.csv"))
cseqssp245 <- read.csv(paste0(setwd_path, "CO2randomforestpredictssp245.csv"))
cseqssp370 <- read.csv(paste0(setwd_path, "CO2randomforestpredictssp370.csv"))
cseqssp585 <- read.csv(paste0(setwd_path, "CO2randomforestpredictssp585.csv"))




# Inner join for paleo data (ch4 and cseq)
paleo <- inner_join(ch4paleo, cseqpaleo, by = c("Wetland.ID", "Year")) %>%
  select(-any_of(c("temporal.x", "temporal.y"))) %>%
  mutate(temporal = "historical")

historical <- inner_join(ch4, cseq, by = c("Wetland.ID", "Year")) %>%
  select(-any_of(c("temporal.x", "temporal.y"))) %>%
  mutate(temporal = "historical")


# Inner joins for SSP data
ssp126 <- inner_join(ch4ssp126, cseqssp126, by = c("Wetland.ID", "Year")) %>%
  select(-any_of(c("temporal.x", "temporal.y"))) %>%
  mutate(temporal = "ssp126")
ssp245 <- inner_join(ch4ssp245, cseqssp245, by = c("Wetland.ID", "Year")) %>%
  select(-any_of(c("temporal.x", "temporal.y"))) %>%
  mutate(temporal = "ssp245")
ssp370 <- inner_join(ch4ssp370, cseqssp370, by = c("Wetland.ID", "Year")) %>%
  select(-any_of(c("temporal.x", "temporal.y"))) %>%
  mutate(temporal = "ssp370")
ssp585 <- inner_join(ch4ssp585, cseqssp585, by = c("Wetland.ID", "Year")) %>%
  select(-any_of(c("temporal.x", "temporal.y"))) %>%
  mutate(temporal = "ssp585")


# Combine all joined datasets into one data frame
inputforcbudget <- bind_rows(paleo, historical, 
                             ssp126,
                             ssp245, ssp370, ssp585) %>%
  arrange(Wetland.ID)

#remove several rewet sites
sitespred <- read.csv(paste0(setwd_path, "cmodelpredictionsinput.csv"))

# Get the AREA info for each wetland
AREA <- sitespred %>% 
  select(Wetland.ID, AREA) %>% 
  distinct()
AREAtotal <- sum(AREA$AREA)
inputforcbudget <- inputforcbudget %>%
  left_join(AREA, by = "Wetland.ID")

#Aggregate Annual Fluxes
aggregated_annual_summary <- inputforcbudget %>%
  group_by(Year, temporal) %>%
  summarize(
    # Sum fluxes weighted by each site's AREA
    total_annualCH4 = sum(mean_annualch4 * AREA, na.rm = TRUE),
    total_cseqmean  = sum(mean_cseq * AREA,      na.rm = TRUE),
    # Combine standard errors using Root-Sum-of-Squares (RSS)
    total_annualSE  = sqrt(sum((sd_annualch4 * AREA)^2, na.rm = TRUE)),
    total_cseqSE    = sqrt(sum((sd_cseq * AREA)^2,       na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    totalarea     = AREAtotal,
    avg_annualCH4 = total_annualCH4 / totalarea,
    avg_cseqmean  = total_cseqmean  / totalarea,
    avg_annualSE  = total_annualSE  / totalarea,
    avg_cseqSE    = total_cseqSE    / totalarea
  )

#Define the Conversion Function with SE Propagation
calc_GWPstar <- function(df) {
  df %>% 
    arrange(Year) %>%
    mutate(
      # Calculate GWPstar_CH4 and its SE
      GWPstar_CH4    = ((4 * avg_annualCH4) - (3.75 * lag(avg_annualCH4, n = 20))) * 27,
      GWPstar_CH4_SE = sqrt((27 * 4)^2 * avg_annualSE^2 + (27 * 3.75)^2 * (lag(avg_annualSE, n = 20))^2),
      
      # Calculate budget and its SE (difference between CH4 conversion and cseq)
      GWPstar_budget = GWPstar_CH4 - avg_cseqmean,
      GWPstar_budget_SE = sqrt(GWPstar_CH4_SE^2 + avg_cseqSE^2),
      
      # Calculate annual RF as the change in the budget
      GWPstar_RF     = GWPstar_budget - lag(GWPstar_budget),
      GWPstar_RF_SE  = sqrt(GWPstar_budget_SE^2 + (lag(GWPstar_budget_SE))^2),
      
      # Cumulative RF is the cumulative sum of annual RF; error is propagated as the square root
      cumGWPstar_RF  = cumsum(replace_na(GWPstar_RF, 0)),
      cumGWPstar_RF_SE = sqrt(cumsum(replace_na(GWPstar_RF_SE^2, 0)))
    ) %>%
    # Propagate NA for cumGWPstar_RF if GWPstar_RF is NA
    mutate(cumGWPstar_RF = ifelse(is.na(GWPstar_RF), NA, cumGWPstar_RF))
}

#Compute GWP* for (1984–2020)
historical_calc <- aggregated_annual_summary %>%
  filter(temporal == "historical") %>%
  calc_GWPstar() %>%
  filter(Year <= 2020)

#Compute GWP* for Each SSP Scenario
ssp_scenarios <- c("ssp126", "ssp245", "ssp370", "ssp585")

ssp_calc <- map_dfr(ssp_scenarios, function(ssp) {
  # Get SSP data (2021–2090) from aggregated summary
  ssp_data <- aggregated_annual_summary %>% filter(temporal == ssp)
  # Combine all historical data with SSP data so that the SSP series starts with historical years.
  combined <- bind_rows(
    aggregated_annual_summary %>% filter(temporal == "historical"),
    ssp_data
  ) %>% arrange(Year)
  
  # Compute GWP* metrics on the combined data
  combined <- calc_GWPstar(combined)
  # Retain data from 2020 onward to ensure the SSP line is connected 
  combined %>% 
    filter(Year >= 2020) %>% 
    mutate(temporal = ssp)
})

#Combine and Arrange the Final Dataset
GWPstar <- bind_rows(
  historical_calc,  # historical period (1984–2020)
  ssp_calc          # SSP series (2020–2090 with historical connection)
) %>%
  # Arrange so that historical (1984-2020) comes first, then SSP scenarios.
  mutate(temporal = factor(temporal, levels = c("historical", ssp_scenarios))) %>%
  arrange(temporal, Year) %>% 
  filter(!is.na(GWPstar_budget))

