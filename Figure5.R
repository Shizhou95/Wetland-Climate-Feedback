# Code for creating Figure 5
# Code by: Shizhou Ma
# Version Date: October 2025

# This script creates Figure 5 using the `GHGpertubmerged` dataset generated
# from the radiative forcing analysis site-level script.
# It generates spatial maps of cumulative radiative forcing outcomes for
# historical and future scenarios at individual wetland sites.
# Note:
#   - Run the radiative forcing analysis site-level script first to generate/load
#     the `GHGpertubmerged` dataset.
# Required packages to run the code:
# dplyr V1.1.4, ggplot2 V3.5.2, sf V1.0-20, ggnewscale V0.5.1

# For historical data: select the last row for each Wetland.ID
end_savings_hist <- GHGpertubmerged %>%
  filter(scenario == "hist") %>%
  group_by(Wetland.ID, scenario) %>%
  slice_tail(n = 1) %>%  # selects the last row per group
  ungroup()



end_savings_ssps <- GHGpertubmerged %>%
  filter(scenario %in% c("ssp126", "ssp245", "ssp370", "ssp585")) %>%  # Keep only SSP scenarios
  group_by(Wetland.ID, scenario) %>%
  arrange(time) %>%   # Ensure rows are sorted by time in each group
  slice(c(max(1, ceiling(n() / 2) - 25), n())) %>%  # Select the row 25 prior to the mid row and the end (last) row
  ungroup()


# Combine the results into one data frame
end_savings <- bind_rows(end_savings_hist, end_savings_ssps)

geocoordinates <- read.csv('/Users/shizhouma/Creed Lab Dropbox/Shizhou Ma/carbon budget analysis/JANallsites_latlong.csv')

end_savings <- end_savings %>%
  left_join(geocoordinates, by = "Wetland.ID")

end_savings <- end_savings %>%
  mutate(status = if_else(savings_cum > 0, "warmer", "cooler")) %>% 
  arrange(Wetland.ID)




end_savings <- end_savings %>%
  group_by(Wetland.ID, scenario) %>% 
  arrange(time) %>% 
  mutate(row_label = case_when(
    scenario == "hist" & row_number() == 1 ~ "historical (1900-2020)",
    scenario %in% c("ssp126", "ssp245", "ssp370", "ssp585") & row_number() == 1 ~ paste(scenario, "mid century"),
    scenario %in% c("ssp126", "ssp245", "ssp370", "ssp585") & row_number() == n() ~ paste(scenario, "end century"),
    TRUE ~ NA_character_
  )) %>%
  ungroup() %>% 
  arrange(Wetland.ID)

library(sf)
library(ggplot2)
library(dplyr)
library(ggnewscale)



# define scale factor
sf <- 1e12
# END CENTURY 
# Subset the data to include only historical and SSP end-century rows
plot_data <- end_savings %>%
  filter(!is.na(row_label)) %>%
  filter(row_label == "historical (1900-2020)" | grepl("end century", row_label)) %>%
  # apply scale factor
  mutate(
    net_RF         = net_RF * sf,
    totalbase_cum  = totalbase_cum * sf,
    savings_cum    = savings_cum * sf
  )

# Jitter the coordinates to avoid overlap
plot_data_jittered <- plot_data %>%
  mutate(
    Longitude = jitter(Longitude, amount = 0.3),
    Latitude  = jitter(Latitude, amount = 0.3)
  )

# Convert to sf and match CRS
points_sf <- st_as_sf(plot_data_jittered, coords = c("Longitude", "Latitude"), crs = 4326)
# Read the shapefile (Prairie Pothole Region boundary) (adjust the path if needed)
shapefile_path <- '/Users/shizhouma/Creed Lab Dropbox/Shizhou Ma/carbon budget analysis/pprshapefile'


# Transform the points to match the coordinate system of the shapefile
points_sf <- st_transform(points_sf, crs = st_crs(shape_data))
shape_data <- st_read(shapefile_path)
points_sf <- st_transform(points_sf, crs = st_crs(shape_data))

# Split into warmer/cooler sites
warmer_points <- points_sf %>% filter(status == "warmer")
cooler_points <- points_sf %>% filter(status == "cooler")

# Custom facet labels
custom_labels_end <- c(
  "historical (1900-2020)" = "Historical (1900–2020)",
  "ssp126 end century"     = "SSP1-2.6 (end century)",
  "ssp245 end century"     = "SSP2-4.5 (end century)",
  "ssp370 end century"     = "SSP3-7.0 (end century)",
  "ssp585 end century"     = "SSP5-8.5 (end century)"
)

# Plot
ggplot() +
  geom_sf(data = shape_data, fill = NA, color = "black") +
  facet_wrap(~ row_label, ncol = 3, labeller = labeller(row_label = custom_labels_end)) +
  geom_sf(data = warmer_points, aes(color = savings_cum), size = 2, alpha = 0.7) +
  scale_color_gradient(
    low    = "lightcoral", high   = "darkred",
    name   = "positive feedback\n(cumRF pW/m²)",
    limits = range(warmer_points$savings_cum, na.rm = TRUE)
  ) +
  new_scale_color() +
  geom_sf(data = cooler_points, aes(color = savings_cum), size = 2) +
  scale_color_gradient(
    low    = "darkblue", high   = "lightblue",
    name   = "negative feedback\n(cumRF pW/m²)",
    limits = range(cooler_points$savings_cum, na.rm = TRUE)
  ) +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  theme(
    legend.position = "right",
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    axis.title      = element_blank(),
    strip.text      = element_text(size = 18),
    legend.text     = element_text(size = 16),
    legend.title    = element_text(size = 18)
  )




