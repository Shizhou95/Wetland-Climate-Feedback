# R code for creating Figure 3
# Code by: Shizhou Ma
# Version Date: October 2025

# This script creates Figure 3 using the `GWPstar` dataset generated from the
# GWPstarbudget script.
# It generates:
#   1. Annual GWP* CO2-equivalent budget trajectories (Figure 3a)
#   2. Cumulative GWP* CO2-equivalent budget trajectories (Figure 3b)
# Note:
#   - Run the GWPstarbudget script first to generate/load the `GWPstar` dataset.
# Required packages to run the code:
# dplyr V1.1.4, ggplot2 V3.5.2, broom V1.0.8, ggnewscale V0.5.1, ggpubr V0.6.0

# Define colors for each temporal group.
data_colors <- c(
  "historical" = "blue", 
  "ssp126"     = "darkgreen", 
  "ssp245"     = "orange", 
  "ssp370"     = "purple", 
  "ssp585"     = "red"
)

scenario_labels <- c(
  historical = "Historical",
  ssp126     = "SSP1-2.6",
  ssp245     = "SSP2-4.5",
  ssp370     = "SSP3-7.0",
  ssp585     = "SSP5-8.5"
)

#plot annual version (Figure 3a)
plot_budget <- ggplot(GWPstar, aes(x = Year, color = temporal, fill = temporal)) +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_ribbon(aes(ymin = GWPstar_budget - GWPstar_budget_SE,
                  ymax = GWPstar_budget + GWPstar_budget_SE),
              alpha = 0.3, color = NA) +
  geom_line(aes(y = GWPstar_budget), size = 1) +
  scale_color_manual(
    name   = "Scenario",
    values = data_colors,
    labels = scenario_labels
  ) +
  scale_fill_manual(
    name   = "Scenario",
    values = data_colors,
    labels = scenario_labels
  ) +
  scale_x_continuous(
    breaks = seq(1920, 2090, by = 10),
    labels = seq(1920, 2090, by = 10),
    limits = c(1921, 2090)
  ) +
  labs(
    y     = expression("GWP* " * CO[2] * "-e.q. budget (kg " * CO[2] * "-e.q./ha/yr)"),
    title = NULL
  ) +
  theme_classic2() +
  theme(
    legend.position = "bottom",
    axis.text       = element_text(size = 16),
    axis.title      = element_text(size = 18),
    strip.text      = element_text(size = 18),
    legend.text     = element_text(size = 16),
    legend.title    = element_text(size = 18)
  )

print(plot_budget)





#cumulative version (Figure 3b)
# Compute raw within-group cumsums (budget + variance)
GWPstar_cum_raw <- GWPstar %>%
  arrange(temporal, Year) %>%
  group_by(temporal) %>%
  mutate(
    cum_raw_budget = cumsum(GWPstar_budget),
    cum_raw_var    = cumsum(GWPstar_budget_SE^2)
  ) %>%
  ungroup()

# Pull out the final historical cumulative values
hist_end <- GWPstar_cum_raw %>%
  filter(temporal == "historical") %>%
  slice_tail(n = 1) %>%
  summarize(
    hist_budget_end = cum_raw_budget,
    hist_var_end    = cum_raw_var
  )

# Build the cumulative series
GWPstar_cum_chained <- GWPstar_cum_raw %>%
  mutate(
    cum_budget = if_else(
      temporal == "historical",
      cum_raw_budget,
      hist_end$hist_budget_end + cum_raw_budget
    ),
    cum_var = if_else(
      temporal == "historical",
      cum_raw_var,
      hist_end$hist_var_end + cum_raw_var
    ),
    cum_SE = sqrt(cum_var)
  )

# Plot
plot_budget_cum_chained <- ggplot(GWPstar_cum_chained, aes(x = Year, color = temporal, fill = temporal)) +
  #geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_ribbon(aes(
    ymin = cum_budget - cum_SE,
    ymax = cum_budget + cum_SE
  ), alpha = 0.3, color = NA) +
  geom_line(aes(y = cum_budget), size = 1) +
  scale_color_manual(values = data_colors) +
  scale_fill_manual(values  = data_colors) +
  scale_x_continuous(
    breaks = seq(1920, 2090, by = 10),
    labels = seq(1920, 2090, by = 10),
    limits = c(1921, 2090)
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(5),
    labels = scales::comma
  ) +
  labs(
    y     = expression("Cumulative GWP* " * CO[2] * "-e.q. budget (kg " * CO[2] * "-e.q./ha)"),
    title = "",
    color = "Scenario",
    fill  = "Scenario"
  ) +
  theme_classic2() +
  theme(
    legend.position = "none",
    axis.text       = element_text(size = 16),
    axis.title      = element_text(size = 16),
    strip.text      = element_text(size = 16),
    legend.text     = element_text(size = 16),
    legend.title    = element_text(size = 16)
  )

plot_budget_cum_chained
